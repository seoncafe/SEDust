module p_sub
   ! Equilibrium-temperature solver and Guhathakurta-Draine (1989)
   ! stochastic-heating temperature-distribution solver.
   !
   ! Includes the following optimizations
   ! relative to the original (justified for use inside a 3D-RT loop
   ! that calls sed_solve millions of times):
   !
   !   - Logarithmic-trapezoid integration weights `dlnlam_int(:)` are
   !     precomputed once via `p_sub_setup(lambda)` instead of being
   !     rebuilt inside calc_Teq on every call. Saves ~3 nlambda
   !     log calls per equilibrium solve.
   !   - `log_lambda(:)` is also cached (used by the highest-bin
   !     correction in calc_P).
   !   - calc_Teq's inner sum is collapsed to a single SUM intrinsic;
   !     gfortran -O2 vectorizes it cleanly.
   !
   ! The numerical algorithm follows Guhathakurta & Draine (1989); outputs match
   ! to roundoff at the regression test points.

   use, intrinsic :: iso_fortran_env, only: real64
   use constants, only: wp
   use sed_mathlib,   only: interp, logadd
   ! The upward-transition rates are driven by PHOTONS, so the band they are
   ! integrated over is the band the field occupies -- not the span of whatever
   ! optics table the model happens to be tabulated on.  hardest_photon_energy
   ! is the library's one answer to "how hard does this field get"; calc_P is
   ! the fourth caller that needs it.
   use radfield,      only: hardest_photon_energy
   implicit none
   private
   public :: p_sub_setup, calc_Teq, calc_P

   real(wp), parameter :: PI     = 3.141592653589793238462643383279502884197d0
   real(wp), parameter :: FOURPI = PI * 4.0d0
   real(wp), parameter :: c     = 2.99792458d10   ! cm/s
   real(wp), parameter :: hp    = 6.62606957d-27  ! cm^2 g s^-1
   real(wp), parameter :: hc    = hp*c*1d4        ! erg um  (lambda is in um)

   ! Sampling of the highest-bin integral, which runs over every photon more
   ! energetic than the top enthalpy bin's gap.
   !
   ! The rule used to be a flat 51 points spanning from lambda(1) -- the short
   ! end of the OPTICS GRID -- to the gap wavelength.  That tied the resolution
   ! to the optics table instead of to the band being integrated.  On a grid
   ! stopping at the Lyman limit it happened to work, giving dln(lambda) =
   ! 0.0018-0.0238 across the bins that carry the correction; on a table
   ! carried to 1e-4 um the same 51 points spread over 912 times the range and
   ! left between 1 and 8 of them inside the band the field actually occupies.
   !
   ! Moving the lower limit to the field's own short end fixes where the band
   ! is, but a fixed count would still make the step follow the band's width,
   ! so the count follows it instead: enough points for a step of
   ! DLNWAV_TARGET, never fewer than NWAV_MIN, never more than NWAV_MAX.
   ! DLNWAV_TARGET is the coarsest step the old 51 points ever gave over the
   ! bands that carry the correction, so a narrow band is sampled at least as
   ! finely as before (NWAV_MIN reproduces it outright), and a band four
   ! decades wide -- which the old rule sampled at dln(lambda) = 0.19 -- is no
   ! longer starved.  The integrand, C_abs(lambda) * lambda^2 * J_lambda, is
   ! smooth over the illuminated band, so a bounded step is all it needs.
   integer,  parameter :: NWAV_MIN      = 51
   integer,  parameter :: NWAV_MAX      = 301
   real(wp), parameter :: DLNWAV_TARGET = 0.02_wp

   real(wp), allocatable :: dlnlam_int(:)
   real(wp), allocatable :: log_lambda(:)
   integer :: NLAM = 0

   ! calc_P matrix workspace, reused per thread across calls to avoid the
   ! nT x nT allocate/zero/deallocate on every call inside the RT loop.
   real(wp), allocatable, save :: Amat_ws(:,:), Bmat_ws(:,:)
   !$omp threadprivate(Amat_ws, Bmat_ws)

contains

   subroutine p_sub_setup(lambda)
      ! Precompute the wavelength integration weights and the cached
      ! log(lambda) array. Call once after the lambda grid is loaded
      ! (e.g., from the Q table) and before any calc_Teq / calc_P call.
      real(wp), intent(in) :: lambda(:)
      integer :: k
      NLAM = size(lambda)
      if (allocated(dlnlam_int)) deallocate(dlnlam_int, log_lambda)
      allocate(dlnlam_int(NLAM), log_lambda(NLAM))
      do k = 1, NLAM
         log_lambda(k) = log(lambda(k))
      end do
      dlnlam_int(1)    = 0.5_wp * (log_lambda(2) - log_lambda(1))
      dlnlam_int(NLAM) = 0.5_wp * (log_lambda(NLAM) - log_lambda(NLAM-1))
      do k = 2, NLAM-1
         dlnlam_int(k) = 0.5_wp * (log_lambda(k+1) - log_lambda(k-1))
      end do
   end subroutine p_sub_setup


   subroutine calc_Teq(lambda, Cabs, Jfield, T, kappB, Teq)
      ! Solve Cabs(lam)*Jfield(lam) integral = kappB(T_eq) for T_eq.
      ! kappB is the precomputed integral of Cabs * B_lam over lambda
      ! at every T grid point (built by sed_init from the same Cabs).
      !
      ! The root is found by log-log interpolation: kappB(T) is a smooth
      ! near-power law (~T^(4+beta)), so interpolating log T against
      ! log kappB is far more accurate than linear interpolation on the
      ! coarse NT grid. Linear interpolation left a ~0.1-0.15% sawtooth
      ! in Teq (a ~0.5% sawtooth in the equilibrium-grain emission)
      ! verified against Draine's single-grain dbdis reference SEDs.
      real(wp), intent(in)  :: lambda(:), Cabs(:), Jfield(:)
      real(wp), intent(in)  :: T(:), kappB(:)
      real(wp), intent(out) :: Teq
      real(wp) :: kappJ, lTeq
      integer  :: nT, j
      kappJ = sum(Cabs * Jfield * lambda * dlnlam_int)
      nT = size(T)
      if (kappJ <= kappB(1)) then
         Teq = T(1)
      else if (kappJ >= kappB(nT)) then
         Teq = T(nT)
      else
         ! local log-log interpolation (kappB is monotonic in T)
         do j = 2, nT
            if (kappB(j) >= kappJ) exit
         end do
         lTeq = log(T(j-1)) + (log(T(j)) - log(T(j-1))) * &
                (log(kappJ) - log(kappB(j-1))) / (log(kappB(j)) - log(kappB(j-1)))
         Teq = exp(lTeq)
      end if
   end subroutine calc_Teq


   subroutine calc_P(lambda, Cabs, Jfield, T, kappB, H, P, lnP, kappCMB)
      ! Guhathakurta & Draine (1989) matrix solver for the temperature
      ! distribution P(T) of a stochastically heated grain.
      ! Numerics follow the Guhathakurta & Draine recursion; only the
      ! `log_lambda` lookup avoids a recomputed log() call inside the
      ! highest-bin correction loop.
      real(wp), intent(in)    :: lambda(:), Cabs(:), Jfield(:)
      real(wp), intent(in)    :: T(:), kappB(:), H(:)
      real(wp), intent(inout) :: P(:), lnP(:)
      real(wp), optional, intent(in) :: kappCMB
      integer  :: nlambda, nT, k, i1, i2, j
      real(wp) :: ener, wavl, kappB_tot
      real(wp) :: Jwavl, delH, Cross
      real(wp) :: sumP
      real(wp) :: sumB
      integer  :: nwav
      real(wp) :: dlnwav, wav
      real(wp) :: u_hard, lam_hard, span

      nlambda = size(lambda)
      nT      = size(T)
      ! Hardest photon the FIELD carries, and the wavelength that goes with it.
      ! Every upward rate below is cut off there: a transition needing more
      ! energy than this has nothing to drive it.  A field that is zero
      ! everywhere returns 0, and the fallback then leaves the old grid-edge
      ! behavior rather than dividing by it.
      u_hard = hardest_photon_energy(lambda, Jfield)
      if (u_hard > 0.0_wp) then
         lam_hard = hc / u_hard
      else
         lam_hard = lambda(1)
         u_hard   = hc / lambda(1)
      end if
      ! Amat_ws/Bmat_ws are reused per thread across calls: (re)allocate only
      ! when unallocated or when nT changes, then zero-fill (the recursion
      ! rebuilds the whole matrix each call).
      if (.not. allocated(Amat_ws)) then
         allocate(Amat_ws(nT, nT), Bmat_ws(nT, nT))
      else if (size(Amat_ws, 1) /= nT) then
         deallocate(Amat_ws, Bmat_ws)
         allocate(Amat_ws(nT, nT), Bmat_ws(nT, nT))
      end if
      Amat_ws = 0d0
      Bmat_ws = 0d0

      ! Downward transitions
      do i1 = 2, nT
         if (present(kappCMB)) then
            kappB_tot = kappB(i1) - kappCMB
         else
            kappB_tot = kappB(i1)
         end if
         Amat_ws(i1-1, i1) = FOURPI / (H(i1) - H(i1-1)) * kappB_tot
      end do

      ! Upward transitions
      do i1 = 1, nT-1
         do i2 = i1+1, nT
            ener = H(i2) - H(i1)
            if (i2 == nT) then
               delH = H(i2) * 0.5_wp * log(H(i2)/H(i2-1))
            else
               delH = H(i2) * 0.5_wp * log(H(i2+1)/H(i2-1))
            end if
            wavl = hc / ener
            call interp(lambda, Jfield, wavl, Jwavl)
            call interp(lambda, Cabs,   wavl, Cross)
            ! No photon this hard exists, so the transition has no driver.
            ! Without the test `interp` clamps instead: asked below lambda(1) it
            ! returns Jfield(1), which on a grid ending at the field's own edge
            ! is the full J there rather than zero, and every gap wider than the
            ! hardest photon is then driven by a photon flux that is not in the
            ! field.  A grid carried past the field makes the clamp harmless,
            ! but the cutoff belongs to the physics, not to the tabulation.
            if (ener > u_hard) Jwavl = 0.0_wp
            Amat_ws(i2, i1) = FOURPI * Cross * hc * delH / ener**3 * Jwavl
            ! Highest-bin: add transitions to all energies above this bin, i.e.
            ! to every photon shorter than wavl.  The field carries none shorter
            ! than lam_hard, so the band is [lam_hard, wavl] -- and it is empty
            ! when the gap already exceeds the hardest photon.  Leaving that
            ! case to the loop would integrate backwards over an unilluminated
            ! band and subtract a rate.
            if (i2 == nT .and. ener < u_hard) then
               span   = log(wavl/lam_hard)
               nwav   = min(NWAV_MAX, max(NWAV_MIN, ceiling(span/DLNWAV_TARGET) + 1))
               dlnwav = span/(nwav-1)
               do k = 1, nwav
                  wav = lam_hard * exp((k-1)*dlnwav)
                  call interp(lambda, Jfield, wav, Jwavl)
                  call interp(lambda, Cabs,   wav, Cross)
                  Amat_ws(i2, i1) = Amat_ws(i2, i1) + &
                     FOURPI * Cross * wav**2 / hc * dlnwav * Jwavl
               end do
            end if
         end do
      end do

      sumB = 0_wp
      i2 = nT
      do i1 = 1, i2-1
         Bmat_ws(i2, i1) = Amat_ws(i2, i1)
         sumB = sumB + Bmat_ws(i2, i1)
      end do
      do i1 = 1, nT-1
         do i2 = nT-1, i1+1, -1
            Bmat_ws(i2, i1) = Bmat_ws(i2+1, i1) + Amat_ws(i2, i1)
            sumB = sumB + Bmat_ws(i2, i1)
         end do
      end do

      P   = 0_wp
      lnP = -1d100
      if (sumB > 0_wp) then
         ! Linear-P method.
         P(1) = 1.0_wp
         do j = 2, nT
            if (Amat_ws(j-1, j) > 0_wp) P(j) = sum(Bmat_ws(j, 1:j-1) * P(1:j-1)) / Amat_ws(j-1, j)
            if (P(j) > 1d50)         P(1:j) = P(1:j) / P(j)
         end do
         sumP = sum(P)
         if (sumP > 0_wp) then
            P   = P / sumP
            !--- P is zero at every temperature the grain does not reach, and
            !--- log(0) raises a divide-by-zero (fatal under -fpe0).  lnP is
            !--- only ever compared against lnP_crit, so the finite sentinel
            !--- set above already carries "unreachable"; leave it in place.
            where (P > 0_wp) lnP = log(P)
         end if
      end if
      ! Amat_ws/Bmat_ws are intentionally not deallocated: they are kept
      ! (per thread) for reuse on the next call.
   end subroutine calc_P

end module p_sub
