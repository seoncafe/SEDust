module p_sub
   ! Equilibrium-temperature solver and Guhathakurta-Draine (1989)
   ! stochastic-heating temperature-distribution solver, with the
   ! finite-bin upward transition rates of Draine & Li (2001).
   !
   ! Absorption and upward transitions use the same sampled field support and
   ! the same interpolant: Cabs * J_lambda * lambda linear in ln(lambda)
   ! between samples, which radiation_absorption integrates exactly (the
   ! logarithmic trapezoid).
   ! The wavelength grid must be positive and strictly increasing. A continuum
   ! field ends at its first and last significant samples; no half interval
   ! toward an exterior zero sample contributes photons beyond those ends.

   use, intrinsic :: iso_fortran_env, only: real64, int64
   use constants, only: wp
   use radfield, only: radiation_field_support, isolated_radiation_support_count, RADIATION_SUPPORT_ISOLATED
   implicit none
   private
   public :: p_sub_setup, calc_Teq, calc_P, radiation_absorption

   real(wp), parameter :: PI     = 3.141592653589793238462643383279502884197d0
   real(wp), parameter :: FOURPI = PI * 4.0d0
   real(wp), parameter :: c     = 2.99792458d10   ! cm/s
   real(wp), parameter :: hp    = 6.62606957d-27  ! cm^2 g s^-1
   real(wp), parameter :: hc    = hp*c*1d4        ! erg um  (lambda is in um)

   ! The absorbed power per unit photon energy on a field's support, built
   ! once for each calc_P call (absorbed_power_spectrum).
   type :: absorbed_power_t
      integer :: m = 0
      real(wp), allocatable :: E(:), y(:), slope(:), F0(:), F1(:)
   end type absorbed_power_t

   ! A photon energy E [erg] on that spectrum: its segment i, y there, and
   ! the absorbed power F0 = int_0^E g dE below it.
   type :: photon_energy_t
      real(wp) :: E = 0.0_wp, y = 0.0_wp, F0 = 0.0_wp
      integer  :: i = 0
   end type photon_energy_t

   ! 1/k, k = 1 .. 15, for the series of exponential_remainders
   real(wp), parameter :: INV_INTEGER(15) = [1.0_wp, 0.5_wp, 1.0_wp/3.0_wp, 0.25_wp, &
      0.2_wp, 1.0_wp/6.0_wp, 1.0_wp/7.0_wp, 0.125_wp, 1.0_wp/9.0_wp, 0.1_wp, &
      1.0_wp/11.0_wp, 1.0_wp/12.0_wp, 1.0_wp/13.0_wp, 1.0_wp/14.0_wp, 1.0_wp/15.0_wp]

   real(wp), allocatable :: cached_lambda(:), log_lambda(:)
   integer :: NLAM = 0

   ! calc_P matrix workspace, reused per thread across calls to avoid the
   ! nT x nT allocate/zero/deallocate on every call inside the RT loop.
   real(wp), allocatable, save :: Amat_ws(:,:), Bmat_ws(:,:)
   !$omp threadprivate(Amat_ws, Bmat_ws)

contains

   subroutine p_sub_setup(lambda)
      ! Cache logarithmic wavelengths after the model grid is loaded.
      ! The public absorption integral also accepts other grids and computes
      ! their logarithms locally without changing this shared read-only cache.
      real(wp), intent(in) :: lambda(:)
      integer :: k
      NLAM = size(lambda)
      if (NLAM < 2) error stop 'p_sub_setup: at least two wavelengths required'
      if (any(lambda <= 0.0_wp)) error stop 'p_sub_setup: nonpositive wavelength'
      if (any(lambda(2:) <= lambda(:NLAM-1))) error stop 'p_sub_setup: unordered wavelengths'
      if (allocated(cached_lambda)) deallocate(cached_lambda, log_lambda)
      allocate(cached_lambda(NLAM), log_lambda(NLAM))
      cached_lambda = lambda
      do k = 1, NLAM
         log_lambda(k) = log(lambda(k))
      end do
   end subroutine p_sub_setup


   subroutine radiation_support(lambda, Jfield, first, last)
      ! The continuum definition is shared with hardest_photon_energy.
      ! Logarithmic integration additionally requires this grid ordering.
      real(wp), intent(in) :: lambda(:), Jfield(:)
      integer, intent(out) :: first, last
      integer :: support_status
      if (size(lambda) /= size(Jfield)) error stop 'radiation_support: size mismatch'
      if (any(lambda <= 0.0_wp)) error stop 'radiation_support: nonpositive wavelength'
      if (size(lambda) > 1) then
         if (any(lambda(2:) <= lambda(:size(lambda)-1))) &
            error stop 'radiation_support: unordered wavelengths'
      end if
      call radiation_field_support(lambda, Jfield, first, last, support_status)
      if (support_status == RADIATION_SUPPORT_ISOLATED) then
         !$omp atomic update
         isolated_radiation_support_count = isolated_radiation_support_count + 1_int64
      end if
   end subroutine radiation_support


   function radiation_absorption(lambda, Cabs, Jfield) result(absorption)
      ! Integral of Cabs * J_lambda d(lambda) by the logarithmic trapezoid
      ! rule, restricted to the sampled continuum support. No 4*pi or unit
      ! conversion is included: lambda and J_lambda must use compatible units.
      ! A host can use this same function for its absorbed-power normalization.
      ! No p_sub_setup call is required for this public integral.
      real(wp), intent(in) :: lambda(:), Cabs(:), Jfield(:)
      real(wp) :: absorption, width, left, right, loglam(size(lambda))
      logical :: cache_matches
      integer :: first, last, k
      if (size(Cabs) /= size(lambda)) error stop 'radiation_absorption: size mismatch'
      call radiation_support(lambda, Jfield, first, last)
      absorption = 0.0_wp
      if (first == 0) return
      cache_matches = .false.
      if (allocated(cached_lambda)) then
         if (size(cached_lambda) == size(lambda)) cache_matches = all(cached_lambda == lambda)
      end if
      if (cache_matches) then
         loglam = log_lambda
      else
         loglam = log(lambda)
      end if
      left = Cabs(first) * Jfield(first) * lambda(first)
      do k = first+1, last
         right = Cabs(k) * Jfield(k) * lambda(k)
         width = loglam(k) - loglam(k-1)
         absorption = absorption + 0.5_wp * (left + right) * width
         left = right
      end do
   end function radiation_absorption


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
      kappJ = radiation_absorption(lambda, Cabs, Jfield)
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
      ! Guhathakurta & Draine (1989, GD89) matrix solver for the temperature
      ! distribution P(T) of a stochastically heated grain, on the enthalpy
      ! grid H(T) (erg), with the GD89 linear recursion for the steady state.
      !
      ! Downward rates are GD89's continuous cooling,
      !    A(u-1, u) = 4 pi [kappB(u) - kappCMB] / (H(u) - H(u-1)).
      ! Upward rates are Draine & Li (2001, ApJ 551, 807; DL01) eqs. (15)-(24)
      ! for l -> u < nT, eq. (16) for the top bin nT, and eq. (28) for the
      ! absorptions that stay inside bin l, added to l -> l+1; see
      ! upward_absorption_rates.  A grain is taken to be spread uniformly over
      ! bin l, and a photon of energy E moves it into bin u in proportion to
      ! the overlap of [H_l,min + E, H_l,max + E] with bin u.  For bins
      ! narrow compared with the photon energies and with the scale on which
      ! Cabs J varies, these reduce to DL01 eq. (26), i.e. GD89 eq. (2.5),
      !    A(u, l) = 4 pi Cabs J_lambda hc Delta H_u / E^3,  E = H_u - H_l,
      ! the point-sampled form used here until 2026-10-21.  That limit fails
      ! wherever a bin is broader than the band the field occupies: sampled
      ! at one photon energy per bin pair, a bin's absorbed power then misses
      ! part of the band or counts part of it twice.
      !
      ! Validity: the rates hold at any bin width.  Summed over u, every bin
      ! l < nT absorbs the power 4 pi radiation_absorption(lambda, Cabs,
      ! Jfield) to rounding; photons that would carry the grain past the
      ! lower edge of the top bin are assigned to the top bin (eq. 16, which
      ! counts their energy, where GD89 eq. (2.6) counted them as photons).
      ! The top bin has no upward transition, so the enthalpy window must
      ! leave P(nT) negligible.  The steady state then emits
      ! 4 pi sum_(u>=2) P(u) [kappB(u) - kappCMB] = (1 - P(nT)) times that
      ! power, an identity of the recursion.  What finite bins still
      ! approximate is the uniform spread within a bin and the cooling
      ! difference quotient, so the emitted spectrum's shape must be checked
      ! by refining the enthalpy grid.  A field that absorbs nothing leaves
      ! P(1) = 1.
      real(wp), intent(in)    :: lambda(:), Cabs(:), Jfield(:)
      real(wp), intent(in)    :: T(:), kappB(:), H(:)
      real(wp), intent(inout) :: P(:), lnP(:)
      real(wp), optional, intent(in) :: kappCMB
      integer  :: nT, i1, i2, j, first, last
      real(wp) :: kappB_tot
      real(wp) :: sumP
      real(wp) :: sumB

      nT      = size(T)
      call radiation_support(lambda, Jfield, first, last)
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

      ! Upward transitions (DL01 eqs. 15-16, 28)
      call upward_absorption_rates(lambda, Cabs, Jfield, first, last, H, nT, Amat_ws)

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
      else
         ! Without absorbed photons the normalized stationary state is the
         ! lowest enthalpy bin, not an identically zero probability vector.
         P(1) = 1.0_wp
         lnP(1) = 0.0_wp
      end if
      ! Amat_ws/Bmat_ws are intentionally not deallocated: they are kept
      ! (per thread) for reuse on the next call.
   end subroutine calc_P

   subroutine upward_absorption_rates(lambda, Cabs, Jfield, first, last, H, nT, Amat)
      ! Photon-absorption rates A(u, l) [s^-1] from enthalpy bin l to every
      ! higher bin u, Draine & Li (2001, ApJ 551, 807) eqs. (15)-(24), (16)
      ! and (28).  Amat(u, l), u > l, must be zero on entry.
      !
      ! Bin u is [E_u,min, E_u,max] = [edge(u-1), edge(u)], with each inner
      ! edge at the geometric mean sqrt(H_u H_(u+1)) of adjacent enthalpies and
      ! the two outer edges one half step beyond the end points, so the width
      ! dH_u agrees with GD89's H_u ln(H_(u+1)/H_(u-1))/2 to second order in
      ! the step.  With g(E) dE = 4 pi Cabs J_lambda dlambda, the power
      ! absorbed per unit photon energy, eq. (15) is
      !    A(u, l) = [int Phi_ul(E) g(E) dE] / (H_u - H_l),
      ! where Phi_ul(E) = dH_u G_ul(E) is the fraction of a grain spread
      ! uniformly over bin l that a photon E moves into bin u (eqs. 17-24):
      ! (E - W1)/dH_l on (W1, W2), min(dH_l, dH_u)/dH_l on (W2, W3) and
      ! (W4 - E)/dH_l on (W3, W4), with W1 = E_u,min - E_l,max,
      ! W2 = W1 + min(dH_l, dH_u), W3 = W1 + max(dH_l, dH_u) and
      ! W4 = W1 + dH_l + dH_u.  In terms of A_k = edge(k) - edge(l) and
      ! B_k = edge(k) - edge(l-1), W1 = A_(u-1), W4 = B_u and {W2, W3} =
      ! {B_(u-1), A_u}, so stepping u adds the two photon energies A_u, B_u.
      !
      ! Eq. (16): the top bin nT takes every photon that carries the grain past
      ! its lower edge, with weight (E - W1)/dH_l up to Wc = W1 + dH_l and 1
      ! beyond.  Eq. (28): photons E < dH_l leave a fraction 1 - E/dH_l of bin
      ! l inside bin l; that power is added to l -> l+1.  Summed over u these
      ! fractions are 1 at every E > 0, so each bin l < nT absorbs the whole
      ! power, 4 pi radiation_absorption(lambda, Cabs, Jfield).
      !
      ! The rising and falling edges are integrated in closed form from their
      ! start over their width (absorption_moments), not as differences of
      ! the cumulants F0 and F1 at their ends: a bin can be many orders of
      ! magnitude narrower than the photon energies, and the moment
      ! int (E - W) g dE taken as F1 - W F0 at E ~ W would lose the digits
      ! of E/dH.  The plateau and the top-bin tail are differences of F0,
      ! whose rounding is a fraction of order 1e-16 of the bin's whole
      ! absorption.
      real(wp), intent(in)    :: lambda(:), Cabs(:), Jfield(:), H(:)
      integer,  intent(in)    :: first, last, nT
      real(wp), intent(inout) :: Amat(:,:)
      type(absorbed_power_t) :: g
      type(photon_energy_t)  :: A_lo, B_lo, A_hi, B_hi, zero, hardest
      real(wp) :: edge(0:max(nT, 1)), dH(max(nT, 1))
      real(wp) :: dH_l, dH_u, rise, power, d0, r1
      integer  :: l, u

      if (first == 0 .or. last <= first .or. nT < 2) return
      call absorbed_power_spectrum(lambda, Cabs, Jfield, first, last, g)

      edge(0) = H(1) * sqrt(H(1) / H(2))
      do l = 1, nT-1
         edge(l) = sqrt(H(l) * H(l+1))
      end do
      edge(nT) = H(nT) * sqrt(H(nT) / H(nT-1))
      ! The widths edge(l) - edge(l-1) from the enthalpy differences, which
      ! are exact, rather than from the difference of two edges, which loses
      ! the digits of H/dH on a narrow window.
      dH(1) = sqrt(H(1) / H(2)) * (H(2) - H(1))
      do l = 2, nT-1
         dH(l) = sqrt(H(l)) * (H(l+1) - H(l-1)) / (sqrt(H(l+1)) + sqrt(H(l-1)))
      end do
      dH(nT) = sqrt(H(nT) / H(nT-1)) * (H(nT) - H(nT-1))

      zero    = photon_energy(g, 0.0_wp, 1)
      hardest = photon_energy(g, g%E(g%m), g%m)
      B_lo%i = 1
      do l = 1, nT-1
         dH_l = dH(l)
         A_lo = zero
         B_lo = photon_energy(g, dH_l, B_lo%i)
         ! eq. (28): weight dH_l - E on (0, dH_l)
         call absorption_moments(g, zero, B_lo, dH_l, d0, r1)
         Amat(l+1, l) = (dH_l * d0 - r1) / (dH_l * (H(l+1) - H(l)))
         do u = l+1, nT
            if (A_lo%i == g%m) exit
            if (u == nT) then
               ! eq. (16): weight E - W1 on (W1, Wc), dH_l on (Wc, E_hard),
               ! with Wc = W1 + dH_l = B_(nT-1)
               call absorption_moments(g, A_lo, B_lo, dH_l, d0, r1)
               power = r1 + dH_l * (hardest%F0 - B_lo%F0)
            else
               dH_u = dH(u)
               rise = min(dH_l, dH_u)
               A_hi = photon_energy(g, edge(u) - edge(l), A_lo%i)
               B_hi = photon_energy(g, edge(u) - edge(l-1), B_lo%i)
               power = 0.0_wp
               if (B_hi%i > 0) then
                  ! rising edge (A_(u-1), W2), plateau (W2, W3), falling
                  ! edge (W3, B_u)
                  if (B_lo%E <= A_hi%E) then
                     call absorption_moments(g, A_lo, B_lo, rise, d0, r1)
                     power = r1 + rise * max(A_hi%F0 - B_lo%F0, 0.0_wp)
                     call absorption_moments(g, A_hi, B_hi, rise, d0, r1)
                  else
                     call absorption_moments(g, A_lo, A_hi, rise, d0, r1)
                     power = r1 + rise * max(B_lo%F0 - A_hi%F0, 0.0_wp)
                     call absorption_moments(g, B_lo, B_hi, rise, d0, r1)
                  end if
                  power = power + (rise * d0 - r1)
               end if
               A_lo = A_hi
               B_lo = B_hi
            end if
            Amat(u, l) = Amat(u, l) + max(power, 0.0_wp) / (dH_l * (H(u) - H(l)))
         end do
      end do
   end subroutine upward_absorption_rates


   subroutine absorbed_power_spectrum(lambda, Cabs, Jfield, first, last, g)
      ! The absorbed power per unit photon energy on the field support, in
      ! ascending photon energy E_n = hc/lambda: g(E) dE = 4 pi y d(ln E),
      ! y = Cabs J_lambda lambda, linear in ln(E) between samples -- the
      ! interpolant of radiation_absorption's logarithmic trapezoid.  The
      ! segment widths tau_n in ln(E) are taken from the sample energies
      ! themselves, ln(1 + relative step), the scale on which every
      ! photon energy inside a segment is placed (photon_energy); the total
      ! power then agrees with radiation_absorption to the rounding of its
      ! widths.  Stores the slopes dy/d(ln E) and the cumulants
      ! F0(E_n) = int_0^E_n g dE and F1(E_n) = int_0^E_n E g dE.
      real(wp), intent(in) :: lambda(:), Cabs(:), Jfield(:)
      integer,  intent(in) :: first, last
      type(absorbed_power_t), intent(out) :: g
      real(wp) :: tau, phi, psi
      integer  :: n, k
      g%m = last - first + 1
      allocate(g%E(g%m), g%y(g%m), g%slope(g%m), g%F0(g%m), g%F1(g%m))
      do n = 1, g%m
         k = last - n + 1
         g%E(n) = hc / lambda(k)
         g%y(n) = Cabs(k) * Jfield(k) * lambda(k)
      end do
      g%F0(1) = 0.0_wp
      g%F1(1) = 0.0_wp
      g%slope(g%m) = 0.0_wp
      do n = 1, g%m - 1
         tau = log_one_plus((g%E(n+1) - g%E(n)) / g%E(n))
         g%slope(n) = (g%y(n+1) - g%y(n)) / tau
         call exponential_remainders(tau, phi, psi)
         g%F0(n+1) = g%F0(n) + FOURPI * (g%y(n) + 0.5_wp * g%slope(n) * tau) * tau
         g%F1(n+1) = g%F1(n) + FOURPI * g%E(n) * &
            (g%y(n) * (phi + tau) + g%slope(n) * (psi + 0.5_wp * tau * tau))
      end do
   end subroutine absorbed_power_spectrum


   pure function photon_energy(g, x, hint) result(pt)
      ! Photon energy x [erg] placed on the absorbed-power spectrum: its
      ! segment i (E(i) <= x < E(i+1); 0 below the softest sample, m at or
      ! above the hardest), y there, and F0(x).  The search steps from hint,
      ! since successive energies lie close together.
      type(absorbed_power_t), intent(in) :: g
      real(wp), intent(in) :: x
      integer,  intent(in) :: hint
      type(photon_energy_t) :: pt
      real(wp) :: t
      integer  :: i
      pt%E = x
      if (x < g%E(1)) then
         pt%i = 0
         pt%y = 0.0_wp
         pt%F0 = 0.0_wp
         return
      end if
      if (x >= g%E(g%m)) then
         pt%i = g%m
         pt%y = g%y(g%m)
         pt%F0 = g%F0(g%m)
         return
      end if
      i = min(max(hint, 1), g%m - 1)
      do while (x < g%E(i))
         i = i - 1
      end do
      do while (x >= g%E(i+1))
         i = i + 1
      end do
      pt%i = i
      t = 0.0_wp
      if (x > g%E(i)) t = log_one_plus((x - g%E(i)) / g%E(i))
      pt%y = g%y(i) + g%slope(i) * t
      pt%F0 = g%F0(i) + FOURPI * (g%y(i) + 0.5_wp * g%slope(i) * t) * t
   end function photon_energy


   subroutine absorption_moments(g, p, q, w, d0, r1)
      ! d0 = int g dE and r1 = int (E - p) g dE over the photon energies
      ! (p, q), q - p = w >= 0, the part outside the field support
      ! contributing nothing.  The interval is split at the samples it
      ! crosses: partial segments are integrated in closed form from their
      ! start and width, whole segments from the cumulants.  q only locates
      ! the end; the width is w, so the accuracy of a narrow interval is that
      ! of w, not that of the difference of two photon energies.
      type(absorbed_power_t), intent(in) :: g
      type(photon_energy_t),  intent(in) :: p, q
      real(wp), intent(in)  :: w
      real(wp), intent(out) :: d0, r1
      real(wp) :: start, y_start, offset, width, piece0, piece1
      integer  :: i, j

      d0 = 0.0_wp
      r1 = 0.0_wp
      if (w <= 0.0_wp .or. p%i == g%m .or. q%i == 0) return
      start = p%E
      y_start = p%y
      i = p%i
      offset = 0.0_wp
      if (i == 0) then
         ! below the softest field photon the integrand is zero
         offset = g%E(1) - p%E
         start = g%E(1)
         y_start = g%y(1)
         i = 1
      end if
      width = w - offset
      j = q%i
      if (j == g%m) then
         ! above the hardest field photon the integrand is zero
         j = g%m - 1
         width = g%E(g%m) - start
      end if

      if (j <= i) then
         call segment_moments(g, i, start, y_start, max(width, 0.0_wp), d0, r1)
      else
         call segment_moments(g, i, start, y_start, g%E(i+1) - start, d0, r1)
         if (j > i + 1) then
            ! whole segments i+1 .. j-1
            piece0 = g%F0(j) - g%F0(i+1)
            r1 = r1 + ((g%F1(j) - g%F1(i+1)) - g%E(i+1) * piece0) + (g%E(i+1) - start) * piece0
            d0 = d0 + piece0
         end if
         call segment_moments(g, j, g%E(j), g%y(j), max(width - (g%E(j) - start), 0.0_wp), &
                              piece0, piece1)
         r1 = r1 + piece1 + (g%E(j) - start) * piece0
         d0 = d0 + piece0
      end if
      r1 = r1 + offset * d0
   end subroutine absorption_moments


   pure subroutine segment_moments(g, i, start, y_start, width, d0, r1)
      ! int g dE and int (E - start) g dE over (start, start + width) inside
      ! segment i, in closed form for y linear in ln(E): with
      ! t = ln(1 + width/start), y_s = y at start and k the segment slope,
      !    d0 = 4 pi (y_s + k t/2) t,
      !    r1 = 4 pi start [y_s (e^t - 1 - t) + k (t e^t - (e^t - 1) - t^2/2)].
      type(absorbed_power_t), intent(in) :: g
      integer,  intent(in)  :: i
      real(wp), intent(in)  :: start, y_start, width
      real(wp), intent(out) :: d0, r1
      real(wp) :: t, phi, psi
      t = log_one_plus(width / start)
      call exponential_remainders(t, phi, psi)
      d0 = FOURPI * (y_start + 0.5_wp * g%slope(i) * t) * t
      r1 = FOURPI * start * (y_start * phi + g%slope(i) * psi)
   end subroutine segment_moments


   pure function log_one_plus(x) result(r)
      ! ln(1 + x) for x >= 0, as ln(u) x / (u - 1) with u = 1 + x, which
      ! cancels the rounding of u (accurate to a few ulps at any x).
      real(wp), intent(in) :: x
      real(wp) :: r, u
      u = 1.0_wp + x
      if (u == 1.0_wp) then
         r = x
      else
         r = log(u) * x / (u - 1.0_wp)
      end if
   end function log_one_plus


   pure subroutine exponential_remainders(t, phi, psi)
      ! phi = e^t - 1 - t = sum_(k>=2) t^k/k! and
      ! psi = t e^t - (e^t - 1) - t^2/2 = sum_(k>=3) (k-1) t^k/k!  for t >= 0.
      ! Below t = 0.3 both are summed as series, truncated where the next
      ! term is under 1e-17 of the sum; above it the closed forms lose to
      ! cancellation at most a factor 40 in relative rounding (psi at 0.3).
      real(wp), intent(in)  :: t
      real(wp), intent(out) :: phi, psi
      real(wp) :: term, e
      integer  :: k, n
      if (t < 0.3_wp) then
         ! highest power kept; psi, whose sum starts at t^3, sets it
         if (t < 1.0e-4_wp) then
            n = 6
         else if (t < 1.0e-3_wp) then
            n = 7
         else if (t < 1.0e-2_wp) then
            n = 9
         else if (t < 0.1_wp) then
            n = 12
         else
            n = 15
         end if
         term = 0.5_wp * t * t
         phi = term
         psi = 0.0_wp
         do k = 3, n
            term = term * t * INV_INTEGER(k)
            phi = phi + term
            psi = psi + real(k - 1, wp) * term
         end do
      else
         e = exp(t)
         phi = e - 1.0_wp - t
         psi = t * e - (e - 1.0_wp) - 0.5_wp * t * t
      end if
   end subroutine exponential_remainders

end module p_sub
