module radfield
   ! Radiation-field constructors and Planck-function utilities for the
   ! astrodust SED solver.
   ! A prepared wavelength grid reuses Planck coefficients across temperatures.

   use, intrinsic :: iso_fortran_env, only: real64, int64
   use constants, only: wp
   implicit none
   private
   public :: J_Mathis, J_diluted_bbody, calc_bbody, bbody
   public :: planck_grid_t, prepare_planck_grid
   public :: hardest_photon_energy, radiation_field_support, J_REL_FLOOR
   public :: isolated_radiation_support_count, RADIATION_SUPPORT_OK, &
             RADIATION_SUPPORT_EMPTY, RADIATION_SUPPORT_ISOLATED
   public :: cmb_temperature
   ! Toggle for Draine's 2008.02.02 correction of the Mathis 4000K
   ! dilution factor (1e-13 -> 1.65e-13) and the modern CMB temperature
   ! (2.9 -> 2.725 K). Default .true. follows Draine
   ! textbook (2011) Table 12.1; .false. recovers the literal Mathis
   ! 1983 values. The corrected field deposits more optical/NIR power,
   ! which raises the equilibrium-grain FIR peak relative to the literal
   ! field (see the memory note fir-15pct-status for the comparison
   ! against HD23 astrodust_irem.dat, and the open questions there).
   public :: use_mathis_corrected

   real(wp), parameter :: c     = 2.99792458d8     ! m/s
   real(wp), parameter :: h     = 6.62606957d-34   ! J.s
   real(wp), parameter :: kB    = 1.3806488d-23    ! J/K
   real(wp), parameter :: hc2   = 2.0d0*h*c**2
   real(wp), parameter :: hc_kB = h*c/kB
   ! h*c in [erg um], written with the same literals and the same order of
   ! operations as HC_ERG_UM in sed_astrodust_mod so that the two agree to the
   ! last bit; the single-photon bound below is compared against a 13.6 eV
   ! constant, and a last-place disagreement there would show up as a spurious
   ! difference in the grains that sit exactly at the switch-over.
   real(wp), parameter :: hc_erg_um = 6.62606957e-27_wp * 2.99792458e10_wp * 1.0e4_wp

   ! Relative intensity threshold for the outer field-support bounds.
   ! Interior samples are retained; see hardest_photon_energy.
   real(wp), parameter :: J_REL_FLOOR = 1.0e-12_wp

   integer, parameter :: RADIATION_SUPPORT_OK = 0
   integer, parameter :: RADIATION_SUPPORT_EMPTY = 1
   integer, parameter :: RADIATION_SUPPORT_ISOLATED = 2
   ! Counts absorption/transition checks rejected for zero bandwidth, not fields
   ! or grains. Read this cumulative diagnostic before/after a solve to detect
   ! an unresolved single-sample field. Pure support/photon queries do not
   ! increment it. Do not reset it during active calls.
   integer(int64), save :: isolated_radiation_support_count = 0_int64

   logical, save :: use_mathis_corrected = .true.

   ! Coefficients shared by all temperatures on one wavelength grid. The
   ! components are private so a prepared grid remains internally consistent.
   type :: planck_grid_t
      private
      real(wp), allocatable :: lambda_m(:), numerator(:)
      logical :: ascending_positive = .false.
   end type planck_grid_t

contains

   pure function cmb_temperature() result(T_cmb)
      ! The CMB temperature this module puts into the radiation field, and the
      ! ONLY place that value is written down.
      !
      ! It has to be a single source of truth because two very different pieces
      ! of the calculation must agree on it. J_Mathis adds a CMB blackbody to
      ! the field the grains absorb, and build_kappCMB integrates C_abs against
      ! a CMB blackbody past 1 mm to form the term calc_P subtracts from the
      ! grain's own emission (kappB - kappCMB is the NET cooling rate, which is
      ! what stops a grain cooling below its surroundings). If the two used
      ! different temperatures the solver would hold grains up against photons
      ! the field never supplied. They did differ -- 2.725 K in the field
      ! against a hard-coded 2.9 K in the cooling term -- once
      ! use_mathis_corrected became the default.
      !
      ! use_mathis_corrected = .true. is the modern value (Mather et al.);
      ! .false. recovers the literal Mathis (1983) 2.9 K together with that
      ! paper's 4000 K dilution factor, so the two move as a pair.
      real(wp) :: T_cmb
      if (use_mathis_corrected) then
         T_cmb = 2.725_wp
      else
         T_cmb = 2.9_wp
      end if
   end function cmb_temperature


   subroutine J_Mathis(U, lambda, J)
      ! Mathis 1983 ISRF, scaled by intensity factor U, plus a CMB
      ! blackbody at long wavelength. With use_mathis_corrected=.true.
      ! the 4000K dilution factor uses Draine's corrected 1.65e-13
      ! (Draine 2008.02.02) and CMB at 2.725 K
      ! (modern Mather et al. value) instead of our historical 1e-13 /
      ! 2.9 K. CMB is added unscaled by U (matches HD23 convention; see
      ! note below).
      real(wp), intent(in)    :: U
      real(wp), intent(in)    :: lambda(:)     ! [um]
      real(wp), intent(inout) :: J(:)
      integer  :: nlambda, i
      real(wp) :: w_4000, T_cmb
      if (use_mathis_corrected) then
         w_4000 = 1.65d-13
      else
         w_4000 = 1.0d-13
      end if
      T_cmb = cmb_temperature()
      nlambda = size(lambda)
      do i = 1, nlambda
         if (lambda(i) < 0.0912d0) then
            J(i) = 0.0d0
         else if (lambda(i) < 0.110d0) then
            J(i) = 3069d0 * lambda(i)**3.4172d0
         else if (lambda(i) < 0.134d0) then
            J(i) = 1.627d0
         else if (lambda(i) < 0.250d0) then
            J(i) = 0.0566d0 * lambda(i)**(-1.6678d0)
         else
            J(i) =   1d-14  * bbody(7500d0, lambda(i)) &
                  + w_4000  * bbody(4000d0, lambda(i)) &
                  + 4d-13   * bbody(3000d0, lambda(i))
         end if
      end do
      J = J * U
      do i = 1, nlambda
         J(i) = J(i) + bbody(T_cmb, lambda(i))
      end do
   end subroutine J_Mathis


   subroutine J_diluted_bbody(itype, lambda, J)
      integer,  intent(in)    :: itype
      real(wp), intent(in)    :: lambda(:)
      real(wp), intent(inout) :: J(:)
      real(wp), parameter :: T(6) = [3000d0, 6000d0, 9000d0, 12000d0, 15000d0, 18000d0]
      real(wp), parameter :: f(6) = [8.28d-12, 2.23d-13, 2.99d-14, 7.23d-15, 2.36d-15, 9.42d-16]
      integer :: i
      do i = 1, size(lambda)
         J(i) = bbody(T(itype), lambda(i)) * f(itype)
      end do
   end subroutine J_diluted_bbody


   subroutine prepare_planck_grid(lambda_um, grid)
      real(wp), intent(in) :: lambda_um(:)
      type(planck_grid_t), intent(out) :: grid
      integer :: i, n

      n = size(lambda_um)
      allocate(grid%lambda_m(n), grid%numerator(n))
      grid%ascending_positive = all(lambda_um > 0.0_wp)
      do i = 2, n
         if (lambda_um(i) < lambda_um(i-1)) grid%ascending_positive = .false.
      end do
      do i = 1, n
         grid%lambda_m(i) = lambda_um(i) * 1.0e-6_wp
         grid%numerator(i) = 0.0_wp
         if (lambda_um(i) > 0.0_wp) &
            grid%numerator(i) = hc2 / grid%lambda_m(i)**5
      end do
   end subroutine prepare_planck_grid


   subroutine calc_bbody(T, lambda_um, spec, grid)
      real(wp), intent(in)  :: T, lambda_um(:)
      real(wp), intent(out) :: spec(:)
      type(planck_grid_t), optional, intent(in) :: grid
      integer :: i, first, low, high, mid, n
      real(wp) :: x

      if (.not. present(grid)) then
         do i = 1, size(lambda_um)
            spec(i) = bbody(T, lambda_um(i))
         end do
         return
      end if
      if (.not. allocated(grid%lambda_m) .or. &
          .not. allocated(grid%numerator)) &
         error stop 'calc_bbody: unprepared Planck grid'
      n = size(lambda_um)
      if (size(grid%lambda_m) /= n .or. size(grid%numerator) /= n .or. &
          size(spec) < n) error stop 'calc_bbody: grid size mismatch'
      if (T <= 0.0_wp) then
         spec(1:n) = 0.0_wp
         return
      end if

      if (grid%ascending_positive) then
         ! On an ascending positive grid x decreases with wavelength. Find the
         ! first x < 700 in logarithmic time, keeping bbody's x arithmetic.
         low = 1
         high = n + 1
         do while (low < high)
            mid = low + (high - low)/2
            x = hc_kB / (T*grid%lambda_m(mid))
            if (x >= 700.0_wp) then
               low = mid + 1
            else
               high = mid
            end if
         end do
         first = low
         if (first > 1) spec(1:first-1) = 0.0_wp
         do i = first, n
            x = hc_kB / (T*grid%lambda_m(i))
            if (x < 1.0e-4_wp) then
               spec(i) = grid%numerator(i) / (2.0_wp * exp(0.5_wp*x) * sinh(0.5_wp*x))
            else
               spec(i) = grid%numerator(i) / (exp(x) - 1.0_wp)
            end if
         end do
      else
         spec(1:n) = 0.0_wp
         do i = 1, n
            if (lambda_um(i) <= 0.0_wp) cycle
            x = hc_kB / (T*grid%lambda_m(i))
            if (x >= 700.0_wp) cycle
            if (x < 1.0e-4_wp) then
               spec(i) = grid%numerator(i) / (2.0_wp * exp(0.5_wp*x) * sinh(0.5_wp*x))
            else
               spec(i) = grid%numerator(i) / (exp(x) - 1.0_wp)
            end if
         end do
      end if
   end subroutine calc_bbody

   pure function bbody(T, lambda_um) result(B)
      ! Planck function B_lambda(T), evaluated stably across the tail. The
      ! three branches keep the production range bit-identical while removing
      ! the over/underflow and the cancellation at the two extremes:
      !   x >= 700  : exp(x) would overflow; the Planck tail is 0 there.
      !   x <  1e-4 : use the exact identity exp(x)-1 = 2 e^{x/2} sinh(x/2),
      !               which has no cancellation as x -> 0.
      !   otherwise : the original exp(x)-1 form. The shipping wavelength/
      !               temperature grid has smallest x ~ 9e-5, so its few
      !               extreme points fall in the sinh branch; the value there
      !               moves by ~1e-12 relative (the removed cancellation
      !               error), below es-format output precision.
      real(wp), intent(in) :: T, lambda_um
      real(wp) :: B, lambda_m, x
      if (T <= 0.0_wp .or. lambda_um <= 0.0_wp) then
         B = 0.0_wp
         return
      end if
      lambda_m = lambda_um * 1.0e-6_wp
      x = hc_kB / (T*lambda_m)
      if (x >= 700.0_wp) then
         B = 0.0_wp
      else if (x < 1.0e-4_wp) then
         B = hc2 / lambda_m**5 / (2.0_wp * exp(0.5_wp*x) * sinh(0.5_wp*x))
      else
         B = hc2 / lambda_m**5 / (exp(x) - 1.0_wp)
      end if
   end function bbody


   pure function hardest_photon_energy(lam_um, J_lam) result(u_photon)
      ! Energy of the hardest photon the radiation field actually carries:
      !   hc / (shortest wavelength at which J_lambda is significantly nonzero).
      !
      ! This is a property of the FIELD, not of the wavelength grid it is
      ! sampled on. The distinction is the whole point of this routine. A dust
      ! model's grid is the grid of its optics tables and can reach far past the
      ! band the transported field occupies, and every shipped model can. The
      ! astrodust and DL07 grid is the T-matrix Q table's, which starts at
      ! 0.0912 um on the default table but at 1.0e-4 um (12398 eV) on the _euv
      ! one -- 912 times harder than a field illuminated to the Lyman limit --
      ! and the Zubko/ZDA optics tables start at 1.0e-3 um (1.24 keV), 91 times
      ! harder, while the field being transported still ends at the Lyman limit.
      ! Measuring the grid instead of the field raises the single-photon bound
      ! by that whole factor with no photon to justify it, which coarsens the
      ! grain-enthalpy bins (umax sets the top of a fixed bin count) and shifts
      ! the emergent SED: 1-2% for Zubko, and more on any grid that reaches the
      ! X-ray -- see docs/EUV_EXTENSION_HOST_REGRESSION.md.
      !
      ! J_lam may be in any units: only ratios within the array are used, so the
      ! caller need not convert. lam_um must be in um. The array order is
      ! irrelevant; every point is examined.
      !
      ! Significance threshold rationale. An exact J > 0 test lets denormal
      ! or roundoff residue at the short-wavelength end set the hardest
      ! photon and unnecessarily enlarge a fixed-count enthalpy grid. The
      ! relative floor 1e-12 rejects such exterior residue while retaining
      ! the continuum samples above it. For example, residue at 1e-30 of
      ! the peak is separated from the floor by eighteen orders of magnitude.
      ! This is a numerical convention, not a universal physical error bound.
      !
      ! The original estimate allowed a factor 1e2 variation in C_abs and
      ! associated a 1e-12 intensity component with a 1e-10 absorption scale.
      ! That estimate also needs comparable wavelength and photon-energy
      ! weights. It does not bound a probability, and 1e-10 exceeds the
      ! 1e-13 temperature-tail criterion. Significant weak hard components
      ! must be resolved; the excluded absorbed power must be checked with
      ! the actual cross sections and wavelength intervals.
      !
      ! The one-sided window argument concerns the INITIAL enthalpy bound:
      ! the narrowing loops expand when the upper state remains populated,
      ! so a low initial estimate can be recovered. Contraction from a high
      ! estimate requires both umax > 1.02*umaxlo and the populated-tail
      ! condition umax > 1.01*U(jcut) (ub in the QM solver). Finite iteration
      ! limits can therefore leave an excessive window after numerical residue.
      ! That argument motivates rejecting residue; it cannot recover photons
      ! excluded from the represented field support in the transition rates.
      !
      ! The shared continuum support needs two distinct significant wavelength
      ! samples. An isolated sample returns zero here and in absorption and
      ! is classified as isolated by radiation_field_support. Absorption and
      ! transition callers increment isolated_radiation_support_count. A line
      ! needs a resolved profile; this pure photon-bound query does not warn.
      real(wp), intent(in) :: lam_um(:)   ! [um]
      real(wp), intent(in) :: J_lam(:)    ! J_lambda, arbitrary units
      real(wp) :: u_photon                ! [erg]

      integer :: first, last

      call radiation_field_support(lam_um, J_lam, first, last)
      u_photon = 0.0_wp
      if (first > 0) u_photon = hc_erg_um / lam_um(first)
   end function hardest_photon_energy


   pure subroutine radiation_field_support(lambda, Jfield, first, last, status)
      ! Significant continuum endpoints, independent of array order. Indices
      ! refer to the shortest and longest accepted wavelengths. Integration
      ! callers additionally require a positive, strictly increasing grid.
      ! Interior zero samples remain part of the represented continuum.
      real(wp), intent(in) :: lambda(:), Jfield(:)
      integer, intent(out) :: first, last
      integer, optional, intent(out) :: status
      real(wp) :: peak, floor
      integer :: k, n

      first = 0
      last = 0
      if (present(status)) status = RADIATION_SUPPORT_EMPTY
      n = min(size(lambda), size(Jfield))
      if (n == 0) return
      peak = maxval(Jfield(1:n))
      if (peak <= 0.0_wp) return
      floor = J_REL_FLOOR * peak
      do k = 1, n
         if (.not. (Jfield(k) > floor .and. lambda(k) > 0.0_wp)) cycle
         if (first == 0) then
            first = k
            last = k
         else
            if (lambda(k) < lambda(first)) first = k
            if (lambda(k) > lambda(last)) last = k
         end if
      end do
      if (first == 0) return
      if (lambda(first) == lambda(last)) then
         first = 0
         last = 0
         if (present(status)) status = RADIATION_SUPPORT_ISOLATED
         return
      end if
      if (present(status)) status = RADIATION_SUPPORT_OK
   end subroutine radiation_field_support

end module radfield
