module radfield
   ! Radiation-field constructors and Planck-function utilities for the
   ! astrodust SED solver.
   ! (J_Mathis, J_diluted_bbody, calc_bbody, bbody) — verbatim except for
   ! moving them out of read_data into a dedicated module and using the
   ! shared `wp` from constants.

   use, intrinsic :: iso_fortran_env, only: real64
   use constants, only: wp
   implicit none
   private
   public :: J_Mathis, J_diluted_bbody, calc_bbody, bbody
   public :: hardest_photon_energy
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

   ! Fraction of peak J_lambda below which a spectral component is treated as
   ! carrying no photons at all; see hardest_photon_energy.
   real(wp), parameter :: J_REL_FLOOR = 1.0e-12_wp

   logical, save :: use_mathis_corrected = .true.

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


   subroutine calc_bbody(T, lambda_um, spec)
      real(wp), intent(in)  :: T, lambda_um(:)
      real(wp), intent(out) :: spec(:)
      integer :: i
      do i = 1, size(lambda_um)
         spec(i) = bbody(T, lambda_um(i))
      end do
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
      ! band the transported field occupies: the astrodust and DL07 Q tables
      ! stop at the Lyman limit (0.0912 um = 13.595 eV) and so happen to agree
      ! with a field that is illuminated to the Lyman limit, but the Zubko/ZDA
      ! DustEM tables start at 1.0e-3 um (1.24 keV), 91 times harder, while the
      ! field being transported still ends at the Lyman limit. Measuring the
      ! grid instead of the field there raises the single-photon bound by a
      ! factor of 91 with no photon to justify it, which coarsens the
      ! grain-enthalpy bins (umax sets the top of a fixed bin count) and shifts
      ! the emergent SED by 1-2% -- see docs/EUV_EXTENSION_HOST_REGRESSION.md.
      !
      ! J_lam may be in any units: only ratios within the array are used, so the
      ! caller need not convert. lam_um must be in um. The array order is
      ! irrelevant; every point is examined.
      !
      ! Significance threshold. Exact zeros must be excluded (J_Mathis returns
      ! exactly 0 below 0.0912 um), but an exact `> 0` test is too fragile: a
      ! host that hands over a field carrying denormal or roundoff-level
      ! residue in its unilluminated bins would push the bound back up by
      ! orders of magnitude. A component a fraction J_REL_FLOOR = 1e-12 below
      ! the peak of J_lambda contributes at most that fraction of the photon
      ! absorption rate; allowing two decades for the wavelength dependence of
      ! C_abs over the illuminated band, the enthalpy states it could populate
      ! carry probability <~ 1e-10 of the peak, at or below the 1e-13 tail
      ! level (PMIN_UP, PMIN_UP_QM) at which both stochastic solvers already
      ! truncate the enthalpy window. Such a component therefore cannot change
      ! a resolved excursion, and 1e-12 sits 18 decades above the residue it is
      ! meant to reject.
      ! The error is deliberately one-sided. Both refinement loops EXPAND the
      ! window when the top bin is still populated, so an underestimate of the
      ! bound costs iterations and is recovered; their contraction is guarded
      ! (umax > 1.02*umaxlo, umax > 1.01*ub(jcut)) and capped at MAX_ITER, so an
      ! overestimate is not. The threshold is set high enough to reject junk for
      ! that reason.
      !
      ! Returns 0 when the field is zero (or non-positive) everywhere. Callers
      ! wrap the result in max() against their own floor, so 0 selects that
      ! floor and no special case is needed at the call site.
      real(wp), intent(in) :: lam_um(:)   ! [um]
      real(wp), intent(in) :: J_lam(:)    ! J_lambda, arbitrary units
      real(wp) :: u_photon                ! [erg]

      real(wp) :: j_floor, lam_short
      integer  :: i, n

      u_photon = 0.0_wp
      n = min(size(lam_um), size(J_lam))
      if (n < 1) return

      j_floor = J_REL_FLOOR * maxval(J_lam(1:n))
      if (j_floor <= 0.0_wp) return       ! field everywhere zero or non-positive

      lam_short = huge(1.0_wp)
      do i = 1, n
         if (J_lam(i) > j_floor .and. lam_um(i) > 0.0_wp) &
            lam_short = min(lam_short, lam_um(i))
      end do
      if (lam_short < huge(1.0_wp)) u_photon = hc_erg_um / lam_short
   end function hardest_photon_energy

end module radfield
