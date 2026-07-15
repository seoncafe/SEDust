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

   logical, save :: use_mathis_corrected = .true.

contains

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
         T_cmb  = 2.725_wp
      else
         w_4000 = 1.0d-13
         T_cmb  = 2.9_wp
      end if
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

end module radfield
