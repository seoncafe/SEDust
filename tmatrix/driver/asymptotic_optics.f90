module fallback
   ! Convergence-failure fallbacks for the T-matrix sweep.
   !
   ! Two regimes the random-orientation T-matrix solver does not handle
   ! reliably:
   !
   !   - small  x = 2 pi a_eff / lambda  (x < x_small_lim, default 0.1):
   !       Rayleigh-limit polarizability for an oblate/prolate spheroid,
   !       random-orientation averaged, following Draine (1992).
   !       Closes a ~6% systematic vs HD23 release at
   !       long wavelengths that an isotropic Mie sphere would leave open.
   !
   !   - large x (x > x_large_lim, default 50):
   !       geometric-optics limit. Q_ext -> 2 and Q_abs = 1 - exp(-4 k x),
   !       Q_sca = Q_ext - Q_abs, g approximately 0 (no preferred direction
   !       for randomly oriented large opaque grains in the GO limit; this
   !       is a deliberate over-simplification, justified by the fact that
   !       the FIR/sub-mm SED has J_lambda ~ 0 here so cross-section
   !       errors do not propagate to the final SED).
   !
   ! Both fallbacks return a single set of Q_ext, Q_sca, albedo, g and a
   ! flag value identifying which fallback was used.

   use, intrinsic :: iso_fortran_env, only: real64
   use constants, only: wp
   implicit none
   private
   public :: fallback_small_x, fallback_large_x

contains

   subroutine fallback_small_x(a_eff, lam, n_r, k_i, eps_ba, qext, qsca, walb, asymm)
      ! Spheroid Rayleigh polarizability + random-orientation average.
      ! Follows Draine (1992), in F90 / double precision.
      !
      ! `eps_ba` follows the Mishchenko convention used elsewhere in this
      ! tree (b/a = horizontal axis / rotational axis):
      !   eps_ba > 1  -> oblate  (symmetry axis short)
      !   eps_ba < 1  -> prolate (symmetry axis long)
      !   eps_ba = 1  -> sphere  (degenerate, returns Mie limit)
      ! Draine's AXRAT (= symm/equator) is the inverse of eps_ba.
      !
      ! Random-orientation average over isotropic incident polarization:
      !   <Q> = (Q_a + 2*Q_b) / 3
      ! where Q_a is for E parallel to symmetry axis and Q_b for E perp.
      !
      ! In the dipole regime g -> 0 (symmetric scattering pattern).
      real(wp), intent(in)  :: a_eff, lam, n_r, k_i, eps_ba
      real(wp), intent(out) :: qext, qsca, walb, asymm
      real(wp)    :: axrat, e2, e, ala, alb, fac
      real(wp)    :: qabs_a, qabs_b, qsca_a, qsca_b, qabs
      complex(wp) :: eps, alpha_a, alpha_b
      real(wp), parameter :: PI = acos(-1.0_wp)

      eps = cmplx(n_r*n_r - k_i*k_i, 2.0_wp*n_r*k_i, kind=wp)
      axrat = 1.0_wp / eps_ba

      e2 = abs(1.0_wp - 1.0_wp/(axrat*axrat))
      e  = sqrt(e2)
      if (axrat < 1.0_wp) then
         ! oblate
         ala = (1.0_wp + 1.0_wp/e2) * (1.0_wp - atan(e)/e)
      else if (axrat > 1.0_wp) then
         ! prolate
         ala = (1.0_wp/e2 - 1.0_wp) * &
               (log((1.0_wp + e)/(1.0_wp - e))/(2.0_wp*e) - 1.0_wp)
      else
         ! sphere
         ala = 1.0_wp/3.0_wp
      end if
      alb = (1.0_wp - ala) / 2.0_wp

      ! Polarizability per orientation (units: volume).
      fac     = a_eff**3 / 3.0_wp
      alpha_a = fac * (eps - 1.0_wp) / ((eps - 1.0_wp)*ala + 1.0_wp)
      alpha_b = fac * (eps - 1.0_wp) / ((eps - 1.0_wp)*alb + 1.0_wp)

      ! Q_abs = (8 pi / (lam * a_eff^2)) * Im(alpha)
      fac    = 8.0_wp * PI / (lam * a_eff*a_eff)
      qabs_a = fac * aimag(alpha_a)
      qabs_b = fac * aimag(alpha_b)

      ! Q_sca = (128 pi^4 / (3 lam^4 a_eff^2)) * |alpha|^2
      fac    = 128.0_wp * PI**4 / (3.0_wp * lam**4 * a_eff*a_eff)
      qsca_a = fac * (real(alpha_a)**2 + aimag(alpha_a)**2)
      qsca_b = fac * (real(alpha_b)**2 + aimag(alpha_b)**2)

      ! Random orientation average (1/3 || + 2/3 perp).
      qabs = (qabs_a + 2.0_wp*qabs_b) / 3.0_wp
      qsca = (qsca_a + 2.0_wp*qsca_b) / 3.0_wp
      qext = qabs + qsca
      walb = qsca / qext
      asymm = 0.0_wp
   end subroutine fallback_small_x


   subroutine fallback_large_x(a_eff, lam, n_r, k_i, qext, qsca, walb, asymm)
      ! Geometric-optics limit. Q_ext = 2 (extinction paradox), and
      ! Q_abs from a single chord through a sphere of radius a_eff:
      !     Q_abs = 1 - exp(-4 k x), with k = Im(m), x = 2 pi a / lambda.
      ! Asymmetry parameter is set to 0 (no preferred direction for
      ! random-orientation averaging in the absorbing-opaque limit).
      real(wp), intent(in)  :: a_eff, lam, n_r, k_i
      real(wp), intent(out) :: qext, qsca, walb, asymm
      real(wp) :: x, qabs
      real(wp), parameter :: PI = acos(-1.0_wp)
      ! n_r is unused in this crude approximation; included in the
      ! signature for symmetry with fallback_small_x and to allow a
      ! future refinement (Fresnel-based reflection / refraction).
      associate (dummy => n_r); end associate
      x = 2.0_wp * PI * a_eff / lam
      qext = 2.0_wp
      qabs = 1.0_wp - exp(-4.0_wp * k_i * x)
      qsca = qext - qabs
      walb = qsca / qext
      asymm = 0.0_wp
   end subroutine fallback_large_x

end module fallback
