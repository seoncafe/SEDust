module asymptotic_optics
   ! Closed-form asymptotic limits of the optical response of a spheroid,
   ! covering the two size-parameter regimes in which the random-
   ! orientation T-matrix solver is either unnecessary or unreliable.
   ! x = 2 pi a_eff / lambda throughout.
   !
   !   - rayleigh_limit, valid for x << 1 (used for x < 0.1):
   !       Rayleigh-limit polarizability for an oblate/prolate spheroid,
   !       random-orientation averaged, following Draine (1992).
   !       Closes a ~6% systematic vs HD23 release at
   !       long wavelengths that an isotropic Mie sphere would leave open.
   !
   !   - geometric_optics_limit, valid for x >> 1 (used for x > 50):
   !       Q_ext -> 2 and Q_abs = 1 - exp(-4 k x),
   !       Q_sca = Q_ext - Q_abs, g approximately 0 (no preferred direction
   !       for randomly oriented large opaque grains in the GO limit; this
   !       is a deliberate over-simplification, justified by the fact that
   !       the FIR/sub-mm SED has J_lambda ~ 0 here so cross-section
   !       errors do not propagate to the final SED).
   !
   ! Both return a single set of Q_ext, Q_sca, albedo, g.
   !
   ! Both also optionally return the six generalized-spherical-function
   ! expansion coefficients of the random-orientation scattering matrix,
   ! in the same convention and normalization as TMD_ONE_SCATMAT (alpha_1(0)
   ! = 1).  The optional arguments are absent in the cross-section-only
   ! sweep, so those calls are unaffected.

   use, intrinsic :: iso_fortran_env, only: real64
   use constants, only: wp
   implicit none
   private
   public :: rayleigh_limit, geometric_optics_limit

contains

   subroutine rayleigh_limit(a_eff, lam, n_r, k_i, eps_ba, qext, qsca, walb, asymm, &
                               al1, al2, al3, al4, be1, be2, lmax)
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
      ! Optional: analytic scattering-matrix expansion coefficients (see
      ! rayleigh_matrix_expansion below).  Arrays must be dimensioned at
      ! least 3.
      real(wp), optional, intent(out) :: al1(:), al2(:), al3(:), al4(:), be1(:), be2(:)
      integer,  optional, intent(out) :: lmax
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

      if (present(al1)) call rayleigh_matrix_expansion(alpha_a, alpha_b, &
                             al1, al2, al3, al4, be1, be2, lmax)
   end subroutine rayleigh_limit


   subroutine rayleigh_matrix_expansion(alpha_a, alpha_b, al1, al2, al3, al4, be1, be2, lmax)
      ! Exact scattering matrix of a randomly oriented axially symmetric
      ! Rayleigh (dipole) scatterer, in expansion-coefficient form.
      !
      ! The body-frame polarizability tensor is diag(alpha_a, alpha_b,
      ! alpha_b).  Writing it as alpha_ij = abar*delta_ij + del*(n_i n_j
      ! - delta_ij/3) with
      !     abar = (alpha_a + 2 alpha_b)/3,   del = alpha_a - alpha_b,
      ! the orientation averages of products of tensor components are
      !     A = <|alpha_xx|^2>       = |abar|^2 + (4/45)|del|^2
      !     B = <alpha_xx alpha_yy*> = |abar|^2 - (2/45)|del|^2
      !     C = <|alpha_xy|^2>       = (1/15)|del|^2       (A = B + 2C)
      ! Feeding these through the amplitude matrix of a dipole gives
      !     F11 = [A(1+u^2) + C(3-u^2)]/2,  F22 = (A-C)(1+u^2)/2,
      ! (continued)
      !     F33 = (B+C)u,  F44 = (B-C)u,  F12 = -(A-C)(1-u^2)/2,  F34 = 0,
      ! with u = cos(scattering angle), all still unnormalized.  Only l =
      ! 0, 1, 2 survive, so LMAX = 2.  Setting R = (A-C)/N with the
      ! normalization N = (2/3)A + (4/3)C (which enforces alpha_1(0) = 1):
      !     alpha_1 = (1, 0, R/3)      alpha_2 = (0, 0, 2R)
      !     alpha_3 = 0                alpha_4 = (0, (B-C)/N, 0)
      !     beta_1  = (0, 0, 2R/sqrt(6))               beta_2 = 0
      ! For an isotropic polarizability (del = 0) this reduces to the
      ! textbook Rayleigh matrix: alpha_1 = (1,0,1/2), alpha_2 = (0,0,3),
      ! alpha_4 = (0,3/2,0), beta_1 = (0,0,sqrt(6)/2), F11 propto
      ! 1 + cos^2(Theta).
      !
      ! The corresponding total cross section, proportional to
      ! (4/3)(A + 2C) = (4/3)(|abar|^2 + (2/9)|del|^2), is identically
      ! (|alpha_a|^2 + 2|alpha_b|^2)/3 -- the same orientation average
      ! that rayleigh_limit uses for Q_sca, so the matrix and the cross
      ! section are mutually consistent by construction.
      !
      ! alpha_1(1) = 0 gives g = 0, matching rayleigh_limit.
      complex(wp), intent(in)  :: alpha_a, alpha_b
      real(wp), optional, intent(out) :: al1(:), al2(:), al3(:), al4(:), be1(:), be2(:)
      integer,  optional, intent(out) :: lmax
      complex(wp) :: abar, del
      real(wp)    :: aa, bb, cc, xnorm, r

      abar = (alpha_a + 2.0_wp*alpha_b) / 3.0_wp
      del  = alpha_a - alpha_b
      aa = abs(abar)**2 + (4.0_wp/45.0_wp) * abs(del)**2
      bb = abs(abar)**2 - (2.0_wp/45.0_wp) * abs(del)**2
      cc = (1.0_wp/15.0_wp) * abs(del)**2
      xnorm = (2.0_wp/3.0_wp)*aa + (4.0_wp/3.0_wp)*cc
      r = (aa - cc) / xnorm

      if (present(al1)) then; al1 = 0.0_wp; al1(1) = 1.0_wp; al1(3) = r/3.0_wp; end if
      if (present(al2)) then; al2 = 0.0_wp; al2(3) = 2.0_wp*r;                  end if
      if (present(al3)) then; al3 = 0.0_wp;                                     end if
      if (present(al4)) then; al4 = 0.0_wp; al4(2) = (bb - cc)/xnorm;           end if
      if (present(be1)) then; be1 = 0.0_wp; be1(3) = 2.0_wp*r/sqrt(6.0_wp);     end if
      if (present(be2)) then; be2 = 0.0_wp;                                     end if
      if (present(lmax)) lmax = 2
   end subroutine rayleigh_matrix_expansion


   subroutine geometric_optics_limit(a_eff, lam, n_r, k_i, qext, qsca, walb, asymm, &
                               al1, al2, al3, al4, be1, be2, lmax)
      ! Geometric-optics limit. Q_ext = 2 (extinction paradox), and
      ! Q_abs from a single chord through a sphere of radius a_eff:
      !     Q_abs = 1 - exp(-4 k x), with k = Im(m), x = 2 pi a / lambda.
      ! Asymmetry parameter is set to 0 (no preferred direction for
      ! random-orientation averaging in the absorbing-opaque limit).
      !
      ! Scattering matrix (optional outputs): isotropic and unpolarizing,
      ! i.e. alpha_1 = (1) and every other coefficient zero, LMAX = 0.
      ! This is the matrix-level counterpart of the g = 0 assumption
      ! above and is equally crude -- the true large-x phase function is
      ! strongly forward-peaked.  Domain of validity: it is used only for
      ! x > 50, where the astrodust size distribution has essentially no
      ! grains (dn/dloga has fallen by many orders of magnitude by
      ! a_eff ~ 4 um), so the contribution to the C_sca-weighted size
      ! integral is negligible at the optical wavelengths this table
      ! targets.  It must NOT be relied on for a size distribution with
      ! significant large-grain weight.
      real(wp), intent(in)  :: a_eff, lam, n_r, k_i
      real(wp), intent(out) :: qext, qsca, walb, asymm
      real(wp), optional, intent(out) :: al1(:), al2(:), al3(:), al4(:), be1(:), be2(:)
      integer,  optional, intent(out) :: lmax
      real(wp) :: x, qabs
      real(wp), parameter :: PI = acos(-1.0_wp)
      ! n_r is unused in this crude approximation; included in the
      ! signature for symmetry with rayleigh_limit and to allow a
      ! future refinement (Fresnel-based reflection / refraction).
      associate (dummy => n_r); end associate
      x = 2.0_wp * PI * a_eff / lam
      qext = 2.0_wp
      qabs = 1.0_wp - exp(-4.0_wp * k_i * x)
      qsca = qext - qabs
      walb = qsca / qext
      asymm = 0.0_wp

      if (present(al1)) then; al1 = 0.0_wp; al1(1) = 1.0_wp; end if
      if (present(al2)) al2 = 0.0_wp
      if (present(al3)) al3 = 0.0_wp
      if (present(al4)) al4 = 0.0_wp
      if (present(be1)) be1 = 0.0_wp
      if (present(be2)) be2 = 0.0_wp
      if (present(lmax)) lmax = 0
   end subroutine geometric_optics_limit

end module asymptotic_optics
