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
   public :: projected_area_extinction, fresnel_opaque_absorption

   ! Default surface-quadrature count per angular direction used by the
   ! orientation-resolved geometric-optics absorption.  The 2D integrand is
   ! smooth and vanishes at the terminator (grazing reflectance -> 1), so a
   ! moderate Gauss-Legendre grid resolves it; the comparison driver confirms
   ! convergence by doubling this.
   integer, parameter :: NQUAD_GO = 64

contains

   subroutine rayleigh_limit(a_eff, lam, n_r, k_i, eps_ba, qext, qsca, walb, asymm, &
                               al1, al2, al3, al4, be1, be2, lmax, &
                               qext_ori, qabs_ori, qsca_ori, qre_ori)
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
      !
      ! Optional orientation-resolved outputs qext_ori/qabs_ori/qsca_ori
      ! (dimensioned 3, indexed by jori) return the per-orientation cross
      ! sections that the random average above is built from, in the
      ! convention documented in sed/src/q_table_jori.f90:
      !   jori=1: k || a           -> E always transverse to the axis -> Q_b
      !   jori=2: k perp a, E || a  -> E along the axis               -> Q_a
      !   jori=3: k perp a, E perp a-> E transverse to the axis       -> Q_b
      ! In the Rayleigh limit jori=1 and jori=3 are identical (both see the
      ! transverse polarizability alpha_b), and (Q1+Q2+Q3)/3 reproduces the
      ! averaged qext/qsca returned above exactly.
      !
      ! Optional qre_ori (length 3) is the birefringence twin of qabs_ori: the
      ! REAL part of the same dipole forward-scattering response whose
      ! imaginary part is the extinction/absorption.  It uses fac*Re(alpha)
      ! where qabs uses fac*Im(alpha), with the same fac, so 0.5*(qre(3)-qre(2))
      ! = 0.5*fac*(Re alpha_b - Re alpha_a) is the closed-form Rayleigh
      ! birefringence and jori=1 = jori=3 identically.
      real(wp), intent(in)  :: a_eff, lam, n_r, k_i, eps_ba
      real(wp), intent(out) :: qext, qsca, walb, asymm
      ! Optional: analytic scattering-matrix expansion coefficients (see
      ! rayleigh_matrix_expansion below).  Arrays must be dimensioned at
      ! least 3.
      real(wp), optional, intent(out) :: al1(:), al2(:), al3(:), al4(:), be1(:), be2(:)
      integer,  optional, intent(out) :: lmax
      ! Optional: orientation-resolved cross sections, arrays of length 3.
      real(wp), optional, intent(out) :: qext_ori(:), qabs_ori(:), qsca_ori(:)
      real(wp), optional, intent(out) :: qre_ori(:)
      real(wp)    :: axrat, e2, e, ala, alb, fac
      real(wp)    :: qabs_a, qabs_b, qsca_a, qsca_b, qabs
      real(wp)    :: qre_a, qre_b
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

      ! Birefringence twin: real part of the same dipole response, same fac.
      qre_a  = fac * real(alpha_a, kind=wp)
      qre_b  = fac * real(alpha_b, kind=wp)

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

      ! Orientation-resolved cross sections (index = jori).  These are the
      ! same qabs_a/qabs_b/qsca_a/qsca_b that the average above collapses:
      ! E along the axis (a) at jori=2, E transverse (b) at jori=1 and 3.
      if (present(qabs_ori)) then
         qabs_ori(1) = qabs_b;  qabs_ori(2) = qabs_a;  qabs_ori(3) = qabs_b
      end if
      if (present(qsca_ori)) then
         qsca_ori(1) = qsca_b;  qsca_ori(2) = qsca_a;  qsca_ori(3) = qsca_b
      end if
      if (present(qext_ori)) then
         qext_ori(1) = qabs_b + qsca_b
         qext_ori(2) = qabs_a + qsca_a
         qext_ori(3) = qabs_b + qsca_b
      end if
      ! Birefringence twin (real part of the forward-amplitude response), same
      ! jori mapping as qabs_ori: transverse alpha_b at jori=1,3, axial
      ! alpha_a at jori=2, so 0.5*(qre(3)-qre(2)) is the Rayleigh birefringence.
      if (present(qre_ori)) then
         qre_ori(1) = qre_b;  qre_ori(2) = qre_a;  qre_ori(3) = qre_b
      end if

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


   subroutine geometric_optics_limit(a_eff, lam, n_r, k_i, eps_ba, qext, qsca, walb, asymm, &
                               al1, al2, al3, al4, be1, be2, lmax, &
                               qext_ori, qabs_ori, qsca_ori, qre_ori)
      ! Geometric-optics limit, x = 2 pi a_eff / lambda >> 1 (used for x > 50).
      !
      ! Averaged (equal-volume sphere) outputs -- UNCHANGED from the original
      ! sphere model, so the cross-section-only sweep in run_tmatrix.f90 stays
      ! byte-identical:
      !     Q_ext = 2 (extinction paradox),
      !     Q_abs = 1 - exp(-4 k_i x)  (a single chord of length 2 a_eff
      !                                 through an absorbing sphere; k_i=Im(m)),
      !     Q_sca = Q_ext - Q_abs, g = 0,
      ! and the optional scattering matrix stays isotropic (alpha_1 = 1, all
      ! other coefficients zero, LMAX = 0).  eps_ba does not enter any of
      ! these; it is consumed only by the orientation-resolved outputs below.
      !
      ! Domain of validity of the averaged g = 0 / isotropic-matrix
      ! assumption: it is used only for x > 50, where the astrodust size
      ! distribution has essentially no grains (dn/dloga has fallen by many
      ! orders of magnitude by a_eff ~ 4 um), so the contribution to the
      ! C_sca-weighted size integral is negligible at the wavelengths this
      ! table targets.  It must NOT be relied on for a size distribution with
      ! significant large-grain weight.
      !
      ! Optional orientation-resolved outputs qext_ori/qabs_ori/qsca_ori
      ! (length 3, jori index of sed/src/q_table_jori.f90) carry the two
      ! pieces of wave-optics physics that survive for a non-spherical opaque
      ! grain in this limit:
      !
      !   * Extinction: the extinction paradox gives Q_ext(jori) = 2 *
      !     A_proj(jori) / (pi a_eff^2), where A_proj is the geometric shadow
      !     area of the spheroid for that incidence, so the orientation split
      !     follows the projected area alone (projected_area_extinction).
      !     Q_ext(2) = Q_ext(3) here, so qpol_ext = 0 at this order.
      !   * Absorption: for an opaque grain (Im(m) x >> 1, refracted ray fully
      !     absorbed) the absorbed fraction at each illuminated surface point
      !     is 1 minus the Fresnel power reflectance for the local angle of
      !     incidence, with the incident polarization decomposed onto the
      !     local s/p directions of the plane of incidence.  Integrating over
      !     the illuminated surface (fresnel_opaque_absorption) makes Q_abs
      !     depend on the incident polarization relative to the symmetry axis:
      !     that dependence is the absorption dichroism by which jori=2 (E||a)
      !     and jori=3 (E perp a) differ.  jori=1 (k||a) is axially symmetric
      !     about k, so its two transverse polarizations absorb equally and
      !     Q_abs(1) is the polarization-averaged value.
      !
      ! Opaque-limit validity: for astrodust in the UV, where x > 50 falls,
      ! Im(m) runs from ~0.06 (near 4 eV) to ~0.7 (>10 eV), so the chord
      ! optical depth 4 Im(m) x stays >~ 11 across the whole x > 50 region and
      ! the fully-absorbed-refracted-ray assumption holds.  It would weaken
      ! only for a weakly absorbing material (Im(m) x ~ 1), where internal
      ! transmission and a second surface crossing would have to be added.
      real(wp), intent(in)  :: a_eff, lam, n_r, k_i, eps_ba
      real(wp), intent(out) :: qext, qsca, walb, asymm
      real(wp), optional, intent(out) :: al1(:), al2(:), al3(:), al4(:), be1(:), be2(:)
      integer,  optional, intent(out) :: lmax
      real(wp), optional, intent(out) :: qext_ori(:), qabs_ori(:), qsca_ori(:)
      ! Birefringence twin: at geometric-optics order the extinction paradox is
      ! polarization-independent (Q_ext(2) = Q_ext(3)), so the birefringence
      ! 0.5*(qre(3)-qre(2)) vanishes.  This is a documented approximation valid
      ! in the x > 50 region, which carries negligible astrodust weight.
      real(wp), optional, intent(out) :: qre_ori(:)
      real(wp) :: x, qabs
      real(wp) :: qe_o(3), qa_o(3)
      complex(wp) :: m
      real(wp), parameter :: PI = acos(-1.0_wp)
      ! Body-frame incidence and polarization for the three orientations.
      ! Body z = spheroid symmetry axis a.
      !   jori=1: k || a  -> k along z, any transverse polarization (use x)
      !   jori=2: k perp a, E || a      -> k along x, E along z
      !   jori=3: k perp a, E perp a    -> k along x, E along y
      real(wp), parameter :: KHAT_PAR(3)  = (/ 0.0_wp, 0.0_wp, 1.0_wp /)
      real(wp), parameter :: EHAT_PAR(3)  = (/ 1.0_wp, 0.0_wp, 0.0_wp /)
      real(wp), parameter :: KHAT_PERP(3) = (/ 1.0_wp, 0.0_wp, 0.0_wp /)
      real(wp), parameter :: EHAT_AXIS(3) = (/ 0.0_wp, 0.0_wp, 1.0_wp /)
      real(wp), parameter :: EHAT_EQ(3)   = (/ 0.0_wp, 1.0_wp, 0.0_wp /)

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

      if (present(qext_ori) .or. present(qabs_ori) .or. present(qsca_ori)) then
         call projected_area_extinction(eps_ba, qe_o)
         m = cmplx(n_r, abs(k_i), kind=wp)
         call fresnel_opaque_absorption(eps_ba, m, KHAT_PAR,  EHAT_PAR,  NQUAD_GO, qa_o(1))
         call fresnel_opaque_absorption(eps_ba, m, KHAT_PERP, EHAT_AXIS, NQUAD_GO, qa_o(2))
         call fresnel_opaque_absorption(eps_ba, m, KHAT_PERP, EHAT_EQ,   NQUAD_GO, qa_o(3))
         if (present(qext_ori)) qext_ori(1:3) = qe_o
         if (present(qabs_ori)) qabs_ori(1:3) = qa_o
         if (present(qsca_ori)) qsca_ori(1:3) = qe_o - qa_o
      end if
      ! Birefringence vanishes at geometric-optics order (see declaration).
      if (present(qre_ori)) qre_ori(1:3) = 0.0_wp
   end subroutine geometric_optics_limit


   subroutine projected_area_extinction(eps_ba, qext_ori)
      ! Orientation-resolved extinction efficiency of a spheroid in the
      ! extinction-paradox limit, Q_ext(jori) = 2 * A_proj(jori)/(pi a_eff^2).
      !
      ! Geometry (same eps_ba convention as rayleigh_limit; eps_ba = b/a with
      ! eps_ba > 1 oblate).  From the equal-volume-sphere radius a_eff the
      ! semi-axes are c = a_eff eps_ba^(-2/3) (symmetry axis) and a_s = a_eff
      ! eps_ba^(1/3) (equatorial), so a_s/c = eps_ba.  The geometric shadow is
      !   jori=1 (k || symmetry axis c):  A_proj = pi a_s^2      = pi a_eff^2 eps_ba^(2/3)
      !   jori=2,3 (k perp c):            A_proj = pi a_s c       = pi a_eff^2 eps_ba^(-1/3)
      ! giving Q_ext(1) = 2 eps_ba^(2/3) and Q_ext(2)=Q_ext(3)=2 eps_ba^(-1/3).
      real(wp), intent(in)  :: eps_ba
      real(wp), intent(out) :: qext_ori(:)
      real(wp) :: two_third, minus_third
      two_third   = eps_ba**( 2.0_wp/3.0_wp)
      minus_third = eps_ba**(-1.0_wp/3.0_wp)
      qext_ori(1) = 2.0_wp * two_third
      qext_ori(2) = 2.0_wp * minus_third
      qext_ori(3) = 2.0_wp * minus_third
   end subroutine projected_area_extinction


   subroutine fresnel_opaque_absorption(eps_ba, m, k_hat, e_hat, nquad, qabs, area_proj)
      ! Absorption efficiency Q_abs = C_abs/(pi a_eff^2) of an opaque spheroid
      ! in the geometric-optics limit, for a plane wave travelling along k_hat
      ! with electric-field unit vector e_hat (both unit vectors in the body
      ! frame whose z-axis is the symmetry axis a; e_hat must be transverse to
      ! k_hat).
      !
      ! For an opaque grain the refracted ray is fully absorbed, so the
      ! absorbed fraction at an illuminated surface point is
      !     A(theta_i) = 1 - ( |e.s|^2 R_s(theta_i) + |e.p|^2 R_p(theta_i) ),
      ! with theta_i the local angle of incidence, s/p the senkrecht/parallel
      ! directions of the local plane of incidence, and R_s, R_p the Fresnel
      ! power reflectances for the complex index m.  The absorbed power is
      !     C_abs = INT_illuminated A(theta_i) dA_proj,
      ! dA_proj = cos(theta_i) dA the projected-area element.  Q_abs = C_abs /
      ! (pi a_eff^2) is independent of a_eff, so the surface is taken at
      ! a_eff = 1 (semi-axes a_s = eps_ba^(1/3), c = eps_ba^(-2/3)).
      !
      ! The polarization dichroism between jori=2 and jori=3 arises entirely
      ! because e_hat (= a_hat vs transverse) projects differently onto the
      ! local s/p basis as it sweeps the surface: R_p < R_s away from normal
      ! incidence, so the polarization that is more nearly p-like over the
      ! illuminated surface absorbs more.
      !
      ! Quadrature: Gauss-Legendre with nquad nodes in u over [0,pi] and nquad
      ! in v over [0,2pi], surface parametrized as
      !   (a_s sin u cos v, a_s sin u sin v, c cos u).
      ! The r_u x r_v cross product carries the area Jacobian, and its dot with
      ! -k_hat is cos(theta_i) dA (positive on the illuminated half), so the
      ! sin(u) weighting and the illumination mask are handled exactly.
      real(wp),    intent(in)  :: eps_ba
      complex(wp), intent(in)  :: m
      real(wp),    intent(in)  :: k_hat(3), e_hat(3)
      integer,     intent(in)  :: nquad
      real(wp),    intent(out) :: qabs
      real(wp), optional, intent(out) :: area_proj

      real(wp), parameter :: PI = acos(-1.0_wp)
      real(wp) :: as, cc
      real(wp), allocatable :: ug(:), uw(:), vg(:), vw(:)
      integer  :: iu, iv
      real(wp) :: su, cu, sv, cv
      real(wp) :: dS(3), dmag, proj, cos_i, sin2, w
      real(wp) :: nhat(3), shat(3), phat(3), smag, es, ep, fs, fp
      real(wp) :: rs2, rp2, reff, absorb
      real(wp) :: qint, aint
      complex(wp) :: m2, cost, rs, rp

      as = eps_ba**( 1.0_wp/3.0_wp)
      cc = eps_ba**(-2.0_wp/3.0_wp)
      m2 = m*m

      allocate(ug(nquad), uw(nquad), vg(nquad), vw(nquad))
      call gauss_legendre(nquad, 0.0_wp,          PI,          ug, uw)
      call gauss_legendre(nquad, 0.0_wp, 2.0_wp*PI,             vg, vw)

      qint = 0.0_wp
      aint = 0.0_wp
      do iu = 1, nquad
         su = sin(ug(iu));  cu = cos(ug(iu))
         do iv = 1, nquad
            cv = cos(vg(iv));  sv = sin(vg(iv))
            ! Outward surface-element vector dS = (r_u x r_v) du dv.
            dS(1) = cc*as * su*su * cv
            dS(2) = cc*as * su*su * sv
            dS(3) = as*as * su    * cu
            ! Illuminated where the outward normal faces the source:
            ! proj = -k.dS = cos(theta_i) dA_element > 0.
            proj = -(k_hat(1)*dS(1) + k_hat(2)*dS(2) + k_hat(3)*dS(3))
            if (proj <= 0.0_wp) cycle
            dmag  = sqrt(dS(1)*dS(1) + dS(2)*dS(2) + dS(3)*dS(3))
            cos_i = proj / dmag
            if (cos_i > 1.0_wp) cos_i = 1.0_wp
            nhat  = dS / dmag

            ! Fresnel power reflectances at this incidence angle.
            sin2 = 1.0_wp - cos_i*cos_i
            cost = sqrt(1.0_wp - sin2/m2)
            if (real(cost, kind=wp) < 0.0_wp) cost = -cost
            rs = (cos_i - m*cost) / (cos_i + m*cost)
            rp = (m*cos_i - cost) / (m*cos_i + cost)
            rs2 = real(rs, kind=wp)**2 + aimag(rs)**2
            rp2 = real(rp, kind=wp)**2 + aimag(rp)**2

            ! Decompose e_hat onto the local s (perp plane of incidence) and
            ! p (in plane, transverse to k) directions.
            shat(1) = k_hat(2)*nhat(3) - k_hat(3)*nhat(2)
            shat(2) = k_hat(3)*nhat(1) - k_hat(1)*nhat(3)
            shat(3) = k_hat(1)*nhat(2) - k_hat(2)*nhat(1)
            smag = sqrt(shat(1)**2 + shat(2)**2 + shat(3)**2)
            if (smag > 1.0e-12_wp) then
               shat = shat / smag
               phat(1) = k_hat(2)*shat(3) - k_hat(3)*shat(2)
               phat(2) = k_hat(3)*shat(1) - k_hat(1)*shat(3)
               phat(3) = k_hat(1)*shat(2) - k_hat(2)*shat(1)
               es = e_hat(1)*shat(1) + e_hat(2)*shat(2) + e_hat(3)*shat(3)
               ep = e_hat(1)*phat(1) + e_hat(2)*phat(2) + e_hat(3)*phat(3)
               fs = es*es;  fp = ep*ep
            else
               ! Normal incidence: plane of incidence undefined but R_s = R_p.
               fs = 0.5_wp;  fp = 0.5_wp
            end if

            reff   = fs*rs2 + fp*rp2
            absorb = 1.0_wp - reff
            w = uw(iu) * vw(iv)
            qint = qint + w * proj * absorb
            aint = aint + w * proj
         end do
      end do
      deallocate(ug, uw, vg, vw)

      qabs = qint / PI
      if (present(area_proj)) area_proj = aint
   end subroutine fresnel_opaque_absorption


   subroutine gauss_legendre(n, a, b, x, w)
      ! Gauss-Legendre nodes x and weights w on [a,b], n-point rule.
      ! Newton iteration on the Legendre polynomial (Numerical Recipes gauleg).
      integer,  intent(in)  :: n
      real(wp), intent(in)  :: a, b
      real(wp), intent(out) :: x(n), w(n)
      real(wp), parameter :: PI = acos(-1.0_wp)
      real(wp), parameter :: TOL = 1.0e-15_wp
      integer  :: i, j, mid
      real(wp) :: p1, p2, p3, pp, z, z1, xm, xl
      mid = (n + 1) / 2
      xm  = 0.5_wp * (b + a)
      xl  = 0.5_wp * (b - a)
      do i = 1, mid
         z = cos(PI * (real(i, wp) - 0.25_wp) / (real(n, wp) + 0.5_wp))
         do
            p1 = 1.0_wp;  p2 = 0.0_wp
            do j = 1, n
               p3 = p2;  p2 = p1
               p1 = (real(2*j-1, wp)*z*p2 - real(j-1, wp)*p3) / real(j, wp)
            end do
            pp = real(n, wp) * (z*p1 - p2) / (z*z - 1.0_wp)
            z1 = z
            z  = z1 - p1/pp
            if (abs(z - z1) <= TOL) exit
         end do
         x(i)     = xm - xl*z
         x(n+1-i) = xm + xl*z
         w(i)     = 2.0_wp * xl / ((1.0_wp - z*z) * pp*pp)
         w(n+1-i) = w(i)
      end do
   end subroutine gauss_legendre

end module asymptotic_optics
