module scattering_matrix_oriented
   ! Fixed-orientation scattering (Mueller) matrix Z of an axially symmetric
   ! spheroid: the angle-resolved counterpart of the orientation-resolved
   ! cross sections in driver/tmatrix_oriented.f90.  Two regimes are covered:
   !
   !   - mueller_matrix_fixed_orientation : T-matrix regime (0.1 < x < 50),
   !     one AMPL (src/ampl_oriented.f) evaluation of the 2x2 complex
   !     amplitude matrix from the converged T-matrix, then Mishchenko's
   !     phase-matrix bilinears.
   !   - rayleigh_mueller_matrix_oriented : analytic electric-dipole limit
   !     (x << 1), amplitude S_pq = k^2 (e_p^sca . alpha . e_q^inc) with the
   !     spheroid polarizability tensor alpha = diag(alpha_b, alpha_b,
   !     alpha_a) of asymptotic_optics.
   !
   ! GEOMETRY AND STOKES CONVENTION
   ! The grain symmetry axis a-hat lies along z (Euler ALPHA = BETA = 0).
   ! Incidence is at polar angle theta_i from the axis, azimuth 0; the
   ! scattered direction at (theta_s, phi).  Because the grain is
   ! axisymmetric only the azimuth difference phi = PL1 - PL0 matters, so
   !     Z = Z(theta_i; theta_s, phi),   16 elements, in um^2 sr^-1.
   ! Stokes vectors use Mishchenko's meridional basis of each propagation
   ! direction: e_1 = theta-hat (V), e_2 = phi-hat (H), real unit vectors,
   ! Q = I_v - I_h.  This is exactly the basis AMPL rotates its amplitudes
   ! into (its AL / AP matrices use theta-hat = (ct cp, ct sp, -st) and
   ! phi-hat = (-sp, cp, 0)), so the analytic dipole and T-matrix matrices
   ! share one convention.
   !
   ! The amplitude matrix S = [[S11,S12],[S21,S22]] = [[VV,VH],[HV,HH]]
   ! (AMPL returns VV, VH, HV, HH, each carrying the dimension of length, so
   ! |S|^2 is an area and Z a differential scattering cross section).  The
   ! 16 phase-matrix elements follow from S by the bilinear combinations of
   ! Mishchenko (Appl. Opt. 39, 1026, 2000), copied verbatim from
   ! src/ampld.lp.f (Z11..Z44) in phase_matrix_from_amplitude below.
   !
   ! SERIAL, COMMON-BASED.  mueller_matrix_fixed_orientation reads the
   ! converged T-matrix from COMMON /TMAT/, which a prior TMD_ONE_SCATMAT
   ! call for the same (a_eff, lam, m, shape) must have left valid (it does,
   ! via its /TMATK/ restore).  It is therefore not thread-safe and must not
   ! be called across a fresh T-matrix solve for a different size.
   !
   ! SYMMETRIES (verified numerically in driver/compare_scatmat_aligned.f90)
   !   - phi mirror: Z(theta_i; theta_s, 360 - phi) equals Z(theta_i;
   !     theta_s, phi) with the sign of the two off-diagonal 2x2 blocks
   !     (elements 13,14,23,24,31,32,41,42) flipped.  Store phi in [0,180].
   !   - theta_i = 0: Z is independent of phi and takes the six-element
   !     block-diagonal form.
   !   - equatorial mirror (oblate z -> -z), used to store theta_i in
   !     [0, 90] only.  The verified mapping (azimuth phi is UNCHANGED) is
   !
   !         Z(180 - theta_i; 180 - theta_s, phi)(i,j)
   !                        = S(i,j) * Z(theta_i; theta_s, phi)(i,j)
   !
   !     with the sign pattern S (rows i = 1..4, cols j = 1..4)
   !
   !         S = [ +  +  -  - ]
   !             [ +  +  -  - ]
   !             [ -  -  +  + ]
   !             [ -  -  +  + ]
   !
   !     i.e. the two off-diagonal 2x2 blocks flip sign, the two diagonal
   !     blocks keep it.  Verified numerically to ~1e-15 relative to Z11
   !     across a grid of (theta_i, theta_s, phi); see anchor E(iii).  This is
   !     the z -> -z reflection of the oblate spheroid: incidence
   !     (theta_i, 0) -> (180 - theta_i, 0) and scattering (theta_s, phi) ->
   !     (180 - theta_s, phi), the improper reflection flipping the handedness
   !     (V, U) off-diagonal blocks.

   use, intrinsic :: iso_fortran_env, only: real64
   use constants, only: wp
   use asymptotic_optics, only: spheroid_dipole_polarizability
   implicit none
   private
   public :: mueller_matrix_fixed_orientation, rayleigh_mueller_matrix_oriented

contains

   subroutine mueller_matrix_fixed_orientation(nmax_tm, lam, theta_i, theta_s, phi, z)
      ! Fixed-orientation Mueller matrix in the T-matrix regime, from the
      ! converged T-matrix left in COMMON /TMAT/ by TMD_ONE_SCATMAT.
      !
      ! INPUT
      !   nmax_tm  multipole truncation of the stored T-matrix (the value
      !            TMD_ONE_SCATMAT returned for this a_eff, lam, m, shape)
      !   lam      wavelength [microns], the same value passed to the solve
      !   theta_i  incidence polar angle from the axis [deg], in [0,180]
      !   theta_s  scattering polar angle [deg], in [0,180]
      !   phi      scattering azimuth relative to the incidence plane [deg],
      !            in [0,360]
      ! OUTPUT
      !   z(4,4)   Mueller matrix [um^2 sr^-1] in the meridional (v,h) basis
      integer,  intent(in)  :: nmax_tm
      real(wp), intent(in)  :: lam, theta_i, theta_s, phi
      real(wp), intent(out) :: z(4,4)

      complex(wp) :: vv, vh, hv, hh
      external :: ampl

      ! Grain axis along z: ALPHA = 0, BETA = 0.  Incidence at (theta_i, 0),
      ! scattering at (theta_s, phi).
      call ampl(nmax_tm, lam, theta_i, theta_s, 0.0_wp, phi, &
                0.0_wp, 0.0_wp, vv, vh, hv, hh)
      call phase_matrix_from_amplitude(vv, vh, hv, hh, z)
   end subroutine mueller_matrix_fixed_orientation


   subroutine rayleigh_mueller_matrix_oriented(a_eff, lam, nr, ki, eps_ba, &
                                               theta_i, theta_s, phi, z)
      ! Analytic electric-dipole Mueller matrix of a spheroid in the same
      ! geometry and (v,h) basis as mueller_matrix_fixed_orientation, valid
      ! for x = 2 pi a_eff / lambda << 1.
      !
      ! The scattered amplitude is the dipole radiation of the induced moment
      ! p = alpha . E_inc, projected onto the scattered meridional basis:
      !     S_pq = k^2 ( e_p^sca . alpha . e_q^inc ),   k = 2 pi / lambda,
      ! with the body-frame polarizability tensor alpha = diag(alpha_b,
      ! alpha_b, alpha_a) (axis along z) in volume units, alpha_a / alpha_b
      ! from asymptotic_optics.  With this normalization the forward optical
      ! theorem (4 pi / k) Im S11(forward) reproduces the dipole absorption
      ! C_abs(jori) of rayleigh_limit exactly (the O(k^4) scattering term of
      ! C_ext is higher order in the bare-polarizability amplitude).
      !
      ! INPUT
      !   a_eff        equivalent-volume-sphere radius [microns]
      !   lam          wavelength [microns]
      !   nr, ki       real and imaginary parts of the refractive index
      !   eps_ba       axis ratio b/a (Mishchenko convention; > 1 oblate)
      !   theta_i, theta_s, phi   angles [deg], as above
      ! OUTPUT
      !   z(4,4)       Mueller matrix [um^2 sr^-1]
      real(wp), intent(in)  :: a_eff, lam, nr, ki, eps_ba
      real(wp), intent(in)  :: theta_i, theta_s, phi
      real(wp), intent(out) :: z(4,4)

      real(wp), parameter :: PI = acos(-1.0_wp)
      real(wp), parameter :: DEG = acos(-1.0_wp) / 180.0_wp
      complex(wp) :: alpha_a, alpha_b, k2
      complex(wp) :: vv, vh, hv, hh
      real(wp)    :: et_i(3), ep_i(3), et_s(3), ep_s(3)
      real(wp)    :: k
      real(wp)    :: ti, ts, ph

      call spheroid_dipole_polarizability(a_eff, nr, ki, eps_ba, alpha_a, alpha_b)
      k  = 2.0_wp * PI / lam
      k2 = cmplx(k*k, 0.0_wp, kind=wp)

      ti = theta_i * DEG
      ts = theta_s * DEG
      ph = phi     * DEG

      ! Meridional unit vectors of the incidence direction (azimuth 0) and
      ! the scattering direction (azimuth phi).
      et_i = (/  cos(ti),        0.0_wp,       -sin(ti) /)     ! theta-hat (V)
      ep_i = (/  0.0_wp,         1.0_wp,        0.0_wp  /)     ! phi-hat   (H)
      et_s = (/  cos(ts)*cos(ph), cos(ts)*sin(ph), -sin(ts) /)
      ep_s = (/ -sin(ph),         cos(ph),         0.0_wp  /)

      ! S_pq = k^2 (e_p^sca . alpha . e_q^inc), mapped to S11=VV, S12=VH,
      ! S21=HV, S22=HH.
      vv = k2 * dipole_bilinear(et_s, et_i, alpha_a, alpha_b)
      vh = k2 * dipole_bilinear(et_s, ep_i, alpha_a, alpha_b)
      hv = k2 * dipole_bilinear(ep_s, et_i, alpha_a, alpha_b)
      hh = k2 * dipole_bilinear(ep_s, ep_i, alpha_a, alpha_b)

      call phase_matrix_from_amplitude(vv, vh, hv, hh, z)
   end subroutine rayleigh_mueller_matrix_oriented


   complex(wp) function dipole_bilinear(e_sca, e_inc, alpha_a, alpha_b) result(s)
      ! e_sca . alpha . e_inc for the axial tensor alpha = diag(alpha_b,
      ! alpha_b, alpha_a) with the symmetry axis along z.  The transverse
      ! (x,y) components carry alpha_b, the axial (z) component alpha_a.
      real(wp),    intent(in) :: e_sca(3), e_inc(3)
      complex(wp), intent(in) :: alpha_a, alpha_b
      s = alpha_b * (e_sca(1)*e_inc(1) + e_sca(2)*e_inc(2)) &
        + alpha_a *  e_sca(3)*e_inc(3)
   end function dipole_bilinear


   subroutine phase_matrix_from_amplitude(s11, s12, s21, s22, z)
      ! Mueller (phase) matrix Z from the 2x2 complex amplitude matrix
      ! S = [[S11,S12],[S21,S22]] = [[VV,VH],[HV,HH]].  The 16 bilinear
      ! combinations are copied verbatim from Mishchenko's src/ampld.lp.f
      ! (its Z11..Z44, Eqs. (13)-(29) of Appl. Opt. 39, 1026, 2000); the real
      ! part is the physical phase-matrix element.
      complex(wp), intent(in)  :: s11, s12, s21, s22
      real(wp),    intent(out) :: z(4,4)
      complex(wp), parameter :: CI = (0.0_wp, 1.0_wp)

      z(1,1) = 0.5_wp*real( s11*conjg(s11)+s12*conjg(s12) &
                           +s21*conjg(s21)+s22*conjg(s22), kind=wp)
      z(1,2) = 0.5_wp*real( s11*conjg(s11)-s12*conjg(s12) &
                           +s21*conjg(s21)-s22*conjg(s22), kind=wp)
      z(1,3) = real(-s11*conjg(s12)-s22*conjg(s21), kind=wp)
      z(1,4) = real(CI*(s11*conjg(s12)-s22*conjg(s21)), kind=wp)

      z(2,1) = 0.5_wp*real( s11*conjg(s11)+s12*conjg(s12) &
                           -s21*conjg(s21)-s22*conjg(s22), kind=wp)
      z(2,2) = 0.5_wp*real( s11*conjg(s11)-s12*conjg(s12) &
                           -s21*conjg(s21)+s22*conjg(s22), kind=wp)
      z(2,3) = real(-s11*conjg(s12)+s22*conjg(s21), kind=wp)
      z(2,4) = real(CI*(s11*conjg(s12)+s22*conjg(s21)), kind=wp)

      z(3,1) = real(-s11*conjg(s21)-s22*conjg(s12), kind=wp)
      z(3,2) = real(-s11*conjg(s21)+s22*conjg(s12), kind=wp)
      z(3,3) = real( s11*conjg(s22)+s12*conjg(s21), kind=wp)
      z(3,4) = real(-CI*(s11*conjg(s22)+s21*conjg(s12)), kind=wp)

      z(4,1) = real(CI*(s21*conjg(s11)+s22*conjg(s12)), kind=wp)
      z(4,2) = real(CI*(s21*conjg(s11)-s22*conjg(s12)), kind=wp)
      z(4,3) = real(-CI*(s22*conjg(s11)-s12*conjg(s21)), kind=wp)
      z(4,4) = real( s22*conjg(s11)-s12*conjg(s21), kind=wp)
   end subroutine phase_matrix_from_amplitude

end module scattering_matrix_oriented
