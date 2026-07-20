module tmatrix_oriented
   ! Orientation-resolved extinction cross section of an axially symmetric
   ! spheroid in the T-matrix regime (0.1 < x < 50, x = 2 pi a_eff / lambda).
   !
   ! The random-orientation solver TMD_ONE_SCATMAT (src/tmd_one.f) leaves the
   ! converged T-matrix in COMMON /TMAT/ and returns its multipole truncation
   ! order NMAX_TM.  Mishchenko's fixed-orientation amplitude routine AMPL
   ! (src/ampl_oriented.f) reads that same common block and returns the 2x2
   ! amplitude matrix (VV, VH, HV, HH) for a chosen incidence/scattering
   ! geometry and particle orientation.  Evaluated in the forward direction
   ! (scattering direction = incidence direction), the co-polar diagonal
   ! amplitudes give the extinction cross section through the optical theorem
   !
   !     C_ext = (4 pi / k) * Im[ S_forward ],   k = 2 pi / lambda,
   !
   ! for the incident polarization whose forward co-polar amplitude is
   ! S_forward.  With k = 2 pi / lambda this is C_ext = 2 * lambda * Im[S].
   !
   ! The three orientations follow the convention of sed/src/q_table_jori.f90
   ! (a = spheroid symmetry axis):
   !   jori=1: k || a
   !   jori=2: k perp a, E || a
   !   jori=3: k perp a, E perp a
   !
   ! They are realized here with incidence fixed along the laboratory z-axis
   ! (polar angle THET0 = 0, so the meridional theta-hat points along x and
   ! phi-hat along y in the forward direction) and the particle symmetry axis
   ! placed by the Euler angles (ALPHA, BETA):
   !   jori=1: BETA = 0   -> axis along z = along k; E is transverse to the
   !                        axis for both polarizations, so VV = HH.  C_ext
   !                        from Im[VV].
   !   jori=2: BETA = 90, ALPHA = 0 -> axis along x.  The V (theta) incident
   !                        polarization has E along x = along the axis, so
   !                        C_ext from Im[VV].
   !   jori=3: BETA = 90, ALPHA = 0 -> axis along x.  The H (phi) incident
   !                        polarization has E along y, transverse to the
   !                        axis, so C_ext from Im[HH].

   ! The orientation-resolved scattering cross section C_sca(jori) follows
   ! from the same fixed-orientation amplitude matrix, integrated over the
   ! full scattering sphere.  For incidence along +z (THET0 = 0) with fixed
   ! polarization p, the differential scattering cross section is |S|^2 (the
   ! AMPL amplitudes already carry the dimension of length), so
   !
   !     C_sca,p = INT_0^{2pi} INT_0^pi ( |S_{theta<-p}|^2 + |S_{phi<-p}|^2 )
   !                                     sin(TL1) dTL1 dPL1,
   !
   ! with (TL1, PL1) the scattering polar/azimuth angles.  For V (theta)
   ! incidence the integrand is |VV|^2 + |HV|^2; for H (phi) incidence it is
   ! |HH|^2 + |VH|^2.  The two sweeps needed are BETA = 0 (jori=1, V) and
   ! BETA = 90 (jori=2 reads V, jori=3 reads H), reusing the stored T-matrix.
   ! Then Q_sca(jori) = C_sca(jori)/(pi a_eff^2) and, by conservation,
   ! Q_abs(jori) = Q_ext(jori) - Q_sca(jori).
   !
   ! The fixed-orientation scattered intensity is a trigonometric polynomial
   ! of degree ~2*NMAX in cos(TL1) and in PL1 (NMAX = T-matrix truncation
   ! order).  Gauss-Legendre nodes in cos(TL1) over [-1,1] with N_theta =
   ! NMAX+2 and a uniform grid in PL1 over [0,2pi) with N_phi = 2*NMAX+2
   ! integrate it essentially exactly.

   use, intrinsic :: iso_fortran_env, only: real64
   use constants, only: wp
   implicit none
   private
   public :: tmatrix_oriented_ext, tmatrix_oriented_cross

contains

   subroutine tmatrix_oriented_ext(a_eff, lam, m, eps_ba, np, ddelt, ndgs, &
                                   qext_ori, ierr)
      ! INPUT
      !   a_eff   equivalent-volume-sphere radius [microns]
      !   lam     wavelength [microns]
      !   m       complex refractive index (Im >= 0 for absorption)
      !   eps_ba  aspect ratio b/a (Mishchenko convention; > 1 oblate)
      !   np      shape flag (-1 spheroid, -2 cylinder, >=0 Chebyshev)
      !   ddelt   convergence tolerance
      !   ndgs    Gauss-quadrature multiplier
      ! OUTPUT
      !   qext_ori(3)  orientation-resolved Q_ext = C_ext/(pi a_eff^2),
      !                indexed by jori = 1,2,3
      !   ierr         status from TMD_ONE_SCATMAT (0 = converged); on a
      !                nonzero return qext_ori is left at 0
      real(wp),    intent(in)  :: a_eff, lam, eps_ba, ddelt
      complex(wp), intent(in)  :: m
      integer,     intent(in)  :: np, ndgs
      real(wp),    intent(out) :: qext_ori(3)
      integer,     intent(out) :: ierr

      integer, parameter :: NPL = 201            ! matches tmd.par.f (NPN2+1)
      real(wp), parameter :: PI = acos(-1.0_wp)
      real(wp) :: mrr, mri
      real(wp) :: qext, qsca, walb, asymm
      real(wp) :: al1(NPL), al2(NPL), al3(NPL), al4(NPL), be1(NPL), be2(NPL)
      integer  :: lmax, nmax_tm
      complex(wp) :: vv, vh, hv, hh
      real(wp) :: cext_fac, cext, area

      external :: tmd_one_scatmat, ampl

      qext_ori = 0.0_wp
      mrr = real(m, kind=wp)
      mri = abs(aimag(m))

      call tmd_one_scatmat(a_eff, lam, mrr, mri, eps_ba, np, ddelt, ndgs, &
                           qext, qsca, walb, asymm, &
                           al1, al2, al3, al4, be1, be2, lmax, ierr, nmax_tm)
      if (ierr /= 0) return

      ! Optical-theorem prefactor 4 pi / k = 2 * lambda.
      cext_fac = 2.0_wp * lam
      area     = PI * a_eff * a_eff

      ! jori=1: k || a  (BETA = 0).  E transverse to the axis; VV = HH.
      call ampl(nmax_tm, lam, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
                0.0_wp, 0.0_wp, vv, vh, hv, hh)
      cext = cext_fac * aimag(vv)
      qext_ori(1) = cext / area

      ! jori=2: k perp a, E || a  (BETA = 90, ALPHA = 0).  V polarization
      ! (theta-hat = x) has E along the axis; C_ext from VV.
      call ampl(nmax_tm, lam, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
                0.0_wp, 90.0_wp, vv, vh, hv, hh)
      cext = cext_fac * aimag(vv)
      qext_ori(2) = cext / area

      ! jori=3: k perp a, E perp a  (BETA = 90, ALPHA = 0).  H polarization
      ! (phi-hat = y) has E transverse to the axis; C_ext from HH.
      cext = cext_fac * aimag(hh)
      qext_ori(3) = cext / area
   end subroutine tmatrix_oriented_ext


   subroutine tmatrix_oriented_cross(a_eff, lam, m, eps_ba, np, ddelt, ndgs, &
                                     qext_ori, qsca_ori, qabs_ori, ierr, qmult)
      ! Orientation-resolved extinction, scattering, and absorption cross
      ! sections in the T-matrix regime.  One T-matrix solve is shared: the
      ! extinction comes from the forward-amplitude optical theorem (exactly
      ! as tmatrix_oriented_ext) and the scattering from the fixed-orientation
      ! phase-matrix integral over the scattering sphere.  Absorption is the
      ! difference C_abs = C_ext - C_sca.
      !
      ! INPUT
      !   a_eff, lam  equivalent-volume-sphere radius and wavelength [microns]
      !   m           complex refractive index (Im >= 0 for absorption)
      !   eps_ba      aspect ratio b/a (Mishchenko convention; > 1 oblate)
      !   np          shape flag (-1 spheroid, -2 cylinder, >=0 Chebyshev)
      !   ddelt       convergence tolerance
      !   ndgs        Gauss-quadrature multiplier for the T-matrix solve
      !   qmult       optional integer (default 1) that scales the angular
      !               quadrature counts N_theta, N_phi; used to confirm the
      !               scattering integral is resolved (qmult = 2 doubles both)
      ! OUTPUT
      !   qext_ori(3), qsca_ori(3), qabs_ori(3)  Q = C/(pi a_eff^2), jori index
      !   ierr        status from TMD_ONE_SCATMAT (0 = converged); on nonzero
      !               return all three arrays are left at 0
      real(wp),    intent(in)  :: a_eff, lam, eps_ba, ddelt
      complex(wp), intent(in)  :: m
      integer,     intent(in)  :: np, ndgs
      real(wp),    intent(out) :: qext_ori(3), qsca_ori(3), qabs_ori(3)
      integer,     intent(out) :: ierr
      integer,     intent(in), optional :: qmult

      integer, parameter :: NPL = 201            ! matches tmd.par.f (NPN2+1)
      real(wp), parameter :: PI = acos(-1.0_wp)
      real(wp) :: mrr, mri
      real(wp) :: qext, qsca, walb, asymm
      real(wp) :: al1(NPL), al2(NPL), al3(NPL), al4(NPL), be1(NPL), be2(NPL)
      integer  :: lmax, nmax_tm, mult
      complex(wp) :: vv, vh, hv, hh
      real(wp) :: cext_fac, area
      real(wp) :: csca_v0, csca_h0, csca_v90, csca_h90

      external :: tmd_one_scatmat, ampl

      qext_ori = 0.0_wp;  qsca_ori = 0.0_wp;  qabs_ori = 0.0_wp
      mult = 1
      if (present(qmult)) mult = max(1, qmult)

      mrr = real(m, kind=wp)
      mri = abs(aimag(m))

      call tmd_one_scatmat(a_eff, lam, mrr, mri, eps_ba, np, ddelt, ndgs, &
                           qext, qsca, walb, asymm, &
                           al1, al2, al3, al4, be1, be2, lmax, ierr, nmax_tm)
      if (ierr /= 0) return

      cext_fac = 2.0_wp * lam                    ! optical-theorem 4 pi / k
      area     = PI * a_eff * a_eff

      ! Extinction: forward-amplitude optical theorem, per jori.
      call ampl(nmax_tm, lam, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
                0.0_wp, 0.0_wp, vv, vh, hv, hh)
      qext_ori(1) = cext_fac * aimag(vv) / area
      call ampl(nmax_tm, lam, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
                0.0_wp, 90.0_wp, vv, vh, hv, hh)
      qext_ori(2) = cext_fac * aimag(vv) / area
      qext_ori(3) = cext_fac * aimag(hh) / area

      ! Scattering: integrate the fixed-orientation phase matrix.  BETA = 0
      ! gives jori=1 (V incidence); BETA = 90 gives jori=2 (V) and jori=3 (H).
      call scatter_sphere_integral(nmax_tm, lam, 0.0_wp,  mult, csca_v0,  csca_h0)
      call scatter_sphere_integral(nmax_tm, lam, 90.0_wp, mult, csca_v90, csca_h90)
      qsca_ori(1) = csca_v0  / area
      qsca_ori(2) = csca_v90 / area
      qsca_ori(3) = csca_h90 / area

      qabs_ori = qext_ori - qsca_ori
   end subroutine tmatrix_oriented_cross


   subroutine scatter_sphere_integral(nmax_tm, lam, beta, mult, csca_v, csca_h)
      ! Integrate the fixed-orientation differential scattering cross section
      ! over the full scattering sphere, at incidence along +z (TL = 0,
      ! PL = 0) and particle orientation (ALPHA = 0, BETA).  Returns
      !   csca_v = INT ( |VV|^2 + |HV|^2 ) dOmega   [V (theta) incidence]
      !   csca_h = INT ( |HH|^2 + |VH|^2 ) dOmega   [H (phi)   incidence]
      ! in the length^2 units of the AMPL amplitudes (i.e. cross-section area).
      !
      ! Quadrature: Gauss-Legendre in x = cos(TL1) over [-1,1] with
      ! N_theta = mult*(NMAX+2) nodes (dOmega = -dx dPL1, so the sin(TL1)
      ! factor is absorbed by the change of variable), and a uniform grid in
      ! PL1 over [0,2pi) with N_phi = mult*(2*NMAX+2) points, weight
      ! 2*pi/N_phi.  mult = 1 is the resolved default; mult = 2 is used only
      ! to confirm convergence.
      integer,  intent(in)  :: nmax_tm, mult
      real(wp), intent(in)  :: lam, beta
      real(wp), intent(out) :: csca_v, csca_h

      real(wp), parameter :: PI = acos(-1.0_wp)
      real(wp), parameter :: RAD2DEG = 180.0_wp / acos(-1.0_wp)
      integer  :: n_theta, n_phi, it, ip
      real(wp), allocatable :: xg(:), wg(:)
      real(wp) :: wphi, tl1_deg, pl1_deg, wt
      complex(wp) :: vv, vh, hv, hh

      external :: ampl, gauss

      n_theta = mult * (nmax_tm + 2)
      n_phi   = mult * (2*nmax_tm + 2)
      allocate(xg(n_theta), wg(n_theta))
      ! GAUSS(N, IND1=0, IND2=0, Z, W): Gauss-Legendre nodes on [-1,1] with
      ! weights summing to 2 (src/tmd_one.f); IND2=0 keeps it silent.
      call gauss(n_theta, 0, 0, xg, wg)
      wphi = 2.0_wp * PI / real(n_phi, kind=wp)

      csca_v = 0.0_wp
      csca_h = 0.0_wp
      do it = 1, n_theta
         tl1_deg = acos(xg(it)) * RAD2DEG          ! [0,180]
         do ip = 1, n_phi
            pl1_deg = real(ip-1, kind=wp) * 360.0_wp / real(n_phi, kind=wp)  ! [0,360)
            call ampl(nmax_tm, lam, 0.0_wp, tl1_deg, 0.0_wp, pl1_deg, &
                      0.0_wp, beta, vv, vh, hv, hh)
            wt = wg(it) * wphi
            csca_v = csca_v + wt * (abs(vv)**2 + abs(hv)**2)
            csca_h = csca_h + wt * (abs(hh)**2 + abs(vh)**2)
         end do
      end do
      deallocate(xg, wg)
   end subroutine scatter_sphere_integral

end module tmatrix_oriented
