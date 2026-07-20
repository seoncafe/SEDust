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

   use, intrinsic :: iso_fortran_env, only: real64
   use constants, only: wp
   implicit none
   private
   public :: tmatrix_oriented_ext

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

end module tmatrix_oriented
