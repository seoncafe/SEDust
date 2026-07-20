program compare_scatmat_aligned
   ! Certification of the fixed-orientation Mueller matrix engine
   ! (driver/scattering_matrix_oriented.f90) and its size-integrated aligned
   ! population (driver/aligned_population_optics.f90).  Anchors A-E work at
   ! single sizes; F-G add the size-integrated products.  Each prints measured
   ! numbers and PASS/FAIL; a nonzero exit status is returned if any fails.
   !
   !   A closure vs the orientation-resolved cross sections (jori table):
   !     INT Z11 dOmega and INT (Z11 +/- Z12) dOmega against C_sca(jori) from
   !     driver/tmatrix_oriented.f90, on the same exact quadrature.
   !   B optical-theorem wiring: (4 pi / k) Im S_forward against C_ext(jori).
   !   C random-orientation average of Z over grain Euler angles against the
   !     GSP-expansion F(Theta) of TMD_ONE_SCATMAT (certifies the bilinear
   !     signs and the (v,h) basis mapping).
   !   D Rayleigh limit: closed-form dipole matrix and optical-theorem
   !     consistency (internal exactness), plus dipole-vs-T-matrix continuity.
   !   E symmetries: phi mirror, theta_i = 0 azimuthal independence and block
   !     structure, and the empirically derived equatorial mapping.
   !   F size-integrated closure of the aligned K and Csca elements against the
   !     same f_align-weighted integrals of the 4-block jori table, at 0.55 um.
   !   G Rayleigh K sin^2 law: dipole forward amplitudes exactly, and the
   !     T-matrix forward amplitudes near x = 0.137 to a few percent.
   !   H meridional rotation at general geometry: 4 pi <Z>/Csca_random from the
   !     orientation-averaged amplitude engine at arbitrary (theta_i, theta_s,
   !     phi) against L(pi - sigma2) F(Theta) L(-sigma1) built with the SAME
   !     rotation logic mueller_matrix_total uses (the definitive sign
   !     certification for the aligned-scattering library).
   !
   ! Runs from tmatrix/.  Dielectric index, shape flags, and tolerances match
   ! run_scatmat.f90 (eps_ba = 1.4, oblate NP = -1, DDELT = 1e-3, NDGS = 2).

   use, intrinsic :: iso_fortran_env, only: error_unit
   use constants,         only: wp
   use read_index,        only: load_index, interp_m
   use tmatrix_oriented,  only: tmatrix_oriented_cross
   use asymptotic_optics, only: rayleigh_limit, spheroid_dipole_polarizability
   use scattering_matrix_oriented, only: mueller_matrix_fixed_orientation, &
                                         rayleigh_mueller_matrix_oriented
   use size_dist_mod,    only: load_size_dist, n_size, a_dist, dn_ad
   use q_table_jori_mod, only: falign_hd23
   use aligned_population_optics, only: accumulate_aligned_population
   implicit none

   character(len=*), parameter :: f_index = &
      '../data/dielectric/index_DH21Ad_P0.20_0.00_1.400'
   character(len=*), parameter :: f_sdist  = '../data/release/size_distribution.dat'
   character(len=*), parameter :: f_qjori  = &
      'output/q_astrodust_jori_P0.20_Fe0.00_1.400.dat.gz'
   character(len=*), parameter :: f_wave   = '../data/dielectric/DH21_wave'
   character(len=*), parameter :: f_aeff   = '../data/dielectric/DH21_aeff'

   integer,  parameter :: NPL     = 201          ! matches src/tmd.par.f
   real(wp), parameter :: EPS_BA  = 1.4_wp
   real(wp), parameter :: DDELT   = 1.0e-3_wp
   integer,  parameter :: NDGS    = 2, NP_OBL = -1
   real(wp), parameter :: X_SMALL = 0.1_wp, X_LARGE = 50.0_wp
   real(wp), parameter :: PI      = acos(-1.0_wp)
   real(wp), parameter :: DEG     = acos(-1.0_wp) / 180.0_wp

   ! T-matrix test points and the Rayleigh continuity pair.
   integer,  parameter :: NTM = 3
   real(wp), parameter :: TM_A(NTM)   = (/ 0.10_wp, 0.30_wp, 0.05_wp  /)
   real(wp), parameter :: TM_L(NTM)   = (/ 0.55_wp, 0.55_wp, 0.365_wp /)

   external :: tmd_one_scatmat, ampl, gauss, scatmat_from_moments

   integer :: nfail_total

   nfail_total = 0

   call load_index(f_index)

   write(*,'(a)') '======================================================================'
   write(*,'(a)') ' Fixed-orientation Mueller matrix -- Stage 1 anchors A-E'
   write(*,'(a)') '======================================================================'

   call anchor_a()
   call anchor_b()
   call anchor_c()
   call anchor_d()
   call anchor_e()
   call anchor_f()
   call anchor_g()
   call anchor_h()

   write(*,'(a)') '======================================================================'
   if (nfail_total == 0) then
      write(*,'(a)') ' ALL ANCHORS PASSED'
   else
      write(*,'(a,i0,a)') ' FAILURES: ', nfail_total, ' check(s) failed'
      stop 1
   end if

contains

   ! ==================================================================
   ! Anchor A: closure of INT Z dOmega against the jori cross sections.
   ! ==================================================================
   subroutine anchor_a()
      integer  :: ip, nmax
      real(wp) :: a, lam, area, maxe
      real(wp) :: qext_o(3), qsca_o(3), qabs_o(3), qsca_rand
      real(wp) :: a1(NPL), a2(NPL), a3(NPL), a4(NPL), b1(NPL), b2(NPL)
      real(wp) :: i11, i11p, i11m
      integer  :: lmax
      logical  :: ok

      write(*,'(a)') '----------------------------------------------------------------------'
      write(*,'(a)') ' Anchor A: INT Z11 dOmega and INT (Z11 +/- Z12) dOmega vs C_sca(jori)'
      write(*,'(a)') '   tol 1e-6 relative (same T-matrix, exact band-limited quadrature)'
      maxe = 0.0_wp
      do ip = 1, NTM
         a = TM_A(ip);  lam = TM_L(ip)
         call solve_full(a, lam, nmax, qext_o, qsca_o, qabs_o, &
                         a1, a2, a3, a4, b1, b2, lmax, qsca_rand, ok)
         if (.not. ok) then
            write(*,'(a,f6.3,a,f6.3)') '   T-matrix did not converge at a=', a, ' lam=', lam
            nfail_total = nfail_total + 1;  cycle
         end if
         area = PI * a * a
         write(*,'(a,f6.3,a,f6.4,a,i0)') '   a_eff=', a, ' um  lam=', lam, ' um   NMAX=', nmax

         ! theta_i = 0: unpolarized closure vs jori=1.
         call sphere_integrals(nmax, lam, 0.0_wp, i11, i11p, i11m)
         maxe = max(maxe, chk_rel('   ti=0  INT Z11      = Csca(j1)', i11, qsca_o(1)*area, 1.0e-6_wp))

         ! theta_i = 90: unpolarized and the two linear-polarization closures.
         call sphere_integrals(nmax, lam, 90.0_wp, i11, i11p, i11m)
         maxe = max(maxe, chk_rel('   ti=90 INT Z11      =(Csca2+3)/2', i11, &
                    0.5_wp*(qsca_o(2)+qsca_o(3))*area, 1.0e-6_wp))
         maxe = max(maxe, chk_rel('   ti=90 INT(Z11+Z12) = Csca(j2)', i11p, qsca_o(2)*area, 1.0e-6_wp))
         maxe = max(maxe, chk_rel('   ti=90 INT(Z11-Z12) = Csca(j3)', i11m, qsca_o(3)*area, 1.0e-6_wp))
      end do
      write(*,'(a,es10.2)') '   Anchor A max relative error = ', maxe
   end subroutine anchor_a


   ! ==================================================================
   ! Anchor B: optical-theorem wiring, (4 pi / k) Im S_fwd vs C_ext(jori).
   ! ==================================================================
   subroutine anchor_b()
      integer  :: ip, nmax, lmax
      real(wp) :: a, lam, area, cext_fac, maxe
      real(wp) :: qext_o(3), qsca_o(3), qabs_o(3), qsca_rand
      real(wp) :: a1(NPL), a2(NPL), a3(NPL), a4(NPL), b1(NPL), b2(NPL)
      complex(wp) :: vv, vh, hv, hh
      logical  :: ok

      write(*,'(a)') '----------------------------------------------------------------------'
      write(*,'(a)') ' Anchor B: (4 pi/k) Im S_forward vs C_ext(jori)   [tol 1e-6 relative]'
      maxe = 0.0_wp
      do ip = 1, NTM
         a = TM_A(ip);  lam = TM_L(ip)
         call solve_full(a, lam, nmax, qext_o, qsca_o, qabs_o, &
                         a1, a2, a3, a4, b1, b2, lmax, qsca_rand, ok)
         if (.not. ok) then
            write(*,'(a,f6.3,a,f6.3)') '   T-matrix did not converge at a=', a, ' lam=', lam
            nfail_total = nfail_total + 1;  cycle
         end if
         area     = PI * a * a
         cext_fac = 2.0_wp * lam                     ! 4 pi / k
         write(*,'(a,f6.3,a,f6.4)') '   a_eff=', a, ' um  lam=', lam, ' um'

         ! theta_i = 0 forward (theta_s = 0, phi = 0): both polarizations -> jori=1.
         call ampl(nmax, lam, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, vv, vh, hv, hh)
         maxe = max(maxe, chk_rel('   ti=0  Im S11 -> Cext(j1)', cext_fac*aimag(vv), qext_o(1)*area, 1.0e-6_wp))
         maxe = max(maxe, chk_rel('   ti=0  Im S22 -> Cext(j1)', cext_fac*aimag(hh), qext_o(1)*area, 1.0e-6_wp))

         ! theta_i = 90 forward: V -> jori=2, H -> jori=3.
         call ampl(nmax, lam, 90.0_wp, 90.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, vv, vh, hv, hh)
         maxe = max(maxe, chk_rel('   ti=90 Im S11 -> Cext(j2)', cext_fac*aimag(vv), qext_o(2)*area, 1.0e-6_wp))
         maxe = max(maxe, chk_rel('   ti=90 Im S22 -> Cext(j3)', cext_fac*aimag(hh), qext_o(3)*area, 1.0e-6_wp))
      end do
      write(*,'(a,es10.2)') '   Anchor B max relative error = ', maxe
   end subroutine anchor_b


   ! ==================================================================
   ! Anchor C: random-orientation average of Z vs the GSP F(Theta).
   ! ==================================================================
   subroutine anchor_c()
      integer,  parameter :: NCP = 2, NTHC = 3
      real(wp), parameter :: CP_A(NCP)   = (/ 0.10_wp, 0.30_wp /)
      real(wp), parameter :: CP_L(NCP)   = (/ 0.55_wp, 0.55_wp /)
      real(wp), parameter :: THC(NTHC)   = (/ 30.0_wp, 90.0_wp, 140.0_wp /)
      integer,  parameter :: IDX(NTHC)   = (/ 31, 91, 141 /)     ! into 181-pt grid
      integer  :: ic, k, nmax, lmax, nb, na
      real(wp) :: a, lam, area, csca_rand, qsca_rand, maxe, f11r
      real(wp) :: qext_o(3), qsca_o(3), qabs_o(3)
      real(wp) :: a1(NPL), a2(NPL), a3(NPL), a4(NPL), b1(NPL), b2(NPL)
      real(wp) :: th181(181), f11(181), f22(181), f33(181), f44(181), f12(181), f34(181)
      real(wp) :: zbar(4,4), fmod(4,4)
      logical  :: ok

      write(*,'(a)') '----------------------------------------------------------------------'
      write(*,'(a)') ' Anchor C: 4 pi <Z>/Csca_random vs GSP F(Theta)   [certifies signs]'
      write(*,'(a)') '   tol 1e-5 rel (|F|>1e-3 F11) / abs 1e-5 F11 (else)'
      maxe = 0.0_wp
      do ic = 1, NCP
         a = CP_A(ic);  lam = CP_L(ic)
         call solve_full(a, lam, nmax, qext_o, qsca_o, qabs_o, &
                         a1, a2, a3, a4, b1, b2, lmax, qsca_rand, ok)
         if (.not. ok) then
            write(*,'(a,f6.3,a,f6.3)') '   T-matrix did not converge at a=', a, ' lam=', lam
            nfail_total = nfail_total + 1;  cycle
         end if
         area      = PI * a * a
         csca_rand = qsca_rand * area
         call scatmat_from_moments(a1, a2, a3, a4, b1, b2, lmax, 181, &
                                   th181, f11, f22, f33, f44, f12, f34)
         nb = 2*nmax + 2                            ! GL nodes in cos(BETA)
         na = 4*nmax + 4                            ! uniform nodes in ALPHA
         write(*,'(a,f6.3,a,f6.4,a,i0,a,i0,a,i0)') '   a_eff=', a, ' um  lam=', lam, &
              ' um   NMAX=', nmax, '   nBETA=', nb, '  nALPHA=', na

         do k = 1, NTHC
            call orient_average_z(nmax, lam, THC(k), nb, na, zbar)
            fmod = 4.0_wp * PI * zbar / csca_rand
            f11r = f11(IDX(k))
            write(*,'(a,f6.1,a)') '     Theta=', THC(k), ' deg'
            maxe = max(maxe, chk_elem('       F11', fmod(1,1), f11(IDX(k)), f11r))
            maxe = max(maxe, chk_elem('       F22', fmod(2,2), f22(IDX(k)), f11r))
            maxe = max(maxe, chk_elem('       F33', fmod(3,3), f33(IDX(k)), f11r))
            maxe = max(maxe, chk_elem('       F44', fmod(4,4), f44(IDX(k)), f11r))
            maxe = max(maxe, chk_elem('       F12', fmod(1,2), f12(IDX(k)), f11r))
            maxe = max(maxe, chk_elem('       F34', fmod(3,4), f34(IDX(k)), f11r))
            ! block structure: off-diagonal 2x2 blocks must vanish.
            maxe = max(maxe, chk_elem('       F13(=0)', fmod(1,3), 0.0_wp, f11r))
            maxe = max(maxe, chk_elem('       F14(=0)', fmod(1,4), 0.0_wp, f11r))
            maxe = max(maxe, chk_elem('       F23(=0)', fmod(2,3), 0.0_wp, f11r))
            maxe = max(maxe, chk_elem('       F24(=0)', fmod(2,4), 0.0_wp, f11r))
            maxe = max(maxe, chk_elem('       F31(=0)', fmod(3,1), 0.0_wp, f11r))
            maxe = max(maxe, chk_elem('       F41(=0)', fmod(4,1), 0.0_wp, f11r))
            maxe = max(maxe, chk_elem('       F32(=0)', fmod(3,2), 0.0_wp, f11r))
            maxe = max(maxe, chk_elem('       F42(=0)', fmod(4,2), 0.0_wp, f11r))
         end do
      end do
      write(*,'(a,es10.2)') '   Anchor C max normalized error = ', maxe
   end subroutine anchor_c


   ! ==================================================================
   ! Anchor D: Rayleigh dipole -- internal exactness and continuity.
   ! ==================================================================
   subroutine anchor_d()
      integer,  parameter :: NDI = 2
      real(wp), parameter :: DA(NDI) = (/ 0.008_wp, 0.012_wp /)   ! continuity pair
      real(wp), parameter :: DL(NDI) = (/ 0.55_wp,  0.55_wp  /)
      integer  :: id, it, is, ipp, nmax, lmax, i, j
      real(wp) :: a, lam, nr, ki, area, k, maxe, maxc
      real(wp) :: u, ts_deg
      real(wp) :: z(4,4), zdip(4,4)
      real(wp) :: qext_o(3), qsca_o(3), qabs_o(3), qsca_rand
      real(wp) :: qext, qsca, walb, asymm
      real(wp) :: qeo(3), qao(3), qso(3)
      complex(wp) :: alpha_a, alpha_b
      real(wp) :: a1(NPL), a2(NPL), a3(NPL), a4(NPL), b1(NPL), b2(NPL)
      real(wp), parameter :: TI3(3) = (/ 0.0_wp, 45.0_wp, 90.0_wp /)
      real(wp), parameter :: TS3(3) = (/ 30.0_wp, 90.0_wp, 140.0_wp /)
      real(wp), parameter :: PH3(3) = (/ 0.0_wp, 60.0_wp, 120.0_wp /)
      real(wp), parameter :: TSC(5) = (/ 30.0_wp, 60.0_wp, 90.0_wp, 120.0_wp, 150.0_wp /)
      logical  :: ok

      write(*,'(a)') '----------------------------------------------------------------------'
      write(*,'(a)') ' Anchor D(i): dipole closed-form identities at theta_i=0  [tol 1e-10]'
      ! Use the small continuity radius; the identities are alpha-independent.
      a = DA(1);  lam = DL(1)
      call interp_m(lam, nr, ki)
      maxe = 0.0_wp
      do is = 1, 5
         ts_deg = TSC(is)
         u = cos(ts_deg*DEG)
         call rayleigh_mueller_matrix_oriented(a, lam, nr, ki, EPS_BA, 0.0_wp, ts_deg, 0.0_wp, zdip)
         write(*,'(a,f6.1,a)') '     theta_s=', ts_deg, ' deg'
         maxe = max(maxe, chk_rel('       Z12/Z11 = -sin^2/(1+cos^2)', zdip(1,2)/zdip(1,1), &
                    -(1.0_wp-u*u)/(1.0_wp+u*u), 1.0e-10_wp))
         maxe = max(maxe, chk_rel('       Z33/Z11 = 2cos/(1+cos^2)  ', zdip(3,3)/zdip(1,1), &
                    2.0_wp*u/(1.0_wp+u*u), 1.0e-10_wp))
         maxe = max(maxe, chk_rel('       Z44 = Z33                 ', zdip(4,4), zdip(3,3), 1.0e-10_wp))
         maxe = max(maxe, chk_elem('       Z34 = 0                  ', zdip(3,4), 0.0_wp, zdip(1,1)))
      end do
      write(*,'(a,es10.2)') '   Anchor D(i) closed-form max error = ', maxe

      write(*,'(a)') ' Anchor D(i): dipole optical theorem (4 pi/k) Im S_fwd vs C_abs(jori)'
      write(*,'(a)') '   [reference is C_abs, not C_ext: the bare-alpha amplitude carries no'
      write(*,'(a)') '    O(k^4) scattering term; the C_ext gap = C_sca is printed too]'
      k = 2.0_wp*PI/lam
      area = PI*a*a
      call spheroid_dipole_polarizability(a, nr, ki, EPS_BA, alpha_a, alpha_b)
      call rayleigh_limit(a, lam, nr, ki, EPS_BA, qext, qsca, walb, asymm, &
                          qext_ori=qeo, qabs_ori=qao, qsca_ori=qso)
      ! ti=0 (V and H both transverse -> alpha_b), ti=90 V -> alpha_a, ti=90 H -> alpha_b.
      maxe = max(maxe, chk_rel('   ti=0  Im S11 -> Cabs(j1)', 2.0_wp*lam*k*k*aimag(alpha_b), qao(1)*area, 1.0e-10_wp))
      maxe = max(maxe, chk_rel('   ti=90 Im S11 -> Cabs(j2)', 2.0_wp*lam*k*k*aimag(alpha_a), qao(2)*area, 1.0e-10_wp))
      maxe = max(maxe, chk_rel('   ti=90 Im S22 -> Cabs(j3)', 2.0_wp*lam*k*k*aimag(alpha_b), qao(3)*area, 1.0e-10_wp))
      write(*,'(a,3es12.4)') '   (info) Csca/Cabs gap per jori = ', qso(1)/qao(1), qso(2)/qao(2), qso(3)/qao(3)

      write(*,'(a)') ' Anchor D(ii): dipole vs T-matrix full matrices near x=0.1 [tol 5%]'
      do id = 1, NDI
         a = DA(id);  lam = DL(id)
         call interp_m(lam, nr, ki)
         call solve_full(a, lam, nmax, qext_o, qsca_o, qabs_o, &
                         a1, a2, a3, a4, b1, b2, lmax, qsca_rand, ok)
         if (.not. ok) then
            write(*,'(a,f6.3,a,f6.3)') '   T-matrix did not converge at a=', a, ' lam=', lam
            nfail_total = nfail_total + 1;  cycle
         end if
         ! Continuity metric: max elementwise |Z_dip - Z_tm| normalized to the
         ! brightest element Z11.  This is the physically meaningful measure of
         ! how far the leading-order dipole matrix is from the full T-matrix
         ! matrix (it -> 0 as x -> 0); a per-element relative error would blow
         ! up on the small elements that sit at a dipole node (Z_dip exactly
         ! zero, Z_tm small but finite) and certify nothing.
         maxc = 0.0_wp
         do it = 1, 3
            do is = 1, 3
               do ipp = 1, 3
                  call mueller_matrix_fixed_orientation(nmax, lam, TI3(it), TS3(is), PH3(ipp), z)
                  call rayleigh_mueller_matrix_oriented(a, lam, nr, ki, EPS_BA, &
                                                        TI3(it), TS3(is), PH3(ipp), zdip)
                  do i = 1, 4
                     do j = 1, 4
                        maxc = max(maxc, abs(zdip(i,j) - z(i,j)) / abs(z(1,1)))
                     end do
                  end do
               end do
            end do
         end do
         write(*,'(a,f6.3,a,f6.4,a,f7.3)') '     a_eff=', a, ' um  lam=', lam, &
              ' um  x=', 2.0_wp*PI*a/lam
         maxe = max(maxe, chk_rel('   dipole vs T-matrix max |dZ|/Z11', maxc, 0.0_wp, 0.05_wp))
      end do
      write(*,'(a,es10.2)') '   Anchor D max error = ', maxe
   end subroutine anchor_d


   ! ==================================================================
   ! Anchor E: symmetries.
   ! ==================================================================
   subroutine anchor_e()
      integer  :: nmax, lmax, is, ipp, i, j, icand, ibest
      real(wp) :: a, lam, maxe, res, bestres
      real(wp) :: qext_o(3), qsca_o(3), qabs_o(3), qsca_rand
      real(wp) :: a1(NPL), a2(NPL), a3(NPL), a4(NPL), b1(NPL), b2(NPL)
      real(wp) :: zp(4,4), zm(4,4), z1(4,4), z2(4,4), z0(4,4), zc(4,4)
      real(wp) :: sgn(4,4), sbest(4,4)
      real(wp), parameter :: TSE(3) = (/ 30.0_wp, 90.0_wp, 140.0_wp /)
      real(wp), parameter :: PHE(2) = (/ 60.0_wp, 120.0_wp /)
      real(wp), parameter :: PHZ(4) = (/ 20.0_wp, 60.0_wp, 110.0_wp, 200.0_wp /)
      real(wp), parameter :: GTI(4) = (/ 20.0_wp, 30.0_wp, 60.0_wp, 70.0_wp /)
      real(wp), parameter :: GTS(3) = (/ 40.0_wp, 80.0_wp, 130.0_wp /)
      real(wp), parameter :: GPH(3) = (/ 30.0_wp, 70.0_wp, 110.0_wp /)
      integer  :: it, jt, kt
      logical  :: ok
      character(len=1) :: srow(4)

      write(*,'(a)') '----------------------------------------------------------------------'
      write(*,'(a)') ' Anchor E: symmetries at T-matrix point (0.10, 0.55)'
      a = 0.10_wp;  lam = 0.55_wp
      call solve_full(a, lam, nmax, qext_o, qsca_o, qabs_o, &
                      a1, a2, a3, a4, b1, b2, lmax, qsca_rand, ok)
      if (.not. ok) then
         write(*,'(a)') '   T-matrix did not converge; anchor E skipped'
         nfail_total = nfail_total + 1;  return
      end if

      ! (i) phi mirror: Z(360-phi) = Z(phi) with off-diagonal blocks flipped.
      call block_flip_signs(sgn)
      maxe = 0.0_wp
      do is = 1, 3
         do ipp = 1, 2
            call mueller_matrix_fixed_orientation(nmax, lam, 45.0_wp, TSE(is), PHE(ipp),           zp)
            call mueller_matrix_fixed_orientation(nmax, lam, 45.0_wp, TSE(is), 360.0_wp-PHE(ipp),  zm)
            res = 0.0_wp
            do i = 1, 4
               do j = 1, 4
                  res = max(res, abs(zm(i,j) - sgn(i,j)*zp(i,j)) / abs(zp(1,1)))
               end do
            end do
            maxe = max(maxe, res)
         end do
      end do
      maxe = chk_rel(' E(i)   phi mirror residual/Z11', maxe, 0.0_wp, 1.0e-5_wp)
      write(*,'(a,es10.2)') '   phi-mirror max residual = ', maxe

      ! (ii) theta_i = 0 axisymmetry.  Incidence along the axis makes the
      ! problem invariant under rotation about z, so:
      !   - Z11 is independent of phi (the unpolarized pattern is axisymmetric);
      !   - at phi = 0 (scattering plane = x-z) Z is the six-element block-
      !     diagonal matrix Zsp(theta_s);
      !   - at general phi, Z(0;theta_s,phi) = Zsp(theta_s) . R(phi), where R is
      !     the Stokes rotation carrying the fixed incident (v,h) basis into the
      !     scattering plane (the scattered meridional basis already lies in it).
      ! The literal "Z(phi1) = Z(phi2)" of the plan holds only for Z11 and only
      ! modulo the 180-degree period of R; the rotation identity is the full,
      ! physically correct statement, verified here.
      maxe = 0.0_wp
      do is = 1, 3
         call mueller_matrix_fixed_orientation(nmax, lam, 0.0_wp, TSE(is), 0.0_wp, z0)  ! Zsp
         ! Z11 phi-independence and block-diagonal Zsp.
         res = 0.0_wp
         res = max(res, abs(z0(1,3))/abs(z0(1,1)));  res = max(res, abs(z0(1,4))/abs(z0(1,1)))
         res = max(res, abs(z0(2,3))/abs(z0(1,1)));  res = max(res, abs(z0(2,4))/abs(z0(1,1)))
         res = max(res, abs(z0(3,1))/abs(z0(1,1)));  res = max(res, abs(z0(3,2))/abs(z0(1,1)))
         res = max(res, abs(z0(4,1))/abs(z0(1,1)));  res = max(res, abs(z0(4,2))/abs(z0(1,1)))
         do ipp = 1, 4
            call mueller_matrix_fixed_orientation(nmax, lam, 0.0_wp, TSE(is), PHZ(ipp), z1)
            res = max(res, abs(z1(1,1) - z0(1,1))/abs(z0(1,1)))     ! Z11 phi-independent
            call stokes_rotate_incident(z0, PHZ(ipp), z2)           ! Zsp . R(phi)
            do i = 1, 4
               do j = 1, 4
                  res = max(res, abs(z1(i,j) - z2(i,j)) / abs(z0(1,1)))
               end do
            end do
         end do
         maxe = max(maxe, res)
      end do
      maxe = chk_rel(' E(ii)  theta_i=0 axisymmetry (Zsp.R)', maxe, 0.0_wp, 1.0e-5_wp)
      write(*,'(a,es10.2)') '   theta_i=0 max residual = ', maxe

      ! (iii) equatorial mirror: derive theta_i -> 180-theta_i mapping.
      write(*,'(a)') ' E(iii) equatorial mirror: derive theta_i -> 180-theta_i mapping'
      call mueller_matrix_fixed_orientation(nmax, lam, 35.0_wp, 75.0_wp, 55.0_wp, z0)
      bestres = huge(1.0_wp);  ibest = 1;  sbest = 1.0_wp
      do icand = 1, 3
         ! candidate azimuth map: 1 -> phi, 2 -> 180-phi, 3 -> 360-phi
         call mueller_matrix_fixed_orientation(nmax, lam, 145.0_wp, 105.0_wp, &
                                               phi_candidate(icand, 55.0_wp), zc)
         ! sign pattern from the generic reference point
         do i = 1, 4
            do j = 1, 4
               if (abs(z0(i,j)) > 1.0e-3_wp*z0(1,1)) then
                  sgn(i,j) = sign(1.0_wp, zc(i,j)/z0(i,j))
               else
                  sgn(i,j) = 1.0_wp
               end if
            end do
         end do
         ! verify across a grid
         res = 0.0_wp
         do it = 1, 4
            do jt = 1, 3
               do kt = 1, 3
                  call mueller_matrix_fixed_orientation(nmax, lam, GTI(it), GTS(jt), GPH(kt), z1)
                  call mueller_matrix_fixed_orientation(nmax, lam, 180.0_wp-GTI(it), &
                       180.0_wp-GTS(jt), phi_candidate(icand, GPH(kt)), z2)
                  do i = 1, 4
                     do j = 1, 4
                        res = max(res, abs(z2(i,j) - sgn(i,j)*z1(i,j)) / abs(z1(1,1)))
                     end do
                  end do
               end do
            end do
         end do
         write(*,'(a,i0,a,es12.4)') '     candidate phi-map ', icand, ' grid residual/Z11 = ', res
         if (res < bestres) then
            bestres = res;  ibest = icand;  sbest = sgn
         end if
      end do

      write(*,'(a)') '     verified mapping: theta_s -> 180-theta_s, phi -> '// &
                     trim(phimap_name(ibest))
      write(*,'(a)') '     sign pattern S(i,j) (rows i=1..4):'
      do i = 1, 4
         do j = 1, 4
            if (sbest(i,j) >= 0.0_wp) then;  srow(j) = '+';  else;  srow(j) = '-';  end if
         end do
         write(*,'(a,4(2x,a1))') '       ', srow(1), srow(2), srow(3), srow(4)
      end do
      maxe = chk_rel(' E(iii) equatorial mapping residual', bestres, 0.0_wp, 1.0e-5_wp)
      write(*,'(a,es10.2)') '   equatorial-mapping residual = ', maxe
   end subroutine anchor_e


   ! ==================================================================
   ! Anchor F: size-integrated aligned K/Csca vs the 4-block jori table.
   ! ==================================================================
   subroutine anchor_f()
      ! At 0.55 um, the aligned extinction-matrix elements (grid-independent,
      ! from the forward amplitudes) and the closure scattering cross sections
      ! from a reduced-grid pass of accumulate_aligned_population are checked
      ! against the SAME f_align-weighted size integrals computed independently
      ! from the orientation-resolved (jori) table:
      !   Cext_al(0)  = SUM dn f Cext(jori=1)
      !   Cext_al(90) = SUM dn f (Cext2 + Cext3)/2
      !   Cpol_al(90) = SUM dn f (Cext3 - Cext2)/2
      !   Cbir_al(90) = SUM dn f qbir area
      !   Csca_al(0)  = SUM dn f Csca(jori=1)
      !   Csca_al(90) = SUM dn f (Csca2 + Csca3)/2
      ! The K elements come from the same forward-amplitude physics the table
      ! stored, so they agree to table precision (~1e-4); the closure C_sca is
      ! quadrature-limited at the test grid.
      !
      ! The q_table_jori module publishes only the qran/qpol/qbir combinations,
      ! not the per-orientation Q the jori = 1 and (2,3) references need, so the
      ! table is read here with a compact local reader in its documented format,
      ! leaving sed/ untouched.
      integer,  parameter :: NTI = 7, NTS = 37, NPH = 13
      real(wp), allocatable :: qext_j(:,:,:), qsca_j(:,:,:), qbir_j(:,:)
      real(wp), allocatable :: lam_j(:), aeff_j(:)
      integer  :: nlam, naeff
      logical  :: has_bir
      real(wp) :: ti(NTI), ts(NTS), ph(NPH)
      real(wp), allocatable :: f_al(:)
      real(wp), allocatable :: z_al(:,:,:,:,:)
      real(wp) :: cext_al(NTI), cpol_al(NTI), cbir_al(NTI), csca_al(NTI)
      real(wp) :: sacc(NPL,6), sref(NPL,6)
      real(wp) :: csca_tot, csca_ref, cext_tot, cext_ref, skipped_weight
      integer  :: lmax_acc, n_small, n_tmat, n_skip, n_fail
      real(wp) :: nr, ki, lam_used, a, area, f, w, dd
      real(wp) :: q1e, q2e, q3e, q1s, q2s, q3s, qb
      real(wp) :: cext0, cext90, cpol90, cbir90, csca0, csca90, maxe
      integer  :: jw, ia, i

      write(*,'(a)') '----------------------------------------------------------------------'
      write(*,'(a)') ' Anchor F: size-integrated aligned K/Csca vs the jori table (0.55 um)'

      call load_size_dist(f_sdist)
      call read_qjori_perori(f_qjori, f_wave, f_aeff, lam_j, aeff_j, nlam, naeff, &
                             qext_j, qsca_j, qbir_j, has_bir)
      if (.not. has_bir) then
         write(*,'(a)') '   jori table has no Q_re (birefringence) block; anchor F skipped'
         nfail_total = nfail_total + 1;  return
      end if

      ! Nearest table wavelength to 0.55 um; the driver pass uses the same lam.
      jw = 1;  dd = huge(1.0_wp)
      do i = 1, nlam
         if (abs(lam_j(i) - 0.55_wp) < dd) then;  dd = abs(lam_j(i) - 0.55_wp);  jw = i;  end if
      end do
      lam_used = lam_j(jw)
      write(*,'(a,f9.5,a)') '   using jori-table lambda = ', lam_used, ' um (nearest to 0.55)'

      ! Independent reference: f_align-weighted size integral from the table,
      ! Q interpolated in log a_eff to each size-distribution radius.
      cext0 = 0.0_wp;  cext90 = 0.0_wp;  cpol90 = 0.0_wp;  cbir90 = 0.0_wp
      csca0 = 0.0_wp;  csca90 = 0.0_wp
      do ia = 1, n_size
         if (dn_ad(ia) <= 0.0_wp) cycle
         a    = a_dist(ia);  area = PI*a*a;  f = falign_hd23(a);  w = dn_ad(ia)*f
         q1e = qinterp_loga(aeff_j, qext_j(jw,:,1), naeff, a)
         q2e = qinterp_loga(aeff_j, qext_j(jw,:,2), naeff, a)
         q3e = qinterp_loga(aeff_j, qext_j(jw,:,3), naeff, a)
         q1s = qinterp_loga(aeff_j, qsca_j(jw,:,1), naeff, a)
         q2s = qinterp_loga(aeff_j, qsca_j(jw,:,2), naeff, a)
         q3s = qinterp_loga(aeff_j, qsca_j(jw,:,3), naeff, a)
         qb  = qinterp_loga(aeff_j, qbir_j(jw,:),   naeff, a)
         cext0  = cext0  + w * q1e * area
         cext90 = cext90 + w * 0.5_wp*(q2e + q3e) * area
         cpol90 = cpol90 + w * 0.5_wp*(q3e - q2e) * area
         cbir90 = cbir90 + w * qb * area
         csca0  = csca0  + w * q1s * area
         csca90 = csca90 + w * 0.5_wp*(q2s + q3s) * area
      end do

      ! Reduced-grid driver pass over the same size distribution and lambda.
      call fill_grid(ti, 0.0_wp, 15.0_wp)
      call fill_grid(ts, 0.0_wp,  5.0_wp)
      call fill_grid(ph, 0.0_wp, 15.0_wp)
      allocate(f_al(n_size), z_al(4,4,NTI,NTS,NPH))
      do ia = 1, n_size
         f_al(ia) = falign_hd23(a_dist(ia))
      end do
      call interp_m(lam_used, nr, ki)
      call accumulate_aligned_population(lam_used, nr, ki, EPS_BA, NP_OBL, DDELT, &
              NDGS, X_SMALL, X_LARGE, a_dist, dn_ad, f_al, ti, ts, ph, &
              z_al, cext_al, cpol_al, cbir_al, csca_al, sacc, sref, lmax_acc, &
              csca_tot, csca_ref, cext_tot, cext_ref, skipped_weight, &
              n_small, n_tmat, n_skip, n_fail)

      ! K elements (grid-independent) at tol 2e-3; closure C_sca at tol 2%.
      maxe = 0.0_wp
      maxe = max(maxe, chk_rel('   Cext_al(0)   vs SUM dn f C1e ', cext_al(1),   cext0,  2.0e-3_wp))
      maxe = max(maxe, chk_rel('   Cext_al(90)  vs SUM dn f(C2e+C3e)/2', cext_al(NTI), cext90, 2.0e-3_wp))
      maxe = max(maxe, chk_rel('   Cpol_al(90)  vs SUM dn f(C3e-C2e)/2', cpol_al(NTI), cpol90, 2.0e-3_wp))
      maxe = max(maxe, chk_rel('   Cbir_al(90)  vs SUM dn f qbir area',  cbir_al(NTI), cbir90, 2.0e-3_wp))
      maxe = max(maxe, chk_rel('   Csca_al(0)   vs SUM dn f C1s  [grid]', csca_al(1),   csca0,  2.0e-2_wp))
      maxe = max(maxe, chk_rel('   Csca_al(90)  vs SUM dn f(C2s+C3s)/2 [grid]', csca_al(NTI), csca90, 2.0e-2_wp))
      write(*,'(a,i0,a,i0,a,i0,a,i0)') '   driver pass: Rayleigh ', n_small, &
         ', T-matrix ', n_tmat, ', GO-skipped ', n_skip, ', redirected ', n_fail
      write(*,'(a,es10.2)') '   Anchor F max error = ', maxe
      deallocate(f_al, z_al, qext_j, qsca_j, qbir_j, lam_j, aeff_j)
   end subroutine anchor_f


   ! ==================================================================
   ! Anchor G: the Rayleigh K sin^2 law from the forward amplitudes.
   ! ==================================================================
   subroutine anchor_g()
      ! The extinction-matrix elements the driver builds from the forward
      ! amplitudes must obey the exact dipole law
      !   C_v(theta_i) = sin^2(theta_i) C(E||a) + cos^2(theta_i) C(E perp a),
      !   C_h(theta_i) = C(E perp a),
      ! for both the imaginary part (absorption level) and the real part
      ! (birefringence level).  In the Rayleigh limit this is exact (the dipole
      ! forward amplitude is the projection e . alpha . e); at x = 0.137 the
      ! T-matrix forward amplitudes follow it to a few percent.
      real(wp), parameter :: TIG(5) = (/ 0.0_wp, 30.0_wp, 45.0_wp, 60.0_wp, 90.0_wp /)
      real(wp) :: a, lam, nr, ki, k, area, fac4, s2, c2, maxe
      real(wp) :: qext, qsca, walb, asymm
      real(wp) :: qao(3), qro(3)
      real(wp) :: cabs_a, cabs_b, cre_a, cre_b
      real(wp) :: cv, ch, rv, rh, ref_cv, ref_ch, ref_rv, ref_rh
      complex(wp) :: alpha_a, alpha_b, k2, svv, shh
      integer  :: it, nmax, lmax
      real(wp) :: qext_o(3), qsca_o(3), qabs_o(3), qsca_rand
      real(wp) :: a1(NPL), a2(NPL), a3(NPL), a4(NPL), b1(NPL), b2(NPL)
      complex(wp) :: vv, vh, hv, hh
      real(wp) :: c_a, c_b, dev
      logical  :: ok

      write(*,'(a)') '----------------------------------------------------------------------'
      write(*,'(a)') ' Anchor G(i): Rayleigh dipole forward amplitudes obey the sin^2 law'
      write(*,'(a)') '   (reference C_abs, Cre from rayleigh_limit qabs_ori/qre_ori) [tol 1e-10]'
      a = 0.008_wp;  lam = 0.55_wp
      call interp_m(lam, nr, ki)
      call spheroid_dipole_polarizability(a, nr, ki, EPS_BA, alpha_a, alpha_b)
      call rayleigh_limit(a, lam, nr, ki, EPS_BA, qext, qsca, walb, asymm, &
                          qabs_ori=qao, qre_ori=qro)
      k    = 2.0_wp*PI/lam
      k2   = cmplx(k*k, 0.0_wp, kind=wp)
      area = PI*a*a
      fac4 = 2.0_wp*lam                          ! 4 pi / k
      ! Axial (E||a) = jori 2; transverse (E perp a) = jori 3.
      cabs_a = qao(2)*area;  cabs_b = qao(3)*area
      cre_a  = qro(2)*area;  cre_b  = qro(3)*area
      maxe = 0.0_wp
      write(*,'(a,f6.4,a,f6.3)') '   a=', a, ' um  x=', 2.0_wp*PI*a/lam
      do it = 1, 5
         s2 = sin(TIG(it)*DEG)**2;  c2 = cos(TIG(it)*DEG)**2
         ! Dipole forward amplitude S_pq = k^2 (e_p . alpha . e_q); at forward
         ! e_v = (cos ti, 0, -sin ti) -> alpha_b cos^2 + alpha_a sin^2, e_h -> alpha_b.
         svv = k2 * (alpha_b*c2 + alpha_a*s2)
         shh = k2 *  alpha_b
         cv = fac4*aimag(svv);  ch = fac4*aimag(shh)
         rv = fac4*real(svv, wp);  rh = fac4*real(shh, wp)
         ref_cv = s2*cabs_a + c2*cabs_b;  ref_ch = cabs_b
         ref_rv = s2*cre_a  + c2*cre_b;   ref_rh = cre_b
         write(*,'(a,f6.1,a)') '     theta_i=', TIG(it), ' deg'
         maxe = max(maxe, chk_rel('       Im C_v = sin^2 Cabs_a + cos^2 Cabs_b', cv, ref_cv, 1.0e-10_wp))
         maxe = max(maxe, chk_rel('       Im C_h = Cabs_b                     ', ch, ref_ch, 1.0e-10_wp))
         maxe = max(maxe, chk_rel('       Re C_v = sin^2 Cre_a  + cos^2 Cre_b ', rv, ref_rv, 1.0e-10_wp))
         maxe = max(maxe, chk_rel('       Re C_h = Cre_b                      ', rh, ref_rh, 1.0e-10_wp))
      end do
      write(*,'(a,es10.2)') '   Anchor G(i) max error = ', maxe

      write(*,'(a)') ' Anchor G(ii): T-matrix forward amplitudes obey the sin^2 law at x=0.137'
      write(*,'(a)') '   C_a = C_v(90) [jori2], C_b = C_h(90) [jori3]                  [tol 5%]'
      a = 0.012_wp;  lam = 0.55_wp
      call solve_full(a, lam, nmax, qext_o, qsca_o, qabs_o, &
                      a1, a2, a3, a4, b1, b2, lmax, qsca_rand, ok)
      if (.not. ok) then
         write(*,'(a)') '   T-matrix did not converge; anchor G(ii) skipped'
         nfail_total = nfail_total + 1;  return
      end if
      fac4 = 2.0_wp*lam
      call ampl(nmax, lam, 90.0_wp, 90.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, vv, vh, hv, hh)
      c_a = fac4*aimag(vv)                        ! E||a at ti=90 (jori 2)
      c_b = fac4*aimag(hh)                        ! E perp a at ti=90 (jori 3)
      write(*,'(a,f6.4,a,f6.3)') '   a=', a, ' um  x=', 2.0_wp*PI*a/lam
      dev = 0.0_wp
      do it = 1, 5
         s2 = sin(TIG(it)*DEG)**2;  c2 = cos(TIG(it)*DEG)**2
         call ampl(nmax, lam, TIG(it), TIG(it), 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, vv, vh, hv, hh)
         cv = fac4*aimag(vv);  ch = fac4*aimag(hh)
         ref_cv = s2*c_a + c2*c_b
         dev = max(dev, abs(cv/ref_cv - 1.0_wp))
         dev = max(dev, abs(ch/c_b   - 1.0_wp))
         write(*,'(a,f6.1,a,es11.4,a,es11.4,a,es11.4)') '     ti=', TIG(it), &
            ' C_v=', cv, ' law=', ref_cv, ' C_h=', ch
      end do
      maxe = chk_rel('   T-matrix sin^2-law max deviation', dev, 0.0_wp, 5.0e-2_wp)
      write(*,'(a,es10.2)') '   Anchor G(ii) measured max deviation = ', dev
   end subroutine anchor_g


   ! ==================================================================
   ! Anchor H: meridional rotation certified at general geometry.
   ! ==================================================================
   subroutine anchor_h()
      ! For each single size (Anchor C list) and a spread of scattering
      ! geometries -- including phi > 180 (the sign flip) and a near-degenerate
      ! phi ~ 1 deg -- the orientation-averaged Mueller matrix, the ground truth
      ! that already contains the rotations,
      !   fmod = 4 pi <Z>(theta_i, theta_s, phi) / Csca_random,
      ! is compared to
      !   L(pi - sigma2) F(Theta) L(-sigma1)
      ! with F(Theta) the block-diagonal GSP scattering matrix
      ! (scatmat_from_moments, interpolated to Theta) and sigma1, sigma2 from the
      ! SAME logic mueller_matrix_total uses (meridional_rotations below is a
      ! sibling copy).  A wrong sign convention fails here: this is the definitive
      ! certification for the aligned-scattering library's random remainder.
      ! NF fixes the F(Theta) evaluation grid (0.005 deg step): the anchor
      ! geometries put Theta off any coarser grid, and linear interpolation of
      ! the oscillatory x ~ 3.4 scattering matrix on a 0.05 deg grid costs up to
      ! ~8e-5 in the normalized measure.  The ~4e-5 residual that remains at
      ! 0.005 deg is insensitive to both further F-grid refinement and
      ! orientation-quadrature doubling (identical to 7 digits), i.e. it is the
      ! common truncation floor of the two representations, not a rotation
      ! residual: a wrong sigma sign fails this anchor at O(1)-O(30).
      integer,  parameter :: NCP = 2, NG = 4, NF = 36001
      real(wp), parameter :: CP_A(NCP) = (/ 0.10_wp, 0.30_wp /)
      real(wp), parameter :: CP_L(NCP) = (/ 0.55_wp, 0.55_wp /)
      real(wp), parameter :: GTI(NG) = (/  40.0_wp,  40.0_wp,  40.0_wp, 40.0_wp /)
      real(wp), parameter :: GTS(NG) = (/  70.0_wp,  70.0_wp,  70.0_wp, 70.0_wp /)
      real(wp), parameter :: GPH(NG) = (/  30.0_wp, 120.0_wp, 250.0_wp,  1.0_wp /)
      real(wp), parameter :: TOL_H   = 5.0e-5_wp
      integer  :: ic, k, i, j, nmax, lmax, nb, na
      real(wp) :: a, lam, area, csca_rand, qsca_rand, maxe, locmax
      real(wp) :: qext_o(3), qsca_o(3), qabs_o(3)
      real(wp) :: a1(NPL), a2(NPL), a3(NPL), a4(NPL), b1(NPL), b2(NPL)
      real(wp) :: thf(NF), f11(NF), f22(NF), f33(NF), f44(NF), f12(NF), f34(NF)
      real(wp) :: zbar(4,4), fmod(4,4), fmat(4,4), rhs(4,4), l1(4,4), l2(4,4)
      real(wp) :: cos_th, big_theta, sigma1, sigma2, ff(6), dummy
      logical  :: ok

      write(*,'(a)') '----------------------------------------------------------------------'
      write(*,'(a)') ' Anchor H: 4 pi <Z>(ti,ts,phi)/Csca_random vs L(pi-s2) F(Theta) L(-s1)'
      write(*,'(a,es8.1,a)') '   general geometry; certifies sigma1, sigma2 signs [tol ', TOL_H, ' norm]'
      maxe = 0.0_wp
      do ic = 1, NCP
         a = CP_A(ic);  lam = CP_L(ic)
         call solve_full(a, lam, nmax, qext_o, qsca_o, qabs_o, &
                         a1, a2, a3, a4, b1, b2, lmax, qsca_rand, ok)
         if (.not. ok) then
            write(*,'(a,f6.3,a,f6.3)') '   T-matrix did not converge at a=', a, ' lam=', lam
            nfail_total = nfail_total + 1;  cycle
         end if
         area      = PI*a*a
         csca_rand = qsca_rand*area
         nb = 2*nmax + 2;  na = 4*nmax + 4
         call scatmat_from_moments(a1, a2, a3, a4, b1, b2, lmax, NF, &
                                   thf, f11, f22, f33, f44, f12, f34)
         write(*,'(a,f6.3,a,f6.4,a,i0,a,i0,a,i0)') '   a_eff=', a, ' um  lam=', lam, &
              ' um   NMAX=', nmax, '   nBETA=', nb, '  nALPHA=', na
         do k = 1, NG
            cos_th = cos(GTI(k)*DEG)*cos(GTS(k)*DEG) &
                   + sin(GTI(k)*DEG)*sin(GTS(k)*DEG)*cos(GPH(k)*DEG)
            cos_th    = max(-1.0_wp, min(1.0_wp, cos_th))
            big_theta = acos(cos_th)/DEG

            ! Ground truth: orientation-averaged Z at the general geometry.
            call orient_average_z_general(nmax, lam, GTI(k), GTS(k), GPH(k), nb, na, zbar)
            fmod = 4.0_wp*PI*zbar/csca_rand

            ! Reference: block-diagonal F(Theta) rotated into the meridional bases.
            call interp_scatmat(thf, f11, f22, f33, f44, f12, f34, NF, big_theta, ff)
            fmat = 0.0_wp
            fmat(1,1) = ff(1);  fmat(2,2) = ff(2);  fmat(3,3) = ff(3);  fmat(4,4) = ff(4)
            fmat(1,2) = ff(5);  fmat(2,1) = ff(5)
            fmat(3,4) = ff(6);  fmat(4,3) = -ff(6)
            call meridional_rotations(GTI(k), GTS(k), GPH(k), sigma1, sigma2)
            call stokes_L(PI - sigma2, l2)
            call stokes_L(-sigma1,     l1)
            rhs = matmul(l2, matmul(fmat, l1))

            locmax = 0.0_wp
            do j = 1, 4
               do i = 1, 4
                  locmax = max(locmax, helem(fmod(i,j), rhs(i,j), ff(1)))
               end do
            end do
            maxe = max(maxe, locmax)
            write(*,'(a,f5.1,a,f5.1,a,f6.1,a,f6.1,a,es10.2)') '     ti=', GTI(k), &
                 ' ts=', GTS(k), ' phi=', GPH(k), ' Theta=', big_theta, &
                 ' deg  max err=', locmax
         end do
      end do
      dummy = chk_rel('   Anchor H max normalized error', maxe, 0.0_wp, TOL_H)
   end subroutine anchor_h


   real(wp) function helem(val, ref, f11r) result(e)
      ! Anchor-H element error: relative where |ref| > 1e-3 F11, else absolute
      ! normalized to F11 -- the same measure as chk_elem, returned silently so a
      ! single verdict covers the whole matrix.
      real(wp), intent(in) :: val, ref, f11r
      if (abs(ref) > 1.0e-3_wp*f11r) then
         e = abs(val/ref - 1.0_wp)
      else
         e = abs(val - ref) / f11r
      end if
   end function helem


   subroutine interp_scatmat(th, f11, f22, f33, f44, f12, f34, n, target, ff)
      ! Linear interpolation of the six GSP scattering-matrix elements to the
      ! scattering angle target [deg] on the uniform grid th (0..180, n points).
      real(wp), intent(in)  :: th(:), f11(:), f22(:), f33(:), f44(:), f12(:), f34(:)
      integer,  intent(in)  :: n
      real(wp), intent(in)  :: target
      real(wp), intent(out) :: ff(6)
      integer  :: lo
      real(wp) :: step, t
      step = th(2) - th(1)
      lo = int(target/step) + 1
      if (lo < 1)   lo = 1
      if (lo > n-1) lo = n-1
      t = (target - th(lo)) / (th(lo+1) - th(lo))
      ff(1) = (1.0_wp-t)*f11(lo) + t*f11(lo+1)
      ff(2) = (1.0_wp-t)*f22(lo) + t*f22(lo+1)
      ff(3) = (1.0_wp-t)*f33(lo) + t*f33(lo+1)
      ff(4) = (1.0_wp-t)*f44(lo) + t*f44(lo+1)
      ff(5) = (1.0_wp-t)*f12(lo) + t*f12(lo+1)
      ff(6) = (1.0_wp-t)*f34(lo) + t*f34(lo+1)
   end subroutine interp_scatmat


   pure subroutine meridional_rotations(theta_i_deg, theta_s_deg, phi_deg, sigma1, sigma2)
      ! Sibling copy of scatmat_aligned_mod's meridional_scattering_angles (this
      ! repo keeps sibling copies rather than sharing across the sed/tmatrix split):
      ! the meridional-to-scattering-plane Stokes rotation angles [rad] for
      ! Z = L(pi - sigma2) F(Theta) L(-sigma1), incidence azimuth 0.  Closed
      ! forms handle the single poles without branches; signs and degenerate
      ! limits documented at the library site, certified by Anchor H here.
      real(wp), intent(in)  :: theta_i_deg, theta_s_deg, phi_deg
      real(wp), intent(out) :: sigma1, sigma2
      real(wp), parameter :: EPS = 1.0e-30_wp
      real(wp) :: ti, ts, ph, sti, cti, sts, cts, sph, cph
      real(wp) :: y1, x1, y2, x2
      ti = theta_i_deg*DEG;  ts = theta_s_deg*DEG;  ph = phi_deg*DEG
      sti = sin(ti);  cti = cos(ti)
      sts = sin(ts);  cts = cos(ts)
      sph = sin(ph);  cph = cos(ph)
      y1 =  sts*sph;   x1 = cti*sts*cph - sti*cts
      y2 = -sti*sph;   x2 = cti*sts - sti*cts*cph
      if (abs(y1) + abs(x1) < EPS) then
         sigma1 = 0.0_wp
      else
         sigma1 = atan2(y1, x1)
      end if
      if (abs(y2) + abs(x2) < EPS) then
         sigma2 = 0.0_wp
      else
         sigma2 = atan2(y2, x2)
      end if
   end subroutine meridional_rotations


   pure subroutine stokes_L(angle, l)
      ! Stokes rotation L(angle) in the Mishchenko convention (sibling copy).
      real(wp), intent(in)  :: angle
      real(wp), intent(out) :: l(4,4)
      real(wp) :: c2, s2
      c2 = cos(2.0_wp*angle);  s2 = sin(2.0_wp*angle)
      l = 0.0_wp
      l(1,1) = 1.0_wp;  l(4,4) = 1.0_wp
      l(2,2) = c2;  l(2,3) =  s2
      l(3,2) = -s2; l(3,3) =  c2
   end subroutine stokes_L


   ! ==================================================================
   ! jori-table reader (per-orientation) and log-a interpolation
   ! ==================================================================
   subroutine read_qjori_perori(qfile, wavefile, aefffile, lam_j, aeff_j, &
                                nlam, naeff, qext_j, qsca_j, qbir_j, has_bir)
      ! Read the orientation-resolved DH21 Q table (its documented format:
      ! 12 header lines, then 3 quantities [ext, abs, sca] x 3 orientations x
      ! nlam records of naeff sizes, then an optional 4th Q_re block) into the
      ! per-orientation Q arrays anchor F needs.  qbir_j = 0.5*(Qre3 - Qre2) is
      ! formed from the Q_re block; has_bir is .false. for an older 3-block file.
      character(len=*),      intent(in)  :: qfile, wavefile, aefffile
      real(wp), allocatable, intent(out) :: lam_j(:), aeff_j(:)
      integer,               intent(out) :: nlam, naeff
      real(wp), allocatable, intent(out) :: qext_j(:,:,:), qsca_j(:,:,:), qbir_j(:,:)
      logical,               intent(out) :: has_bir
      integer, parameter :: NW = 1129, NA = 169, NHEAD = 12
      real(wp), allocatable :: qre_j(:,:,:), row(:)
      character(len=256) :: scratch
      character(len=512) :: line
      integer :: u, ios, iq, jori, jw, i, estat, cstat

      nlam = NW;  naeff = NA
      allocate(lam_j(NW), aeff_j(NA), row(NA))
      allocate(qext_j(NW,NA,3), qsca_j(NW,NA,3), qre_j(NW,NA,3))

      call read_axis(wavefile, NW, lam_j)
      call read_axis(aefffile, NA, aeff_j)

      scratch = 'q_jori_anchor_scratch.dat'
      call execute_command_line('gzip -dc "'//trim(qfile)//'" > "'//trim(scratch)//'"', &
                                exitstat=estat, cmdstat=cstat)
      if (cstat /= 0 .or. estat /= 0) then
         write(*,'(a,a)') '   ERROR: gzip -dc failed on ', trim(qfile)
         nfail_total = nfail_total + 1;  has_bir = .false.;  return
      end if

      open(newunit=u, file=trim(scratch), status='old', action='read')
      do i = 1, NHEAD
         read(u,'(a)') line
      end do
      ! ext, abs, sca (abs discarded), each 3 orientations x NW records.
      do iq = 1, 3
         do jori = 1, 3
            do jw = 1, NW
               read(u,*) row(1:NA)
               if (iq == 1) qext_j(jw,:,jori) = row(1:NA)
               if (iq == 3) qsca_j(jw,:,jori) = row(1:NA)
            end do
         end do
      end do
      ! Optional 4th (Q_re) block.
      read(u,*,iostat=ios) row(1:NA)
      if (is_iostat_end(ios)) then
         has_bir = .false.
      else
         qre_j(1,:,1) = row(1:NA)
         do jori = 1, 3
            do jw = 1, NW
               if (jori == 1 .and. jw == 1) cycle
               read(u,*) row(1:NA)
               qre_j(jw,:,jori) = row(1:NA)
            end do
         end do
         has_bir = .true.
      end if
      close(u, status='delete')

      allocate(qbir_j(NW,NA))
      if (has_bir) then
         qbir_j = 0.5_wp * (qre_j(:,:,3) - qre_j(:,:,2))
      else
         qbir_j = 0.0_wp
      end if
      deallocate(qre_j, row)
   end subroutine read_qjori_perori


   subroutine read_axis(fname, n, arr)
      ! DH21_wave / DH21_aeff: two title lines then n free-format values.
      character(len=*), intent(in)  :: fname
      integer,          intent(in)  :: n
      real(wp),         intent(out) :: arr(n)
      integer :: u
      character(len=512) :: line
      open(newunit=u, file=fname, status='old', action='read')
      read(u,'(a)') line;  read(u,'(a)') line
      read(u,*) arr(1:n)
      close(u)
   end subroutine read_axis


   real(wp) function qinterp_loga(aeff, q, n, a) result(qi)
      ! Linear interpolation of Q(aeff) at a, in log10(a_eff), clamped at ends.
      real(wp), intent(in) :: aeff(:), q(:), a
      integer,  intent(in) :: n
      real(wp) :: t, x
      integer  :: lo, hi, mid
      if (a <= aeff(1)) then
         qi = q(1);  return
      else if (a >= aeff(n)) then
         qi = q(n);  return
      end if
      lo = 1;  hi = n
      do while (hi - lo > 1)
         mid = (lo + hi)/2
         if (aeff(mid) <= a) then;  lo = mid;  else;  hi = mid;  end if
      end do
      x = log10(a)
      t = (x - log10(aeff(lo))) / (log10(aeff(hi)) - log10(aeff(lo)))
      qi = (1.0_wp - t)*q(lo) + t*q(hi)
   end function qinterp_loga


   subroutine fill_grid(x, a0, step)
      ! Fill x(1:size(x)) with a0, a0+step, ...
      real(wp), intent(out) :: x(:)
      real(wp), intent(in)  :: a0, step
      integer :: i
      do i = 1, size(x)
         x(i) = a0 + real(i-1, wp)*step
      end do
   end subroutine fill_grid


   ! ==================================================================
   ! shared solvers and quadratures
   ! ==================================================================
   subroutine solve_full(a, lam, nmax_tm, qext_ori, qsca_ori, qabs_ori, &
                         a1, a2, a3, a4, b1, b2, lmax, qsca_rand, ok)
      ! One (a_eff, lam) solve serving every anchor: the oriented cross
      ! sections (reference for A/B) from tmatrix_oriented_cross, and the
      ! truncation order NMAX_TM, random-orientation Q_sca, and GSP expansion
      ! coefficients from TMD_ONE_SCATMAT.  Both are the identical T-matrix
      ! solve, so the second leaves COMMON /TMAT/ valid with the returned
      ! NMAX_TM for the subsequent AMPL evaluations.
      real(wp), intent(in)  :: a, lam
      integer,  intent(out) :: nmax_tm, lmax
      real(wp), intent(out) :: qext_ori(3), qsca_ori(3), qabs_ori(3)
      real(wp), intent(out) :: a1(NPL), a2(NPL), a3(NPL), a4(NPL), b1(NPL), b2(NPL)
      real(wp), intent(out) :: qsca_rand
      logical,  intent(out) :: ok
      real(wp) :: nr, ki, qext, walb, asymm
      complex(wp) :: m
      integer  :: ierr, ierr2

      call interp_m(lam, nr, ki)
      m = cmplx(nr, ki, kind=wp)
      call tmatrix_oriented_cross(a, lam, m, EPS_BA, NP_OBL, DDELT, NDGS, &
                                  qext_ori, qsca_ori, qabs_ori, ierr)
      call tmd_one_scatmat(a, lam, nr, ki, EPS_BA, NP_OBL, DDELT, NDGS, &
                           qext, qsca_rand, walb, asymm, &
                           a1, a2, a3, a4, b1, b2, lmax, ierr2, nmax_tm)
      ok = (ierr == 0 .and. ierr2 == 0)
   end subroutine solve_full


   subroutine sphere_integrals(nmax, lam, theta_i, i11, i11p12, i11m12)
      ! INT Z11, INT (Z11+Z12), INT (Z11-Z12) over the scattering sphere,
      ! on the SAME exact quadrature as scatter_sphere_integral in
      ! driver/tmatrix_oriented.f90: Gauss-Legendre in cos(theta_s) with
      ! N_theta = NMAX+2 nodes, uniform phi with N_phi = 2*NMAX+2 points.
      integer,  intent(in)  :: nmax
      real(wp), intent(in)  :: lam, theta_i
      real(wp), intent(out) :: i11, i11p12, i11m12
      integer  :: nth, nph, it, ip
      real(wp), allocatable :: xg(:), wg(:)
      real(wp) :: wphi, ts_deg, phi_deg, wt, z(4,4)

      nth = nmax + 2;  nph = 2*nmax + 2
      allocate(xg(nth), wg(nth))
      call gauss(nth, 0, 0, xg, wg)
      wphi = 2.0_wp*PI / real(nph, kind=wp)
      i11 = 0.0_wp;  i11p12 = 0.0_wp;  i11m12 = 0.0_wp
      do it = 1, nth
         ts_deg = acos(xg(it)) / DEG
         do ip = 1, nph
            phi_deg = real(ip-1, kind=wp) * 360.0_wp / real(nph, kind=wp)
            call mueller_matrix_fixed_orientation(nmax, lam, theta_i, ts_deg, phi_deg, z)
            wt = wg(it) * wphi
            i11    = i11    + wt *  z(1,1)
            i11p12 = i11p12 + wt * (z(1,1) + z(1,2))
            i11m12 = i11m12 + wt * (z(1,1) - z(1,2))
         end do
      end do
      deallocate(xg, wg)
   end subroutine sphere_integrals


   subroutine orient_average_z(nmax, lam, theta_deg, nb, na, zbar)
      ! Orientation average of the fixed-orientation Mueller matrix, incidence
      ! along z (TL = 0, PL = 0) and scattering at (theta_deg, 0), over the
      ! grain axis direction: Gauss-Legendre in cos(BETA) with nb nodes and
      ! uniform ALPHA with na nodes.  The averaged matrix depends only on the
      ! incidence-to-scattering angle theta_deg and, divided by Csca/4pi,
      ! reproduces the random-orientation F(Theta).
      integer,  intent(in)  :: nmax, nb, na
      real(wp), intent(in)  :: lam, theta_deg
      real(wp), intent(out) :: zbar(4,4)
      real(wp), allocatable :: bx(:), bw(:)
      integer  :: ib, ia
      real(wp) :: beta_deg, alpha_deg, wb, z(4,4)
      complex(wp) :: vv, vh, hv, hh

      allocate(bx(nb), bw(nb))
      call gauss(nb, 0, 0, bx, bw)              ! cos(BETA) in [-1,1], weights sum to 2
      zbar = 0.0_wp
      do ib = 1, nb
         beta_deg = acos(bx(ib)) / DEG
         wb = 0.5_wp * bw(ib)                    ! normalize the cos(BETA) average
         do ia = 1, na
            alpha_deg = real(ia-1, kind=wp) * 360.0_wp / real(na, kind=wp)
            call ampl(nmax, lam, 0.0_wp, theta_deg, 0.0_wp, 0.0_wp, &
                      alpha_deg, beta_deg, vv, vh, hv, hh)
            call mueller_from_amplitude(vv, vh, hv, hh, z)
            zbar = zbar + (wb / real(na, kind=wp)) * z
         end do
      end do
      deallocate(bx, bw)
   end subroutine orient_average_z


   subroutine orient_average_z_general(nmax, lam, theta_i, theta_s, phi, nb, na, zbar)
      ! Orientation average of the fixed-orientation Mueller matrix at ARBITRARY
      ! incidence (theta_i, azimuth 0) and scattering (theta_s, phi), over the
      ! grain axis direction: Gauss-Legendre in cos(BETA) with nb nodes and
      ! uniform ALPHA with na nodes (the same scheme as orient_average_z, which
      ! is the theta_i = 0, phi = 0 special case).  The ensemble average obeys
      ! the rotation identity Z = L(pi - sigma2) F(Theta) L(-sigma1), which
      ! Anchor H checks.
      integer,  intent(in)  :: nmax, nb, na
      real(wp), intent(in)  :: lam, theta_i, theta_s, phi
      real(wp), intent(out) :: zbar(4,4)
      real(wp), allocatable :: bx(:), bw(:)
      integer  :: ib, ia
      real(wp) :: beta_deg, alpha_deg, wb, z(4,4)
      complex(wp) :: vv, vh, hv, hh

      allocate(bx(nb), bw(nb))
      call gauss(nb, 0, 0, bx, bw)
      zbar = 0.0_wp
      do ib = 1, nb
         beta_deg = acos(bx(ib)) / DEG
         wb = 0.5_wp * bw(ib)
         do ia = 1, na
            alpha_deg = real(ia-1, kind=wp) * 360.0_wp / real(na, kind=wp)
            call ampl(nmax, lam, theta_i, theta_s, 0.0_wp, phi, &
                      alpha_deg, beta_deg, vv, vh, hv, hh)
            call mueller_from_amplitude(vv, vh, hv, hh, z)
            zbar = zbar + (wb / real(na, kind=wp)) * z
         end do
      end do
      deallocate(bx, bw)
   end subroutine orient_average_z_general


   subroutine mueller_from_amplitude(s11, s12, s21, s22, z)
      ! Mueller matrix from the amplitude matrix, the same verbatim Mishchenko
      ! bilinears as phase_matrix_from_amplitude in the engine module; kept
      ! local so the orientation average can form Z at arbitrary grain Euler
      ! angles (ALPHA, BETA), which the two-angle engine interface fixes to 0.
      complex(wp), intent(in)  :: s11, s12, s21, s22
      real(wp),    intent(out) :: z(4,4)
      complex(wp), parameter :: CI = (0.0_wp, 1.0_wp)
      z(1,1) = 0.5_wp*real( s11*conjg(s11)+s12*conjg(s12)+s21*conjg(s21)+s22*conjg(s22), kind=wp)
      z(1,2) = 0.5_wp*real( s11*conjg(s11)-s12*conjg(s12)+s21*conjg(s21)-s22*conjg(s22), kind=wp)
      z(1,3) = real(-s11*conjg(s12)-s22*conjg(s21), kind=wp)
      z(1,4) = real(CI*(s11*conjg(s12)-s22*conjg(s21)), kind=wp)
      z(2,1) = 0.5_wp*real( s11*conjg(s11)+s12*conjg(s12)-s21*conjg(s21)-s22*conjg(s22), kind=wp)
      z(2,2) = 0.5_wp*real( s11*conjg(s11)-s12*conjg(s12)-s21*conjg(s21)+s22*conjg(s22), kind=wp)
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
   end subroutine mueller_from_amplitude


   subroutine stokes_rotate_incident(zsp, phi_deg, z)
      ! Apply the Stokes rotation R(phi) to the incident side of the
      ! scattering-plane matrix Zsp: z = Zsp . R(phi).  R rotates the fixed
      ! incident (v,h) meridional basis into the scattering plane, which for
      ! incidence along the axis makes the angle phi with it:
      !   R(phi) = [ 1     0        0      0 ]
      !            [ 0   cos2phi  -sin2phi 0 ]
      !            [ 0   sin2phi   cos2phi 0 ]
      !            [ 0     0        0      1 ]
      real(wp), intent(in)  :: zsp(4,4), phi_deg
      real(wp), intent(out) :: z(4,4)
      real(wp) :: r(4,4), c, s
      c = cos(2.0_wp*phi_deg*DEG)
      s = sin(2.0_wp*phi_deg*DEG)
      r = 0.0_wp
      r(1,1) = 1.0_wp;  r(4,4) = 1.0_wp
      r(2,2) = c;  r(2,3) = -s
      r(3,2) = s;  r(3,3) =  c
      z = matmul(zsp, r)
   end subroutine stokes_rotate_incident


   subroutine block_flip_signs(sgn)
      ! Sign pattern that flips the two off-diagonal 2x2 blocks (elements
      ! 13,14,23,24,31,32,41,42) and keeps the two diagonal blocks.
      real(wp), intent(out) :: sgn(4,4)
      integer :: i, j
      do i = 1, 4
         do j = 1, 4
            if ((i <= 2 .and. j >= 3) .or. (i >= 3 .and. j <= 2)) then
               sgn(i,j) = -1.0_wp
            else
               sgn(i,j) = 1.0_wp
            end if
         end do
      end do
   end subroutine block_flip_signs


   real(wp) function phi_candidate(icand, phi) result(pm)
      integer,  intent(in) :: icand
      real(wp), intent(in) :: phi
      select case (icand)
      case (1);  pm = phi
      case (2);  pm = 180.0_wp - phi
      case default;  pm = 360.0_wp - phi
      end select
   end function phi_candidate


   function phimap_name(icand) result(nm)
      integer, intent(in) :: icand
      character(len=16) :: nm
      select case (icand)
      case (1);  nm = 'phi'
      case (2);  nm = '180-phi'
      case default;  nm = '360-phi'
      end select
   end function phimap_name


   ! ==================================================================
   ! check reporters
   ! ==================================================================
   real(wp) function chk_rel(label, val, ref, tol) result(err)
      ! Relative check |val/ref - 1| <= tol (absolute |val-ref| when ref = 0).
      character(len=*), intent(in) :: label
      real(wp),         intent(in) :: val, ref, tol
      character(len=4) :: verdict
      if (ref /= 0.0_wp) then
         err = abs(val/ref - 1.0_wp)
      else
         err = abs(val - ref)
      end if
      if (err <= tol) then
         verdict = 'PASS'
      else
         verdict = 'FAIL';  nfail_total = nfail_total + 1
      end if
      write(*,'(a,a,es14.6,a,es14.6,a,es9.2,2x,a4)') label, ' : val=', val, &
           ' ref=', ref, ' err=', err, verdict
   end function chk_rel


   real(wp) function chk_elem(label, val, ref, f11ref) result(err)
      ! Anchor-C element check: relative tol 1e-5 where |ref| > 1e-3*F11,
      ! else absolute tol 1e-5*F11.  Returns the normalized error.
      character(len=*), intent(in) :: label
      real(wp),         intent(in) :: val, ref, f11ref
      real(wp) :: tol
      character(len=4) :: verdict
      if (abs(ref) > 1.0e-3_wp*f11ref) then
         err = abs(val/ref - 1.0_wp);  tol = 1.0e-5_wp
      else
         err = abs(val - ref) / f11ref;  tol = 1.0e-5_wp
      end if
      if (err <= tol) then
         verdict = 'PASS'
      else
         verdict = 'FAIL';  nfail_total = nfail_total + 1
      end if
      write(*,'(a,a,es14.6,a,es14.6,a,es9.2,2x,a4)') label, ' : val=', val, &
           ' ref=', ref, ' err=', err, verdict
   end function chk_elem

end program compare_scatmat_aligned
