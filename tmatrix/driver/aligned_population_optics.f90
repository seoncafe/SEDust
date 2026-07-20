module aligned_population_optics
   ! Size-integrated scattering and extinction optics of a partially aligned
   ! astrodust population, at one wavelength.  This is the SHARED PHYSICS of
   ! the aligned-grain products: the CLI/I/O driver run_scatmat_aligned.f90 and
   ! the verification program compare_scatmat_aligned.f90 both call it, so the
   ! numbers a verification anchor checks are produced by the exact code the
   ! production file is written from.
   !
   ! PARTIAL-ALIGNMENT DECOMPOSITION
   ! A fraction f_align(a) of the grains of radius a is treated as perfectly
   ! aligned (symmetry axis along z = the local alignment axis) and the
   ! remainder 1 - f_align(a) as randomly oriented.  For a size distribution
   ! n(a) da (already bin-integrated in dn) the products accumulated here are,
   ! all per H:
   !
   !   Z_al(theta_i; theta_s, phi) = SUM dn f(a) Z(a; theta_i, theta_s, phi)
   !                                 the aligned phase matrix, 16 elements,
   !                                 um^2 sr^-1;
   !   K elements on the theta_i grid, from the forward amplitudes:
   !     C_ext_al(theta_i) = SUM dn f(a) (C_ext_v + C_ext_h)/2
   !     C_pol_al(theta_i) = SUM dn f(a) (C_ext_h - C_ext_v)/2   (dichroism)
   !     C_bir_al(theta_i) = SUM dn f(a) (C_re_h  - C_re_v )/2   (birefringence)
   !   F_tot(Theta) : GSP coefficients weighted dn C_sca over EVERY size bin,
   !                  the randomly oriented matrix identical to run_scatmat.f90;
   !   F_ref(Theta) : GSP coefficients weighted dn f(a) C_sca, the aligned
   !                  population's random-orientation twin;
   !   C_sca_al(theta_i) : INT Z_al11 dOmega, by quadrature over the stored
   !                  angular grid (the closure scattering cross section).
   !
   ! ETA CONTRACT (used by the radiative-transfer host).  With a cell-level
   ! alignment scale eta the aligned matrix scales linearly, Z_al,cell =
   ! eta Z_al, and the unaligned remainder is F_tot - eta F_ref.  The K matrix
   ! scales the same way, the unaligned population adding the isotropic C_ext.
   ! The linearity is exact in f_align, so the three stored integrals give the
   ! optics of any eta with no re-integration.
   !
   ! GEOMETRY AND STOKES BASIS
   ! theta_i is the polar angle of incidence from the alignment axis, theta_s
   ! and phi the scattering polar/azimuth angles, all in degrees.  Stokes
   ! vectors use Mishchenko's meridional (v, h) = (theta-hat, phi-hat) basis,
   ! Q = I_v - I_h, exactly as scattering_matrix_oriented.f90.  v-polarization
   ! at theta_i = 90 is jori = 2 (E parallel to axis), h is jori = 3 (E
   ! transverse), so at theta_i = 90 the K elements land on the tree's jori
   ! convention C_pol = 0.5 (C3 - C2), C_bir = 0.5 (Cre3 - Cre2).
   !
   ! SIZE-PARAMETER REGIMES, x = 2 pi a / lambda, matching run_scatmat.f90:
   !   x < 0.1        analytic electric-dipole matrix (Rayleigh); the K
   !                  elements follow the exact dipole law
   !                    C_v(theta_i) = sin^2(theta_i) C(E||a)
   !                                 + cos^2(theta_i) C(E perp a),
   !                    C_h(theta_i) = C(E perp a),
   !                  built from rayleigh_limit's oriented C_ext(jori) and
   !                  their birefringence twin qre_ori.
   !   0.1 <= x <= 50 one TMD_ONE_SCATMAT solve serves all products: its GSP
   !                  coefficients feed F_tot/F_ref, the /TMAT/ it leaves valid
   !                  feeds the AMPL Z loop, and its forward amplitudes give K.
   !   x > 50         no oriented matrix exists (geometric optics).  The bin is
   !                  omitted from the aligned products and enters F_tot as its
   !                  random-orientation geometric-optics matrix (so it is
   !                  scattered as unaligned).  Its C_sca is tallied in
   !                  skipped_weight; the caller stops if the skipped fraction
   !                  of the total scattering weight is not negligible.
   ! A T-matrix non-convergence redirects exactly as run_scatmat.f90: x < 1 to
   ! the Rayleigh path (fully aligned + random), x >= 1 to geometric optics
   ! (F_tot only, tallied in skipped_weight).
   !
   ! THREADING.  The only parallel region is the per-size Z-node loop in
   ! oriented_mueller_grid.  AMPL and VIGAMPL (src/ampl_oriented.f) only READ
   ! COMMON /TMAT/ and write local automatic arrays, so distinct nodes are
   ! independent; each grid slot is written once, making the result identical
   ! for any thread count.  The TMD_ONE_SCATMAT solve stays serial.

   use, intrinsic :: iso_fortran_env, only: real64, int64
   use constants, only: wp
   use asymptotic_optics, only: rayleigh_limit, geometric_optics_limit
   use scattering_matrix_oriented, only: mueller_matrix_fixed_orientation, &
                                         rayleigh_mueller_matrix_oriented
   implicit none
   private
   public :: accumulate_aligned_population, oriented_mueller_grid
   public :: scattering_cross_section_from_grid

   real(wp), parameter :: PI  = acos(-1.0_wp)
   real(wp), parameter :: DEG = acos(-1.0_wp) / 180.0_wp

contains

   subroutine accumulate_aligned_population(lam, nr, ki, eps_ba, np, ddelt, ndgs, &
                       x_small, x_large, a_arr, dn_arr, f_arr, ti, ts, ph, &
                       z_al, cext_al, cpol_al, cbir_al, csca_al_grid, &
                       sacc, sref, lmax_acc, &
                       csca_tot, csca_ref, cext_tot, cext_ref, skipped_weight, &
                       n_small, n_tmat, n_skip, n_fail, t_tmat, t_node)
      ! One wavelength: sweep the size distribution and accumulate every
      ! aligned-population product listed in the module header.  All cross
      ! sections are in um^2 (a_arr in microns), all sums per H (dn_arr already
      ! bin-integrated).
      !
      ! INPUT
      !   lam, nr, ki      wavelength [um] and refractive index at lam
      !   eps_ba           axis ratio b/a (Mishchenko; > 1 oblate)
      !   np, ddelt, ndgs  T-matrix shape flag, tolerance, quadrature multiplier
      !   x_small, x_large size-parameter regime bounds (0.1, 50)
      !   a_arr(:)         effective radii [um]
      !   dn_arr(:)        number per H in each size bin
      !   f_arr(:)         alignment fraction f_align(a) per size bin
      !   ti(:), ts(:), ph(:)   incidence, scattering, azimuth grids [deg];
      !                    ti in [0,90], ts in [0,180], ph in [0,180]
      ! OUTPUT
      !   z_al(4,4,nti,nts,nph)         aligned phase matrix [um^2 sr^-1 per H]
      !   cext_al(nti), cpol_al(nti), cbir_al(nti)   K elements [um^2 per H]
      !   csca_al_grid(nti)             INT Z_al11 dOmega, grid closure [um^2/H]
      !   sacc(npl,6), sref(npl,6)      un-normalized GSP accumulators for F_tot
      !                                 (dn C_sca) and F_ref (dn f C_sca); the
      !                                 six columns are alpha1..4, beta1, beta2
      !   lmax_acc                      1 + highest expansion order reached
      !   csca_tot, cext_tot            SUM dn C_sca, SUM dn C_ext [um^2 per H]
      !   csca_ref, cext_ref            SUM dn f C_sca, SUM dn f C_ext
      !   skipped_weight                SUM dn C_sca over bins omitted from the
      !                                 aligned products (x > 50 and T-matrix
      !                                 failures at x >= 1)
      !   n_small, n_tmat, n_skip, n_fail   regime counts
      !   t_tmat, t_node   optional wall seconds spent in T-matrix solves and in
      !                    the Z-node loop; used by the driver's run-time estimate
      real(wp), intent(in)  :: lam, nr, ki, eps_ba, ddelt, x_small, x_large
      integer,  intent(in)  :: np, ndgs
      real(wp), intent(in)  :: a_arr(:), dn_arr(:), f_arr(:)
      real(wp), intent(in)  :: ti(:), ts(:), ph(:)
      real(wp), intent(out) :: z_al(:,:,:,:,:)
      real(wp), intent(out) :: cext_al(:), cpol_al(:), cbir_al(:), csca_al_grid(:)
      real(wp), intent(out) :: sacc(:,:), sref(:,:)
      integer,  intent(out) :: lmax_acc
      real(wp), intent(out) :: csca_tot, csca_ref, cext_tot, cext_ref, skipped_weight
      integer,  intent(out) :: n_small, n_tmat, n_skip, n_fail
      real(wp), intent(out), optional :: t_tmat, t_node

      integer  :: nsz, nti, nts, nph, npl
      integer  :: ia, lmax, ierr, nmax_tm
      real(wp) :: a, area, x, f, dn, qext, qsca, walb, asymm, csca, cext
      real(wp) :: qeo(3), qro(3)
      real(wp), allocatable :: al1(:), al2(:), al3(:), al4(:), be1(:), be2(:)
      real(wp), allocatable :: z_grid(:,:,:,:,:)
      real(wp) :: t_tm, t_nd
      integer(int64) :: c0, c1, crate
      external :: tmd_one_scatmat

      nsz = size(a_arr);  nti = size(ti);  nts = size(ts);  nph = size(ph)
      npl = size(sacc, 1)

      z_al = 0.0_wp
      cext_al = 0.0_wp;  cpol_al = 0.0_wp;  cbir_al = 0.0_wp;  csca_al_grid = 0.0_wp
      sacc = 0.0_wp;  sref = 0.0_wp
      lmax_acc = 0
      csca_tot = 0.0_wp;  csca_ref = 0.0_wp
      cext_tot = 0.0_wp;  cext_ref = 0.0_wp
      skipped_weight = 0.0_wp
      n_small = 0;  n_tmat = 0;  n_skip = 0;  n_fail = 0
      t_tm = 0.0_wp;  t_nd = 0.0_wp
      qeo = 0.0_wp;  qro = 0.0_wp

      allocate(al1(npl), al2(npl), al3(npl), al4(npl), be1(npl), be2(npl))
      allocate(z_grid(4,4,nti,nts,nph))
      call system_clock(count_rate=crate)

      do ia = 1, nsz
         if (dn_arr(ia) <= 0.0_wp) cycle
         a    = a_arr(ia)
         area = PI * a * a
         x    = 2.0_wp * PI * a / lam
         f    = f_arr(ia)
         dn   = dn_arr(ia)

         if (x < x_small) then
            ! --- Rayleigh dipole: aligned Z, K, and both random matrices -----
            call rayleigh_limit(a, lam, nr, ki, eps_ba, qext, qsca, walb, asymm, &
                                al1, al2, al3, al4, be1, be2, lmax, &
                                qext_ori=qeo, qre_ori=qro)
            call add_aligned_size(a, lam, nr, ki, eps_ba, .true., 0, dn, f, area, &
                                  qext, qsca, qeo, qro, al1, al2, al3, al4, be1, be2, lmax, &
                                  ti, ts, ph, z_grid, z_al, cext_al, cpol_al, cbir_al, &
                                  sacc, sref, lmax_acc, csca_tot, csca_ref, cext_tot, &
                                  cext_ref, t_nd, crate)
            n_small = n_small + 1

         else if (x <= x_large) then
            ! --- T-matrix: one solve leaves /TMAT/ valid for the AMPL loop ---
            call system_clock(c0)
            call tmd_one_scatmat(a, lam, nr, ki, eps_ba, np, ddelt, ndgs, &
                                 qext, qsca, walb, asymm, &
                                 al1, al2, al3, al4, be1, be2, lmax, ierr, nmax_tm)
            call system_clock(c1)
            t_tm = t_tm + real(c1 - c0, wp) / real(crate, wp)

            if (ierr == 0) then
               call add_aligned_size(a, lam, nr, ki, eps_ba, .false., nmax_tm, dn, f, area, &
                                     qext, qsca, qeo, qro, al1, al2, al3, al4, be1, be2, lmax, &
                                     ti, ts, ph, z_grid, z_al, cext_al, cpol_al, cbir_al, &
                                     sacc, sref, lmax_acc, csca_tot, csca_ref, cext_tot, &
                                     cext_ref, t_nd, crate)
               n_tmat = n_tmat + 1
            else if (x < 1.0_wp) then
               ! Redirect to Rayleigh (run_scatmat rule): fully processed.
               call rayleigh_limit(a, lam, nr, ki, eps_ba, qext, qsca, walb, asymm, &
                                   al1, al2, al3, al4, be1, be2, lmax, &
                                   qext_ori=qeo, qre_ori=qro)
               call add_aligned_size(a, lam, nr, ki, eps_ba, .true., 0, dn, f, area, &
                                     qext, qsca, qeo, qro, al1, al2, al3, al4, be1, be2, lmax, &
                                     ti, ts, ph, z_grid, z_al, cext_al, cpol_al, cbir_al, &
                                     sacc, sref, lmax_acc, csca_tot, csca_ref, cext_tot, &
                                     cext_ref, t_nd, crate)
               n_fail = n_fail + 1
            else
               ! No oriented geometric-optics matrix: unaligned only.
               call geometric_optics_limit(a, lam, nr, ki, eps_ba, qext, qsca, walb, asymm, &
                                           al1, al2, al3, al4, be1, be2, lmax)
               csca = qsca * area;  cext = qext * area
               call add_random_coeffs(sacc, dn*csca, al1, al2, al3, al4, be1, be2, lmax, lmax_acc)
               csca_tot = csca_tot + dn*csca;  cext_tot = cext_tot + dn*cext
               skipped_weight = skipped_weight + dn*csca
               n_fail = n_fail + 1
            end if

         else
            ! --- x > 50: geometric optics, unaligned only -------------------
            call geometric_optics_limit(a, lam, nr, ki, eps_ba, qext, qsca, walb, asymm, &
                                        al1, al2, al3, al4, be1, be2, lmax)
            csca = qsca * area;  cext = qext * area
            call add_random_coeffs(sacc, dn*csca, al1, al2, al3, al4, be1, be2, lmax, lmax_acc)
            csca_tot = csca_tot + dn*csca;  cext_tot = cext_tot + dn*cext
            skipped_weight = skipped_weight + dn*csca
            n_skip = n_skip + 1
         end if
      end do

      ! Closure scattering cross section per theta_i from the finished Z_al.
      call scattering_cross_section_from_grid(z_al, ts, ph, csca_al_grid)

      deallocate(al1, al2, al3, al4, be1, be2, z_grid)
      if (present(t_tmat)) t_tmat = t_tm
      if (present(t_node)) t_node = t_nd
   end subroutine accumulate_aligned_population


   subroutine add_aligned_size(a, lam, nr, ki, eps_ba, is_rayleigh, nmax_tm, &
                               dn, f, area, qext, qsca, qeo, qro, &
                               al1, al2, al3, al4, be1, be2, lmax, ti, ts, ph, &
                               z_grid, z_al, cext_al, cpol_al, cbir_al, &
                               sacc, sref, lmax_acc, csca_tot, csca_ref, &
                               cext_tot, cext_ref, t_nd, crate)
      ! Add one fully-treatable size (Rayleigh or converged T-matrix) to every
      ! aligned product: the random-orientation matrices F_tot and F_ref, the
      ! forward-amplitude K elements, and the aligned phase matrix Z_al.  For a
      ! Rayleigh size qeo/qro carry the oriented C_ext/birefringence twins used
      ! by the sin^2 dipole law; for a T-matrix size K comes from AMPL forward
      ! amplitudes and qeo/qro are ignored.
      real(wp), intent(in)    :: a, lam, nr, ki, eps_ba, dn, f, area, qext, qsca
      real(wp), intent(in)    :: qeo(3), qro(3)
      logical,  intent(in)    :: is_rayleigh
      integer,  intent(in)    :: nmax_tm, lmax
      real(wp), intent(in)    :: al1(:), al2(:), al3(:), al4(:), be1(:), be2(:)
      real(wp), intent(in)    :: ti(:), ts(:), ph(:)
      real(wp), intent(inout) :: z_grid(:,:,:,:,:)
      real(wp), intent(inout) :: z_al(:,:,:,:,:)
      real(wp), intent(inout) :: cext_al(:), cpol_al(:), cbir_al(:)
      real(wp), intent(inout) :: sacc(:,:), sref(:,:)
      integer,  intent(inout) :: lmax_acc
      real(wp), intent(inout) :: csca_tot, csca_ref, cext_tot, cext_ref, t_nd
      integer(int64), intent(in) :: crate

      integer  :: it, i, j, is, ip
      real(wp) :: csca, cext, cv, ch, rv, rh, s2, c2
      integer(int64) :: c0, c1
      complex(wp) :: vv, vh, hv, hh
      external :: ampl

      csca = qsca * area;  cext = qext * area

      ! Random-orientation matrices: total (dn C_sca) and aligned twin (dn f C_sca).
      call add_random_coeffs(sacc, dn*csca,   al1, al2, al3, al4, be1, be2, lmax, lmax_acc)
      call add_random_coeffs(sref, dn*f*csca, al1, al2, al3, al4, be1, be2, lmax, lmax_acc)
      csca_tot = csca_tot + dn*csca;    cext_tot = cext_tot + dn*cext
      csca_ref = csca_ref + dn*f*csca;  cext_ref = cext_ref + dn*f*cext

      ! K elements from the forward amplitudes, on the theta_i grid.
      do it = 1, size(ti)
         if (is_rayleigh) then
            ! Exact dipole law: axial (E||a) = jori 2 = qeo(2); transverse
            ! (E perp a) = jori 3 = qeo(3).  C_v mixes them by sin^2/cos^2 of
            ! the incidence angle; C_h is purely transverse.
            s2 = sin(ti(it)*DEG)**2;  c2 = cos(ti(it)*DEG)**2
            cv = (s2*qeo(2) + c2*qeo(3)) * area
            ch =  qeo(3) * area
            rv = (s2*qro(2) + c2*qro(3)) * area
            rh =  qro(3) * area
         else
            ! Optical theorem on the forward amplitude, C = (4 pi / k) S =
            ! 2 lambda S; birefringence twin from the real part.
            call ampl(nmax_tm, lam, ti(it), ti(it), 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
                      vv, vh, hv, hh)
            cv = 2.0_wp * lam * aimag(vv);  ch = 2.0_wp * lam * aimag(hh)
            rv = 2.0_wp * lam * real(vv, wp);  rh = 2.0_wp * lam * real(hh, wp)
         end if
         cext_al(it) = cext_al(it) + dn*f * 0.5_wp*(cv + ch)
         cpol_al(it) = cpol_al(it) + dn*f * 0.5_wp*(ch - cv)
         cbir_al(it) = cbir_al(it) + dn*f * 0.5_wp*(rh - rv)
      end do

      ! Aligned phase matrix on the (theta_i, theta_s, phi) grid.
      call system_clock(c0)
      call oriented_mueller_grid(is_rayleigh, nmax_tm, a, lam, nr, ki, eps_ba, &
                                 ti, ts, ph, z_grid)
      call system_clock(c1)
      t_nd = t_nd + real(c1 - c0, wp) / real(crate, wp)
      do ip = 1, size(ph)
         do is = 1, size(ts)
            do it = 1, size(ti)
               do i = 1, 4
                  do j = 1, 4
                     z_al(i,j,it,is,ip) = z_al(i,j,it,is,ip) + dn*f * z_grid(i,j,it,is,ip)
                  end do
               end do
            end do
         end do
      end do
   end subroutine add_aligned_size


   subroutine oriented_mueller_grid(is_rayleigh, nmax_tm, a_eff, lam, nr, ki, eps_ba, &
                                    ti, ts, ph, z_grid)
      ! Fixed-orientation Mueller matrix at every (theta_i, theta_s, phi) node
      ! for ONE size, in um^2 sr^-1.  T-matrix regime reads the converged
      ! T-matrix from COMMON /TMAT/ (valid from a prior TMD_ONE_SCATMAT for this
      ! size) through AMPL; Rayleigh regime uses the analytic dipole matrix.
      !
      ! The node loop only reads shared state and writes each grid slot once, so
      ! it is parallelized with OpenMP and is identical for any thread count.
      logical,  intent(in)  :: is_rayleigh
      integer,  intent(in)  :: nmax_tm
      real(wp), intent(in)  :: a_eff, lam, nr, ki, eps_ba
      real(wp), intent(in)  :: ti(:), ts(:), ph(:)
      real(wp), intent(out) :: z_grid(:,:,:,:,:)
      integer  :: it, is, ip
      real(wp) :: ztmp(4,4)

      !$omp parallel do collapse(3) default(shared) private(it, is, ip, ztmp) &
      !$omp&            schedule(static)
      do ip = 1, size(ph)
         do is = 1, size(ts)
            do it = 1, size(ti)
               if (is_rayleigh) then
                  call rayleigh_mueller_matrix_oriented(a_eff, lam, nr, ki, eps_ba, &
                                     ti(it), ts(is), ph(ip), ztmp)
               else
                  call mueller_matrix_fixed_orientation(nmax_tm, lam, &
                                     ti(it), ts(is), ph(ip), ztmp)
               end if
               z_grid(:,:,it,is,ip) = ztmp
            end do
         end do
      end do
      !$omp end parallel do
   end subroutine oriented_mueller_grid


   subroutine scattering_cross_section_from_grid(z_al, ts, ph, csca)
      ! C_sca_al(theta_i) = INT Z_al11 dOmega over the full scattering sphere,
      ! by trapezoid quadrature over the STORED grid.  The azimuth is stored on
      ! [0,180] only; the mirror symmetry Z11(360 - phi) = Z11(phi) doubles the
      ! [0,180] azimuth integral to cover [180,360].  Trapezoid in cos(theta_s)
      ! over ts in [0,180] absorbs the sin(theta_s) area factor.
      real(wp), intent(in)  :: z_al(:,:,:,:,:), ts(:), ph(:)
      real(wp), intent(out) :: csca(:)
      integer  :: it, is, nti, nts
      real(wp) :: pint_lo, pint_hi, u_lo, u_hi, acc

      nti = size(z_al, 3);  nts = size(ts)
      do it = 1, nti
         acc = 0.0_wp
         ! Trapezoid in u = cos(theta_s); u decreases from +1 (ts=0) to -1.
         call azimuth_integral(z_al(1,1,it,1,:), ph, pint_lo)
         u_lo = cos(ts(1)*DEG)
         do is = 1, nts - 1
            call azimuth_integral(z_al(1,1,it,is+1,:), ph, pint_hi)
            u_hi = cos(ts(is+1)*DEG)
            acc  = acc + 0.5_wp*(pint_lo + pint_hi) * (u_lo - u_hi)
            pint_lo = pint_hi;  u_lo = u_hi
         end do
         csca(it) = 2.0_wp * acc      ! phi mirror covers [180,360]
      end do
   end subroutine scattering_cross_section_from_grid


   subroutine azimuth_integral(z11, ph, val)
      ! Trapezoid of z11(phi) over the stored azimuth grid ph [deg], returning
      ! INT_0^pi z11 dphi in radians.
      real(wp), intent(in)  :: z11(:), ph(:)
      real(wp), intent(out) :: val
      integer :: ip
      val = 0.0_wp
      do ip = 1, size(ph) - 1
         val = val + 0.5_wp*(z11(ip) + z11(ip+1)) * (ph(ip+1) - ph(ip)) * DEG
      end do
   end subroutine azimuth_integral


   subroutine add_random_coeffs(s, w, al1, al2, al3, al4, be1, be2, lmax, lmax_acc)
      ! Accumulate one size's GSP coefficients into the six-column store s with
      ! weight w, as tmd.lp.f averages over its size quadrature (weight =
      ! bin population times C_sca).  Columns are alpha1..4, beta1, beta2.
      real(wp), intent(inout) :: s(:,:)
      real(wp), intent(in)    :: w
      real(wp), intent(in)    :: al1(:), al2(:), al3(:), al4(:), be1(:), be2(:)
      integer,  intent(in)    :: lmax
      integer,  intent(inout) :: lmax_acc
      integer :: i, l1m
      l1m = lmax + 1
      lmax_acc = max(lmax_acc, l1m)
      do i = 1, l1m
         s(i,1) = s(i,1) + al1(i)*w
         s(i,2) = s(i,2) + al2(i)*w
         s(i,3) = s(i,3) + al3(i)*w
         s(i,4) = s(i,4) + al4(i)*w
         s(i,5) = s(i,5) + be1(i)*w
         s(i,6) = s(i,6) + be2(i)*w
      end do
   end subroutine add_random_coeffs

end module aligned_population_optics
