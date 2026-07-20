program run_scatmat_aligned
   ! Size-integrated fixed-orientation scattering (Mueller) matrix of the
   ! aligned astrodust population, for polarized radiative transfer.
   !
   ! For each requested wavelength this sweeps the HD23 astrodust size
   ! distribution and, through aligned_population_optics, accumulates the
   ! aligned phase matrix Z_al(theta_i; theta_s, phi), the extinction-matrix
   ! elements C_ext/C_pol/C_bir on the theta_i grid (from the forward
   ! amplitudes), the closure scattering cross section C_sca_al(theta_i), and
   ! the two random-orientation matrices F_tot (every grain) and F_ref (the
   ! aligned population, weighted by f_align).  See the module header for the
   ! partial-alignment decomposition and the eta contract the radiative-
   ! transfer host uses.
   !
   ! Usage:
   !   ./run_scatmat_aligned.x                 ! default UBVRI bands
   !   ./run_scatmat_aligned.x 0.55            ! one wavelength [um]
   !   ./run_scatmat_aligned.x 0.44 0.55 0.79  ! several wavelengths
   !   ./run_scatmat_aligned.x test            ! one band (0.55), reduced grid,
   !                                           ! run-time estimate + OpenMP check
   !   ...  profile=FILE                       ! optional last argument: read a
   !                                           ! two-column (a_eff[um], f_align)
   !                                           ! profile instead of falign_hd23
   !
   ! Output (text, ASCII):
   !   output/scatmat_aligned_astrodust_P0.20_Fe0.00_1.400.dat        (default)
   !   output/scatmat_aligned_astrodust_P0.20_Fe0.00_1.400.test.dat   (test)

   use, intrinsic :: iso_fortran_env, only: real64, error_unit
   use constants,     only: wp
   use read_index,    only: load_index, interp_m
   use size_dist_mod, only: load_size_dist, n_size, a_dist, dn_ad
   use q_table_jori_mod, only: falign_hd23, A_ALIGN, ALPHA_ALIGN, FMAX_ALIGN
   use aligned_population_optics, only: accumulate_aligned_population, &
                                        oriented_mueller_grid
   !$ use omp_lib
   implicit none

   character(len=*), parameter :: f_index = &
      '../data/dielectric/index_DH21Ad_P0.20_0.00_1.400'
   character(len=*), parameter :: f_sdist = '../data/release/size_distribution.dat'
   character(len=*), parameter :: f_stem  = &
      'output/scatmat_aligned_astrodust_P0.20_Fe0.00_1.400'

   ! NPL must match the PARAMETER of the same name in src/tmd.par.f
   ! (NPL = 2*NPN1 + 1 = 201); repeated here because the fixed-form include
   ! cannot be consumed by free-form Fortran.  A run-time guard in the
   ! accumulation catches a mismatch through the GSP truncation order.
   integer, parameter :: NPL = 201
   integer, parameter :: NF  = 181            ! scattering-angle rows in F block

   real(wp), parameter :: EPS_BA  = 1.4_wp
   real(wp), parameter :: DDELT   = 1.0e-3_wp
   integer,  parameter :: NDGS    = 2
   integer,  parameter :: NP_OBL  = -1
   real(wp), parameter :: X_SMALL = 0.1_wp
   real(wp), parameter :: X_LARGE = 50.0_wp
   real(wp), parameter :: SKIP_TOL = 1.0e-6_wp
   real(wp), parameter :: PI      = acos(-1.0_wp)

   ! Angular grids.
   real(wp), allocatable :: ti(:), ts(:), ph(:)
   integer  :: nti, nts, nph

   ! Requested wavelengths and alignment profile.
   real(wp), allocatable :: wl_req(:)
   integer  :: nwl
   character(len=256) :: f_out, prof_file
   logical  :: test_mode, use_profile
   real(wp), allocatable :: prof_a(:), prof_f(:)
   integer  :: n_prof

   real(wp), allocatable :: f_al(:)

   ! Accumulators returned per wavelength.
   real(wp), allocatable :: z_al(:,:,:,:,:)
   real(wp), allocatable :: cext_al(:), cpol_al(:), cbir_al(:), csca_al_grid(:)
   real(wp) :: sacc(NPL,6), sref(NPL,6)
   real(wp) :: theta(NF)
   real(wp) :: f11t(NF), f22t(NF), f33t(NF), f44t(NF), f12t(NF), f34t(NF)
   real(wp) :: f11r(NF), f22r(NF), f33r(NF), f44r(NF), f12r(NF), f34r(NF)

   real(wp) :: nr, ki
   real(wp) :: csca_tot, csca_ref, cext_tot, cext_ref, skipped_weight
   real(wp) :: t_tmat, t_node
   integer  :: lmax_acc, n_small, n_tmat, n_skip, n_fail
   integer  :: iw, u_out, ia, c

   external :: scatmat_from_moments

   call parse_cli(wl_req, nwl, f_out, test_mode, use_profile, prof_file)

   call load_index(f_index)
   call load_size_dist(f_sdist)
   if (use_profile) call read_profile(prof_file, prof_a, prof_f, n_prof)

   ! Alignment fraction per size bin: the HD23 power law, or a tabulated
   ! profile interpolated in log(a_eff) and clamped at the ends.
   allocate(f_al(n_size))
   do ia = 1, n_size
      if (use_profile) then
         f_al(ia) = falign_interp_loga(a_dist(ia), prof_a, prof_f, n_prof)
      else
         f_al(ia) = falign_hd23(a_dist(ia))
      end if
   end do

   call build_grids(test_mode, ti, ts, ph)
   nti = size(ti);  nts = size(ts);  nph = size(ph)
   allocate(z_al(4,4,nti,nts,nph))
   allocate(cext_al(nti), cpol_al(nti), cbir_al(nti), csca_al_grid(nti))

   write(*,'(a,i0,a,es11.4,a,es11.4,a)') ' size distribution: ', n_size, &
      ' bins, a_eff = ', a_dist(1), ' .. ', a_dist(n_size), ' um'
   write(*,'(a,i0,a,i0,a,i0,a,i0,a)') ' grid: theta_i x theta_s x phi = ', &
      nti, ' x ', nts, ' x ', nph, ' = ', nti*nts*nph, ' nodes/band'
   write(*,'(a,i0)') ' wavelengths requested: ', nwl
   write(*,'(a,a)')  ' output = ', trim(f_out)

   open(newunit=u_out, file=trim(f_out), status='replace', action='write')
   call write_global_header(u_out)

   do iw = 1, nwl
      call interp_m(wl_req(iw), nr, ki)
      call accumulate_aligned_population(wl_req(iw), nr, ki, EPS_BA, NP_OBL, DDELT, &
              NDGS, X_SMALL, X_LARGE, a_dist, dn_ad, f_al, ti, ts, ph, &
              z_al, cext_al, cpol_al, cbir_al, csca_al_grid, sacc, sref, lmax_acc, &
              csca_tot, csca_ref, cext_tot, cext_ref, skipped_weight, &
              n_small, n_tmat, n_skip, n_fail, t_tmat, t_node)

      if (csca_tot <= 0.0_wp) then
         write(error_unit,'(a,es12.5)') ' ERROR: zero total C_sca at lambda = ', wl_req(iw)
         stop 1
      end if
      if (skipped_weight / csca_tot > SKIP_TOL) then
         write(error_unit,'(a,es11.4,a,es11.4)') &
            ' ERROR: aligned products omit a non-negligible scattering fraction ', &
            skipped_weight/csca_tot, ' > ', SKIP_TOL
         stop 1
      end if

      ! Normalize the GSP accumulators exactly as run_scatmat (alpha1(0) = 1)
      ! and expand F_tot (dn C_sca weight) and F_ref (dn f C_sca weight).
      do c = 1, 6
         sacc(:,c) = sacc(:,c) / csca_tot
         if (csca_ref > 0.0_wp) sref(:,c) = sref(:,c) / csca_ref
      end do
      call scatmat_from_moments(sacc(:,1), sacc(:,2), sacc(:,3), sacc(:,4), &
                                sacc(:,5), sacc(:,6), lmax_acc-1, NF, theta, &
                                f11t, f22t, f33t, f44t, f12t, f34t)
      call scatmat_from_moments(sref(:,1), sref(:,2), sref(:,3), sref(:,4), &
                                sref(:,5), sref(:,6), lmax_acc-1, NF, theta, &
                                f11r, f22r, f33r, f44r, f12r, f34r)

      call write_band_block(u_out, wl_req(iw), nr, ki)

      write(*,'(a,i0,a,i0,a,es11.4,a,i0,a,i0,a,i0,a,i0,a,es9.2)') &
         ' [', iw, '/', nwl, ']  lambda=', wl_req(iw), ' um  Ray=', n_small, &
         ' Tmat=', n_tmat, ' GOskip=', n_skip, ' fail=', n_fail, &
         '  skip/Csca=', skipped_weight/csca_tot
   end do

   close(u_out)
   write(*,'(a,a)') ' wrote ', trim(f_out)

   if (test_mode) call test_diagnostics(t_tmat, t_node)

contains

   ! ==================================================================
   ! command line, profile, and angular grids
   ! ==================================================================
   subroutine parse_cli(wl, nwl, fout, testm, useprof, proffile)
      real(wp), allocatable, intent(out) :: wl(:)
      integer,               intent(out) :: nwl
      character(len=*),      intent(out) :: fout, proffile
      logical,               intent(out) :: testm, useprof
      real(wp), parameter :: BANDS(5) = (/ 0.36_wp, 0.44_wp, 0.55_wp, 0.64_wp, 0.79_wp /)
      integer  :: narg, nkeep, k, ios
      character(len=256) :: a1
      character(len=64), allocatable :: keep(:)

      narg = command_argument_count()
      testm = .false.;  useprof = .false.;  proffile = ''

      ! A trailing "profile=FILE" argument is peeled off first.
      allocate(keep(max(narg,1)))
      nkeep = 0
      do k = 1, narg
         call get_command_argument(k, a1)
         if (a1(1:8) == 'profile=') then
            proffile = trim(a1(9:))
            useprof  = .true.
         else
            nkeep = nkeep + 1
            keep(nkeep) = a1
         end if
      end do

      if (nkeep == 0) then
         nwl = size(BANDS)
         allocate(wl(nwl));  wl = BANDS
         fout = f_stem//'.dat'
      else if (trim(keep(1)) == 'test') then
         testm = .true.
         nwl = 1
         allocate(wl(1));  wl = (/ 0.55_wp /)
         fout = f_stem//'.test.dat'
      else
         nwl = nkeep
         allocate(wl(nwl))
         do k = 1, nwl
            read(keep(k),*,iostat=ios) wl(k)
            if (ios /= 0 .or. wl(k) <= 0.0_wp) then
               write(error_unit,'(a,a,a)') ' ERROR: "', trim(keep(k)), &
                  '" is not a positive wavelength'
               stop 1
            end if
         end do
         fout = f_stem//'.dat'
      end if
      deallocate(keep)
   end subroutine parse_cli


   subroutine read_profile(fname, pa, pf, np)
      ! Two-column alignment profile: a_eff [um], f_align.  Skips '#' and blank
      ! lines; requires strictly ascending a_eff for the log interpolation.
      character(len=*),      intent(in)  :: fname
      real(wp), allocatable, intent(out) :: pa(:), pf(:)
      integer,               intent(out) :: np
      integer :: u, ios, i
      character(len=512) :: line

      open(newunit=u, file=trim(fname), status='old', action='read', iostat=ios)
      if (ios /= 0) then
         write(error_unit,'(a,a)') ' ERROR: cannot open profile ', trim(fname)
         stop 1
      end if
      np = 0
      do
         read(u,'(a)',iostat=ios) line
         if (ios /= 0) exit
         line = adjustl(line)
         if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
         np = np + 1
      end do
      if (np < 2) then
         write(error_unit,'(a)') ' ERROR: profile needs >= 2 rows'
         stop 1
      end if
      allocate(pa(np), pf(np))
      rewind(u)
      i = 0
      do
         read(u,'(a)',iostat=ios) line
         if (ios /= 0) exit
         line = adjustl(line)
         if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
         i = i + 1
         read(line,*) pa(i), pf(i)
      end do
      close(u)
      do i = 2, np
         if (pa(i) <= pa(i-1)) then
            write(error_unit,'(a)') ' ERROR: profile a_eff not strictly ascending'
            stop 1
         end if
      end do
   end subroutine read_profile


   pure function falign_interp_loga(a, pa, pf, np) result(f)
      ! Linear interpolation of f_align in log10(a_eff), clamped to the end
      ! values outside the tabulated range.
      real(wp), intent(in) :: a, pa(:), pf(:)
      integer,  intent(in) :: np
      real(wp) :: f, t, x
      integer  :: lo, hi, mid
      if (a <= pa(1)) then
         f = pf(1);  return
      else if (a >= pa(np)) then
         f = pf(np);  return
      end if
      lo = 1;  hi = np
      do while (hi - lo > 1)
         mid = (lo + hi) / 2
         if (pa(mid) <= a) then;  lo = mid;  else;  hi = mid;  end if
      end do
      x = log10(a)
      t = (x - log10(pa(lo))) / (log10(pa(hi)) - log10(pa(lo)))
      f = (1.0_wp - t)*pf(lo) + t*pf(hi)
   end function falign_interp_loga


   subroutine build_grids(testm, ti, ts, ph)
      ! Production grid: theta_i 0(5)90, theta_s 0(1)180, phi 0(5)180.
      ! Test grid (reduced): theta_i 0(15)90, theta_s 0(5)180, phi 0(15)180.
      logical, intent(in) :: testm
      real(wp), allocatable, intent(out) :: ti(:), ts(:), ph(:)
      if (testm) then
         call linspace_step(0.0_wp, 90.0_wp,  15.0_wp, ti)
         call linspace_step(0.0_wp, 180.0_wp,  5.0_wp, ts)
         call linspace_step(0.0_wp, 180.0_wp, 15.0_wp, ph)
      else
         call linspace_step(0.0_wp, 90.0_wp,   5.0_wp, ti)
         call linspace_step(0.0_wp, 180.0_wp,  1.0_wp, ts)
         call linspace_step(0.0_wp, 180.0_wp,  5.0_wp, ph)
      end if
   end subroutine build_grids


   subroutine linspace_step(a, b, step, x)
      real(wp), intent(in)  :: a, b, step
      real(wp), allocatable, intent(out) :: x(:)
      integer :: n, i
      n = nint((b - a)/step) + 1
      allocate(x(n))
      do i = 1, n
         x(i) = a + real(i-1, wp)*step
      end do
   end subroutine linspace_step


   ! ==================================================================
   ! output
   ! ==================================================================
   subroutine write_global_header(u)
      integer, intent(in) :: u
      write(u,'(a)') '# SEDust -- size-integrated fixed-orientation (aligned) scattering'
      write(u,'(a)') '# matrix of DH21 astrodust spheroids, for polarized radiative transfer.'
      write(u,'(a)') '#   P = 0.20, f_Fe = 0.00, axis ratio b/a = 1.4 (oblate, NP = -1)'
      write(u,'(a,es9.2,a,i0)') '#   T-matrix: DDELT = ', DDELT, ', NDGS = ', NDGS
      write(u,'(a,a)') '#   dielectric index : ', trim(f_index)
      write(u,'(a,a)') '#   size distribution: ', trim(f_sdist)
      if (use_profile) then
         write(u,'(a,a)') '#   alignment profile: tabulated file ', trim(prof_file)
      else
         write(u,'(a,f5.3,a,f6.4,a,f4.2,a)') &
            '#   alignment profile: falign_hd23  f_max = ', FMAX_ALIGN, &
            ', a_align = ', A_ALIGN, ' um, alpha = ', ALPHA_ALIGN, ' (HD23)'
      end if
      write(u,'(a,i0,a,i0,a,i0)') '#   grid: theta_i(0..90) = ', nti, &
         ', theta_s(0..180) = ', nts, ', phi(0..180) = ', nph
      write(u,'(a)') '#'
      write(u,'(a)') '# PRODUCTS (all per H).  Z_al is the aligned phase matrix, 16 elements,'
      write(u,'(a)') '# um^2 sr^-1.  F_tot and F_ref are the six-element random-orientation'
      write(u,'(a)') '# matrices of, respectively, EVERY grain and the aligned population'
      write(u,'(a)') '# (weighted by f_align), each alpha1-normalized so (1/2) INT F11 dcos = 1;'
      write(u,'(a)') '# multiply by the Csca_tot / Csca_ref given per band to restore um^2/H.'
      write(u,'(a)') '# The K elements Cext_al, Cpol_al, Cbir_al [um^2/H] and Csca_al [um^2/H]'
      write(u,'(a)') '# are listed on the theta_i grid.'
      write(u,'(a)') '#'
      write(u,'(a)') '# ETA CONTRACT.  For a cell alignment scale eta the aligned matrix scales'
      write(u,'(a)') '# linearly, Z_al,cell = eta Z_al, and the unaligned remainder scattering'
      write(u,'(a)') '# matrix is F_unal = F_tot - eta F_ref (in absolute units, i.e. after'
      write(u,'(a)') '# restoring Csca_tot and Csca_ref).  The extinction matrix scales the same'
      write(u,'(a)') '# way: K_al(theta_i) -> eta K_al(theta_i), and the unaligned population'
      write(u,'(a)') '# adds the isotropic Cext = Cext_tot - Cext_ref (its Cpol = Cbir = 0).'
      write(u,'(a)') '# The linearity in f_align is exact, so these stored integrals give any eta.'
      write(u,'(a)') '#'
      write(u,'(a)') '# STOKES BASIS.  Mishchenko meridional (v,h) = (theta-hat, phi-hat) of each'
      write(u,'(a)') '# propagation direction in the grain frame (z = alignment axis), Q = Iv - Ih.'
      write(u,'(a)') '# theta_i is the incidence polar angle from the axis; (theta_s, phi) the'
      write(u,'(a)') '# scattering polar/azimuth angles.  At theta_i = 90 the v-polarization is'
      write(u,'(a)') '# jori = 2 (E||axis) and h is jori = 3 (E perp axis), so Cpol = 0.5(C3-C2)'
      write(u,'(a)') '# and Cbir = 0.5(Cre3-Cre2).'
      write(u,'(a)') '#'
      write(u,'(a)') '# SYMMETRIES (to reconstruct the unstored ranges; verified in Stage 1):'
      write(u,'(a)') '#  phi -> 360-phi : Z unchanged except the two off-diagonal 2x2 blocks'
      write(u,'(a)') '#     (elements 13,14,23,24,31,32,41,42) flip sign.  Covers phi in [180,360].'
      write(u,'(a)') '#  theta_i -> 180-theta_i : maps to theta_s -> 180-theta_s (phi unchanged)'
      write(u,'(a)') '#     with the same off-diagonal-block sign flip.  Covers theta_i in (90,180].'
      write(u,'(a)') '#  theta_i = 0 : Z(0;theta_s,phi) = Z(0;theta_s,0) R(phi), the phi = 0'
      write(u,'(a)') '#     matrix being six-element block-diagonal.'
      write(u,'(a)') '#'
      write(u,'(a)') '# Each band is a "# lambda =" block with a K block (theta_i rows:'
      write(u,'(a)') '#   theta_i  Cext_al  Cpol_al  Cbir_al  Csca_al[grid closure]),'
      write(u,'(a)') '# an F block (181 rows: Theta  F_tot(11,22,33,44,12,34)  F_ref(...)),'
      write(u,'(a)') '# and a Z block (theta_i theta_s phi  Z11 Z12 ... Z44, um^2 sr^-1 per H).'
      write(u,'(a)') '#'
   end subroutine write_global_header


   subroutine write_band_block(u, lam, nr_, ki_)
      integer,  intent(in) :: u
      real(wp), intent(in) :: lam, nr_, ki_
      integer :: it, is, ip, i, j
      write(u,'(a)') '#'
      write(u,'(a,es15.7)')          '# lambda [um]   = ', lam
      write(u,'(a,es15.7,a,es15.7)') '# m = ', nr_, '  + i ', ki_
      write(u,'(a,es15.7)')          '# Cext_tot/H [um^2] = ', cext_tot
      write(u,'(a,es15.7)')          '# Csca_tot/H [um^2] = ', csca_tot
      write(u,'(a,es15.7)')          '# Cext_ref/H [um^2] = ', cext_ref
      write(u,'(a,es15.7)')          '# Csca_ref/H [um^2] = ', csca_ref
      write(u,'(a,i0,a,i0,a,i0,a,i0)') '# size bins: Rayleigh ', n_small, &
         ', T-matrix ', n_tmat, ', GO-skipped ', n_skip, ', redirected ', n_fail
      write(u,'(a,es11.4)')          '# aligned-omitted Csca fraction = ', &
                                       skipped_weight/csca_tot
      write(u,'(a,i0)')              '# Lmax = ', lmax_acc-1
      ! K block.
      write(u,'(a)') '# K block: theta_i[deg]  Cext_al  Cpol_al  Cbir_al  Csca_al(grid closure)  [um^2/H]'
      do it = 1, nti
         write(u,'(f8.2,4es16.7e2)') ti(it), cext_al(it), cpol_al(it), cbir_al(it), &
                                     csca_al_grid(it)
      end do
      ! F block.
      write(u,'(a)') '# F block: Theta[deg]  F_tot(11 22 33 44 12 34)  F_ref(11 22 33 44 12 34)'
      do is = 1, NF
         write(u,'(f8.2,12es14.6e2)') theta(is), &
            f11t(is), f22t(is), f33t(is), f44t(is), f12t(is), f34t(is), &
            f11r(is), f22r(is), f33r(is), f44r(is), f12r(is), f34r(is)
      end do
      ! Z block.
      write(u,'(a)') '# Z block: theta_i theta_s phi  Z11 Z12 Z13 Z14 Z21 Z22 Z23 Z24 Z31 Z32 Z33 Z34 Z41 Z42 Z43 Z44  [um^2 sr^-1 /H]'
      do it = 1, nti
         do is = 1, nts
            do ip = 1, nph
               write(u,'(3f8.2,16es12.4e2)') ti(it), ts(is), ph(ip), &
                  ((z_al(i,j,it,is,ip), j=1,4), i=1,4)
            end do
         end do
      end do
   end subroutine write_band_block


   ! ==================================================================
   ! test-mode diagnostics: OpenMP identity and full-run time estimate
   ! ==================================================================
   subroutine test_diagnostics(t_tm, t_nd)
      real(wp), intent(in) :: t_tm, t_nd
      integer  :: ia_ref, nmax_tm, ierr, lmax, nthreads
      real(wp) :: a, lam, nr_, ki_, x, dnbest
      real(wp) :: qext, qsca, walb, asymm
      real(wp) :: a1(NPL), a2(NPL), a3(NPL), a4(NPL), b1(NPL), b2(NPL)
      real(wp), allocatable :: zg1(:,:,:,:,:), zgn(:,:,:,:,:)
      real(wp) :: dmax
      integer  :: nti_p, nts_p, nph_p
      real(wp) :: nodes_prod, nodes_test, t_full
      external :: tmd_one_scatmat

      lam = wl_req(1)
      nthreads = 1
      !$ nthreads = omp_get_max_threads()

      write(*,'(a)') '----------------------------------------------------------------------'
      write(*,'(a)') ' test diagnostics'

      ! Representative T-matrix size: the most populated bin with 0.1 <= x <= 50.
      ia_ref = 0;  dnbest = 0.0_wp
      do ia = 1, n_size
         if (dn_ad(ia) <= 0.0_wp) cycle
         x = 2.0_wp*PI*a_dist(ia)/lam
         if (x >= X_SMALL .and. x <= X_LARGE .and. dn_ad(ia) > dnbest) then
            dnbest = dn_ad(ia);  ia_ref = ia
         end if
      end do

      if (ia_ref > 0) then
         a = a_dist(ia_ref)
         call interp_m(lam, nr_, ki_)
         call tmd_one_scatmat(a, lam, nr_, ki_, EPS_BA, NP_OBL, DDELT, NDGS, &
                              qext, qsca, walb, asymm, a1, a2, a3, a4, b1, b2, &
                              lmax, ierr, nmax_tm)
         allocate(zg1(4,4,nti,nts,nph), zgn(4,4,nti,nts,nph))
         !$ call omp_set_num_threads(1)
         call oriented_mueller_grid(.false., nmax_tm, a, lam, nr_, ki_, EPS_BA, &
                                    ti, ts, ph, zg1)
         !$ call omp_set_num_threads(nthreads)
         call oriented_mueller_grid(.false., nmax_tm, a, lam, nr_, ki_, EPS_BA, &
                                    ti, ts, ph, zgn)
         dmax = maxval(abs(zg1 - zgn))
         write(*,'(a,f6.4,a,f6.3,a,i0)') '   OpenMP identity size a=', a, &
            ' um  x=', 2.0_wp*PI*a/lam, '  threads=', nthreads
         if (dmax == 0.0_wp) then
            write(*,'(a,es12.4,a)') '   max |Z(1 thread) - Z(N threads)| = ', dmax, &
               '   (bitwise identical)'
         else
            write(*,'(a,es12.4,a)') '   max |Z(1 thread) - Z(N threads)| = ', dmax, &
               '   (DIFFERS -- drop OpenMP)'
         end if
         deallocate(zg1, zgn)
      else
         write(*,'(a)') '   no T-matrix size found for the OpenMP check'
      end if

      ! Full-run time estimate: the T-matrix time is per size (grid-independent),
      ! the node time scales with the node count.  Production is 5 bands on the
      ! full grid.
      nti_p = nint(90.0_wp/5.0_wp)  + 1
      nts_p = nint(180.0_wp/1.0_wp) + 1
      nph_p = nint(180.0_wp/5.0_wp) + 1
      nodes_prod = real(nti_p, wp)*real(nts_p, wp)*real(nph_p, wp)
      nodes_test = real(nti, wp)*real(nts, wp)*real(nph, wp)
      t_full = 5.0_wp * (t_tm + t_nd * nodes_prod/nodes_test)
      write(*,'(a,f8.2,a,f8.2,a)') '   measured: T-matrix ', t_tm, ' s + node loop ', &
         t_nd, ' s (this test band)'
      write(*,'(a,f8.1,a,i0,a,i0,a)') '   estimated full 5-band run: ', t_full, &
         ' s with ', nthreads, ' threads (', int(nodes_prod), ' nodes/band)'
   end subroutine test_diagnostics

end program run_scatmat_aligned
