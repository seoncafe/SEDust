program run_scatmat
   ! Size-integrated scattering (Mueller) matrix of randomly oriented DH21
   ! astrodust spheroids, for use by polarized radiative transfer.
   !
   ! For each requested wavelength this loops over the HD23 astrodust size
   ! distribution, obtains the six generalized-spherical-function expansion
   ! coefficients of the single-size random-orientation scattering matrix
   ! from TMD_ONE_SCATMAT, accumulates them with C_sca * dn/nH weighting exactly
   ! as Mishchenko's tmd.lp.f does over its own size quadrature, normalizes,
   ! and expands the result onto a scattering-angle grid.
   !
   ! Randomly oriented particles with a plane of symmetry have only six
   ! independent scattering-matrix elements, not sixteen:
   !
   !     [ F11  F12   0    0  ]
   !     [ F12  F22   0    0  ]
   !     [  0    0   F33  F34 ]
   !     [  0    0  -F34  F44 ]
   !
   ! Usage:
   !   ./run_scatmat.x 0.55                 ! one wavelength [um]
   !   ./run_scatmat.x 0.44 0.55 0.75       ! several wavelengths
   !   ./run_scatmat.x list <file>          ! one wavelength [um] per line
   !   ./run_scatmat.x all                  ! every lambda on the DH21 grid
   !                                        ! (1129 points; hours -- opt-in)
   !   ./run_scatmat.x test                 ! 3-wavelength smoke test
   !
   ! Output (text, ASCII):
   !   output/scatmat_astrodust_P0.20_Fe0.00_1.400.dat        (default)
   !   output/scatmat_astrodust_P0.20_Fe0.00_1.400.test.dat   (test)
   !   output/scatmat_astrodust_P0.20_Fe0.00_1.400.all.dat    (all)
   !
   ! This program does not touch the cross-section table written by
   ! run_tmatrix.x; the two executables are independent.

   use, intrinsic :: iso_fortran_env, only: real64, error_unit
   use constants,    only: wp
   use read_index,   only: load_index, interp_m
   use asymptotic_optics, only: rayleigh_limit, geometric_optics_limit
   use size_dist_mod, only: load_size_dist, n_size, a_dist, dn_ad
   implicit none

   character(len=*), parameter :: f_wave  = '../data/dielectric/DH21_wave'
   character(len=*), parameter :: f_index = &
      '../data/dielectric/index_DH21Ad_P0.20_0.00_1.400'
   character(len=*), parameter :: f_sdist = '../data/release/size_distribution.dat'
   character(len=*), parameter :: f_stem  = &
      'output/scatmat_astrodust_P0.20_Fe0.00_1.400'

   ! NPL must match the PARAMETER of the same name in src/tmd.par.f
   ! (NPL = 2*NPN1 + 1 = 201).  The fixed-form include file cannot be
   ! consumed by free-form Fortran, so the value is repeated here; a
   ! run-time guard below catches a mismatch via LMAX overflow.
   integer, parameter :: NPL = 201

   integer,  parameter :: NW_GRID = 1129        ! DH21 wavelength grid size
   integer,  parameter :: NTHETA  = 181         ! 0..180 deg, 1 deg steps
   real(wp), parameter :: EPS_BA  = 1.4_wp
   real(wp), parameter :: DDELT   = 1.0e-3_wp
   integer,  parameter :: NDGS    = 2
   integer,  parameter :: NP_OBL  = -1
   real(wp), parameter :: X_SMALL = 0.1_wp
   real(wp), parameter :: X_LARGE = 50.0_wp
   real(wp), parameter :: PI      = acos(-1.0_wp)

   ! Single-size coefficients returned by TMD_ONE_SCATMAT or by the
   ! asymptotic limits.
   real(wp) :: al1(NPL), al2(NPL), al3(NPL), al4(NPL), be1(NPL), be2(NPL)
   ! C_sca-weighted accumulators over the size distribution.
   real(wp) :: sa1(NPL), sa2(NPL), sa3(NPL), sa4(NPL), sb1(NPL), sb2(NPL)
   real(wp) :: theta(NTHETA)
   real(wp) :: f11(NTHETA), f22(NTHETA), f33(NTHETA), f44(NTHETA)
   real(wp) :: f12(NTHETA), f34(NTHETA)

   real(wp), allocatable :: wl_req(:)
   real(wp) :: nr, ki, x, qext, qsca, walb, asymm, csca, cext
   real(wp) :: csca_tot, cext_tot, wgi, g_ref, g_coeff, alb_tot
   real(wp) :: f11_int
   integer  :: lmax, lmax_acc, l1m, ierr_t, kontr, lviol, nmax_tm
   integer  :: i, ia, iw, nwl, u_out, n_small, n_large, n_fail
   character(len=256) :: f_out
   character(len=64)  :: arg

   call wavelength_list_from_cli(wl_req, nwl, f_out)

   call load_index(f_index)
   call load_size_dist(f_sdist)
   write(*,'(a,i0,a,es11.4,a,es11.4,a)') ' size distribution: ', n_size, &
      ' bins, a_eff = ', a_dist(1), ' .. ', a_dist(n_size), ' um'
   write(*,'(a,i0)') ' wavelengths requested: ', nwl
   write(*,'(a,a)')  ' output = ', trim(f_out)

   open(newunit=u_out, file=trim(f_out), status='replace', action='write')
   call write_scatmat_header(u_out)

   do iw = 1, nwl
      call interp_m(wl_req(iw), nr, ki)

      sa1 = 0.0_wp;  sa2 = 0.0_wp;  sa3 = 0.0_wp
      sa4 = 0.0_wp;  sb1 = 0.0_wp;  sb2 = 0.0_wp
      csca_tot = 0.0_wp
      cext_tot = 0.0_wp
      g_ref    = 0.0_wp
      lmax_acc = 0
      n_small  = 0;  n_large = 0;  n_fail = 0

      do ia = 1, n_size
         if (dn_ad(ia) <= 0.0_wp) cycle       ! empty bin: nothing to add
         x = 2.0_wp * PI * a_dist(ia) / wl_req(iw)

         if (x < X_SMALL) then
            call rayleigh_limit(a_dist(ia), wl_req(iw), nr, ki, EPS_BA, &
                                  qext, qsca, walb, asymm, &
                                  al1, al2, al3, al4, be1, be2, lmax)
            n_small = n_small + 1
         else if (x > X_LARGE) then
            call geometric_optics_limit(a_dist(ia), wl_req(iw), nr, ki, &
                                  qext, qsca, walb, asymm, &
                                  al1, al2, al3, al4, be1, be2, lmax)
            n_large = n_large + 1
         else
            call tmd_one_scatmat(a_dist(ia), wl_req(iw), nr, ki, EPS_BA, NP_OBL, &
                            DDELT, NDGS, qext, qsca, walb, asymm, &
                            al1, al2, al3, al4, be1, be2, lmax, ierr_t, nmax_tm)
            if (ierr_t /= 0) then
               ! Same redirection rule as run_tmatrix.x: take the result
               ! from whichever asymptotic limit x is closer to.
               n_fail = n_fail + 1
               if (x < 1.0_wp) then
                  call rayleigh_limit(a_dist(ia), wl_req(iw), nr, ki, EPS_BA, &
                                        qext, qsca, walb, asymm, &
                                        al1, al2, al3, al4, be1, be2, lmax)
               else
                  call geometric_optics_limit(a_dist(ia), wl_req(iw), nr, ki, &
                                        qext, qsca, walb, asymm, &
                                        al1, al2, al3, al4, be1, be2, lmax)
               end if
            end if
         end if

         ! Cross sections in um^2 (Q convention of TMD_ONE: Q = C/(pi a^2)).
         csca = qsca * PI * a_dist(ia)**2
         cext = qext * PI * a_dist(ia)**2

         ! Accumulation follows tmd.lp.f lines 660-690: the expansion
         ! coefficients are averaged with weight (bin population) x C_sca,
         ! because each grain contributes to the emergent scattered
         ! radiation in proportion to how much it scatters.  dn_ad is
         ! already integrated over its size bin, so it plays the role of
         ! Mishchenko's quadrature weight WG1(I) directly.
         wgi = dn_ad(ia) * csca
         l1m = lmax + 1
         if (l1m > NPL) then
            ! Guards the NPL value repeated above against a change to
            ! src/tmd.par.f: TMD_ONE_SCATMAT can never return LMAX+1
            ! beyond its own NPL, so this can only fire on a mismatch.
            write(error_unit,'(a,i0,a,i0)') &
               ' ERROR: LMAX+1 = ', l1m, ' exceeds NPL = ', NPL
            stop 1
         end if
         lmax_acc = max(lmax_acc, l1m)
         do i = 1, l1m
            sa1(i) = sa1(i) + al1(i)*wgi
            sa2(i) = sa2(i) + al2(i)*wgi
            sa3(i) = sa3(i) + al3(i)*wgi
            sa4(i) = sa4(i) + al4(i)*wgi
            sb1(i) = sb1(i) + be1(i)*wgi
            sb2(i) = sb2(i) + be2(i)*wgi
         end do
         csca_tot = csca_tot + wgi
         cext_tot = cext_tot + cext * dn_ad(ia)
         ! Independent reference for validation item 1: the C_sca-weighted
         ! mean of the single-size <cos theta> returned by TMD_ONE_SCATMAT.
         g_ref = g_ref + asymm * wgi
      end do

      if (csca_tot <= 0.0_wp) then
         write(error_unit,'(a,es12.5)') &
            ' ERROR: zero total C_sca at lambda = ', wl_req(iw)
         stop 1
      end if

      sa1(1:lmax_acc) = sa1(1:lmax_acc) / csca_tot
      sa2(1:lmax_acc) = sa2(1:lmax_acc) / csca_tot
      sa3(1:lmax_acc) = sa3(1:lmax_acc) / csca_tot
      sa4(1:lmax_acc) = sa4(1:lmax_acc) / csca_tot
      sb1(1:lmax_acc) = sb1(1:lmax_acc) / csca_tot
      sb2(1:lmax_acc) = sb2(1:lmax_acc) / csca_tot
      g_ref   = g_ref / csca_tot
      g_coeff = sa1(2) / 3.0_wp
      ! csca_tot is sum(dn * C_sca); cext_tot is sum(dn * C_ext).
      alb_tot = csca_tot / cext_tot

      call vdm_hovenier_test(lmax_acc, sa1, sa2, sa3, sa4, sb1, sb2, kontr, lviol)

      call scatmat_from_moments(sa1, sa2, sa3, sa4, sb1, sb2, lmax_acc-1, NTHETA, &
                    theta, f11, f22, f33, f44, f12, f34)

      ! Normalization check: (1/2) int_{-1}^{1} F11 d(cos) should be 1.
      f11_int = phase_function_norm(theta, f11)

      call write_scatmat_block(u_out, wl_req(iw), nr, ki, cext_tot, csca_tot, &
                       alb_tot, g_coeff, g_ref, lmax_acc-1, kontr, lviol, &
                       sa1(1), f11_int, n_small, n_large, n_fail, &
                       theta, f11, f22, f33, f44, f12, f34)

      write(*,'(a,i0,a,i0,a,es11.4,a,f9.6,a,es9.2,a,i0,a,f10.7,a,a)') &
         ' [', iw, '/', nwl, ']  lambda=', wl_req(iw), &
         ' um  g=', g_coeff, '  dg=', abs(g_coeff-g_ref), &
         '  Lmax=', lmax_acc-1, '  <F11>=', f11_int, &
         '  hovenr=', trim(merge('OK  ', 'FAIL', kontr == 1))
      if (n_fail > 0) write(*,'(a,i0,a)') &
         '        note: ', n_fail, ' size bins fell back on an asymptotic limit after T-matrix failure'
   end do

   close(u_out)
   write(*,'(a,a)') ' wrote ', trim(f_out)

contains

   subroutine wavelength_list_from_cli(wl, nwl, fout)
      ! Parses the command line into an explicit list of wavelengths.
      real(wp), allocatable, intent(out) :: wl(:)
      integer,               intent(out) :: nwl
      character(len=*),      intent(out) :: fout
      real(wp) :: grid(NW_GRID)
      integer  :: narg, k, ios, u, n
      character(len=256) :: fname
      character(len=64)  :: a1

      narg = command_argument_count()
      if (narg < 1) then
         write(*,'(a)') ' usage:'
         write(*,'(a)') '   run_scatmat.x LAMBDA [LAMBDA ...]   wavelengths in um'
         write(*,'(a)') '   run_scatmat.x list FILE             one wavelength per line'
         write(*,'(a)') '   run_scatmat.x all                   full DH21 grid (slow)'
         write(*,'(a)') '   run_scatmat.x test                  3-wavelength smoke test'
         stop 1
      end if

      call get_command_argument(1, a1)
      select case (trim(a1))
      case ('test')
         nwl = 3
         allocate(wl(nwl))
         wl = [0.44_wp, 0.55_wp, 2.20_wp]
         fout = f_stem//'.test.dat'
      case ('all')
         call read_wave_grid(f_wave, NW_GRID, grid)
         nwl = NW_GRID
         allocate(wl(nwl))
         wl = grid
         fout = f_stem//'.all.dat'
      case ('list')
         if (narg < 2) then
            write(error_unit,'(a)') ' usage: run_scatmat.x list FILE'
            stop 1
         end if
         call get_command_argument(2, fname)
         open(newunit=u, file=trim(fname), status='old', action='read', iostat=ios)
         if (ios /= 0) then
            write(error_unit,'(a,a)') ' ERROR: cannot open ', trim(fname)
            stop 1
         end if
         n = 0
         do
            read(u,*,iostat=ios) x
            if (ios /= 0) exit
            n = n + 1
         end do
         rewind(u)
         if (n < 1) then
            write(error_unit,'(a,a)') ' ERROR: no wavelengths in ', trim(fname)
            stop 1
         end if
         nwl = n
         allocate(wl(nwl))
         do k = 1, nwl
            read(u,*) wl(k)
         end do
         close(u)
         fout = f_stem//'.dat'
      case default
         nwl = narg
         allocate(wl(nwl))
         do k = 1, nwl
            call get_command_argument(k, arg)
            read(arg,*,iostat=ios) wl(k)
            if (ios /= 0) then
               write(error_unit,'(a,a,a)') ' ERROR: "', trim(arg), &
                  '" is not a wavelength; see run_scatmat.x with no arguments'
               stop 1
            end if
            if (wl(k) <= 0.0_wp) then
               write(error_unit,'(a)') ' ERROR: wavelength must be positive'
               stop 1
            end if
         end do
         fout = f_stem//'.dat'
      end select
   end subroutine wavelength_list_from_cli


   subroutine read_wave_grid(filename, n, xw)
      ! DH21_wave: 2 header lines, then all n values on one long record.
      character(len=*), intent(in)  :: filename
      integer,          intent(in)  :: n
      real(wp),         intent(out) :: xw(n)
      integer :: u, ios
      character(len=512) :: header
      open(newunit=u, file=filename, status='old', action='read', iostat=ios)
      if (ios /= 0) then
         write(error_unit,'(a,a)') ' ERROR: cannot open ', trim(filename)
         stop 1
      end if
      read(u,'(a)') header
      read(u,'(a)') header
      read(u,*) xw(1:n)
      close(u)
   end subroutine read_wave_grid


   function phase_function_norm(th, f) result(s)
      ! (1/2) * integral_{-1}^{1} f d(cos theta), by the trapezoid rule on
      ! the equally spaced theta grid.  For the normalization convention
      ! used here this equals 1 when f = F11.
      real(wp), intent(in) :: th(:), f(:)
      real(wp) :: s, u1, u2
      integer  :: k
      s = 0.0_wp
      do k = 1, size(th) - 1
         u1 = cos(th(k)   * PI / 180.0_wp)
         u2 = cos(th(k+1) * PI / 180.0_wp)
         s = s + 0.5_wp * (f(k) + f(k+1)) * (u1 - u2)
      end do
      s = 0.5_wp * s
   end function phase_function_norm


   subroutine write_scatmat_header(u)
      integer, intent(in) :: u
      write(u,'(a)') '# SEDust -- size-integrated scattering matrix,'
      write(u,'(a)') '# randomly oriented DH21 astrodust spheroids.'
      write(u,'(a)') '#   P = 0.20, f_Fe = 0.00, axis ratio b/a = 1.4 (oblate)'
      write(u,'(a)') '#   size distribution: data/release/size_distribution.dat,'
      write(u,'(a)') '#     astrodust column, C_sca-weighted average over sizes.'
      write(u,'(a)') '#'
      write(u,'(a)') '# Randomly oriented particles with a plane of symmetry have'
      write(u,'(a)') '# six independent scattering-matrix elements:'
      write(u,'(a)') '#     [ F11  F12   0    0  ]'
      write(u,'(a)') '#     [ F12  F22   0    0  ]'
      write(u,'(a)') '#     [  0    0   F33  F34 ]'
      write(u,'(a)') '#     [  0    0  -F34  F44 ]'
      write(u,'(a)') '# Theta is the scattering angle in degrees, 0 = forward.'
      write(u,'(a)') '#'
      write(u,'(a)') '# Normalization (Mishchenko/Hovenier, as produced by MATR):'
      write(u,'(a)') '#   (1/2) * integral_{-1}^{+1} F11 d(cos Theta) = 1,'
      write(u,'(a)') '#   so F11 = 1 everywhere for isotropic scattering.'
      write(u,'(a)') '#   F22, F33, F44, F12, F34 are NOT divided by F11.'
      write(u,'(a)') '#   Degree of linear polarization for unpolarized incident'
      write(u,'(a)') '#   light is -F12/F11.  Multiply F_ij by C_sca (given per'
      write(u,'(a)') '#   block, in um^2 per H) to recover absolute units.'
      write(u,'(a)') '#'
      write(u,'(a)') '# Blocks are separated by a "# lambda =" line; each block'
      write(u,'(a)') '# has 181 rows at Theta = 0, 1, ..., 180 degrees.'
      write(u,'(a)') '#'
   end subroutine write_scatmat_header


   subroutine write_scatmat_block(u, lam, nr_, ki_, ce, cs, alb, gg, gg_ref, lmx, &
                          kon, lv, a1_zero, f11int, ns, nl, nf, &
                          th, a11, a22, a33, a44, a12, a34)
      integer,  intent(in) :: u, lmx, kon, lv, ns, nl, nf
      real(wp), intent(in) :: lam, nr_, ki_, ce, cs, alb, gg, gg_ref, a1_zero, f11int
      real(wp), intent(in) :: th(:), a11(:), a22(:), a33(:), a44(:), a12(:), a34(:)
      integer :: k
      write(u,'(a)') '#'
      write(u,'(a,es15.7)')   '# lambda [um]   = ', lam
      write(u,'(a,es15.7,a,es15.7)') '# m = ', nr_, '  + i ', ki_
      write(u,'(a,es15.7)')   '# Cext/H [um^2] = ', ce
      write(u,'(a,es15.7)')   '# Csca/H [um^2] = ', cs
      write(u,'(a,es15.7)')   '# albedo        = ', alb
      write(u,'(a,es15.7)')   '# g = a1(2)/3   = ', gg
      write(u,'(a,es15.7)')   '# g (Csca-wtd single-size, cross-check) = ', gg_ref
      write(u,'(a,i0)')       '# Lmax          = ', lmx
      write(u,'(a,es22.14)')  '# alpha_1(0) (exactly 1 by construction) = ', a1_zero
      write(u,'(a,es15.7)')   '# (1/2)int F11 dcos on this 1-deg grid   = ', f11int
      write(u,'(a)')          '#   (differs from 1 only by trapezoid error on the'
      write(u,'(a)')          '#    forward peak; alpha_1(0) is the exact statement)'
      if (kon == 1) then
         write(u,'(a)')       '# van der Mee & Hovenier test: SATISFIED'
      else
         write(u,'(a,i0)')    '# van der Mee & Hovenier test: VIOLATED first at L = ', lv
      end if
      write(u,'(a,i0,a,i0,a,i0)') '# size bins: small-x limit ', ns, &
         ', large-x limit ', nl, ', T-matrix failures redirected ', nf
      write(u,'(a)') '#   Theta[deg]          F11            F22            F33            F44            F12            F34'
      do k = 1, size(th)
         write(u,'(f10.2,6es15.6)') th(k), a11(k), a22(k), a33(k), a44(k), a12(k), a34(k)
      end do
   end subroutine write_scatmat_block

end program run_scatmat
