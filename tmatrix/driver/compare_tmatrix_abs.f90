program compare_tmatrix_abs
   ! Stage-C check of the orientation-resolved T-matrix absorption.
   !
   ! Over a strided subset of the T-matrix regime (0.1 < x < 50, with
   ! x = 2 pi a_eff / lambda) of the DH21 astrodust grid, this computes the
   ! oriented cross sections Q_ext(jori), Q_sca(jori), Q_abs(jori) from
   ! tmatrix_oriented_cross (optical theorem for extinction, phase-matrix
   ! integral for scattering, absorption by difference) and compares the
   ! derived combinations, in the convention of sed/src/q_table_jori.f90:
   !
   !   qpol_abs = 0.5 * (Q_abs(jori=3) - Q_abs(jori=2))   polarized absorption
   !   qran_abs = (Q_abs(1) + Q_abs(2) + Q_abs(3)) / 3    mean absorption
   !
   ! against the Q_abs block of the Hensley & Draine (2021/2023)
   ! orientation-resolved table.
   !
   ! Three de-risking anchors are enforced/reported:
   !   Anchor 1  physical bounds: Q_sca(jori) >= 0, Q_abs(jori) >= 0,
   !             albedo in [0,1] at every node.  Violations are counted; the
   !             run STOPs if more than a rounding-level handful occur.
   !   Anchor 2  Rayleigh overlap: at x ~ 0.1 the T-matrix Q_sca(jori) must
   !             match the closed-form Rayleigh Q_sca(jori) (rayleigh_limit).
   !   Anchor 3  quadrature convergence: at a few nodes the scattering
   !             integral is recomputed with the angular counts doubled and
   !             the relative change reported.
   !
   ! The refractive index m(lambda), the spheroid shape (oblate b/a = 1.4),
   ! and the x-window match run_tmatrix.f90.  A T-matrix solve plus two
   ! angular sweeps is done per node, so the grid is subsampled with a fixed
   ! stride (printed in the report), never silently capped.

   use, intrinsic :: iso_fortran_env, only: real64, error_unit
   use constants,        only: wp
   use read_index,       only: load_index, interp_m
   use asymptotic_optics, only: rayleigh_limit
   use tmatrix_oriented, only: tmatrix_oriented_cross
   implicit none

   character(len=*), parameter :: f_aeff  = '../data/dielectric/DH21_aeff'
   character(len=*), parameter :: f_wave  = '../data/dielectric/DH21_wave'
   character(len=*), parameter :: f_index = &
      '../data/dielectric/index_DH21Ad_P0.20_0.00_1.400'
   character(len=*), parameter :: f_qjori = &
      '../data/dielectric/q_DH21Ad_P0.20_Fe0.00_1.400.dat.gz'
   character(len=*), parameter :: f_scratch = 'q_jori_abs_cmp_scratch.dat'

   integer,  parameter :: NA = 169, NW = 1129, NORI = 3, NHEAD = 12
   real(wp), parameter :: EPS_BA  = 1.4_wp
   real(wp), parameter :: DDELT   = 1.0e-3_wp
   integer,  parameter :: NDGS    = 2, NP_OBL = -1
   real(wp), parameter :: X_SMALL = 0.1_wp, X_LARGE = 50.0_wp
   real(wp), parameter :: PI      = acos(-1.0_wp)

   ! The two angular sweeps make each node heavier than the Stage-B (ext)
   ! check, so the grid is strided more aggressively and large x is capped.
   integer,  parameter :: JW_STEP = 16, JA_STEP = 8
   real(wp), parameter :: X_CAP   = 30.0_wp   ! skip x > X_CAP (stated, not silent)
   integer,  parameter :: MAX_VIOL = 3        ! Anchor-1 tolerance (rounding-level)

   real(wp) :: a_eff(NA), lambda(NW)
   real(wp) :: nr_cache(NW), ki_cache(NW)
   real(wp), allocatable :: qabs_hd(:,:,:)          ! (jw, ja, jori)
   real(wp), allocatable :: dev_pol(:), dev_ran(:)
   real(wp) :: qext(NORI), qsca(NORI), qabs(NORI)
   complex(wp) :: m
   real(wp) :: x, qpol_ours, qran_ours, qpol_hd, qran_hd, floor_pol
   real(wp) :: walb
   integer  :: jw, ja, jo, ierr, n_win, n_solved, n_fail, n_pol, n_ran, n_viol
   real(wp) :: t0, t1

   call cpu_time(t0)

   call read_one_col(f_aeff, NA, a_eff)
   call read_one_col(f_wave, NW, lambda)
   call load_index(f_index)
   do jw = 1, NW
      call interp_m(lambda(jw), nr_cache(jw), ki_cache(jw))
   end do

   allocate(qabs_hd(NW, NA, NORI))
   call read_qabs_jori(f_qjori, f_scratch, qabs_hd)

   ! ---- Anchor 2 first: Rayleigh overlap at x ~ 0.1 --------------------
   call anchor_rayleigh_overlap()

   ! ---- Anchor 3: quadrature convergence at small/mid/large x ----------
   call anchor_quadrature_convergence()

   ! ---- main comparison sweep ------------------------------------------
   allocate(dev_pol(NW*NA), dev_ran(NW*NA))
   n_win = 0;  n_solved = 0;  n_fail = 0;  n_pol = 0;  n_ran = 0;  n_viol = 0

   do jw = 1, NW, JW_STEP
      do ja = 1, NA, JA_STEP
         x = 2.0_wp * PI * a_eff(ja) / lambda(jw)
         if (x <= X_SMALL .or. x >= X_LARGE) cycle
         if (x > X_CAP) cycle
         n_win = n_win + 1

         m = cmplx(nr_cache(jw), ki_cache(jw), kind=wp)
         call tmatrix_oriented_cross(a_eff(ja), lambda(jw), m, EPS_BA, NP_OBL, &
                                     DDELT, NDGS, qext, qsca, qabs, ierr)
         if (ierr /= 0) then
            n_fail = n_fail + 1
            cycle
         end if
         n_solved = n_solved + 1

         ! Anchor 1: physical bounds per orientation.
         do jo = 1, NORI
            walb = 0.0_wp
            if (qext(jo) > 0.0_wp) walb = qsca(jo) / qext(jo)
            if (qsca(jo) < 0.0_wp .or. qabs(jo) < 0.0_wp .or. &
                walb < 0.0_wp .or. walb > 1.0_wp) then
               n_viol = n_viol + 1
               write(error_unit,'(a,i0,a,3(1x,es12.4))') &
                  ' ANCHOR-1 VIOLATION jori=', jo, &
                  '  Qext,Qsca,Qabs =', qext(jo), qsca(jo), qabs(jo)
               write(error_unit,'(a,es12.4,a,es12.4,a,es12.4)') &
                  '   at lambda=', lambda(jw), '  a_eff=', a_eff(ja), '  x=', x
            end if
         end do

         qpol_ours = 0.5_wp * (qabs(3) - qabs(2))
         qran_ours = (qabs(1) + qabs(2) + qabs(3)) / 3.0_wp

         qpol_hd = 0.5_wp * (qabs_hd(jw,ja,3) - qabs_hd(jw,ja,2))
         qran_hd = (qabs_hd(jw,ja,1) + qabs_hd(jw,ja,2) + qabs_hd(jw,ja,3)) / 3.0_wp

         if (abs(qran_hd) > 0.0_wp) then
            n_ran = n_ran + 1
            dev_ran(n_ran) = abs(qran_ours/qran_hd - 1.0_wp)
         end if

         ! Polarized absorption can pass through zero, so form the ratio only
         ! where the HD23 signal is meaningful relative to the local mean.
         floor_pol = 1.0e-4_wp * abs(qran_hd)
         if (abs(qpol_hd) > floor_pol) then
            n_pol = n_pol + 1
            dev_pol(n_pol) = abs(qpol_ours/qpol_hd - 1.0_wp)
         end if
      end do
   end do

   if (n_viol > MAX_VIOL) then
      write(error_unit,'(a,i0,a,i0,a)') ' STOP: Anchor-1 violations (', &
         n_viol, ') exceed the rounding-level tolerance (', MAX_VIOL, ').'
      stop 1
   end if

   call cpu_time(t1)

   write(*,'(a)') '======================================================================'
   write(*,'(a)') ' Stage-C T-matrix oriented Q_abs vs HD23 (0.1 < x < 50)'
   write(*,'(a)') '======================================================================'
   write(*,'(a,i0,a,i0)') ' lambda stride : every ', JW_STEP, ' of ', NW
   write(*,'(a,i0,a,i0)') ' a_eff  stride : every ', JA_STEP, ' of ', NA
   write(*,'(a,f0.1)')    ' x cap (skip x >): ', X_CAP
   write(*,'(a,i0)')      ' strided nodes in 0.1 < x < min(50,cap) : ', n_win
   write(*,'(a,i0)')      '   T-matrix converged           : ', n_solved
   write(*,'(a,i0)')      '   T-matrix IERR /= 0 (skipped)  : ', n_fail
   write(*,'(a,i0)')      ' Anchor-1 physical-bound violations : ', n_viol
   write(*,'(a)') '----------------------------------------------------------------------'
   write(*,'(a)') ' qran_abs = (Q1+Q2+Q3)/3   [mean absorption]'
   write(*,'(a,i0)')     '   nodes compared         : ', n_ran
   write(*,'(a,es12.4)') '   median |ours/HD23 - 1| : ', median(dev_ran, n_ran)
   write(*,'(a,es12.4)') '   max    |ours/HD23 - 1| : ', maxval_n(dev_ran, n_ran)
   write(*,'(a)') '----------------------------------------------------------------------'
   write(*,'(a)') ' qpol_abs = 0.5*(Q3-Q2)    [polarized absorption]'
   write(*,'(a,i0,a)')   '   nodes compared         : ', n_pol, &
        '   (|qpol_HD| > 1e-4 * qran_HD)'
   write(*,'(a,es12.4)') '   median |ours/HD23 - 1| : ', median(dev_pol, n_pol)
   write(*,'(a,es12.4)') '   max    |ours/HD23 - 1| : ', maxval_n(dev_pol, n_pol)
   write(*,'(a)') '----------------------------------------------------------------------'
   write(*,'(a,f8.1,a)') ' wall (cpu) time : ', t1 - t0, ' s'
   write(*,'(a)') '======================================================================'

contains

   subroutine anchor_rayleigh_overlap()
      ! Find grid nodes with x just below the Rayleigh/T-matrix boundary
      ! (0.08 <= x <= 0.10) and compare the T-matrix oriented Q_sca(jori)
      ! against the closed-form Rayleigh Q_sca(jori).  Both code paths must
      ! agree there; this pins the absolute normalization of the scattering
      ! integral, as Stage B's Anchor 2 did for extinction.
      integer  :: jwl, jr, nshow
      real(wp) :: xx, qe_r, qs_r, wa_r, as_r
      real(wp) :: qext_r(NORI), qabs_r(NORI), qsca_r(NORI)
      real(wp) :: qe_t(NORI), qs_t(NORI), qa_t(NORI)
      complex(wp) :: mm
      integer :: kerr

      write(*,'(a)') '======================================================================'
      write(*,'(a)') ' Anchor 2: T-matrix vs Rayleigh Q_sca(jori) at x ~ 0.08-0.10'
      write(*,'(a)') '======================================================================'
      write(*,'(a)') '   node            x     jori    Q_sca(Tmat)   Q_sca(Rayl)   rel.diff'
      nshow = 0
      do jwl = 1, NW, 1
         do jr = 1, NA, 1
            xx = 2.0_wp * PI * a_eff(jr) / lambda(jwl)
            if (xx < 0.08_wp .or. xx > 0.10_wp) cycle
            mm = cmplx(nr_cache(jwl), ki_cache(jwl), kind=wp)

            call tmatrix_oriented_cross(a_eff(jr), lambda(jwl), mm, EPS_BA, &
                                        NP_OBL, DDELT, NDGS, qe_t, qs_t, qa_t, kerr)
            if (kerr /= 0) cycle

            call rayleigh_limit(a_eff(jr), lambda(jwl), nr_cache(jwl), &
                                ki_cache(jwl), EPS_BA, qe_r, qs_r, wa_r, as_r, &
                                qext_ori=qext_r, qabs_ori=qabs_r, qsca_ori=qsca_r)

            call jori_block(jr, jwl, xx, qs_t, qsca_r)
            nshow = nshow + 1
            if (nshow >= 3) then
               write(*,'(a)') '======================================================================'
               return
            end if
         end do
      end do
      if (nshow == 0) write(*,'(a)') '   (no grid node fell in 0.08 <= x <= 0.10)'
      write(*,'(a)') '======================================================================'
   end subroutine anchor_rayleigh_overlap

   subroutine jori_block(jr, jwl, xx, qs_t, qsca_r)
      integer,  intent(in) :: jr, jwl
      real(wp), intent(in) :: xx, qs_t(NORI), qsca_r(NORI)
      integer  :: jo
      real(wp) :: rd
      do jo = 1, NORI
         rd = 0.0_wp
         if (abs(qsca_r(jo)) > 0.0_wp) rd = qs_t(jo)/qsca_r(jo) - 1.0_wp
         write(*,'(i5,i5,f10.4,i6,3x,es13.5,1x,es13.5,1x,es11.3)') &
            jwl, jr, xx, jo, qs_t(jo), qsca_r(jo), rd
      end do
   end subroutine jori_block

   subroutine anchor_quadrature_convergence()
      ! Recompute Q_sca(jori) with the angular quadrature doubled (qmult=2)
      ! at representative small/mid/large-x nodes and report the relative
      ! change against the default (qmult=1).  A negligible change proves the
      ! scattering integral is resolved by the default counts.
      integer :: k, jwl, jr, jo, kerr
      real(wp) :: xx
      real(wp) :: qe1(NORI), qs1(NORI), qa1(NORI)
      real(wp) :: qe2(NORI), qs2(NORI), qa2(NORI)
      real(wp) :: rd, rdmax
      complex(wp) :: mm
      real(wp) :: x_targets(3)

      x_targets = (/ 0.3_wp, 3.0_wp, 20.0_wp /)
      write(*,'(a)') '======================================================================'
      write(*,'(a)') ' Anchor 3: scattering-integral convergence (double N_theta, N_phi)'
      write(*,'(a)') '======================================================================'
      write(*,'(a)') '   target_x    node        x     max_jori |Qsca(2x)/Qsca(1x) - 1|'
      do k = 1, 3
         call nearest_node(x_targets(k), jwl, jr, xx)
         mm = cmplx(nr_cache(jwl), ki_cache(jwl), kind=wp)
         call tmatrix_oriented_cross(a_eff(jr), lambda(jwl), mm, EPS_BA, NP_OBL, &
                                     DDELT, NDGS, qe1, qs1, qa1, kerr, qmult=1)
         if (kerr /= 0) cycle
         call tmatrix_oriented_cross(a_eff(jr), lambda(jwl), mm, EPS_BA, NP_OBL, &
                                     DDELT, NDGS, qe2, qs2, qa2, kerr, qmult=2)
         if (kerr /= 0) cycle
         rdmax = 0.0_wp
         do jo = 1, NORI
            if (abs(qs1(jo)) > 0.0_wp) then
               rd = abs(qs2(jo)/qs1(jo) - 1.0_wp)
               rdmax = max(rdmax, rd)
            end if
         end do
         write(*,'(f10.3,i8,i5,f10.4,6x,es13.5)') x_targets(k), jwl, jr, xx, rdmax
      end do
      write(*,'(a)') '======================================================================'
   end subroutine anchor_quadrature_convergence

   subroutine nearest_node(xtarget, jwl_out, jr_out, x_out)
      ! Grid node whose x = 2 pi a_eff / lambda is closest to xtarget.
      real(wp), intent(in)  :: xtarget
      integer,  intent(out) :: jwl_out, jr_out
      real(wp), intent(out) :: x_out
      integer  :: jwl, jr
      real(wp) :: xx, best
      best = huge(1.0_wp)
      jwl_out = 1;  jr_out = 1;  x_out = 0.0_wp
      do jwl = 1, NW
         do jr = 1, NA
            xx = 2.0_wp * PI * a_eff(jr) / lambda(jwl)
            if (abs(xx - xtarget) < best) then
               best = abs(xx - xtarget)
               jwl_out = jwl;  jr_out = jr;  x_out = xx
            end if
         end do
      end do
   end subroutine nearest_node

   subroutine read_one_col(filename, n, arr)
      character(len=*), intent(in)  :: filename
      integer,          intent(in)  :: n
      real(wp),         intent(out) :: arr(n)
      integer :: u, ios
      character(len=512) :: header
      open(newunit=u, file=filename, status='old', action='read', iostat=ios)
      if (ios /= 0) then
         write(error_unit,'(a,a)') ' ERROR: cannot open ', trim(filename)
         stop 1
      end if
      read(u,'(a)') header
      read(u,'(a)') header
      read(u,*) arr(1:n)
      close(u)
   end subroutine read_one_col

   subroutine read_qabs_jori(gz_file, scratch, qabs_tab)
      ! Decompress the gzip'd HD23 table once and read only its Q_abs block.
      ! On-disk stream (see sed/src/q_table_jori.f90): 12 header lines, then
      ! for each quantity (Q_ext, Q_abs, Q_sca), for each orientation
      ! jori = 1..3, for each wavelength (1129 records), one record of the
      ! 169 radii.  Only the second quantity block (Q_abs) is kept.
      character(len=*), intent(in)  :: gz_file, scratch
      real(wp),         intent(out) :: qabs_tab(NW, NA, NORI)
      integer :: u, ios, estat, cstat, iq, jori, jwl, i
      real(wp) :: row(NA)
      character(len=512) :: line

      call execute_command_line('gzip -dc "'//trim(gz_file)//'" > "'// &
                                trim(scratch)//'"', exitstat=estat, cmdstat=cstat)
      if (cstat /= 0 .or. estat /= 0) then
         write(error_unit,'(a,a)') ' ERROR: gzip -dc failed on ', trim(gz_file)
         call delete_scratch(scratch);  stop 1
      end if

      open(newunit=u, file=trim(scratch), status='old', action='read', iostat=ios)
      if (ios /= 0) then
         write(error_unit,'(a,a)') ' ERROR: cannot open ', trim(scratch)
         call delete_scratch(scratch);  stop 1
      end if
      do i = 1, NHEAD
         read(u,'(a)',iostat=ios) line
         if (ios /= 0) then
            write(error_unit,'(a)') ' ERROR: unexpected EOF in header'
            close(u);  call delete_scratch(scratch);  stop 1
         end if
      end do

      do iq = 1, 3
         do jori = 1, NORI
            do jwl = 1, NW
               read(u,*,iostat=ios) row(1:NA)
               if (ios /= 0) then
                  write(error_unit,'(a)') ' ERROR: read error in Q block'
                  close(u);  call delete_scratch(scratch);  stop 1
               end if
               if (iq == 2) qabs_tab(jwl, :, jori) = row(1:NA)
            end do
         end do
      end do
      close(u)
      call delete_scratch(scratch)
   end subroutine read_qabs_jori

   subroutine delete_scratch(path)
      character(len=*), intent(in) :: path
      integer :: u, ios
      open(newunit=u, file=path, status='old', iostat=ios)
      if (ios == 0) close(u, status='delete')
   end subroutine delete_scratch

   real(wp) function median(v, n)
      real(wp), intent(in) :: v(:)
      integer,  intent(in) :: n
      real(wp), allocatable :: s(:)
      integer :: i, j
      real(wp) :: t
      if (n <= 0) then
         median = 0.0_wp;  return
      end if
      allocate(s(n));  s = v(1:n)
      do i = 2, n
         t = s(i);  j = i - 1
         do while (j >= 1)
            if (s(j) <= t) exit
            s(j+1) = s(j);  j = j - 1
         end do
         s(j+1) = t
      end do
      if (mod(n,2) == 1) then
         median = s((n+1)/2)
      else
         median = 0.5_wp * (s(n/2) + s(n/2 + 1))
      end if
      deallocate(s)
   end function median

   real(wp) function maxval_n(v, n)
      real(wp), intent(in) :: v(:)
      integer,  intent(in) :: n
      if (n <= 0) then
         maxval_n = 0.0_wp
      else
         maxval_n = maxval(v(1:n))
      end if
   end function maxval_n

end program compare_tmatrix_abs
