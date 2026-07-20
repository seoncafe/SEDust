program compare_tmatrix_jori
   ! Stage-B check of the orientation-resolved T-matrix extinction.
   !
   ! Over a strided subset of the T-matrix regime (0.1 < x < 50, with
   ! x = 2 pi a_eff / lambda) of the DH21 astrodust grid, this compares the
   ! oriented Q_ext(jori) produced by tmatrix_oriented_ext (optical theorem
   ! on Mishchenko's fixed-orientation amplitude matrix) against the
   ! Hensley & Draine (2021/2023) orientation-resolved table.  Two
   ! combinations are compared, in the convention of sed/src/q_table_jori.f90:
   !
   !   qpol_ext = 0.5 * (Q_ext(jori=3) - Q_ext(jori=2))   polarized extinction
   !   qran_ext = (Q_ext(1) + Q_ext(2) + Q_ext(3)) / 3    mean extinction
   !
   ! The refractive index m(lambda), the spheroid shape (oblate b/a = 1.4),
   ! and the x-window cutoffs match run_tmatrix.f90.
   !
   ! A T-matrix solve is done per node, so the grid is subsampled with a
   ! fixed stride (printed in the report).  The stride is stated, never
   ! silently capped.

   use, intrinsic :: iso_fortran_env, only: real64, error_unit
   use constants,        only: wp
   use read_index,       only: load_index, interp_m
   use tmatrix_oriented, only: tmatrix_oriented_ext
   implicit none

   character(len=*), parameter :: f_aeff  = '../data/dielectric/DH21_aeff'
   character(len=*), parameter :: f_wave  = '../data/dielectric/DH21_wave'
   character(len=*), parameter :: f_index = &
      '../data/dielectric/index_DH21Ad_P0.20_0.00_1.400'
   character(len=*), parameter :: f_qjori = &
      '../data/dielectric/q_DH21Ad_P0.20_Fe0.00_1.400.dat.gz'
   character(len=*), parameter :: f_scratch = 'q_jori_ext_cmp_scratch.dat'

   integer,  parameter :: NA = 169, NW = 1129, NORI = 3, NHEAD = 12
   real(wp), parameter :: EPS_BA  = 1.4_wp
   real(wp), parameter :: DDELT   = 1.0e-3_wp
   integer,  parameter :: NDGS    = 2, NP_OBL = -1
   real(wp), parameter :: X_SMALL = 0.1_wp, X_LARGE = 50.0_wp
   real(wp), parameter :: PI      = acos(-1.0_wp)

   ! Grid strides: subsample so a T-matrix solve per node stays affordable.
   integer, parameter :: JW_STEP = 8, JA_STEP = 4

   real(wp) :: a_eff(NA), lambda(NW)
   real(wp) :: nr_cache(NW), ki_cache(NW)
   real(wp), allocatable :: qext_hd(:,:,:)          ! (jw, ja, jori)
   real(wp), allocatable :: dev_pol(:), dev_ran(:)
   real(wp) :: qt(NORI)
   complex(wp) :: m
   real(wp) :: x, qpol_ours, qran_ours, qpol_hd, qran_hd, floor_pol
   integer  :: jw, ja, ierr, n_win, n_solved, n_fail, n_pol, n_ran

   call read_one_col(f_aeff, NA, a_eff)
   call read_one_col(f_wave, NW, lambda)
   call load_index(f_index)
   do jw = 1, NW
      call interp_m(lambda(jw), nr_cache(jw), ki_cache(jw))
   end do

   allocate(qext_hd(NW, NA, NORI))
   call read_qext_jori(f_qjori, f_scratch, qext_hd)

   allocate(dev_pol(NW*NA), dev_ran(NW*NA))
   n_win = 0;  n_solved = 0;  n_fail = 0;  n_pol = 0;  n_ran = 0

   do jw = 1, NW, JW_STEP
      do ja = 1, NA, JA_STEP
         x = 2.0_wp * PI * a_eff(ja) / lambda(jw)
         if (x <= X_SMALL .or. x >= X_LARGE) cycle
         n_win = n_win + 1

         m = cmplx(nr_cache(jw), ki_cache(jw), kind=wp)
         call tmatrix_oriented_ext(a_eff(ja), lambda(jw), m, EPS_BA, NP_OBL, &
                                   DDELT, NDGS, qt, ierr)
         if (ierr /= 0) then
            n_fail = n_fail + 1
            cycle
         end if
         n_solved = n_solved + 1

         qpol_ours = 0.5_wp * (qt(3) - qt(2))
         qran_ours = (qt(1) + qt(2) + qt(3)) / 3.0_wp

         qpol_hd = 0.5_wp * (qext_hd(jw,ja,3) - qext_hd(jw,ja,2))
         qran_hd = (qext_hd(jw,ja,1) + qext_hd(jw,ja,2) + qext_hd(jw,ja,3)) / 3.0_wp

         if (abs(qran_hd) > 0.0_wp) then
            n_ran = n_ran + 1
            dev_ran(n_ran) = abs(qran_ours/qran_hd - 1.0_wp)
         end if

         ! Polarized extinction can pass through zero, so form the ratio only
         ! where the HD23 signal is meaningful relative to the local mean.
         floor_pol = 1.0e-4_wp * abs(qran_hd)
         if (abs(qpol_hd) > floor_pol) then
            n_pol = n_pol + 1
            dev_pol(n_pol) = abs(qpol_ours/qpol_hd - 1.0_wp)
         end if
      end do
   end do

   write(*,'(a)') '======================================================================'
   write(*,'(a)') ' Stage-B T-matrix oriented Q_ext vs HD23 (0.1 < x < 50)'
   write(*,'(a)') '======================================================================'
   write(*,'(a,i0,a,i0)') ' lambda stride : every ', JW_STEP, ' of ', NW
   write(*,'(a,i0,a,i0)') ' a_eff  stride : every ', JA_STEP, ' of ', NA
   write(*,'(a,i0)')      ' strided nodes in 0.1 < x < 50 : ', n_win
   write(*,'(a,i0)')      '   T-matrix converged           : ', n_solved
   write(*,'(a,i0)')      '   T-matrix IERR /= 0 (skipped)  : ', n_fail
   write(*,'(a)') '----------------------------------------------------------------------'
   write(*,'(a)') ' qran_ext = (Q1+Q2+Q3)/3   [mean extinction]'
   write(*,'(a,i0)')     '   nodes compared         : ', n_ran
   write(*,'(a,es12.4)') '   median |ours/HD23 - 1| : ', median(dev_ran, n_ran)
   write(*,'(a,es12.4)') '   max    |ours/HD23 - 1| : ', maxval_n(dev_ran, n_ran)
   write(*,'(a)') '----------------------------------------------------------------------'
   write(*,'(a)') ' qpol_ext = 0.5*(Q3-Q2)    [polarized extinction]'
   write(*,'(a,i0,a)')   '   nodes compared         : ', n_pol, &
        '   (|qpol_HD| > 1e-4 * qran_HD)'
   write(*,'(a,es12.4)') '   median |ours/HD23 - 1| : ', median(dev_pol, n_pol)
   write(*,'(a,es12.4)') '   max    |ours/HD23 - 1| : ', maxval_n(dev_pol, n_pol)
   write(*,'(a)') '======================================================================'

contains

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

   subroutine read_qext_jori(gz_file, scratch, qext)
      ! Decompress the gzip'd HD23 table once and read only its Q_ext block.
      ! On-disk stream (see sed/src/q_table_jori.f90): 12 header lines, then
      ! for each quantity (Q_ext, Q_abs, Q_sca), for each orientation
      ! jori = 1..3, for each wavelength (1129 records), one record of the
      ! 169 radii.  Only the first quantity block (Q_ext) is kept.
      character(len=*), intent(in)  :: gz_file, scratch
      real(wp),         intent(out) :: qext(NW, NA, NORI)
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
               if (iq == 1) qext(jwl, :, jori) = row(1:NA)
            end do
         end do
      end do
      close(u)
      call delete_scratch(scratch)
   end subroutine read_qext_jori

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

end program compare_tmatrix_jori
