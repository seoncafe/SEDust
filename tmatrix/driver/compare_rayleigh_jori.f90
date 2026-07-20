program compare_rayleigh_jori
   ! Stage-A check of the orientation-resolved Rayleigh cross sections.
   !
   ! Over the small-x subset (x = 2 pi a_eff / lambda < 0.1) of the DH21
   ! astrodust grid, this compares the oriented Q_abs that the extended
   ! rayleigh_limit produces against the Hensley & Draine (2021/2023)
   ! orientation-resolved table.  Two emission-relevant combinations are
   ! compared, in the convention of sed/src/q_table_jori.f90:
   !
   !   qpol_abs = 0.5 * (Q_abs(jori=3) - Q_abs(jori=2))   polarized absorption
   !   qran_abs = (Q_abs(1) + Q_abs(2) + Q_abs(3)) / 3    mean absorption
   !
   ! The refractive index m(lambda) is read the same way run_tmatrix.f90
   ! reads it (read_index::load_index / interp_m), the spheroid shape is the
   ! same oblate b/a = 1.4, and the Rayleigh branch is entered under the same
   ! x < 0.1 cutoff.
   !
   ! Orientation convention (jori index): the Rayleigh limit predicts
   ! Q(jori=1) = Q(jori=3) (both transverse to the axis), so the HD23 table's
   ! own |Q(1)-Q(3)|/Q(1) is reported as a check that the shipped file carries
   ! that property in this regime.

   use, intrinsic :: iso_fortran_env, only: real64, error_unit
   use constants,       only: wp
   use read_index,      only: load_index, interp_m
   use asymptotic_optics, only: rayleigh_limit
   implicit none

   ! Reference parameters and paths (relative to tmatrix/, where the
   ! executable is launched from), matching run_tmatrix.f90.
   character(len=*), parameter :: f_aeff  = '../data/dielectric/DH21_aeff'
   character(len=*), parameter :: f_wave  = '../data/dielectric/DH21_wave'
   character(len=*), parameter :: f_index = &
      '../data/dielectric/index_DH21Ad_P0.20_0.00_1.400'
   character(len=*), parameter :: f_qjori = &
      '../data/dielectric/q_DH21Ad_P0.20_Fe0.00_1.400.dat.gz'
   character(len=*), parameter :: f_scratch = 'q_jori_cmp_scratch.dat'

   integer,  parameter :: NA = 169, NW = 1129, NORI = 3, NHEAD = 12
   real(wp), parameter :: EPS_BA  = 1.4_wp
   real(wp), parameter :: X_SMALL = 0.1_wp
   real(wp), parameter :: PI      = acos(-1.0_wp)

   real(wp) :: a_eff(NA), lambda(NW)
   real(wp) :: nr_cache(NW), ki_cache(NW)
   ! HD23 orientation-resolved table, (jw, ja, jori).
   real(wp), allocatable :: qabs_hd(:,:,:)
   ! Deviation samples collected over the x < 0.1 subset.
   real(wp), allocatable :: dev_pol(:), dev_ran(:), sym_hd(:)
   real(wp) :: qext, qsca, walb, asymm
   real(wp) :: qext_ori(NORI), qabs_ori(NORI), qsca_ori(NORI)
   real(wp) :: x, qpol_ours, qran_ours, qpol_hd, qran_hd, qabs1, qabs3
   real(wp) :: floor_pol
   integer  :: jw, ja, n_small, n_pol, n_ran, n_sym

   ! ---- grids and index --------------------------------------------------
   call read_one_col(f_aeff, NA, a_eff)
   call read_one_col(f_wave, NW, lambda)
   call load_index(f_index)
   do jw = 1, NW
      call interp_m(lambda(jw), nr_cache(jw), ki_cache(jw))
   end do

   ! ---- HD23 orientation-resolved Q_abs table ----------------------------
   allocate(qabs_hd(NW, NA, NORI))
   call read_qabs_jori(f_qjori, f_scratch, qabs_hd)

   ! ---- sweep the small-x subset -----------------------------------------
   allocate(dev_pol(NW*NA), dev_ran(NW*NA), sym_hd(NW*NA))
   n_small = 0;  n_pol = 0;  n_ran = 0;  n_sym = 0

   do jw = 1, NW
      do ja = 1, NA
         x = 2.0_wp * PI * a_eff(ja) / lambda(jw)
         if (x >= X_SMALL) cycle
         n_small = n_small + 1

         call rayleigh_limit(a_eff(ja), lambda(jw), nr_cache(jw), ki_cache(jw), &
                             EPS_BA, qext, qsca, walb, asymm, &
                             qext_ori=qext_ori, qabs_ori=qabs_ori, qsca_ori=qsca_ori)

         qpol_ours = 0.5_wp * (qabs_ori(3) - qabs_ori(2))
         qran_ours = (qabs_ori(1) + qabs_ori(2) + qabs_ori(3)) / 3.0_wp

         qpol_hd = 0.5_wp * (qabs_hd(jw,ja,3) - qabs_hd(jw,ja,2))
         qran_hd = (qabs_hd(jw,ja,1) + qabs_hd(jw,ja,2) + qabs_hd(jw,ja,3)) / 3.0_wp

         ! Mean absorption: qran_hd is strictly positive here, so the ratio
         ! is always well defined.
         if (abs(qran_hd) > 0.0_wp) then
            n_ran = n_ran + 1
            dev_ran(n_ran) = abs(qran_ours/qran_hd - 1.0_wp)
         end if

         ! Polarized absorption: qpol can pass through zero, so form the
         ! ratio only where the HD23 value carries meaningful signal
         ! (measured against the local mean absorption).  This guards both
         ! divide-by-zero and near-zero sign flips.
         floor_pol = 1.0e-4_wp * abs(qran_hd)
         if (abs(qpol_hd) > floor_pol) then
            n_pol = n_pol + 1
            dev_pol(n_pol) = abs(qpol_ours/qpol_hd - 1.0_wp)
         end if

         ! HD23 self-check: jori=1 and jori=3 should coincide in the Rayleigh
         ! limit.  Reported as |Q(1)-Q(3)|/Q(1) on Q_abs.
         qabs1 = qabs_hd(jw,ja,1)
         qabs3 = qabs_hd(jw,ja,3)
         if (abs(qabs1) > 0.0_wp) then
            n_sym = n_sym + 1
            sym_hd(n_sym) = abs(qabs1 - qabs3) / abs(qabs1)
         end if
      end do
   end do

   ! ---- report -----------------------------------------------------------
   write(*,'(a)') '======================================================================'
   write(*,'(a)') ' Stage-A Rayleigh oriented Q_abs vs HD23 (x = 2 pi a_eff/lambda < 0.1)'
   write(*,'(a)') '======================================================================'
   write(*,'(a,i0,a,i0,a)') ' grid nodes total          : ', NW*NA, &
        '   (', NW, ' lambda x 169 a_eff)'
   write(*,'(a,i0)')        ' nodes with x < 0.1        : ', n_small
   write(*,'(a)') '----------------------------------------------------------------------'
   write(*,'(a)') ' qran_abs = (Q1+Q2+Q3)/3   [mean absorption]'
   write(*,'(a,i0)')          '   nodes compared          : ', n_ran
   write(*,'(a,es12.4)')      '   median |ours/HD23 - 1|  : ', median(dev_ran, n_ran)
   write(*,'(a,es12.4)')      '   max    |ours/HD23 - 1|  : ', maxval_n(dev_ran, n_ran)
   write(*,'(a)') '----------------------------------------------------------------------'
   write(*,'(a)') ' qpol_abs = 0.5*(Q3-Q2)    [polarized absorption]'
   write(*,'(a,i0,a)')        '   nodes compared          : ', n_pol, &
        '   (|qpol_HD| > 1e-4 * qran_HD)'
   write(*,'(a,es12.4)')      '   median |ours/HD23 - 1|  : ', median(dev_pol, n_pol)
   write(*,'(a,es12.4)')      '   max    |ours/HD23 - 1|  : ', maxval_n(dev_pol, n_pol)
   write(*,'(a)') '----------------------------------------------------------------------'
   write(*,'(a)') ' HD23 self-check: |Q_abs(jori=1) - Q_abs(jori=3)| / Q_abs(jori=1)'
   write(*,'(a,i0)')          '   nodes                   : ', n_sym
   write(*,'(a,es12.4)')      '   median                  : ', median(sym_hd, n_sym)
   write(*,'(a,es12.4)')      '   max                     : ', maxval_n(sym_hd, n_sym)
   write(*,'(a)') '======================================================================'

contains

   subroutine read_one_col(filename, n, arr)
      ! DH21_aeff / DH21_wave: 2 header lines, then n free-format values.
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

   subroutine read_qabs_jori(gz_file, scratch, qabs)
      ! Decompress the gzip'd HD23 table once and read only its Q_abs block.
      !
      ! On-disk stream (see sed/src/q_table_jori.f90): 12 header lines, then
      ! for each quantity (Q_ext, Q_abs, Q_sca) in turn, for each orientation
      ! jori = 1..3, for each wavelength jw = 1..1129, one record of the 169
      ! effective radii.  Only the second quantity block (Q_abs) is kept.
      character(len=*), intent(in)  :: gz_file, scratch
      real(wp),         intent(out) :: qabs(NW, NA, NORI)
      integer :: u, ios, estat, cstat, iq, jori, jwl, i
      real(wp) :: row(NA)
      character(len=512) :: line

      call execute_command_line('gzip -dc "'//trim(gz_file)//'" > "'// &
                                trim(scratch)//'"', exitstat=estat, cmdstat=cstat)
      if (cstat /= 0 .or. estat /= 0) then
         write(error_unit,'(a,a)') ' ERROR: gzip -dc failed on ', trim(gz_file)
         call delete_scratch(scratch)
         stop 1
      end if

      open(newunit=u, file=trim(scratch), status='old', action='read', iostat=ios)
      if (ios /= 0) then
         write(error_unit,'(a,a)') ' ERROR: cannot open ', trim(scratch)
         call delete_scratch(scratch)
         stop 1
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
               if (iq == 2) qabs(jwl, :, jori) = row(1:NA)
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
      ! Median of v(1:n); returns 0 for n = 0.  Sorts a local copy.
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

end program compare_rayleigh_jori
