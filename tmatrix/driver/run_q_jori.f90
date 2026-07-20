program run_q_jori
   ! Precomputes the orientation-resolved optics (Q_ext, Q_abs, Q_sca for
   ! jori = 1, 2, 3) of the DH21 astrodust spheroid on its native
   ! (a_eff, lambda) grid, and writes them in the exact format that
   ! sed/src/q_table_jori.f90::load_q_table_jori reads.
   !
   ! Each (a_eff, lambda) node is handed to oriented_cross_sections
   ! (tmatrix/driver/tmatrix_oriented.f90), the shared first-principles core
   ! that selects the size-parameter regime and returns the three orientations.
   !
   ! Usage:
   !   ./run_q_jori.x                    ! full sweep, 169 x 1129 points
   !   ./run_q_jori.x test               ! ~7 x 7 sample + full-sweep time estimate
   !   ./run_q_jori.x range JW1 JW2      ! full a range, jw in [JW1, JW2]
   !   ./run_q_jori.x merge FILE ...     ! assemble range outputs into the full file
   !
   ! Output (text, ASCII), in tmatrix/output/:
   !   q_astrodust_jori_P0.20_Fe0.00_1.400.dat              (full / merge)
   !   q_astrodust_jori_P0.20_Fe0.00_1.400.test.dat         (test, diagnostic columns)
   !   q_astrodust_jori_P0.20_Fe0.00_1.400.jwJW1-JW2.dat    (range)
   !
   ! Stream format of the full and range files (see the 12-line header written
   ! below and q_table_jori.f90:198-223): 12 header lines, then free-format
   !   do iq = 1, 3          ! iq = 1 Q_ext, 2 Q_abs, 3 Q_sca
   !     do jori = 1, 3
   !       do jw = 1, NW
   !         one record of NA = 169 a_eff values
   ! A range file writes only the jw in [JW1, JW2] records but keeps the same
   ! (iq, jori, jw) nesting.  Because jw is the innermost of the three loops,
   ! jw-window files are NOT concatenable into the full stream with a plain
   ! cat: the full stream needs all 1129 jw records contiguous inside each
   ! (iq, jori) block, whereas cat would interleave the blocks.  Use the
   ! `merge` mode, which reads the range files and re-emits them in full
   ! (iq, jori, jw) order.

   use, intrinsic :: iso_fortran_env, only: real64, int64, error_unit
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use constants, only: wp
   use read_index, only: load_index, interp_m
   use tmatrix_oriented, only: oriented_cross_sections
   implicit none

   ! Reference parameters (HD23 best fit).  Paths are relative to tmatrix/,
   ! where the Makefile drops run_q_jori.x.
   character(len=*), parameter :: f_aeff  = '../data/dielectric/DH21_aeff'
   character(len=*), parameter :: f_wave  = '../data/dielectric/DH21_wave'
   character(len=*), parameter :: f_index = &
      '../data/dielectric/index_DH21Ad_P0.20_0.00_1.400'
   character(len=*), parameter :: f_out_full = &
      'output/q_astrodust_jori_P0.20_Fe0.00_1.400.dat'
   character(len=*), parameter :: f_out_test = &
      'output/q_astrodust_jori_P0.20_Fe0.00_1.400.test.dat'

   integer, parameter :: NA = 169, NW = 1129, NORI = 3
   real(wp), parameter :: EPS_BA  = 1.4_wp
   real(wp), parameter :: DDELT   = 1.0e-3_wp
   integer,  parameter :: NDGS    = 2
   integer,  parameter :: NP_OBL  = -1            ! oblate spheroid, Mishchenko convention
   real(wp), parameter :: X_SMALL = 0.1_wp
   real(wp), parameter :: X_LARGE = 50.0_wp
   real(wp), parameter :: PI      = acos(-1.0_wp)

   real(wp) :: a_eff(NA), lambda(NW)
   real(wp) :: nr_cache(NW), ki_cache(NW)
   complex(wp) :: m_cache(NW)
   real(wp) :: nr, ki

   integer, parameter :: MODE_FULL=0, MODE_TEST=1, MODE_RANGE=2, MODE_MERGE=3
   integer  :: mode
   integer  :: jw_lo, jw_hi
   character(len=256) :: f_out
   character(len=64)  :: arg, arg2, arg3
   integer  :: ios, jw

   ! ---- CLI ---------------------------------------------------------------
   mode  = MODE_FULL
   jw_lo = 1
   jw_hi = NW
   if (command_argument_count() >= 1) then
      call get_command_argument(1, arg)
      select case (trim(arg))
      case ('test')
         mode = MODE_TEST
      case ('range')
         if (command_argument_count() < 3) then
            write(*,'(a)') ' usage: run_q_jori.x range JW1 JW2'
            stop 1
         end if
         mode = MODE_RANGE
         call get_command_argument(2, arg2)
         call get_command_argument(3, arg3)
         read(arg2,*,iostat=ios) jw_lo
         if (ios /= 0) then; write(*,'(a)') ' bad JW1'; stop 1; end if
         read(arg3,*,iostat=ios) jw_hi
         if (ios /= 0) then; write(*,'(a)') ' bad JW2'; stop 1; end if
         if (jw_lo < 1 .or. jw_hi > NW .or. jw_lo > jw_hi) then
            write(*,'(a,i0,a,i0,a,i0,a)') ' JW range out of bounds: [', &
               jw_lo, ', ', jw_hi, '] not in [1, ', NW, ']'
            stop 1
         end if
      case ('merge')
         mode = MODE_MERGE
      case default
         write(*,'(a,a,a)') ' unknown mode "', trim(arg), &
              '" -- expected one of: (none), test, range, merge'
         stop 1
      end select
   end if

   ! merge does not need the optics core; it only reassembles range files.
   if (mode == MODE_MERGE) then
      call merge_range_files()
      stop 0
   end if

   call read_one_col(f_aeff, NA, a_eff)
   call read_one_col(f_wave, NW, lambda)
   call load_index(f_index)
   do jw = 1, NW
      call interp_m(lambda(jw), nr, ki)
      nr_cache(jw) = nr
      ki_cache(jw) = ki
      m_cache(jw)  = cmplx(nr, ki, kind=wp)
   end do

   select case (mode)
   case (MODE_TEST)
      call run_test()
   case (MODE_RANGE)
      write(f_out,'(a,i0,a,i0,a)') &
         'output/q_astrodust_jori_P0.20_Fe0.00_1.400.jw', jw_lo, '-', jw_hi, '.dat'
      call sweep_and_write(jw_lo, jw_hi, trim(f_out))
   case default
      call sweep_and_write(1, NW, f_out_full)
   end select

contains

   subroutine sweep_and_write(jlo, jhi, fname)
      ! Full/range sweep: compute the three orientations at every (jw in
      ! [jlo,jhi], ja) node and write them in (iq, jori, jw) stream order.
      integer,          intent(in) :: jlo, jhi
      character(len=*), intent(in) :: fname
      real(wp), allocatable :: qe(:,:,:), qa(:,:,:), qs(:,:,:)  ! (NA, NORI, jlo:jhi)
      real(wp) :: qext_ori(3), qabs_ori(3), qsca_ori(3)
      integer  :: ja, jjw, flag, u, n_total, n_done

      allocate(qe(NA, NORI, jlo:jhi), qa(NA, NORI, jlo:jhi), qs(NA, NORI, jlo:jhi))

      n_total = (jhi - jlo + 1) * NA
      n_done  = 0
      write(*,'(a)')        ' mode = sweep'
      write(*,'(a,i0,a,i0,a)') ' jw range = [', jlo, ', ', jhi, ']'
      write(*,'(a,i0)')     ' total points = ', n_total
      write(*,'(a,a)')      ' output = ', trim(fname)

      do jjw = jlo, jhi
         do ja = 1, NA
            call oriented_cross_sections(a_eff(ja), lambda(jjw), m_cache(jjw), &
                     EPS_BA, NP_OBL, DDELT, NDGS, qext_ori, qabs_ori, qsca_ori, flag)
            qe(ja, :, jjw) = qext_ori
            qa(ja, :, jjw) = qabs_ori
            qs(ja, :, jjw) = qsca_ori
            n_done = n_done + 1
         end do
         if (mod(jjw - jlo + 1, max((jhi-jlo+1)/20, 1)) == 0) then
            write(*,'(a,i0,a,i0,a,f6.1,a)') ' progress: ', n_done, '/', &
               n_total, '  (', 100.0_wp*real(n_done,wp)/real(n_total,wp), '%)'
         end if
      end do

      open(newunit=u, file=trim(fname), status='replace', action='write')
      call write_header(u, jlo, jhi)
      call write_block(u, qe, jlo, jhi)   ! iq = 1  Q_ext
      call write_block(u, qa, jlo, jhi)   ! iq = 2  Q_abs
      call write_block(u, qs, jlo, jhi)   ! iq = 3  Q_sca
      close(u)
      write(*,'(a,a)') ' wrote ', trim(fname)

      deallocate(qe, qa, qs)
   end subroutine sweep_and_write


   subroutine write_block(u, q, jlo, jhi)
      ! One quantity: for jori = 1..3, for jw = jlo..jhi, one record of NA
      ! sizes.  Matches the inner two loops of q_table_jori.f90:201-221.
      integer,  intent(in) :: u, jlo, jhi
      real(wp), intent(in) :: q(NA, NORI, jlo:jhi)
      integer :: jori, jjw
      do jori = 1, NORI
         do jjw = jlo, jhi
            write(u,'(*(es13.5))') q(1:NA, jori, jjw)
         end do
      end do
   end subroutine write_block


   subroutine write_header(u, jlo, jhi)
      ! Exactly 12 header lines, so load_q_table_jori (NHEAD = 12) skips them.
      ! Line 12 carries a machine-readable jw window that merge mode parses.
      integer, intent(in) :: u, jlo, jhi
      write(u,'(a)') '# DH21 astrodust orientation-resolved Q, P = 0.20, fFe = 0.00, b/a = 1.4'
      write(u,'(a)') '# Q = C/(pi a_eff^2); a_eff from DH21_aeff (169), lambda from DH21_wave (1129).'
      write(u,'(a)') '# jori convention (a = spheroid symmetry axis):'
      write(u,'(a)') '#   jori=1: k || a'
      write(u,'(a)') '#   jori=2: k perp a, E || a'
      write(u,'(a)') '#   jori=3: k perp a, E perp a'
      write(u,'(a)') '# Stream order after this header (free format, one record = 169 a_eff values):'
      write(u,'(a)') '#   do iq=1,3 ; do jori=1,3 ; do jw=1,1129 ; write one record ; end'
      write(u,'(a)') '#   iq = 1 Q_ext, 2 Q_abs, 3 Q_sca'
      write(u,'(a)') '# Regime: x<0.1 Rayleigh, 0.1<=x<=50 T-matrix, x>50 geometric optics.'
      write(u,'(a)') '# Generated by tmatrix/driver/run_q_jori.f90.'
      write(u,'(a,i0,1x,i0)') '# JW_WINDOW ', jlo, jhi
   end subroutine write_header


   subroutine run_test()
      ! Strided ~7x7 sample.  Times each node, classifies it as cheap
      ! (Rayleigh / geometric optics) or expensive (T-matrix attempted), and
      ! extrapolates the two averages onto the full grid for a wall-time
      ! estimate.  Writes a diagnostic columnar file (not the stream format).
      integer, parameter :: JA_STEP = 28, JW_STEP = 188
      integer  :: ja, jjw, flag, u
      real(wp) :: qext_ori(3), qabs_ori(3), qsca_ori(3), x
      integer(int64) :: c0, c1, crate
      real(wp) :: dt
      integer  :: n_cheap_s, n_exp_s
      real(wp) :: t_cheap_s, t_exp_s
      integer  :: n_cheap_f, n_exp_f
      real(wp) :: avg_cheap, avg_exp, est
      logical  :: expensive

      call system_clock(count_rate=crate)

      ! Full-grid regime census (arithmetic only, no solve).
      n_cheap_f = 0;  n_exp_f = 0
      do jjw = 1, NW
         do ja = 1, NA
            x = 2.0_wp*PI*a_eff(ja)/lambda(jjw)
            if (x < X_SMALL .or. x > X_LARGE) then
               n_cheap_f = n_cheap_f + 1
            else
               n_exp_f = n_exp_f + 1
            end if
         end do
      end do

      open(newunit=u, file=f_out_test, status='replace', action='write')
      write(u,'(a)') '# run_q_jori test sample (strided).  Columns:'
      write(u,'(a)') '#  lambda[um]  a_eff[um]  x  flag  Qext(1:3)  Qabs(1:3)  Qsca(1:3)'

      n_cheap_s = 0;  n_exp_s = 0
      t_cheap_s = 0.0_wp;  t_exp_s = 0.0_wp
      write(*,'(a)') ' mode = test (strided sample)'
      do jjw = 1, NW, JW_STEP
         do ja = 1, NA, JA_STEP
            x = 2.0_wp*PI*a_eff(ja)/lambda(jjw)
            expensive = (x >= X_SMALL .and. x <= X_LARGE)
            call system_clock(c0)
            call oriented_cross_sections(a_eff(ja), lambda(jjw), m_cache(jjw), &
                     EPS_BA, NP_OBL, DDELT, NDGS, qext_ori, qabs_ori, qsca_ori, flag)
            call system_clock(c1)
            dt = real(c1 - c0, kind=wp) / real(crate, kind=wp)
            if (expensive) then
               n_exp_s = n_exp_s + 1;  t_exp_s = t_exp_s + dt
            else
               n_cheap_s = n_cheap_s + 1;  t_cheap_s = t_cheap_s + dt
            end if
            write(u,'(2es13.5,es12.4,i5,9es13.5)') lambda(jjw), a_eff(ja), x, flag, &
               qext_ori, qabs_ori, qsca_ori
            write(*,'(a,es11.4,a,es11.4,a,es10.3,a,i4)') &
               '  lam=', lambda(jjw), ' a=', a_eff(ja), ' x=', x, ' flag=', flag
         end do
      end do
      close(u)

      avg_cheap = 0.0_wp;  if (n_cheap_s > 0) avg_cheap = t_cheap_s/real(n_cheap_s,wp)
      avg_exp   = 0.0_wp;  if (n_exp_s   > 0) avg_exp   = t_exp_s  /real(n_exp_s,wp)
      est = real(n_cheap_f,wp)*avg_cheap + real(n_exp_f,wp)*avg_exp

      write(*,'(a)')            ' --- full-sweep wall-time estimate ---'
      write(*,'(a,a)')          ' wrote ', f_out_test
      write(*,'(a,i0,a,es10.3,a)') ' cheap nodes sampled : ', n_cheap_s, &
         '   avg ', avg_cheap, ' s'
      write(*,'(a,i0,a,es10.3,a)') ' T-matrix nodes samp.: ', n_exp_s, &
         '   avg ', avg_exp, ' s'
      write(*,'(a,i0,a,i0)')    ' full grid: cheap = ', n_cheap_f, '   T-matrix = ', n_exp_f
      write(*,'(a,f10.1,a,f8.2,a,f7.3,a)') ' estimated full sweep: ', est, &
         ' s  = ', est/60.0_wp, ' min  = ', est/3600.0_wp, ' h'
   end subroutine run_test


   subroutine merge_range_files()
      ! Reads the range files named on the command line (argv 2 onward),
      ! places each jw record into the full (iq, jori, jw) arrays, verifies the
      ! union covers jw = 1..NW with no gaps or overlaps, and writes the full
      ! file in q_table_jori.f90 stream order.
      real(wp), allocatable :: qe(:,:,:), qa(:,:,:), qs(:,:,:)  ! (NA, NORI, NW)
      logical  :: covered(NW)
      integer  :: nfiles, k, u, ios, i, iq, jori, jjw, jlo, jhi
      character(len=256) :: path
      character(len=512) :: line
      real(wp) :: row(NA)

      nfiles = command_argument_count() - 1
      if (nfiles < 1) then
         write(*,'(a)') ' usage: run_q_jori.x merge FILE [FILE ...]'
         stop 1
      end if

      allocate(qe(NA, NORI, NW), qa(NA, NORI, NW), qs(NA, NORI, NW))
      qe = 0.0_wp;  qa = 0.0_wp;  qs = 0.0_wp
      covered = .false.

      do k = 1, nfiles
         call get_command_argument(k+1, path)
         open(newunit=u, file=trim(path), status='old', action='read', iostat=ios)
         if (ios /= 0) then
            write(error_unit,'(a,a)') ' merge: cannot open ', trim(path)
            stop 1
         end if
         ! Read 12 header lines; recover the jw window from the JW_WINDOW line.
         jlo = 0;  jhi = -1
         do i = 1, 12
            read(u,'(a)',iostat=ios) line
            if (ios /= 0) then
               write(error_unit,'(a,a)') ' merge: short header in ', trim(path)
               stop 1
            end if
            if (index(line, 'JW_WINDOW') > 0) &
               read(line(index(line,'JW_WINDOW')+9:), *) jlo, jhi
         end do
         if (jlo < 1 .or. jhi > NW .or. jlo > jhi) then
            write(error_unit,'(a,a)') ' merge: bad JW_WINDOW in ', trim(path)
            stop 1
         end if
         ! Body: (iq, jori, jw in [jlo,jhi]) records of NA values.
         do iq = 1, 3
            do jori = 1, NORI
               do jjw = jlo, jhi
                  read(u,*,iostat=ios) row
                  if (ios /= 0) then
                     write(error_unit,'(a,a)') ' merge: short body in ', trim(path)
                     stop 1
                  end if
                  select case (iq)
                  case (1);  qe(:, jori, jjw) = row
                  case (2);  qa(:, jori, jjw) = row
                  case (3);  qs(:, jori, jjw) = row
                  end select
               end do
            end do
         end do
         close(u)
         do jjw = jlo, jhi
            if (covered(jjw)) then
               write(error_unit,'(a,i0)') ' merge: jw covered twice at jw=', jjw
               stop 1
            end if
            covered(jjw) = .true.
         end do
         write(*,'(a,a,a,i0,a,i0,a)') ' read ', trim(path), '  (jw ', jlo, '..', jhi, ')'
      end do

      do jjw = 1, NW
         if (.not. covered(jjw)) then
            write(error_unit,'(a,i0)') ' merge: jw not covered at jw=', jjw
            stop 1
         end if
      end do

      open(newunit=u, file=f_out_full, status='replace', action='write')
      call write_header(u, 1, NW)
      call write_block(u, qe, 1, NW)
      call write_block(u, qa, 1, NW)
      call write_block(u, qs, 1, NW)
      close(u)
      write(*,'(a,a)') ' wrote ', f_out_full
      deallocate(qe, qa, qs)
   end subroutine merge_range_files


   subroutine read_one_col(filename, n, x)
      ! DH21_aeff / DH21_wave: 2 header lines, then all n values on one long
      ! whitespace-separated record.
      character(len=*), intent(in)  :: filename
      integer,          intent(in)  :: n
      real(wp),         intent(out) :: x(n)
      integer :: u, ios
      character(len=512) :: header
      open(newunit=u, file=filename, status='old', action='read', iostat=ios)
      if (ios /= 0) then
         write(*,'(a,a)') ' ERROR: cannot open ', trim(filename)
         stop 1
      end if
      read(u,'(a)') header
      read(u,'(a)') header
      read(u,*) x(1:n)
      close(u)
   end subroutine read_one_col

end program run_q_jori
