program run_tmatrix
   ! Precomputes (Q_ext, Q_abs, Q_sca, albedo, g) for the DH21 astrodust
   ! grain on its native (a_eff, lambda) grid.
   ! Every point goes through `spheroid_q`, which selects the size-parameter
   ! regime and, in the T-matrix regime, calls the default full-direct
   ! reentrant library backend.
   !
   ! Usage:
   !   ./run_tmatrix.x                       ! full sweep, NA x NW points
   !   ./run_tmatrix.x test                  ! subset:  7 x 7  for smoke test
   !   ./run_tmatrix.x range JW1 JW2         ! partial sweep, jw in [JW1, JW2]
   !                                         ! (used by run_parallel.sh)
   !
   ! NA and NW are the lengths of the two grid files, counted at run time; no
   ! grid size is written down here.
   !
   ! Output (text, ASCII):
   !   tmatrix/output/q_astrodust_P0.20_Fe0.00_1.400.dat              (full)
   !   tmatrix/output/q_astrodust_P0.20_Fe0.00_1.400.test.dat         (subset)
   !   tmatrix/output/q_astrodust_P0.20_Fe0.00_1.400.jwJW1-JW2.dat    (range)
   !
   ! Columns:
   !   lambda[um]  a_eff[um]  Q_ext  Q_abs  Q_sca  albedo  g  flag
   !
   ! Regime selection, the flag convention below, and the physical-consistency
   ! bounds all live in driver/spheroid_optics.f90; this program only supplies
   ! the astrodust grid and refractive index and writes the table.
   !
   ! flag legend:
   !    0    T-matrix converged
   !   10    Rayleigh dipole limit (x < X_RAYLEIGH_MAX)
   !   20    geometric optics limit (x > X_GEOMETRIC_MIN)
   !   1..9  T-matrix returned IERR=1..9 (IERR=5 is the Gaussian-quadrature
   !         refinement loop failing to converge and IERR=6..9 are
   !         internal/LAPACK failures),
   !         then the result was taken from whichever limit x is closer to:
   !           IERR in 1..9 with x <  X_LIMIT_SPLIT : Rayleigh dipole limit
   !                                                  (flag = IERR + 10)
   !           IERR in 1..9 with x >= X_LIMIT_SPLIT : geometric optics limit
   !                                                  (flag = IERR + 20)
   !  100+  : failed physical-consistency check (see the warning on stdout).
   !         The base flag is preserved in the low two digits; +100 marks a
   !         written row whose Q / albedo / g fell outside the finiteness,
   !         non-negativity, or range bounds.

   use, intrinsic :: iso_fortran_env, only: iostat_eor
   use tmatrix_api, only: wp, tmatrix_options_t, &
                          tmatrix_workspace_t, tmatrix_workspace_init, &
                          tmatrix_workspace_finalize
   use tmatrix_status, only: TMATRIX_SUCCESS
   use read_index, only: load_index, interp_m
   use spheroid_optics, only: spheroid_q, check_physical_bounds
   implicit none

   ! Reference parameters (HD23 best fit)
   ! Paths are relative to the directory the executable is launched from
   ! (i.e. tmatrix/, where the Makefile drops run_tmatrix.x).
   character(len=*), parameter :: f_aeff  = &
      '../data/dielectric/DH21_aeff'
   ! Wavelength axis.  DH21_wave stops at 0.0912 um (13.6 eV); DH21_wave_to_12keV
   ! is that grid with the dielectric-function's own energy nodes below it
   ! prepended, carrying the table to 1.0e-4 um (12398 eV).  Those nodes are
   ! reused rather than resampled because that band is full of absorption edges
   ! tabulated as steps across a close pair of energies: counting them takes a
   ! threshold, and there are 23 places where k jumps by more than 0.5% across a
   ! pair closer than 5e-4 in relative energy, 14 of them closer than 1e-4.  A
   ! uniform 200-per-decade axis would average each edge into a ramp and lose,
   ! for instance, the +211% jump of k at Fe K (7124 eV).  On the tabulated
   ! nodes interp_m below evaluates at zero interpolation weight and reproduces
   ! the tabulated index exactly on both sides of every edge.
   !
   ! One point of the axis is an exception: 0.0912*(1 - 1e-4), just below the
   ! Lyman limit, is NOT a dielectric-table node, so interp_m interpolates it.
   ! That is harmless for the grain -- the astrodust dielectric function is
   ! smooth across 13.6 eV and has no edge there -- because the point is placed
   ! for the RADIATION FIELD.  An ISRF is zero below the Lyman limit and finite
   ! above it, and the dielectric axis has no node between 13.595 and 14.000 eV,
   ! so without this point a consumer's wavelength integral has one 2.98%-wide
   ! cell straddling that step and manufactures absorption where the field
   ! carries no photon (1.74% too much heating for the smallest grains; 0.006%
   ! with the point in place).  Resolving a step with a close pair of nodes is
   ! what the extension already does for the absorption edges; this applies it
   ! to the one step that belongs to the field rather than to the grain.
   character(len=*), parameter :: f_wave  = &
      '../data/dielectric/DH21_wave_to_12keV'
   character(len=*), parameter :: f_index = &
      '../data/dielectric/index_DH21Ad_P0.20_0.00_1.400'
   ! This driver sweeps the whole f_wave axis, so what it writes is the EUV
   ! product.  The 1129-wavelength companion that ships beside it,
   ! output/q_astrodust_P0.20_Fe0.00_1.400.dat, is this file with every
   ! wavelength shortward of the Lyman limit dropped -- row selection only, no
   ! value recomputed -- and `make lyman_cut` in tmatrix/ is what makes it.
   ! Regenerating the table therefore takes both steps.
   character(len=*), parameter :: f_out_full = &
      'output/q_astrodust_P0.20_Fe0.00_1.400_euv.dat'
   character(len=*), parameter :: f_out_test = &
      'output/q_astrodust_P0.20_Fe0.00_1.400_euv.test.dat'

   real(wp), parameter :: EPS_BA  = 1.4_wp
   real(wp), parameter :: DDELT   = 1.0e-3_wp
   integer,  parameter :: NDGS    = 2
   integer,  parameter :: NP_OBL  = -1            ! oblate spheroid, Mishchenko convention
   ! Number of samples the smoke subset takes along each axis.
   integer,  parameter :: N_SMOKE = 7

   integer  :: NA, NW
   real(wp), allocatable :: a_eff(:), lambda(:)
   real(wp), allocatable :: nr_cache(:), ki_cache(:)
   real(wp) :: nr, ki, qext, qabs, qsca, walb, asymm
   type(tmatrix_workspace_t) :: work
   type(tmatrix_options_t) :: tm_options
   integer  :: ja, jw, flag, tm_status
   integer  :: ja_step, jw_step, jw_lo, jw_hi, n_total, n_done
   integer  :: u_out, ios
   integer  :: n_viol
   logical  :: phys_ok
   character(len=160) :: tm_message
   character(len=96)  :: viol_reason
   character(len=256) :: f_out
   character(len=256) :: output_override
   character(len=32)  :: arg, arg2, arg3
   integer, parameter :: MODE_FULL=0, MODE_TEST=1, MODE_RANGE=2
   integer  :: mode

   ! ---- grids -------------------------------------------------------------
   ! Read before the CLI is parsed: the `range` bounds are checked against NW.
   call read_grid(f_aeff, NA, a_eff)
   call read_grid(f_wave, NW, lambda)
   call load_index(f_index)

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
            write(*,'(a)') ' usage: run_tmatrix.x range JW1 JW2'
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
      case default
         write(*,'(a,a,a)') ' unknown mode "', trim(arg), &
              '" — expected one of: (none), test, range'
         stop 1
      end select
   end if

   tm_options%aspect_ratio = EPS_BA
   tm_options%tolerance    = DDELT
   tm_options%shape        = NP_OBL
   tm_options%ndgs         = NDGS
   ! A workspace that cannot be allocated is a setup failure, not a point
   ! failure: report it and stop before the sweep rather than letting every
   ! point report the same thing.
   call tmatrix_workspace_init(work, tm_status, tm_message)
   if (tm_status /= TMATRIX_SUCCESS) then
      write(*,'(a,a)') ' ERROR: libtmatrix: ', trim(tm_message)
      stop 2
   end if

   ! Cache m(lambda) once per wavelength (lambda-loop outer)
   allocate(nr_cache(NW), ki_cache(NW))
   do jw = 1, NW
      call interp_m(lambda(jw), nr_cache(jw), ki_cache(jw))
   end do

   select case (mode)
   case (MODE_TEST)
      ! Stride to span the full ranges with N_SMOKE x N_SMOKE sample points,
      ! exercising small-x, mid-x (T-matrix), and large-x regimes.  Derived
      ! from the grid lengths so the subset keeps that coverage whatever
      ! wavelength axis is in use (on the 169 x 1129 grid this reproduces the
      ! strides 28 and 188 the subset has always used).
      ja_step = max((NA - 1)/(N_SMOKE - 1), 1)
      jw_step = max((NW - 1)/(N_SMOKE - 1), 1)
      f_out   = f_out_test
   case (MODE_RANGE)
      ! Partial sweep over jw in [jw_lo, jw_hi], full a range.
      ja_step = 1
      jw_step = 1
      ! Filename pattern: q_astrodust_..._.jwJW1-JW2.dat
      write(f_out,'(a,a,i0,a,i0,a)') &
         'output/q_astrodust_P0.20_Fe0.00_1.400', &
         '.jw', jw_lo, '-', jw_hi, '.dat'
   case default
      ja_step = 1
      jw_step = 1
      f_out   = f_out_full
   end select
   ! Test automation can redirect the generated table outside the source tree.
   ! Ordinary standalone operation deliberately keeps its historical output/
   ! paths when TMATRIX_OUTPUT_FILE is not set.
   call get_environment_variable('TMATRIX_OUTPUT_FILE', output_override, status=ios)
   if (ios == 0 .and. len_trim(output_override) > 0) f_out = trim(output_override)
   n_total = ((NA - 1)/ja_step + 1) * ((jw_hi - jw_lo)/jw_step + 1)
   n_done  = 0
   n_viol  = 0

   open(newunit=u_out, file=trim(f_out), status='replace', action='write')
   write(u_out,'(a)')  '# DH21 astrodust, P = 0.20, fFe = 0.00, b/a = 1.4'
   write(u_out,'(a,i0,a,i0)')  '# a_eff stride: every ', ja_step, ' of ', NA
   write(u_out,'(a,i0,a,i0)')  '# lambda stride: every ', jw_step, ' of ', NW
   if (mode == MODE_RANGE) then
      write(u_out,'(a,i0,a,i0)') '# lambda range: jw in ', jw_lo, ' .. ', jw_hi
   end if
   write(u_out,'(a)') '#   lambda[um]   a_eff[um]      Q_ext         Q_abs          Q_sca         albedo         g         flag'

   write(*,'(a,i0)')           ' mode (0=full, 1=test, 2=range) = ', mode
   if (mode == MODE_RANGE) write(*,'(a,i0,a,i0,a)') ' jw range = [', jw_lo, ', ', jw_hi, ']'
   write(*,'(a,i0)')           ' total points = ', n_total
   write(*,'(a,a)')            ' output = ', trim(f_out)

   ! Outer loop on lambda (m is cached); inner loop on a_eff.
   do jw = jw_lo, jw_hi, jw_step
      nr = nr_cache(jw)
      ki = ki_cache(jw)
      do ja = 1, NA, ja_step
         call spheroid_q(work, a_eff(ja), lambda(jw), nr, ki, tm_options, &
                         qext, qabs, qsca, walb, asymm, flag, tm_status, tm_message)

         ! A nonzero status with flag = 0 is a library failure the IERR
         ! redirection does not cover.  Keep the row and let the bounds check
         ! below mark it, rather than abandoning a multi-hour sweep over one
         ! point: the reader rebuilds the grid from the row count.
         if (tm_status /= TMATRIX_SUCCESS .and. flag == 0) then
            write(*,'(a,es13.6,a,es13.6,a,a)') &
               ' WARNING libtmatrix: lambda=', lambda(jw), &
               ' um  a_eff=', a_eff(ja), ' um  ', trim(tm_message)
         end if

         ! Physical-consistency check.  All three regimes (Rayleigh dipole,
         ! T-matrix, geometric optics) converge here, so validating just
         ! before the write covers every result path.  A failure keeps the
         ! row but marks its flag with +100 and emits a one-line warning.
         call check_physical_bounds(qext, qabs, qsca, walb, asymm, phys_ok, viol_reason)
         if (.not. phys_ok) then
            n_viol = n_viol + 1
            write(*,'(a,es13.6,a,es13.6,a,i0,a,a)') &
               ' WARNING physical-consistency: lambda=', lambda(jw), &
               ' um  a_eff=', a_eff(ja), ' um  flag=', flag, &
               '  broken:', trim(viol_reason)
            flag = flag + 100
         end if

         write(u_out,'(2es15.6,5es15.6,i6)') &
            lambda(jw), a_eff(ja), qext, qabs, qsca, walb, asymm, flag

         n_done = n_done + 1
         if (mod(n_done, max(n_total/20, 1)) == 0) then
            write(*,'(a,i0,a,i0,a,f6.1,a)') ' progress: ', n_done, '/', &
                  n_total, '  (', 100.0_wp * real(n_done,wp) / real(n_total,wp), '%)'
         end if
      end do
   end do
   close(u_out)
   call tmatrix_workspace_finalize(work)
   write(*,'(a,a)') ' wrote ', trim(f_out)
   write(*,'(a,i0,a,i0,a)') ' physical-consistency violations: ', n_viol, &
      ' of ', n_done, ' rows'

contains

   subroutine read_grid(filename, n, x)
      ! DH21_aeff, DH21_wave, and DH21_wave_to_12keV have 2 header lines, then
      ! ALL values on a single very long whitespace-separated line.
      !
      ! The count is taken from the record itself rather than from a parameter
      ! here, so a grid file and this program cannot fall out of step: swapping
      ! the wavelength axis for a longer one needs no edit below.  The record is
      ! accumulated with non-advancing reads because it is tens of kilobytes
      ! long and no fixed buffer should have to bound it.
      character(len=*),      intent(in)  :: filename
      integer,               intent(out) :: n
      real(wp), allocatable, intent(out) :: x(:)
      integer :: u, ios, nch, i
      logical :: in_token
      character(len=4096) :: chunk
      character(len=512)  :: header
      character(len=:), allocatable :: record

      open(newunit=u, file=filename, status='old', action='read', iostat=ios)
      if (ios /= 0) then
         write(*,'(a,a)') ' ERROR: cannot open ', trim(filename)
         stop 1
      end if
      read(u,'(a)') header
      read(u,'(a)') header
      record = ''
      reading: do
         do
            read(u,'(a)',advance='no',size=nch,iostat=ios) chunk
            if (nch > 0) record = record // chunk(1:nch)
            if (ios /= 0) exit
         end do
         record = record // ' '
         if (ios /= iostat_eor) exit reading    ! end of file, or a read error
      end do reading
      close(u)

      n = 0
      in_token = .false.
      do i = 1, len(record)
         if (record(i:i) == ' ' .or. record(i:i) == achar(9)) then
            in_token = .false.
         else if (.not. in_token) then
            in_token = .true.
            n = n + 1
         end if
      end do
      if (n < 2) then
         write(*,'(a,a)') ' ERROR: fewer than 2 grid values in ', trim(filename)
         stop 1
      end if
      allocate(x(n))
      read(record,*) x(1:n)

      do i = 2, n
         if (x(i) <= x(i-1)) then
            write(*,'(a,a,a,i0)') ' ERROR: ', trim(filename), &
               ' is not strictly increasing at i = ', i
            stop 1
         end if
      end do
   end subroutine read_grid

end program run_tmatrix
