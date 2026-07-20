program run_tmatrix
   ! Precomputes (Q_ext, Q_abs, Q_sca, albedo, g) for the DH21 astrodust
   ! grain on its native (a_eff, lambda) grid.
   !
   ! Usage:
   !   ./run_tmatrix.x                       ! full sweep, 169 x 1129 points
   !   ./run_tmatrix.x test                  ! subset:  7 x 7  for smoke test
   !   ./run_tmatrix.x range JW1 JW2         ! partial sweep, jw in [JW1, JW2]
   !                                         ! (used by run_parallel.sh)
   !
   ! Output (text, ASCII):
   !   tmatrix/output/q_astrodust_P0.20_Fe0.00_1.400.dat              (full)
   !   tmatrix/output/q_astrodust_P0.20_Fe0.00_1.400.test.dat         (subset)
   !   tmatrix/output/q_astrodust_P0.20_Fe0.00_1.400.jwJW1-JW2.dat    (range)
   !
   ! Columns:
   !   lambda[um]  a_eff[um]  Q_ext  Q_abs  Q_sca  albedo  g  flag
   !
   ! flag legend:
   !    0    T-matrix converged
   !   10    small-x limit (Rayleigh dipole)
   !   20    large-x limit (geometric optics)
   !   1..5  T-matrix returned IERR=1..5 (see tmd_one.f header; IERR=5 is
   !         the Gaussian-quadrature refinement loop failing to converge),
   !         then the result was taken from whichever limit x is closer to:
   !           IERR in 1..5 with x < 1.0    : redirected to the small-x limit
   !                                          (flag = IERR + 10)
   !           IERR in 1..5 with x >= 1.0   : redirected to the large-x limit
   !                                          (flag = IERR + 20)
   !  100+  : failed physical-consistency check (see stderr).  The base flag
   !         is preserved in the low two digits; +100 marks a written row
   !         whose Q / albedo / g fell outside the finiteness, non-negativity,
   !         or range bounds.

   use, intrinsic :: iso_fortran_env, only: real64
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use constants, only: wp
   use read_index, only: load_index, interp_m
   use asymptotic_optics, only: rayleigh_limit, geometric_optics_limit
   implicit none

   ! Reference parameters (HD23 best fit)
   ! Paths are relative to the directory the executable is launched from
   ! (i.e. tmatrix/, where the Makefile drops run_tmatrix.x).
   character(len=*), parameter :: f_aeff  = &
      '../data/dielectric/DH21_aeff'
   character(len=*), parameter :: f_wave  = &
      '../data/dielectric/DH21_wave'
   character(len=*), parameter :: f_index = &
      '../data/dielectric/index_DH21Ad_P0.20_0.00_1.400'
   character(len=*), parameter :: f_out_full = &
      'output/q_astrodust_P0.20_Fe0.00_1.400.dat'
   character(len=*), parameter :: f_out_test = &
      'output/q_astrodust_P0.20_Fe0.00_1.400.test.dat'

   integer, parameter :: NA = 169, NW = 1129
   real(wp), parameter :: EPS_BA  = 1.4_wp
   real(wp), parameter :: DDELT   = 1.0e-3_wp
   integer,  parameter :: NDGS    = 2
   integer,  parameter :: NP_OBL  = -1            ! oblate spheroid, Mishchenko convention
   real(wp), parameter :: X_SMALL = 0.1_wp        ! small-x cutoff
   real(wp), parameter :: X_LARGE = 50.0_wp       ! large-x cutoff
   real(wp), parameter :: PI      = acos(-1.0_wp)

   real(wp) :: a_eff(NA), lambda(NW)
   real(wp) :: nr_cache(NW), ki_cache(NW)
   real(wp) :: nr, ki, x, qext, qabs, qsca, walb, asymm
   integer  :: ja, jw, ierr_t, flag
   integer  :: ja_step, jw_step, jw_lo, jw_hi, n_total, n_done
   integer  :: u_out, ios
   integer  :: n_viol
   logical  :: phys_ok
   character(len=96)  :: viol_reason
   character(len=256) :: f_out
   character(len=32)  :: arg, arg2, arg3
   integer, parameter :: MODE_FULL=0, MODE_TEST=1, MODE_RANGE=2
   integer  :: mode

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

   call read_one_col(f_aeff, NA, a_eff)
   call read_one_col(f_wave, NW, lambda)
   call load_index(f_index)

   ! Cache m(lambda) once per wavelength (lambda-loop outer)
   do jw = 1, NW
      call interp_m(lambda(jw), nr_cache(jw), ki_cache(jw))
   end do

   select case (mode)
   case (MODE_TEST)
      ! Stride to span the full ranges with ~7 x 7 ~ 49 sample points,
      ! exercising small-x, mid-x (T-matrix), and large-x regimes.
      ja_step = 28
      jw_step = 188
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
         x = 2.0_wp * PI * a_eff(ja) / lambda(jw)

         if (x < X_SMALL) then
            call rayleigh_limit(a_eff(ja), lambda(jw), nr, ki, EPS_BA, qext, qsca, walb, asymm)
            qabs = qext - qsca
            flag = 10
         else if (x > X_LARGE) then
            call geometric_optics_limit(a_eff(ja), lambda(jw), nr, ki, EPS_BA, qext, qsca, walb, asymm)
            qabs = qext - qsca
            flag = 20
         else
            call tmd_one(a_eff(ja), lambda(jw), nr, ki, EPS_BA, NP_OBL, &
                         DDELT, NDGS, qext, qsca, walb, asymm, ierr_t)
            if (ierr_t /= 0) then
               if (x < 1.0_wp) then
                  call rayleigh_limit(a_eff(ja), lambda(jw), nr, ki, EPS_BA, qext, qsca, walb, asymm)
                  flag = ierr_t + 10
               else
                  call geometric_optics_limit(a_eff(ja), lambda(jw), nr, ki, EPS_BA, qext, qsca, walb, asymm)
                  flag = ierr_t + 20
               end if
            else
               flag = 0
            end if
            qabs = qext - qsca
         end if

         ! Physical-consistency check.  All three branches (T-matrix,
         ! small-x, large-x) converge here, so validating just before the
         ! write covers every result path.  A failure keeps the row (the
         ! downstream reader counts rows to rebuild the grid) but marks its
         ! flag with +100 and emits a one-line warning.
         call check_physics(qext, qabs, qsca, walb, asymm, phys_ok, viol_reason)
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
   write(*,'(a,a)') ' wrote ', trim(f_out)
   write(*,'(a,i0,a,i0,a)') ' physical-consistency violations: ', n_viol, &
      ' of ', n_done, ' rows'

contains

   subroutine read_one_col(filename, n, x)
      ! DH21_aeff and DH21_wave have 2 header lines, then ALL n values
      ! on a single very long whitespace-separated line. List-directed
      ! read with implied DO across one record handles this correctly.
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

   subroutine check_physics(qe, qa, qs, wa, gg, ok, reason)
      ! Physical-consistency test for one (lambda, a_eff) result.
      ! Sets ok = .false. and appends a short tag to reason for every bound
      ! that is violated; reason is empty when everything passes.
      real(wp),         intent(in)  :: qe, qa, qs, wa, gg   ! qext,qabs,qsca,albedo,g
      logical,          intent(out) :: ok
      character(len=*), intent(out) :: reason
      real(wp), parameter :: tol = 1.0e-9_wp

      reason = ''
      ! Finiteness (NaN / Inf).
      if (.not. ieee_is_finite(qe)) reason = trim(reason)//' qext-nonfinite'
      if (.not. ieee_is_finite(qs)) reason = trim(reason)//' qsca-nonfinite'
      if (.not. ieee_is_finite(qa)) reason = trim(reason)//' qabs-nonfinite'
      if (.not. ieee_is_finite(wa)) reason = trim(reason)//' albedo-nonfinite'
      if (.not. ieee_is_finite(gg)) reason = trim(reason)//' g-nonfinite'
      ! Non-negativity (small negative roundoff tolerated) and qext >= qsca.
      if (ieee_is_finite(qe) .and. qe < -tol) reason = trim(reason)//' qext<0'
      if (ieee_is_finite(qs) .and. qs < -tol) reason = trim(reason)//' qsca<0'
      if (ieee_is_finite(qa) .and. qa < -tol) reason = trim(reason)//' qabs<0'
      if (ieee_is_finite(qe) .and. ieee_is_finite(qs) .and. qe < qs - tol) &
         reason = trim(reason)//' qext<qsca'
      ! Albedo in [0,1] and asymmetry g in [-1,1].
      if (ieee_is_finite(wa) .and. (wa < -tol .or. wa > 1.0_wp + tol)) &
         reason = trim(reason)//' albedo-range'
      if (ieee_is_finite(gg) .and. (gg < -1.0_wp - tol .or. gg > 1.0_wp + tol)) &
         reason = trim(reason)//' g-range'

      ok = (len_trim(reason) == 0)
   end subroutine check_physics

end program run_tmatrix
