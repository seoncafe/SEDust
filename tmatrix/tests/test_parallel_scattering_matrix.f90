program test_parallel_scattering_matrix
   !! The scattering-matrix expansion and the fixed-orientation amplitude
   !! matrix must reproduce their serial values when distinct caller-owned
   !! workspaces are used concurrently.  Run this executable with 1, 2, 4,
   !! and 8 OpenMP threads.
   !!
   !! Both entry points are deterministic and touch only their own
   !! workspace, so the comparison is bitwise: any difference means state
   !! leaked between workspaces.
   use, intrinsic :: iso_fortran_env, only: real64, int64
   use omp_lib, only: omp_get_thread_num, omp_get_max_threads
   use tmatrix_api, only: tmatrix_options_t, tmatrix_result_t, tmatrix_workspace_t, &
                          tmatrix_scatmat_t, tmatrix_workspace_init, tmatrix_workspace_finalize, &
                          tmatrix_eval_scattering_matrix, ampl
   use tmatrix_test_cases, only: n_cases, tolerance, ndgs, get_case
   implicit none

   integer, parameter :: n_angles = 4
   real(real64), parameter :: angle_tl(n_angles)    = [0.0d0, 90.0d0, 30.0d0,  45.0d0]
   real(real64), parameter :: angle_tl1(n_angles)   = [0.0d0, 90.0d0, 60.0d0, 135.0d0]
   real(real64), parameter :: angle_pl(n_angles)    = [0.0d0,  0.0d0,  0.0d0,  20.0d0]
   real(real64), parameter :: angle_pl1(n_angles)   = [0.0d0,  0.0d0, 90.0d0, 200.0d0]
   real(real64), parameter :: angle_alpha(n_angles) = [0.0d0,  0.0d0,  0.0d0,  30.0d0]
   real(real64), parameter :: angle_beta(n_angles)  = [0.0d0,  0.0d0,  0.0d0,  60.0d0]

   type(tmatrix_workspace_t), allocatable :: work(:)
   type(tmatrix_options_t) :: options
   type(tmatrix_result_t) :: serial_result(n_cases), parallel_result(n_cases)
   type(tmatrix_scatmat_t) :: serial_scatmat(n_cases), parallel_scatmat(n_cases)
   complex(real64) :: serial_amp(4,n_angles,n_cases), parallel_amp(4,n_angles,n_cases)
   integer :: serial_amp_status(n_angles,n_cases), parallel_amp_status(n_angles,n_cases)
   integer :: i, nwork, failures
   real(real64) :: a, lambda, nr, ki, aspect_ratio
   integer :: shape

   options%tolerance = tolerance
   options%ndgs = ndgs
   nwork = max(1, omp_get_max_threads())
   allocate(work(nwork))
   do i = 1, nwork
      call tmatrix_workspace_init(work(i))
      if (.not. work(i)%initialized) error stop 'scattering-matrix workspace allocation failed'
   end do

   do i = 1, n_cases
      call get_case(i, a, lambda, nr, ki, aspect_ratio, shape)
      options%aspect_ratio = aspect_ratio
      options%shape = shape
      call evaluate(work(1), a, lambda, nr, ki, options, serial_result(i), serial_scatmat(i), &
                    serial_amp(:,:,i), serial_amp_status(:,i))
   end do

   !$omp parallel do default(none) private(i,a,lambda,nr,ki,aspect_ratio,shape,options) &
   !$omp shared(work,parallel_result,parallel_scatmat,parallel_amp,parallel_amp_status) schedule(static)
   do i = 1, n_cases
      call get_case(i, a, lambda, nr, ki, aspect_ratio, shape)
      options%tolerance = tolerance
      options%ndgs = ndgs
      options%aspect_ratio = aspect_ratio
      options%shape = shape
      call evaluate(work(omp_get_thread_num()+1), a, lambda, nr, ki, options, parallel_result(i), &
                    parallel_scatmat(i), parallel_amp(:,:,i), parallel_amp_status(:,i))
   end do
   !$omp end parallel do

   failures = 0
   do i = 1, n_cases
      if (parallel_result(i)%status /= serial_result(i)%status .or. &
          parallel_result(i)%legacy_ierr /= serial_result(i)%legacy_ierr .or. &
          parallel_scatmat(i)%lmax /= serial_scatmat(i)%lmax .or. &
          parallel_scatmat(i)%nmax_tm /= serial_scatmat(i)%nmax_tm .or. &
          bits_differ(parallel_scatmat(i)%al1, serial_scatmat(i)%al1) .or. &
          bits_differ(parallel_scatmat(i)%al2, serial_scatmat(i)%al2) .or. &
          bits_differ(parallel_scatmat(i)%al3, serial_scatmat(i)%al3) .or. &
          bits_differ(parallel_scatmat(i)%al4, serial_scatmat(i)%al4) .or. &
          bits_differ(parallel_scatmat(i)%be1, serial_scatmat(i)%be1) .or. &
          bits_differ(parallel_scatmat(i)%be2, serial_scatmat(i)%be2) .or. &
          any(parallel_amp_status(:,i) /= serial_amp_status(:,i)) .or. &
          bits_differ_complex(parallel_amp(:,:,i), serial_amp(:,:,i))) failures = failures + 1
   end do
   do i = 1, nwork
      call tmatrix_workspace_finalize(work(i))
   end do
   if (failures /= 0) error stop 'parallel scattering-matrix or amplitude mismatch'
   write(*,'(a,i0,a)') 'PASS: OpenMP scattering-matrix and amplitude regression (', nwork, &
      ' independent workspaces)'

contains

   subroutine evaluate(work, a, lambda, nr, ki, options, result, scatmat, amp, amp_status)
      !! One random-orientation solve followed by the fixed-orientation
      !! amplitude matrix of the same particle, read from the same workspace.
      type(tmatrix_workspace_t), intent(inout) :: work
      real(real64), intent(in) :: a, lambda, nr, ki
      type(tmatrix_options_t), intent(in) :: options
      type(tmatrix_result_t), intent(out) :: result
      type(tmatrix_scatmat_t), intent(out) :: scatmat
      complex(real64), intent(out) :: amp(4,n_angles)
      integer, intent(out) :: amp_status(n_angles)
      integer :: ia

      amp = (0.0d0, 0.0d0)
      amp_status = 0
      call tmatrix_eval_scattering_matrix(work, a, lambda, nr, ki, options, result, scatmat)
      if (result%status /= 0) return
      do ia = 1, n_angles
         call ampl(work, scatmat, angle_tl(ia), angle_tl1(ia), angle_pl(ia), &
                   angle_pl1(ia), angle_alpha(ia), angle_beta(ia), &
                   amp(1,ia), amp(2,ia), amp(3,ia), amp(4,ia), amp_status(ia))
      end do
   end subroutine evaluate

   logical function bits_differ(actual, expected)
      !! Compare IEEE bit patterns.  The two runs execute the same code on
      !! the same inputs, so anything short of an identical pattern is a
      !! leak between workspaces, not a rounding difference.
      real(real64), intent(in) :: actual(:), expected(:)
      bits_differ = any(transfer(actual, 0_int64, size(actual)) /= &
                        transfer(expected, 0_int64, size(expected)))
   end function bits_differ

   logical function bits_differ_complex(actual, expected)
      complex(real64), intent(in) :: actual(:,:), expected(:,:)
      bits_differ_complex = any(transfer(actual, 0_int64, 2*size(actual)) /= &
                                transfer(expected, 0_int64, 2*size(expected)))
   end function bits_differ_complex

end program test_parallel_scattering_matrix
