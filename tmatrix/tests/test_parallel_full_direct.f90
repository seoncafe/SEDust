program test_parallel_full_direct
   !! Different caller-owned workspaces must produce the serial full-direct
   !! reference concurrently.  Run this executable with 1, 2, 4, and 8
   !! OpenMP threads; it intentionally contains no serialized backend call.
   use, intrinsic :: iso_fortran_env, only: real64
   use omp_lib, only: omp_get_thread_num, omp_get_max_threads
   use tmatrix_api, only: tmatrix_options_t, tmatrix_result_t, tmatrix_workspace_t, &
                          tmatrix_workspace_init, tmatrix_workspace_finalize, tmatrix_eval
   use tmatrix_test_cases, only: n_cases, tolerance, ndgs, get_case
   implicit none

   type(tmatrix_workspace_t), allocatable :: work(:)
   type(tmatrix_options_t) :: options
   type(tmatrix_result_t) :: serial_result(n_cases), parallel_result(n_cases)
   integer :: i, nwork, failures
   real(real64) :: a, lambda, nr, ki, aspect_ratio
   integer :: shape

   options%tolerance = tolerance
   options%ndgs = ndgs
   nwork = max(1, omp_get_max_threads())
   allocate(work(nwork))
   do i = 1, nwork
      call tmatrix_workspace_init(work(i))
      if (.not. work(i)%initialized) error stop 'full-direct workspace allocation failed'
   end do

   do i = 1, n_cases
      call get_case(i, a, lambda, nr, ki, aspect_ratio, shape)
      options%aspect_ratio = aspect_ratio
      options%shape = shape
      call tmatrix_eval(work(1), a, lambda, nr, ki, options, serial_result(i))
   end do

   !$omp parallel do default(none) private(i,a,lambda,nr,ki,aspect_ratio,shape,options) &
   !$omp shared(work,parallel_result) schedule(static)
   do i = 1, n_cases
      call get_case(i, a, lambda, nr, ki, aspect_ratio, shape)
      options%tolerance = tolerance
      options%ndgs = ndgs
      options%aspect_ratio = aspect_ratio
      options%shape = shape
      call tmatrix_eval(work(omp_get_thread_num()+1), a, lambda, nr, ki, options, parallel_result(i))
   end do
   !$omp end parallel do

   failures = 0
   do i = 1, n_cases
      if (parallel_result(i)%status /= serial_result(i)%status .or. &
          parallel_result(i)%legacy_ierr /= serial_result(i)%legacy_ierr .or. &
          parallel_result(i)%lapack_info /= serial_result(i)%lapack_info .or. &
          .not. same(parallel_result(i)%qext, serial_result(i)%qext) .or. &
          .not. same(parallel_result(i)%qsca, serial_result(i)%qsca) .or. &
          .not. same(parallel_result(i)%albedo, serial_result(i)%albedo) .or. &
          .not. same(parallel_result(i)%asymmetry, serial_result(i)%asymmetry)) failures = failures + 1
   end do
   do i = 1, nwork
      call tmatrix_workspace_finalize(work(i))
   end do
   if (failures /= 0) error stop 'parallel full-direct API mismatch'
   write(*,'(a,i0,a)') 'PASS: OpenMP full-direct regression (', nwork, ' independent workspaces)'

contains

   logical function same(actual, expected)
      real(real64), intent(in) :: actual, expected
      same = abs(actual-expected) <= max(1.0e-12_real64, 1.0e-10_real64*abs(expected))
   end function same

end program test_parallel_full_direct
