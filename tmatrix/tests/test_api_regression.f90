program test_api_regression
   use, intrinsic :: iso_fortran_env, only: int64, real64, error_unit
   use tmatrix_api, only: tmatrix_options_t, tmatrix_result_t, tmatrix_workspace_t, &
                          tmatrix_workspace_init, tmatrix_workspace_finalize, &
                          tmatrix_workspace_storage_bytes, tmatrix_eval
   use tmatrix_test_cases, only: n_cases, tolerance, ndgs, get_case
   implicit none

   character(len=*), parameter :: golden_path = 'tests/golden/legacy_cases.dat'
   integer(int64), parameter :: expected_default_workspace_bytes = 83926096_int64
   type(tmatrix_workspace_t) :: work
   type(tmatrix_options_t) :: options
   type(tmatrix_result_t) :: result
   integer :: i, ios, expected_ierr, expected_lapack_info, case_id, failed
   real(real64) :: a, lambda, nr, ki, aspect_ratio
   integer :: shape
   real(real64) :: expected_qext, expected_qsca, expected_walb, expected_asymm

   options%tolerance = tolerance
   options%ndgs = ndgs
   call tmatrix_workspace_init(work)
   open(unit=10, file=golden_path, status='old', action='read', iostat=ios)
   if (ios /= 0) error stop 'cannot open '//golden_path
   failed = 0
   do i = 1, n_cases
      read(10,*,iostat=ios) case_id, expected_qext, expected_qsca, expected_walb, expected_asymm, expected_ierr, expected_lapack_info
      if (ios /= 0 .or. case_id /= i) error stop 'malformed legacy golden file'
      call get_case(i, a, lambda, nr, ki, aspect_ratio, shape)
      options%aspect_ratio = aspect_ratio
      options%shape = shape
      call tmatrix_eval(work, a, lambda, nr, ki, options, result)
      if (i == 1 .and. tmatrix_workspace_storage_bytes(work) /= expected_default_workspace_bytes) then
         failed = failed + 1
         write(error_unit,'(a,i0,a,i0)') 'default workspace bytes=', tmatrix_workspace_storage_bytes(work), &
            ', expected=', expected_default_workspace_bytes
      end if
      if (result%status /= 0 .or. result%legacy_ierr /= expected_ierr .or. &
          result%lapack_info /= expected_lapack_info .or. &
          .not. same(result%qext, expected_qext) .or. .not. same(result%qsca, expected_qsca) .or. &
          .not. same(result%albedo, expected_walb) .or. .not. same(result%asymmetry, expected_asymm)) then
         failed = failed + 1
         write(error_unit,'(a,i0,a,i0,a,a)') 'API regression failed for case ', i, &
            ', status=', result%status, ': ', trim(result%message)
      end if
   end do
   close(10)
   call tmatrix_workspace_finalize(work)
   if (failed /= 0) error stop 'API regression failures'
   write(*,'(a,i0,a)') 'PASS: public API golden regression (', n_cases, ' cases)'

contains

   logical function same(actual, expected)
      real(real64), intent(in) :: actual, expected
      same = abs(actual-expected) <= max(1.0e-12_real64, 1.0e-10_real64*abs(expected))
   end function same

end program test_api_regression
