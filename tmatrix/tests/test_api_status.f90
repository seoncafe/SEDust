program test_api_status
   use tmatrix_api, only: wp, tmatrix_options_t, tmatrix_result_t, tmatrix_workspace_t, &
                          tmatrix_workspace_init, tmatrix_workspace_finalize, tmatrix_eval
   use tmatrix_status, only: TMATRIX_ERR_INVALID_INPUT, TMATRIX_ERR_NOT_INITIALIZED
   implicit none

   type(tmatrix_workspace_t) :: work
   type(tmatrix_options_t) :: options
   type(tmatrix_result_t) :: result

   call tmatrix_eval(work, 0.1_wp, 0.55_wp, 1.7_wp, 0.02_wp, options, result)
   if (result%status /= TMATRIX_ERR_NOT_INITIALIZED) error stop 'uninitialized workspace status failed'

   call tmatrix_workspace_init(work)
   call tmatrix_eval(work, 0.1_wp, 0.0_wp, 1.7_wp, 0.02_wp, options, result)
   if (result%status /= TMATRIX_ERR_INVALID_INPUT) error stop 'zero wavelength status failed'
   call tmatrix_eval(work, 0.1_wp, 0.55_wp, 1.7_wp, -0.02_wp, options, result)
   if (result%status /= TMATRIX_ERR_INVALID_INPUT) error stop 'negative absorption status failed'
   options%ndgs = 0
   call tmatrix_eval(work, 0.1_wp, 0.55_wp, 1.7_wp, 0.02_wp, options, result)
   if (result%status /= TMATRIX_ERR_INVALID_INPUT) error stop 'invalid options status failed'
   call tmatrix_workspace_finalize(work)
   write(*,'(a)') 'PASS: public API invalid-input and lifecycle status checks'
end program test_api_status
