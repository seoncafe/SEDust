program evaluate_one
   !! Minimal external-consumer example: the default reentrant full-direct API
   !! needs no SEDust module or input file.
   use tmatrix_api, only: wp, tmatrix_options_t, tmatrix_result_t, &
                          tmatrix_workspace_t, tmatrix_workspace_init, &
                          tmatrix_workspace_finalize, tmatrix_eval
   implicit none

   type(tmatrix_workspace_t) :: work
   type(tmatrix_options_t) :: options
   type(tmatrix_result_t) :: result

   call tmatrix_workspace_init(work)
   call tmatrix_eval(work, 0.10_wp, 0.55_wp, 1.70_wp, 0.02_wp, options, result)
   if (result%status /= 0) error stop trim(result%message)
   write(*,'(a,4(1x,es14.6))') 'Qext Qsca Qabs g:', result%qext, result%qsca, &
                                result%qabs, result%asymmetry
   call tmatrix_workspace_finalize(work)
end program evaluate_one
