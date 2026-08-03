module tmatrix_full_direct_bridge
   !! Final direct-storage bridge.  Every mutable array in the full-direct
   !! fixed-form call chain belongs to the caller's workspace; this module
   !! intentionally contains no serializing OpenMP directive.
   use tmatrix_kinds, only: wp
   use tmatrix_core, only: tmd_one_full_direct
   use tmatrix_types, only: tmatrix_workspace_t, tmatrix_expansion_size
   implicit none
   private
   public :: full_direct_eval, full_direct_eval_scattering_matrix

contains

   subroutine full_direct_eval(work, a_um, lambda_um, nr, ki, aspect_ratio, shape, tolerance, ndgs, &
                               qext, qsca, albedo, asymmetry, legacy_ierr, lapack_info, nmax_tm)
      !! Cross-section-only evaluation.  It returns the multipole truncation
      !! order of the converged T-matrix left in work%full_tstore as well, so
      !! that the workspace can record what it holds; the core computes it
      !! either way and no numerical statement depends on the argument.
      type(tmatrix_workspace_t), intent(inout) :: work
      real(wp), intent(in) :: a_um, lambda_um, nr, ki, aspect_ratio, tolerance
      integer, intent(in) :: shape, ndgs
      real(wp), intent(out) :: qext, qsca, albedo, asymmetry
      integer, intent(out) :: legacy_ierr, lapack_info, nmax_tm

      call tmd_one_full_direct(a_um, lambda_um, nr, ki, aspect_ratio, shape, tolerance, ndgs, &
                               qext, qsca, albedo, asymmetry, legacy_ierr, lapack_info, &
                               work%ct_tr, work%ct_ti, work%ctt_qr, work%ctt_qi, &
                               work%ctt_rgqr, work%ctt_rgqi, work%cbess_j, work%cbess_y, &
                               work%cbess_jr, work%cbess_ji, work%cbess_dj, work%cbess_dy, &
                               work%cbess_djr, work%cbess_dji, work%cbess_gsp_b1r, &
                               work%cbess_gsp_b1i, work%cbess_gsp_b2r, work%cbess_gsp_b2i, &
                               work%full_tstore, work%full_tmatr_work, work%full_gsp_d1, &
                               work%full_gsp_d2, work%full_gsp_d3, work%full_gsp_d4, &
                               work%full_gsp_d5r, work%full_gsp_d5i, work%fac, work%ssign, &
                               work, NMAX_TM_OUT=nmax_tm)
   end subroutine full_direct_eval

   subroutine full_direct_eval_scattering_matrix(work, a_um, lambda_um, nr, ki, aspect_ratio, shape, &
                                                 tolerance, ndgs, qext, qsca, albedo, asymmetry, &
                                                 legacy_ierr, lapack_info, &
                                                 al1, al2, al3, al4, be1, be2, lmax, nmax_tm)
      !! Same evaluation as full_direct_eval, additionally returning the
      !! generalized-spherical-function expansion of the random-orientation
      !! scattering matrix and the multipole truncation order of the
      !! converged T-matrix left in work%full_tstore.
      type(tmatrix_workspace_t), intent(inout) :: work
      real(wp), intent(in) :: a_um, lambda_um, nr, ki, aspect_ratio, tolerance
      integer, intent(in) :: shape, ndgs
      real(wp), intent(out) :: qext, qsca, albedo, asymmetry
      integer, intent(out) :: legacy_ierr, lapack_info
      !! Explicit shape: the core declares these dummies as NPL-element
      !! arrays, so matching the extent here avoids a copy-in/copy-out
      !! temporary.
      real(wp), intent(out) :: al1(tmatrix_expansion_size), al2(tmatrix_expansion_size), &
                               al3(tmatrix_expansion_size), al4(tmatrix_expansion_size), &
                               be1(tmatrix_expansion_size), be2(tmatrix_expansion_size)
      integer, intent(out) :: lmax, nmax_tm

      call tmd_one_full_direct(a_um, lambda_um, nr, ki, aspect_ratio, shape, tolerance, ndgs, &
                               qext, qsca, albedo, asymmetry, legacy_ierr, lapack_info, &
                               work%ct_tr, work%ct_ti, work%ctt_qr, work%ctt_qi, &
                               work%ctt_rgqr, work%ctt_rgqi, work%cbess_j, work%cbess_y, &
                               work%cbess_jr, work%cbess_ji, work%cbess_dj, work%cbess_dy, &
                               work%cbess_djr, work%cbess_dji, work%cbess_gsp_b1r, &
                               work%cbess_gsp_b1i, work%cbess_gsp_b2r, work%cbess_gsp_b2i, &
                               work%full_tstore, work%full_tmatr_work, work%full_gsp_d1, &
                               work%full_gsp_d2, work%full_gsp_d3, work%full_gsp_d4, &
                               work%full_gsp_d5r, work%full_gsp_d5i, work%fac, work%ssign, &
                               work, al1, al2, al3, al4, be1, be2, lmax, nmax_tm)
   end subroutine full_direct_eval_scattering_matrix

end module tmatrix_full_direct_bridge
