module tmatrix_api
   !! SEDust-independent public API for one randomly oriented particle.
   !! `tmatrix_eval` is the reentrant full-direct backend: different caller-
   !! owned workspaces may be evaluated concurrently.  The shared-state
   !! compatibility paths deliberately live in tmatrix_regression_api and are
   !! not part of this production module or archive.
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use, intrinsic :: iso_fortran_env, only: int64
   use tmatrix_kinds, only: wp
   use tmatrix_status, only: TMATRIX_SUCCESS, TMATRIX_ERR_INVALID_INPUT, &
                             TMATRIX_ERR_NOT_INITIALIZED, TMATRIX_ERR_LEGACY_CORE, &
                             TMATRIX_ERR_ALLOCATION, TMATRIX_ERR_LAPACK_FACTOR, TMATRIX_ERR_LAPACK_INVERT, &
                             TMATRIX_ERR_NUMERICAL_NONFINITE
   use tmatrix_types, only: tmatrix_options_t, tmatrix_result_t, tmatrix_workspace_t, &
                            tmatrix_scatmat_t, tmatrix_expansion_size
   use tmatrix_full_direct_bridge, only: full_direct_eval, full_direct_eval_scattering_matrix
   use fixed_orientation_amplitude, only: ampl
   use fixed_orientation_phase_matrix, only: phase_matrix_from_amplitude
   use tmatrix_core, only: gauss => gauss_full_direct
   implicit none
   private

   public :: wp
   public :: tmatrix_options_t, tmatrix_result_t, tmatrix_workspace_t
   public :: tmatrix_scatmat_t, tmatrix_expansion_size
   public :: tmatrix_workspace_init, tmatrix_workspace_finalize, tmatrix_workspace_storage_bytes
   public :: tmatrix_eval, tmatrix_eval_scattering_matrix
   !! Fixed-orientation amplitude matrix.  Re-exported here because its
   !! workspace contract is part of this API: it reads the T-matrix that
   !! tmatrix_eval_scattering_matrix left in the same workspace.
   public :: ampl
   !! Phase (Mueller) matrix built from that amplitude matrix, Mishchenko's
   !! Z11...Z44 bilinears.  Re-exported next to AMPL because the two are used
   !! as a pair; it reads no workspace of its own.
   public :: phase_matrix_from_amplitude
   !! Gauss-Legendre nodes and weights, GAUSS(N, IND1, IND2, Z, W): nodes on
   !! (-1,1) with weights summing to 2 for IND1 = 0, on (0,1) for IND1 = 1.
   !! This is the same rule the T-matrix core integrates its own surface and
   !! angular quadratures with, so a caller that integrates the scattering
   !! matrix this API returns lands on identical nodes.  Mishchenko's original
   !! name is kept; internally the core calls it GAUSS_FULL_DIRECT, after the
   !! full-direct workspace variant it belongs to.
   public :: gauss

contains

   subroutine tmatrix_workspace_init(work, status, message)
      type(tmatrix_workspace_t), intent(inout) :: work
      integer, intent(out), optional :: status
      character(len=*), intent(out), optional :: message
      integer :: alloc_status

      call tmatrix_workspace_finalize(work)
      allocate(work%ct_tr(200,200), work%ct_ti(200,200), work%ctt_qr(200,200), &
               work%ctt_qi(200,200), work%ctt_rgqr(200,200), work%ctt_rgqi(200,200), stat=alloc_status)
      if (alloc_status /= 0) then
         call tmatrix_workspace_finalize(work)
         if (present(status)) status = TMATRIX_ERR_ALLOCATION
         if (present(message)) message = 'unable to allocate CT/CTT workspace components'
         return
      end if
      allocate(work%cbess_j(600,100), work%cbess_y(600,100), work%cbess_jr(600,100), &
               work%cbess_ji(600,100), work%cbess_dj(600,100), work%cbess_dy(600,100), &
               work%cbess_djr(600,100), work%cbess_dji(600,100), &
               work%cbess_gsp_b1r(161,161,80), work%cbess_gsp_b1i(161,161,80), &
               work%cbess_gsp_b2r(161,161,80), work%cbess_gsp_b2i(161,161,80), stat=alloc_status)
      if (alloc_status /= 0) then
         call tmatrix_workspace_finalize(work)
         if (present(status)) status = TMATRIX_ERR_ALLOCATION
         if (present(message)) message = 'unable to allocate CBESS workspace components'
         return
      end if
      work%initialized = .true.
      call forget_stored_tmatrix(work)
      work%tm_tag = 0_int64
      if (present(status)) status = TMATRIX_SUCCESS
      if (present(message)) message = ''
   end subroutine tmatrix_workspace_init

   subroutine tmatrix_workspace_finalize(work)
      type(tmatrix_workspace_t), intent(inout) :: work
      if (allocated(work%ct_tr)) deallocate(work%ct_tr)
      if (allocated(work%ct_ti)) deallocate(work%ct_ti)
      if (allocated(work%ctt_qr)) deallocate(work%ctt_qr)
      if (allocated(work%ctt_qi)) deallocate(work%ctt_qi)
      if (allocated(work%ctt_rgqr)) deallocate(work%ctt_rgqr)
      if (allocated(work%ctt_rgqi)) deallocate(work%ctt_rgqi)
      if (allocated(work%cbess_j)) deallocate(work%cbess_j)
      if (allocated(work%cbess_y)) deallocate(work%cbess_y)
      if (allocated(work%cbess_jr)) deallocate(work%cbess_jr)
      if (allocated(work%cbess_ji)) deallocate(work%cbess_ji)
      if (allocated(work%cbess_dj)) deallocate(work%cbess_dj)
      if (allocated(work%cbess_dy)) deallocate(work%cbess_dy)
      if (allocated(work%cbess_djr)) deallocate(work%cbess_djr)
      if (allocated(work%cbess_dji)) deallocate(work%cbess_dji)
      if (allocated(work%cbess_gsp_b1r)) deallocate(work%cbess_gsp_b1r)
      if (allocated(work%cbess_gsp_b1i)) deallocate(work%cbess_gsp_b1i)
      if (allocated(work%cbess_gsp_b2r)) deallocate(work%cbess_gsp_b2r)
      if (allocated(work%cbess_gsp_b2i)) deallocate(work%cbess_gsp_b2i)
      if (allocated(work%full_tstore)) deallocate(work%full_tstore)
      if (allocated(work%full_tmatr_work)) deallocate(work%full_tmatr_work)
      if (allocated(work%full_gsp_d1)) deallocate(work%full_gsp_d1)
      if (allocated(work%full_gsp_d2)) deallocate(work%full_gsp_d2)
      if (allocated(work%full_gsp_d3)) deallocate(work%full_gsp_d3)
      if (allocated(work%full_gsp_d4)) deallocate(work%full_gsp_d4)
      if (allocated(work%full_gsp_d5r)) deallocate(work%full_gsp_d5r)
      if (allocated(work%full_gsp_d5i)) deallocate(work%full_gsp_d5i)
      call deallocate_tt_tmatr_gsp_scratch(work)
      work%initialized = .false.
      call forget_stored_tmatrix(work)
      work%tm_tag = 0_int64
   end subroutine tmatrix_workspace_finalize

   subroutine forget_stored_tmatrix(work)
      !! Mark the workspace as holding no usable T-matrix.  The generation tag
      !! is deliberately untouched here: it must keep advancing so that a
      !! scatmat obtained earlier can never match again.
      type(tmatrix_workspace_t), intent(inout) :: work
      work%tm_filled    = .false.
      work%tm_lambda_um = 0.0_wp
      work%tm_nmax      = 0
   end subroutine forget_stored_tmatrix

   subroutine record_stored_tmatrix(work, status, lambda_um, nmax_tm)
      !! Stamp what an evaluation just left in work%full_tstore.  The tag
      !! advances on every evaluation, successful or not, because a failed one
      !! may already have overwritten part of the storage: any scatmat obtained
      !! before this call therefore refers to a T-matrix that is gone.
      type(tmatrix_workspace_t), intent(inout) :: work
      integer, intent(in) :: status, nmax_tm
      real(wp), intent(in) :: lambda_um

      work%tm_tag = work%tm_tag + 1_int64
      if (status == TMATRIX_SUCCESS) then
         work%tm_filled    = .true.
         work%tm_lambda_um = lambda_um
         work%tm_nmax      = nmax_tm
      else
         call forget_stored_tmatrix(work)
      end if
   end subroutine record_stored_tmatrix

   subroutine deallocate_tt_tmatr_gsp_scratch(work)
      !! Deallocate the TT/TMATR/GSP scratch arrays that used to be automatic
      !! arrays on the stack of the fixed-form routines.
      type(tmatrix_workspace_t), intent(inout) :: work
      if (allocated(work%tt_zq)) deallocate(work%tt_zq)
      if (allocated(work%tmatr_d1)) deallocate(work%tmatr_d1)
      if (allocated(work%tmatr_d2)) deallocate(work%tmatr_d2)
      if (allocated(work%gsp_tr1)) deallocate(work%gsp_tr1)
      if (allocated(work%gsp_tr2)) deallocate(work%gsp_tr2)
      if (allocated(work%gsp_ti1)) deallocate(work%gsp_ti1)
      if (allocated(work%gsp_ti2)) deallocate(work%gsp_ti2)
      if (allocated(work%gsp_g1)) deallocate(work%gsp_g1)
      if (allocated(work%gsp_g2)) deallocate(work%gsp_g2)
      if (allocated(work%gsp_fr)) deallocate(work%gsp_fr)
      if (allocated(work%gsp_fi)) deallocate(work%gsp_fi)
      if (allocated(work%gsp_ff)) deallocate(work%gsp_ff)
   end subroutine deallocate_tt_tmatr_gsp_scratch

   subroutine tmatrix_workspace_prepare_full_direct(work, status, message)
      !! Allocate final TMAT/GSP storage lazily for the production full-direct
      !! backend.
      type(tmatrix_workspace_t), intent(inout) :: work
      integer, intent(out), optional :: status
      character(len=*), intent(out), optional :: message
      integer :: alloc_status

      if (.not. work%initialized) then
         if (present(status)) status = TMATRIX_ERR_NOT_INITIALIZED
         if (present(message)) message = 'workspace has not been initialized'
         return
      end if
      if (allocated(work%full_tstore)) then
         if (present(status)) status = TMATRIX_SUCCESS
         if (present(message)) message = ''
         return
      end if
      allocate(work%full_tstore(81,80,80,8), work%full_tmatr_work(100,100,16), &
               work%full_gsp_d1(161,80,80), work%full_gsp_d2(161,80,80), &
               work%full_gsp_d3(161,80,80), work%full_gsp_d4(161,80,80), &
               work%full_gsp_d5r(161,80,80), work%full_gsp_d5i(161,80,80), stat=alloc_status)
      if (alloc_status /= 0) then
         if (allocated(work%full_tstore)) deallocate(work%full_tstore)
         if (allocated(work%full_tmatr_work)) deallocate(work%full_tmatr_work)
         if (allocated(work%full_gsp_d1)) deallocate(work%full_gsp_d1)
         if (allocated(work%full_gsp_d2)) deallocate(work%full_gsp_d2)
         if (allocated(work%full_gsp_d3)) deallocate(work%full_gsp_d3)
         if (allocated(work%full_gsp_d4)) deallocate(work%full_gsp_d4)
         if (allocated(work%full_gsp_d5r)) deallocate(work%full_gsp_d5r)
         if (allocated(work%full_gsp_d5i)) deallocate(work%full_gsp_d5i)
         if (present(status)) status = TMATRIX_ERR_ALLOCATION
         if (present(message)) message = 'unable to allocate full-direct TMAT/GSP workspace components'
         return
      end if
      ! Extents repeat the fixed-form parameters NPN2=200, NPNG2=600,
      ! NPN1=100, NPL1=161, NPN4=80, and NPN6=81 of the numerical core.
      allocate(work%tt_zq(200,200), &
               work%tmatr_d1(600,100), work%tmatr_d2(600,100), &
               work%gsp_tr1(161,80), work%gsp_tr2(161,80), &
               work%gsp_ti1(161,80), work%gsp_ti2(161,80), &
               work%gsp_g1(161,81), work%gsp_g2(161,81), &
               work%gsp_fr(80,80), work%gsp_fi(80,80), work%gsp_ff(80,80), stat=alloc_status)
      if (alloc_status /= 0) then
         call deallocate_tt_tmatr_gsp_scratch(work)
         ! full_tstore is the allocation flag read on entry; drop the whole
         ! full-direct set so a later call retries instead of reporting success.
         deallocate(work%full_tstore, work%full_tmatr_work, work%full_gsp_d1, &
                    work%full_gsp_d2, work%full_gsp_d3, work%full_gsp_d4, &
                    work%full_gsp_d5r, work%full_gsp_d5i)
         if (present(status)) status = TMATRIX_ERR_ALLOCATION
         if (present(message)) message = 'unable to allocate full-direct TT/TMATR/GSP workspace components'
         return
      end if
      if (present(status)) status = TMATRIX_SUCCESS
      if (present(message)) message = ''
   end subroutine tmatrix_workspace_prepare_full_direct

   function tmatrix_workspace_storage_bytes(work) result(bytes)
      type(tmatrix_workspace_t), intent(in) :: work
      integer(int64) :: bytes

      bytes = int(size(work%fac), int64) * storage_size(work%fac)/8_int64 + &
              int(size(work%ssign), int64) * storage_size(work%ssign)/8_int64
      if (allocated(work%ct_tr)) bytes = bytes + int(size(work%ct_tr), int64) * storage_size(work%ct_tr)/8_int64
      if (allocated(work%ct_ti)) bytes = bytes + int(size(work%ct_ti), int64) * storage_size(work%ct_ti)/8_int64
      if (allocated(work%ctt_qr)) bytes = bytes + int(size(work%ctt_qr), int64) * storage_size(work%ctt_qr)/8_int64
      if (allocated(work%ctt_qi)) bytes = bytes + int(size(work%ctt_qi), int64) * storage_size(work%ctt_qi)/8_int64
      if (allocated(work%ctt_rgqr)) bytes = bytes + int(size(work%ctt_rgqr), int64) * storage_size(work%ctt_rgqr)/8_int64
      if (allocated(work%ctt_rgqi)) bytes = bytes + int(size(work%ctt_rgqi), int64) * storage_size(work%ctt_rgqi)/8_int64
      if (allocated(work%cbess_j)) bytes = bytes + int(size(work%cbess_j), int64) * storage_size(work%cbess_j)/8_int64
      if (allocated(work%cbess_y)) bytes = bytes + int(size(work%cbess_y), int64) * storage_size(work%cbess_y)/8_int64
      if (allocated(work%cbess_jr)) bytes = bytes + int(size(work%cbess_jr), int64) * storage_size(work%cbess_jr)/8_int64
      if (allocated(work%cbess_ji)) bytes = bytes + int(size(work%cbess_ji), int64) * storage_size(work%cbess_ji)/8_int64
      if (allocated(work%cbess_dj)) bytes = bytes + int(size(work%cbess_dj), int64) * storage_size(work%cbess_dj)/8_int64
      if (allocated(work%cbess_dy)) bytes = bytes + int(size(work%cbess_dy), int64) * storage_size(work%cbess_dy)/8_int64
      if (allocated(work%cbess_djr)) bytes = bytes + int(size(work%cbess_djr), int64) * storage_size(work%cbess_djr)/8_int64
      if (allocated(work%cbess_dji)) bytes = bytes + int(size(work%cbess_dji), int64) * storage_size(work%cbess_dji)/8_int64
      if (allocated(work%cbess_gsp_b1r)) bytes = bytes + int(size(work%cbess_gsp_b1r), int64) * storage_size(work%cbess_gsp_b1r)/8_int64
      if (allocated(work%cbess_gsp_b1i)) bytes = bytes + int(size(work%cbess_gsp_b1i), int64) * storage_size(work%cbess_gsp_b1i)/8_int64
      if (allocated(work%cbess_gsp_b2r)) bytes = bytes + int(size(work%cbess_gsp_b2r), int64) * storage_size(work%cbess_gsp_b2r)/8_int64
      if (allocated(work%cbess_gsp_b2i)) bytes = bytes + int(size(work%cbess_gsp_b2i), int64) * storage_size(work%cbess_gsp_b2i)/8_int64
      if (allocated(work%full_tstore)) bytes = bytes + int(size(work%full_tstore), int64) * storage_size(work%full_tstore)/8_int64
      if (allocated(work%full_tmatr_work)) bytes = bytes + int(size(work%full_tmatr_work), int64) * storage_size(work%full_tmatr_work)/8_int64
      if (allocated(work%full_gsp_d1)) bytes = bytes + int(size(work%full_gsp_d1), int64) * storage_size(work%full_gsp_d1)/8_int64
      if (allocated(work%full_gsp_d2)) bytes = bytes + int(size(work%full_gsp_d2), int64) * storage_size(work%full_gsp_d2)/8_int64
      if (allocated(work%full_gsp_d3)) bytes = bytes + int(size(work%full_gsp_d3), int64) * storage_size(work%full_gsp_d3)/8_int64
      if (allocated(work%full_gsp_d4)) bytes = bytes + int(size(work%full_gsp_d4), int64) * storage_size(work%full_gsp_d4)/8_int64
      if (allocated(work%full_gsp_d5r)) bytes = bytes + int(size(work%full_gsp_d5r), int64) * storage_size(work%full_gsp_d5r)/8_int64
      if (allocated(work%full_gsp_d5i)) bytes = bytes + int(size(work%full_gsp_d5i), int64) * storage_size(work%full_gsp_d5i)/8_int64
      if (allocated(work%tt_zq)) bytes = bytes + int(size(work%tt_zq), int64) * storage_size(work%tt_zq)/8_int64
      if (allocated(work%tmatr_d1)) bytes = bytes + int(size(work%tmatr_d1), int64) * storage_size(work%tmatr_d1)/8_int64
      if (allocated(work%tmatr_d2)) bytes = bytes + int(size(work%tmatr_d2), int64) * storage_size(work%tmatr_d2)/8_int64
      if (allocated(work%gsp_tr1)) bytes = bytes + int(size(work%gsp_tr1), int64) * storage_size(work%gsp_tr1)/8_int64
      if (allocated(work%gsp_tr2)) bytes = bytes + int(size(work%gsp_tr2), int64) * storage_size(work%gsp_tr2)/8_int64
      if (allocated(work%gsp_ti1)) bytes = bytes + int(size(work%gsp_ti1), int64) * storage_size(work%gsp_ti1)/8_int64
      if (allocated(work%gsp_ti2)) bytes = bytes + int(size(work%gsp_ti2), int64) * storage_size(work%gsp_ti2)/8_int64
      if (allocated(work%gsp_g1)) bytes = bytes + int(size(work%gsp_g1), int64) * storage_size(work%gsp_g1)/8_int64
      if (allocated(work%gsp_g2)) bytes = bytes + int(size(work%gsp_g2), int64) * storage_size(work%gsp_g2)/8_int64
      if (allocated(work%gsp_fr)) bytes = bytes + int(size(work%gsp_fr), int64) * storage_size(work%gsp_fr)/8_int64
      if (allocated(work%gsp_fi)) bytes = bytes + int(size(work%gsp_fi), int64) * storage_size(work%gsp_fi)/8_int64
      if (allocated(work%gsp_ff)) bytes = bytes + int(size(work%gsp_ff), int64) * storage_size(work%gsp_ff)/8_int64
   end function tmatrix_workspace_storage_bytes

   subroutine tmatrix_eval(work, a_um, lambda_um, nr, ki, options, result)
      !! Default reentrant backend.  A workspace belongs to exactly one active
      !! evaluation; distinct initialized workspaces may be evaluated in
      !! parallel without serialization.
      type(tmatrix_workspace_t), intent(inout) :: work
      real(wp), intent(in) :: a_um, lambda_um, nr, ki
      type(tmatrix_options_t), intent(in) :: options
      type(tmatrix_result_t), intent(out) :: result
      integer :: nmax_tm

      nmax_tm = 0
      if (workspace_is_ready(work, a_um, lambda_um, nr, ki, options, result)) then
         call full_direct_eval(work, a_um, lambda_um, nr, ki, options%aspect_ratio, options%shape, &
                               options%tolerance, options%ndgs, result%qext, result%qsca, &
                               result%albedo, result%asymmetry, result%legacy_ierr, &
                               result%lapack_info, nmax_tm)
         call derive_qabs_and_status(result)
      end if
      call record_stored_tmatrix(work, result%status, lambda_um, nmax_tm)
   end subroutine tmatrix_eval

   subroutine tmatrix_eval_scattering_matrix(work, a_um, lambda_um, nr, ki, options, result, scatmat)
      !! Same evaluation and same concurrency contract as `tmatrix_eval`,
      !! additionally returning the generalized-spherical-function expansion
      !! of the random-orientation scattering matrix and the multipole
      !! truncation order of the converged T-matrix.
      !!
      !! On success the workspace holds that T-matrix and the returned
      !! `scatmat` identifies it, so `ampl` may be called on the same
      !! workspace with that `scatmat` to get the fixed-orientation amplitude
      !! matrix of the same particle.  Any further evaluation on that
      !! workspace replaces the T-matrix and retires the `scatmat` with it.
      type(tmatrix_workspace_t), intent(inout) :: work
      real(wp), intent(in) :: a_um, lambda_um, nr, ki
      type(tmatrix_options_t), intent(in) :: options
      type(tmatrix_result_t), intent(out) :: result
      type(tmatrix_scatmat_t), intent(out) :: scatmat

      scatmat = tmatrix_scatmat_t()
      if (workspace_is_ready(work, a_um, lambda_um, nr, ki, options, result)) then
         call full_direct_eval_scattering_matrix(work, a_um, lambda_um, nr, ki, options%aspect_ratio, &
                            options%shape, options%tolerance, options%ndgs, result%qext, result%qsca, &
                            result%albedo, result%asymmetry, result%legacy_ierr, result%lapack_info, &
                            scatmat%al1, scatmat%al2, scatmat%al3, scatmat%al4, &
                            scatmat%be1, scatmat%be2, scatmat%lmax, scatmat%nmax_tm)
         call derive_qabs_and_status(result)
      end if
      call record_stored_tmatrix(work, result%status, lambda_um, scatmat%nmax_tm)
      if (result%status /= TMATRIX_SUCCESS) then
         scatmat = tmatrix_scatmat_t()
      else
         scatmat%state_tag = work%tm_tag
      end if
   end subroutine tmatrix_eval_scattering_matrix

   logical function workspace_is_ready(work, a_um, lambda_um, nr, ki, options, result)
      !! Initialize the result, screen the inputs, and make sure the
      !! full-direct components of the workspace exist.  A false return
      !! leaves the failure status and message in `result`.
      type(tmatrix_workspace_t), intent(inout) :: work
      real(wp), intent(in) :: a_um, lambda_um, nr, ki
      type(tmatrix_options_t), intent(in) :: options
      type(tmatrix_result_t), intent(out) :: result
      integer :: alloc_status
      character(len=160) :: alloc_message

      result = tmatrix_result_t()
      workspace_is_ready = .false.
      if (.not. work%initialized) then
         call fail(result, TMATRIX_ERR_NOT_INITIALIZED, 'workspace has not been initialized')
         return
      end if
      if (.not. valid_input(a_um, lambda_um, nr, ki, options)) then
         call fail(result, TMATRIX_ERR_INVALID_INPUT, 'invalid size, wavelength, refractive index, or options')
         return
      end if
      call tmatrix_workspace_prepare_full_direct(work, alloc_status, alloc_message)
      if (alloc_status /= TMATRIX_SUCCESS) then
         call fail(result, alloc_status, alloc_message)
         return
      end if
      workspace_is_ready = .true.
   end function workspace_is_ready

   subroutine derive_qabs_and_status(result)
      !! Derive Q_abs and map the numerical core's IERR and LAPACK status
      !! onto the public status code.
      type(tmatrix_result_t), intent(inout) :: result

      result%qabs = result%qext - result%qsca
      if (result%legacy_ierr /= 0) then
         if (result%legacy_ierr == 8) then
            call fail(result, TMATRIX_ERR_LAPACK_FACTOR, 'full-direct LAPACK factorization failed')
         else if (result%legacy_ierr == 9) then
            call fail(result, TMATRIX_ERR_LAPACK_INVERT, 'full-direct LAPACK inversion failed')
         else if (result%lapack_info /= 0) then
            call fail(result, TMATRIX_ERR_LEGACY_CORE, 'full-direct LAPACK call failed')
         else
            call fail(result, TMATRIX_ERR_LEGACY_CORE, 'full-direct backend returned a nonzero IERR')
         end if
         return
      end if
      if (.not. ieee_is_finite(result%qext) .or. .not. ieee_is_finite(result%qsca) .or. &
          .not. ieee_is_finite(result%qabs) .or. .not. ieee_is_finite(result%albedo) .or. &
          .not. ieee_is_finite(result%asymmetry)) then
         call fail(result, TMATRIX_ERR_NUMERICAL_NONFINITE, 'full-direct backend returned NaN or Inf')
         return
      end if
      result%status = TMATRIX_SUCCESS
      result%message = ''
   end subroutine derive_qabs_and_status

   logical function valid_input(a_um, lambda_um, nr, ki, options)
      real(wp), intent(in) :: a_um, lambda_um, nr, ki
      type(tmatrix_options_t), intent(in) :: options

      valid_input = ieee_is_finite(a_um) .and. ieee_is_finite(lambda_um) .and. &
                    ieee_is_finite(nr) .and. ieee_is_finite(ki) .and. &
                    ieee_is_finite(options%aspect_ratio) .and. ieee_is_finite(options%tolerance)
      valid_input = valid_input .and. a_um > 0.0_wp .and. lambda_um > 0.0_wp .and. &
                    nr > 0.0_wp .and. ki >= 0.0_wp .and. options%aspect_ratio > 0.0_wp .and. &
                    options%tolerance > 0.0_wp .and. options%ndgs > 0
      valid_input = valid_input .and. (options%shape == -1 .or. options%shape == -2 .or. options%shape >= 0)
   end function valid_input

   subroutine fail(result, status, message)
      type(tmatrix_result_t), intent(inout) :: result
      integer, intent(in) :: status
      character(len=*), intent(in) :: message

      result%status = status
      result%message = message
   end subroutine fail

end module tmatrix_api
