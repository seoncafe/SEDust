module dust_lib
   ! RT-facing entry point for the model-agnostic dust thermal-emission
   ! library. A Fortran 3D radiative-transfer code links this module and:
   !
   !   use dust_lib
   !   type(dust_model_t) :: m
   !   call build_astrodust(m, qtab, sizedist, NT, T_lo, T_hi)   ! once
   !   ...
   !   do icell = 1, ncells
   !      ! ... assemble local mean intensity J_lam(:) on m%lam ...
   !      call dust_emission(m, J_lam, lamI_total [, lamI_chan])  ! per cell
   !   end do
   !
   ! Two usage modes:
   !   (a) single-cell EXACT solve: dust_emission(m, J_lam, ...)  -- arbitrary J(lambda).
   !   (b) precomputed TABLE + interpolation, for when the cell field is a
   !       fixed reference SHAPE scaled by an intensity U:
   !         call dust_build_table(m, J_ref, U_grid, tab)        ! once
   !         call dust_emission_interp(tab, U, lamI_total [, lamI_chan]) ! per cell
   !
   ! NONLINEARITY CAVEAT: stochastic heating makes the emission a NONLINEAR
   ! functional of the full J(lambda), not of a scalar. The table is valid
   ! ONLY when the cell field is  U * (fixed J_ref shape).  For cells whose
   ! field SHAPE departs from J_ref (e.g. hardened spectra near hot stars),
   ! use the single-cell exact dust_emission instead. The table is built from
   ! exact solves, so it is exact AT the U grid points; interpolation between
   ! them is the usual smooth-in-U approximation.
   !
   ! dust_emission takes an optional final argument, status (integer):
   !   call dust_emission(m, J_lam, lamI_total [, lamI_chan] [, status])
   ! On return status = 0 means success, 1 an unknown m%stoch_method, and 2 a
   ! 'qm' model whose populations are missing their radii. When status is
   ! present a bad model is reported through it instead of stopping the
   ! process; when it is omitted such a model stops the run, as before.
   !
   ! The validated solver core (sed_grain_loop & helpers in sed_astrodust_mod)
   ! is untouched; this module only re-exports the model API and adds the
   ! table/interpolation layer.
   use constants,         only: wp
   use mathlib,           only: interp
   use sed_astrodust_mod, only: dust_model_t, &
                                build_astrodust, build_dl07, build_zubko, build_from_files, &
                                dust_emission, dust_emission_single_teq
   implicit none
   private

   ! Re-exported model API
   public :: dust_model_t, build_astrodust, build_dl07, build_zubko, build_from_files
   public :: dust_emission, dust_emission_single_teq
   ! Table API
   public :: dust_emis_table_t, dust_build_table, dust_emission_interp, dust_free_table
   ! Convenience accessors
   public :: dust_nlam, dust_lambda, dust_n_channel, dust_channel_name

   type :: dust_emis_table_t
      integer               :: NLAM = 0, n_channel = 0, NU = 0
      real(wp), allocatable :: U(:)            ! (NU)   intensity-scaling grid
      real(wp), allocatable :: lam(:)          ! (NLAM) [um]
      real(wp), allocatable :: J_ref(:)        ! (NLAM) reference field shape (U=1)
      real(wp), allocatable :: total(:,:)      ! (NLAM, NU)
      real(wp), allocatable :: chan(:,:,:)     ! (NLAM, n_channel, NU)
   end type dust_emis_table_t

contains

   ! --- accessors -------------------------------------------------------
   pure integer function dust_nlam(m)
      type(dust_model_t), intent(in) :: m
      dust_nlam = m%NLAM
   end function dust_nlam

   pure integer function dust_n_channel(m)
      type(dust_model_t), intent(in) :: m
      dust_n_channel = m%n_channel
   end function dust_n_channel

   function dust_lambda(m) result(lam)
      type(dust_model_t), intent(in) :: m
      real(wp), allocatable :: lam(:)
      lam = m%lam
   end function dust_lambda

   function dust_channel_name(m, ic) result(name)
      type(dust_model_t), intent(in) :: m
      integer,            intent(in) :: ic
      character(len=16) :: name
      name = m%channel_name(ic)
   end function dust_channel_name

   ! --- emission table over an intensity grid ---------------------------
   subroutine dust_build_table(m, J_ref, U_grid, tab)
      ! Precompute lamI(lambda) for J = U*J_ref at each U in U_grid (which
      ! should be ascending). m must be the active model.
      type(dust_model_t),      intent(in)  :: m
      real(wp),                intent(in)  :: J_ref(:)    ! (NLAM) shape at U=1
      real(wp),                intent(in)  :: U_grid(:)   ! (NU)
      type(dust_emis_table_t), intent(out) :: tab
      real(wp), allocatable :: total(:), chan(:,:)
      integer :: iu

      call dust_free_table(tab)
      tab%NLAM = m%NLAM;  tab%n_channel = m%n_channel;  tab%NU = size(U_grid)
      allocate(tab%U(tab%NU), tab%lam(tab%NLAM), tab%J_ref(tab%NLAM))
      allocate(tab%total(tab%NLAM, tab%NU))
      allocate(tab%chan(tab%NLAM, tab%n_channel, tab%NU))
      tab%U = U_grid;  tab%lam = m%lam;  tab%J_ref = J_ref

      allocate(total(m%NLAM), chan(m%NLAM, m%n_channel))
      do iu = 1, tab%NU
         call dust_emission(m, U_grid(iu)*J_ref, total, chan)
         tab%total(:, iu)  = total
         tab%chan(:, :, iu) = chan
      end do
      deallocate(total, chan)
   end subroutine dust_build_table

   subroutine dust_emission_interp(tab, U, lamI_total, lamI_chan)
      ! Log-log interpolate the table at intensity U (per wavelength &
      ! channel). Exact at the U grid points (modulo exp(log) round-off);
      ! clamps to the grid ends outside [U(1), U(NU)].
      type(dust_emis_table_t), intent(in)  :: tab
      real(wp),                intent(in)  :: U
      real(wp),                intent(out) :: lamI_total(:)      ! (NLAM)
      real(wp), optional,      intent(out) :: lamI_chan(:,:)     ! (NLAM, n_channel)
      real(wp), allocatable :: lU(:)
      real(wp) :: lr, lUq
      integer  :: k, c

      allocate(lU(tab%NU));  lU = log(tab%U)
      lUq = log(U)
      do k = 1, tab%NLAM
         call interp(lU, log(max(tab%total(k, :), 1.0e-300_wp)), lUq, lr)
         lamI_total(k) = exp(lr)
      end do
      if (present(lamI_chan)) then
         do c = 1, tab%n_channel
            do k = 1, tab%NLAM
               call interp(lU, log(max(tab%chan(k, c, :), 1.0e-300_wp)), lUq, lr)
               lamI_chan(k, c) = exp(lr)
            end do
         end do
      end if
      deallocate(lU)
   end subroutine dust_emission_interp

   subroutine dust_free_table(tab)
      type(dust_emis_table_t), intent(inout) :: tab
      if (allocated(tab%U))     deallocate(tab%U)
      if (allocated(tab%lam))   deallocate(tab%lam)
      if (allocated(tab%J_ref)) deallocate(tab%J_ref)
      if (allocated(tab%total)) deallocate(tab%total)
      if (allocated(tab%chan))  deallocate(tab%chan)
      tab%NLAM = 0;  tab%n_channel = 0;  tab%NU = 0
   end subroutine dust_free_table

end module dust_lib
