module dust_model_mod
   ! Derived types for the model-agnostic dust thermal-emission library.
   !
   ! A dust *model* (DL07, astrodust/HD23, Zubko/ZDA, or a file-defined model)
   ! is represented as a `dust_model_t`: shared wavelength / size / temperature
   ! grids plus an array of grain *populations* (`grain_pop_t`). Each population
   ! is one stochastically-heated species/charge-state, carrying its own optics,
   ! enthalpy, and Planck-integral tables. The generic solver
   ! (`dust_solver_mod`) and builders (`dust_builders_mod`) operate on these
   ! types; nothing here knows about a specific model.
   !
   ! Conventions (matching the existing sed_astrodust_mod globals):
   !   lam      [um]      wavelength grid (NLAM)
   !   aeff     [um]      effective-radius grid (NA)
   !   T_first  [K]       temperature grid (NT)
   !   dn       [1/H]     number of grains per H atom in each size bin (NA)
   !   Cabs/Csca[cm^2]    absorption/scattering cross section (NLAM, NA)
   !   gsca               scattering asymmetry <cos> (NLAM, NA), read only by
   !                      dust_extinction
   !   kappB    [..]      Planck-integral table used by calc_P/calc_Teq (NT, NA)
   !   H        [erg]     grain enthalpy (NT, NA)
   !   kappCMB  [..]      CMB-pumped term for calc_P (NA)
   !
   ! Types plus trivial (de)allocation helpers.
   use constants, only: wp
   implicit none
   private
   public :: grain_pop_t, dust_model_t, free_dust_model

   ! One stochastically-heated population (= one species/charge state).
   type :: grain_pop_t
      character(len=8)      :: grain_type = 'sil'   ! 'sil' | 'pah' | 'gra'
      integer               :: out_channel = 1      ! index into dust_model_t channels
      real(wp), allocatable :: aeff(:)              ! (NA) [um] effective-radius grid
      real(wp), allocatable :: dn(:)                ! (NA)
      real(wp), allocatable :: Cabs(:,:), Csca(:,:) ! (NLAM, NA)
      ! Scattering asymmetry <cos>, read only by dust_extinction. Left
      ! unallocated for a population that does not scatter (the PAHs), which is
      ! how the extinction size integral recognizes a zero contribution.
      real(wp), allocatable :: gsca(:,:)            ! (NLAM, NA) scattering asymmetry <cos>
      real(wp), allocatable :: kappB(:,:), log_kappB(:,:)   ! (NT, NA)
      real(wp), allocatable :: H(:,:),     log_H(:,:)       ! (NT, NA)
      real(wp), allocatable :: kappCMB(:)           ! (NA)
   end type grain_pop_t

   ! A full dust model = shared grids + a set of populations grouped into
   ! named output channels.
   type :: dust_model_t
      character(len=32)     :: name = ''
      integer               :: NA = 0, NLAM = 0, NT = 0
      real(wp), allocatable :: lam(:)          ! (NLAM) [um]
      real(wp), allocatable :: aeff(:)         ! (NA)   [um]
      real(wp), allocatable :: T_first(:)      ! (NT)   [K]
      real(wp), allocatable :: log_T_first(:)  ! (NT)
      type(grain_pop_t), allocatable :: pops(:)
      integer               :: n_channel = 0
      character(len=16), allocatable :: channel_name(:)   ! (n_channel)
      logical               :: use_induced_emission = .false.
      character(len=16)     :: stoch_method = 'heuristic'
      ! When .false. (default), the library solve path stays silent; when
      ! .true. it emits the same solver diagnostics as the CLI drivers.
      logical               :: verbose = .false.
   end type dust_model_t

contains

   ! Deallocate everything held by a model (safe on a zero/partly-filled model).
   subroutine free_dust_model(m)
      type(dust_model_t), intent(inout) :: m
      integer :: i
      if (allocated(m%pops)) then
         do i = 1, size(m%pops)
            call free_pop(m%pops(i))
         end do
         deallocate(m%pops)
      end if
      if (allocated(m%lam))          deallocate(m%lam)
      if (allocated(m%aeff))         deallocate(m%aeff)
      if (allocated(m%T_first))      deallocate(m%T_first)
      if (allocated(m%log_T_first))  deallocate(m%log_T_first)
      if (allocated(m%channel_name)) deallocate(m%channel_name)
      m%NA = 0; m%NLAM = 0; m%NT = 0; m%n_channel = 0
   end subroutine free_dust_model

   subroutine free_pop(p)
      type(grain_pop_t), intent(inout) :: p
      if (allocated(p%aeff))      deallocate(p%aeff)
      if (allocated(p%dn))        deallocate(p%dn)
      if (allocated(p%Cabs))      deallocate(p%Cabs)
      if (allocated(p%Csca))      deallocate(p%Csca)
      if (allocated(p%gsca))      deallocate(p%gsca)
      if (allocated(p%kappB))     deallocate(p%kappB)
      if (allocated(p%log_kappB)) deallocate(p%log_kappB)
      if (allocated(p%H))         deallocate(p%H)
      if (allocated(p%log_H))     deallocate(p%log_H)
      if (allocated(p%kappCMB))   deallocate(p%kappCMB)
   end subroutine free_pop

end module dust_model_mod
