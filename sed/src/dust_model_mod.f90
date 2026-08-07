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
   !   rho_bulk [g/cm^3]  solid mass density of the grain material, which turns
   !                      dn into the population's mass per H
   !   Cabs/Csca[cm^2]    absorption/scattering cross section (NLAM, NA)
   !   gsca               scattering asymmetry <cos> (NLAM, NA), read only by
   !                      size_integrated_extinction
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
      ! Scattering asymmetry <cos>, read only by size_integrated_extinction.
      ! Left unallocated for a population that does not scatter (the PAHs),
      ! which is how the extinction size integral recognizes a zero contribution.
      real(wp), allocatable :: gsca(:,:)            ! (NLAM, NA) scattering asymmetry <cos>
      ! Solid mass density of the grain material [g/cm^3], as the model defines
      ! it -- for a porous grain the density already reduced by the porosity.
      ! It converts the binned number per H into a mass per H,
      !   M/H = rho_bulk * sum_a (4/3) pi a_cm^3 dn(a),
      ! which is what dust_mass_per_H sums over the populations.  A model that
      ! states no density leaves it at 0, and that population then contributes
      ! nothing to the model's dust mass.
      real(wp)              :: rho_bulk = 0.0_wp    ! [g/cm^3]
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
      ! Size-integrated extinction curve this model serves to an RT host
      ! (dust_extinction), read from a precomputed data/kext_*.dat table.
      ! It is loaded by the builder, not by dust_extinction, because the table
      ! path is relative to the sed/ directory and a host is free to change
      ! directory once the model is built.  kext_n = 0 means no table was
      ! loaded and dust_extinction has nothing to return.
      character(len=512)    :: kext_path = ''
      integer               :: kext_n    = 0
      real(wp), allocatable :: kext_lam(:)                ! (kext_n) [um] ascending
      real(wp), allocatable :: kext_Cext(:), kext_Cabs(:) ! (kext_n) [cm^2/H]
      real(wp), allocatable :: kext_Csca(:)               ! (kext_n) [cm^2/H]
      real(wp), allocatable :: kext_gbar(:)               ! (kext_n) <cos>
      ! Which build produced this model.  The emission side is solved on the
      ! module-level grids of sed_astrodust_mod, which the LAST build filled,
      ! so dust_emission and size_integrated_extinction can only answer for the
      ! model built most recently.  A host keeping two models alive and
      ! querying them alternately used to get the wrong one's numbers in
      ! silence; the builders stamp this, those routines compare it against the
      ! stamp of the active build, and a mismatch is a status instead.  0 marks
      ! a model that was never built.  dust_extinction is exempt: it serves
      ! kext_* off the model argument and reads no module grid.
      integer               :: build_id  = 0
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
      if (allocated(m%kext_lam))     deallocate(m%kext_lam)
      if (allocated(m%kext_Cext))    deallocate(m%kext_Cext)
      if (allocated(m%kext_Cabs))    deallocate(m%kext_Cabs)
      if (allocated(m%kext_Csca))    deallocate(m%kext_Csca)
      if (allocated(m%kext_gbar))    deallocate(m%kext_gbar)
      m%NA = 0; m%NLAM = 0; m%NT = 0; m%n_channel = 0
      m%kext_n = 0; m%kext_path = ''
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
