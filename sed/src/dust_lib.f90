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
   ! exact solves, so it reproduces the U grid points to round-off, except that
   ! dust_emission_interp floors the stored emissivities at 1e-300 before the
   ! log, so a node whose true value is 0 comes back as 1e-300, not 0;
   ! interpolation between grid points is the usual smooth-in-U approximation.
   !
   ! dust_emission takes an optional final argument, status (integer):
   !   call dust_emission(m, J_lam, lamI_total [, lamI_chan] [, status])
   ! On return status = 0 means success, 1 an unknown m%stoch_method, and 2 a
   ! 'qm' model whose populations are missing their radii. When status is
   ! present a bad model is reported through it instead of stopping the
   ! process; when it is omitted such a model stops the run, as before.
   !
   ! dust_emission also takes an optional lamI_pol (same shape and units as
   ! lamI_total):
   !   call dust_emission(m, J_lam, lamI_total [, lamI_chan] [, status] [, lamI_pol])
   ! It returns the INTRINSIC polarized emission -- what a population of
   ! perfectly aligned grains seen with their symmetry axes in the plane of
   ! the sky would radiate -- summed over the populations that carry both
   ! Cpol and falign. Populations without polarized optics (PAHs, and any
   ! model built without them) contribute zero. The geometric sin^2(gamma)
   ! projection and any turbulent depolarization are left to the radiative
   ! transfer and are NOT applied here. Being optional, it leaves every
   ! existing caller valid.
   !
   ! EXTINCTION. dust_extinction is the extinction counterpart of
   ! dust_emission, so an RT host takes its opacity from the same model object
   ! and on the same wavelength grid (m%lam) as its emission, rather than
   ! parsing data/kext_astrodust_MW.dat and interpolating off that file's grid:
   !   call dust_extinction(m, Cext, Cabs, Csca [, gbar] [, Cpol_ext] [, status])
   ! All three required outputs are (m%NLAM) cross sections per H atom
   ! [cm^2/H], integrated over the size distribution of every population:
   ! Cext = Cabs + Csca. Optional gbar is the scattering-weighted asymmetry
   ! <cos>, sum dn*Csca*g / sum dn*Csca, and is 0 at wavelengths where nothing
   ! scatters. Optional Cpol_ext is the dichroic (polarized) extinction,
   ! sum dn*Cpol_ext*f_align [cm^2/H]; the size integral and the alignment
   ! weight are done here, but the sin^2(gamma) geometry factor and any
   ! turbulent depolarization are left to the radiative transfer, exactly as
   ! for lamI_pol. Populations without scattering or polarized optics (the
   ! PAHs) contribute zero to those terms and enter through absorption only.
   ! status is 0 on success and 1 if an output array is not of size m%NLAM;
   ! when it is omitted such a call stops the run.
   !
   ! GRAIN ALIGNMENT. Both the polarized emission (lamI_pol) and the dichroic
   ! extinction (Cpol_ext) are weighted by an alignment efficiency
   ! f_align(a_eff), which build_astrodust initializes to the Hensley & Draine
   ! (2023) Table 1 fit. A host that wants a cell-dependent alignment state
   ! overrides it, either as that power law with its own parameters
   !     call dust_set_alignment(m, f_max, a_align, alpha_align [, status])
   !     f_align(a) = f_max / (1 + (a_align/a)**alpha_align)
   ! or as an arbitrary tabulated profile, interpolated in log(a) onto each
   ! population's radius grid and clamped at the ends of the table
   !     call dust_set_alignment_profile(m, aeff_in, falign_in [, status])
   ! for prescriptions that power law cannot express (a RAT-derived Rayleigh
   ! reduction factor, the GRADE-POL exponential). Both leave every population
   ! without polarized optics -- the PAHs, which HD23 take to be unaligned --
   ! untouched and contributing zero.
   !
   ! Calling neither leaves the HD23 alignment in place, so an existing host
   ! sees no change. The current state is readable off the model as
   ! m%align_fmax, m%align_a [um], m%align_alpha and m%align_tabulated; the
   ! three scalars describe the loaded efficiency only while align_tabulated
   ! is .false. They are reported, not applied -- assigning to them does not
   ! re-fill f_align; only the two setters do.
   !
   ! Alignment is a size WEIGHT and enters nowhere in the energy balance, so
   ! resetting it does not invalidate any P(T) solution and does not require a
   ! re-solve: it is one function evaluation on the size grid, and the total
   ! unpolarized emission is bit-for-bit unchanged. A single setter feeds both
   ! dust_emission's lamI_pol and dust_extinction's Cpol_ext, which therefore
   ! cannot fall out of step.
   !
   ! Division of labor, as for lamI_pol and Cpol_ext throughout: SEDust does
   ! the size-distribution integral and the alignment weight; the sin^2(gamma)
   ! projection onto the plane of the sky, any turbulent depolarization
   ! F_turb, and the position angle are the host's job.
   !
   ! status codes (0 = success; when omitted a bad argument stops the run):
   !   dust_set_alignment:          1 a_align <= 0
   !                                2 alpha_align <= 0
   !                                3 f_max outside [0, 1]
   !   dust_set_alignment_profile:  1 aeff_in/falign_in size mismatch or < 2 points
   !                                2 aeff_in not positive and strictly increasing
   !                                3 a falign_in value outside [-1, 1]
   ! A tabulated efficiency outside [-1, 1] is rejected rather than clamped:
   ! |f| <= 1 is the physical bound, so exceeding it is a caller error.
   ! NEGATIVE values are accepted on purpose -- a grain in the wrong internal
   ! alignment state has a negative Rayleigh reduction factor, flipping the
   ! polarization direction by 90 degrees, and clamping at zero would delete
   ! that effect.
   !
   ! dust_build_table and dust_emission_interp take the same optional final
   ! status argument (0 = success); when present a bad argument is reported
   ! through it instead of stopping the process; when omitted such a call stops
   ! the run.
   !   dust_build_table:      1 size(J_ref) /= m%NLAM
   !                          2 size(U_grid) < 2
   !                          3 U_grid not positive-and-strictly-increasing
   !   dust_emission_interp:  1 U <= 0
   !                          2 size(lamI_total) /= tab%NLAM
   !                          3 lamI_chan present but not (tab%NLAM, tab%n_channel)
   !
   ! The model builders take the same optional final argument, status (integer,
   ! 0 = success). When present, a missing or malformed input file is reported
   ! through it and the model is NOT built (so an RT host can recover); when
   ! omitted such a failure stops the process, which the CLI drivers rely on.
   ! The exact non-zero code only distinguishes the failing stage; the contract
   ! is simply "0 = built, non-zero = build failed". Codes per builder:
   !   build_astrodust / build_dl07:  1 Q-table load failed
   !                                  2 size-distribution load failed
   !   build_zubko:   1 config read failed        2 fewer than 3 components
   !                  3 a component's optics read  4 grid inconsistent
   !                  5 a component's calorimetry read failed
   !   build_from_files: 1 descriptor open   2 too many pop: lines
   !                     3 invalid channel   4 no pop: lines
   !                     5 optics read       6 grid inconsistent
   !                     7 size-dist read    8 calorimetry read failed
   !
   ! The validated solver core (sed_grain_loop & helpers in sed_astrodust_mod)
   ! is untouched; this module only re-exports the model API and adds the
   ! table/interpolation layer.
   use constants,         only: wp
   use mathlib,           only: locate
   use sed_astrodust_mod, only: dust_model_t, &
                                build_astrodust, build_dl07, build_zubko, build_from_files, &
                                dust_emission, dust_emission_single_teq, &
                                dust_extinction, &
                                dust_set_alignment, dust_set_alignment_profile
   implicit none
   private

   ! Re-exported model API
   public :: dust_model_t, build_astrodust, build_dl07, build_zubko, build_from_files
   public :: dust_emission, dust_emission_single_teq, dust_extinction
   public :: dust_set_alignment, dust_set_alignment_profile
   ! Table API
   public :: dust_emis_table_t, dust_build_table, dust_emission_interp, dust_free_table
   ! Convenience accessors
   public :: dust_nlam, dust_lambda, dust_n_channel, dust_channel_name

   type :: dust_emis_table_t
      integer               :: NLAM = 0, n_channel = 0, NU = 0
      real(wp), allocatable :: U(:)            ! (NU)   intensity-scaling grid
      real(wp), allocatable :: logU(:)         ! (NU)   log(U), cached for interp
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
   subroutine dust_build_table(m, J_ref, U_grid, tab, status)
      ! Precompute lamI(lambda) for J = U*J_ref at each U in U_grid (which must
      ! be positive and strictly ascending). m must be the active model.
      type(dust_model_t),      intent(in)  :: m
      real(wp),                intent(in)  :: J_ref(:)    ! (NLAM) shape at U=1
      real(wp),                intent(in)  :: U_grid(:)   ! (NU)
      type(dust_emis_table_t), intent(out) :: tab
      ! Optional error report (0 = success); see the module header for codes.
      integer, optional,       intent(out) :: status
      real(wp), allocatable :: total(:), chan(:,:)
      integer :: iu, nu

      if (present(status)) status = 0
      nu = size(U_grid)

      if (size(J_ref) /= m%NLAM) then
         if (present(status)) then
            status = 1;  return
         else
            write(*,'(a,i0,a,i0)') 'dust_build_table: size(J_ref)=', size(J_ref), &
                                    ' /= m%NLAM=', m%NLAM
            stop 1
         end if
      end if
      if (nu < 2) then
         if (present(status)) then
            status = 2;  return
         else
            write(*,'(a,i0)') 'dust_build_table: need size(U_grid) >= 2, got ', nu
            stop 1
         end if
      end if
      if (U_grid(1) <= 0.0_wp .or. any(U_grid(2:nu) <= U_grid(1:nu-1))) then
         if (present(status)) then
            status = 3;  return
         else
            write(*,'(a)') 'dust_build_table: U_grid must be positive and strictly increasing'
            stop 1
         end if
      end if

      call dust_free_table(tab)
      tab%NLAM = m%NLAM;  tab%n_channel = m%n_channel;  tab%NU = size(U_grid)
      allocate(tab%U(tab%NU), tab%logU(tab%NU), tab%lam(tab%NLAM), tab%J_ref(tab%NLAM))
      allocate(tab%total(tab%NLAM, tab%NU))
      allocate(tab%chan(tab%NLAM, tab%n_channel, tab%NU))
      tab%U = U_grid;  tab%lam = m%lam;  tab%J_ref = J_ref
      tab%logU = log(tab%U)   ! cached once; the U bracket is grid-fixed at interp

      allocate(total(m%NLAM), chan(m%NLAM, m%n_channel))
      do iu = 1, tab%NU
         call dust_emission(m, U_grid(iu)*J_ref, total, chan)
         tab%total(:, iu)  = total
         tab%chan(:, :, iu) = chan
      end do
      deallocate(total, chan)
   end subroutine dust_build_table

   subroutine dust_emission_interp(tab, U, lamI_total, lamI_chan, status)
      ! Log-log interpolate the table at intensity U (per wavelength and
      ! channel). Reproduces the U grid points to round-off, except that the
      ! stored emissivities are floored at 1e-300 before the log, so a node
      ! whose true value is 0 comes back as 1e-300, not 0; clamps to the grid
      ! ends outside [U(1), U(NU)].
      type(dust_emis_table_t), intent(in)  :: tab
      real(wp),                intent(in)  :: U
      real(wp),                intent(out) :: lamI_total(:)      ! (NLAM)
      real(wp), optional,      intent(out) :: lamI_chan(:,:)     ! (NLAM, n_channel)
      ! Optional error report (0 = success); see the module header for codes.
      ! Validation is scalar/size-only to keep this on the per-cell hot path.
      integer, optional,       intent(out) :: status
      real(wp) :: ly(tab%NU)
      real(wp) :: lr, lUq
      integer  :: k, c, jlo

      if (present(status)) status = 0
      if (U <= 0.0_wp) then
         if (present(status)) then
            status = 1;  return
         else
            write(*,'(a,es12.4)') 'dust_emission_interp: need U > 0, got ', U
            stop 1
         end if
      end if
      if (size(lamI_total) /= tab%NLAM) then
         if (present(status)) then
            status = 2;  return
         else
            write(*,'(a,i0,a,i0)') 'dust_emission_interp: size(lamI_total)=', &
                                    size(lamI_total), ' /= tab%NLAM=', tab%NLAM
            stop 1
         end if
      end if
      if (present(lamI_chan)) then
         if (size(lamI_chan,1) /= tab%NLAM .or. size(lamI_chan,2) /= tab%n_channel) then
            if (present(status)) then
               status = 3;  return
            else
               write(*,'(a)') 'dust_emission_interp: lamI_chan must be (tab%NLAM, tab%n_channel)'
               stop 1
            end if
         end if
      end if

      ! tab%logU (= log(tab%U)) is cached at build time and is strictly
      ! ascending because U_grid is positive and strictly increasing. The U
      ! bracket is the same for every wavelength and channel, so locate it once
      ! and reuse jlo below. The interpolation arithmetic is the ascending
      ! branch of interp1, kept term-for-term (log-log, clamped to the ends).
      lUq = log(U)
      call locate(tab%logU, lUq, jlo)
      do k = 1, tab%NLAM
         ly = log(max(tab%total(k, :), 1.0e-300_wp))
         if (jlo == 0) then
            lr = ly(1)
         else if (jlo == tab%NU) then
            lr = ly(tab%NU)
         else
            lr = ly(jlo) + (ly(jlo+1)-ly(jlo))*(lUq-tab%logU(jlo))/(tab%logU(jlo+1)-tab%logU(jlo))
         end if
         lamI_total(k) = exp(lr)
      end do
      if (present(lamI_chan)) then
         do c = 1, tab%n_channel
            do k = 1, tab%NLAM
               ly = log(max(tab%chan(k, c, :), 1.0e-300_wp))
               if (jlo == 0) then
                  lr = ly(1)
               else if (jlo == tab%NU) then
                  lr = ly(tab%NU)
               else
                  lr = ly(jlo) + (ly(jlo+1)-ly(jlo))*(lUq-tab%logU(jlo))/(tab%logU(jlo+1)-tab%logU(jlo))
               end if
               lamI_chan(k, c) = exp(lr)
            end do
         end do
      end if
   end subroutine dust_emission_interp

   subroutine dust_free_table(tab)
      type(dust_emis_table_t), intent(inout) :: tab
      if (allocated(tab%U))     deallocate(tab%U)
      if (allocated(tab%logU))  deallocate(tab%logU)
      if (allocated(tab%lam))   deallocate(tab%lam)
      if (allocated(tab%J_ref)) deallocate(tab%J_ref)
      if (allocated(tab%total)) deallocate(tab%total)
      if (allocated(tab%chan))  deallocate(tab%chan)
      tab%NLAM = 0;  tab%n_channel = 0;  tab%NU = 0
   end subroutine dust_free_table

end module dust_lib
