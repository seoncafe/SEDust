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
   !   Cpol     [cm^2]    polarized absorption cross section (NLAM, NA)
   !   Cpol_ext [cm^2]    polarized extinction cross section (NLAM, NA)
   !   gsca               scattering asymmetry <cos> (NLAM, NA)
   !   falign             alignment efficiency, 0 for a population that does
   !                      not contribute to polarization (NA)
   !   kappB    [..]      Planck-integral table used by calc_P/calc_Teq (NT, NA)
   !   H        [erg]     grain enthalpy (NT, NA)
   !   kappCMB  [..]      CMB-pumped term for calc_P (NA)
   !
   ! Types plus trivial (de)allocation helpers.
   use constants,        only: wp
   use q_table_jori_mod, only: falign_powerlaw, &
                               A_ALIGN, ALPHA_ALIGN, FMAX_ALIGN
   implicit none
   private
   public :: grain_pop_t, dust_model_t, free_dust_model
   public :: dust_set_alignment, dust_set_alignment_profile

   ! One stochastically-heated population (= one species/charge state).
   type :: grain_pop_t
      character(len=8)      :: grain_type = 'sil'   ! 'sil' | 'pah' | 'gra'
      integer               :: out_channel = 1      ! index into dust_model_t channels
      real(wp), allocatable :: aeff(:)              ! (NA) [um] effective-radius grid
      real(wp), allocatable :: dn(:)                ! (NA)
      real(wp), allocatable :: Cabs(:,:), Csca(:,:) ! (NLAM, NA)
      real(wp), allocatable :: Cpol(:,:)            ! (NLAM, NA) [cm^2] polarized absorption
      ! Extinction-side optics, used by dust_extinction. Left unallocated for a
      ! population that neither scatters nor polarizes (the PAHs), which is how
      ! the size integral recognizes a zero contribution.
      real(wp), allocatable :: Cpol_ext(:,:)        ! (NLAM, NA) [cm^2] polarized extinction
      real(wp), allocatable :: gsca(:,:)            ! (NLAM, NA) scattering asymmetry <cos>
      real(wp), allocatable :: falign(:)            ! (NA) alignment efficiency
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
      ! Grain-alignment state currently loaded into the populations' falign.
      ! The defaults are the Hensley & Draine (2023) Table 1 fit, which is
      ! what build_astrodust installs; dust_set_alignment overwrites them.
      ! They are read-only bookkeeping -- editing them does NOT re-fill
      ! falign. Meaningful only while align_tabulated is .false.; after
      ! dust_set_alignment_profile the loaded efficiency is an arbitrary
      ! function that no power law describes, and the three scalars still
      ! hold whatever was set before. Models without polarized optics (DL07,
      ! Zubko) carry the defaults but never apply them, since no population
      ! of theirs has a Cpol to align.
      real(wp)              :: align_fmax  = FMAX_ALIGN
      real(wp)              :: align_a     = A_ALIGN      ! [um]
      real(wp)              :: align_alpha = ALPHA_ALIGN
      logical               :: align_tabulated = .false.
      ! When .false. (default), the library solve path stays silent; when
      ! .true. it emits the same solver diagnostics as the CLI drivers.
      logical               :: verbose = .false.
   end type dust_model_t

contains

   ! ---- grain alignment -------------------------------------------------
   !
   ! The alignment efficiency f_align(a) enters ONLY as a size weight on the
   ! polarized optics: the grain loop forms dn(a)*f_align(a)*Cpol(lam,a) for
   ! the polarized emission and dust_extinction forms dn(a)*f_align(a)*
   ! Cpol_ext(lam,a) for the dichroic extinction. It appears nowhere in the
   ! energy balance, so P(T) and the total (unpolarized) emission are
   ! independent of it: changing the alignment does NOT invalidate any
   ! temperature solution and costs one function evaluation on the size grid.
   ! An RT host can therefore give every cell its own alignment state at
   ! essentially no cost.
   !
   ! Both setters refill falign only for populations that carry polarized
   ! optics (Cpol for emission or Cpol_ext for extinction). A population
   ! without either -- the PAHs, which HD23 take to be unaligned -- is left
   ! untouched and keeps contributing nothing to polarization.

   ! Install the HD23 power-law rolloff
   !     f_align(a) = f_max / (1 + (a_align/a)**alpha_align)
   ! evaluated on each population's own radius grid (the grids may differ
   ! between populations).
   !
   ! status (optional): 0 = success; when absent a bad argument stops the run.
   !   1  a_align <= 0
   !   2  alpha_align <= 0
   !   3  f_max outside [0, 1]
   subroutine dust_set_alignment(m, f_max, a_align, alpha_align, status)
      type(dust_model_t), intent(inout) :: m
      real(wp),           intent(in)    :: f_max, a_align, alpha_align
      integer, optional,  intent(out)   :: status
      integer :: ip, ia, na_p

      if (present(status)) status = 0

      if (a_align <= 0.0_wp) then
         call alignment_error(1, 'a_align must be > 0')
         return
      end if
      if (alpha_align <= 0.0_wp) then
         call alignment_error(2, 'alpha_align must be > 0')
         return
      end if
      if (f_max < 0.0_wp .or. f_max > 1.0_wp) then
         call alignment_error(3, 'f_max must lie in [0, 1]')
         return
      end if

      if (allocated(m%pops)) then
         do ip = 1, size(m%pops)
            if (.not. polarizable(m%pops(ip))) cycle
            na_p = size(m%pops(ip)%aeff)
            if (.not. allocated(m%pops(ip)%falign)) allocate(m%pops(ip)%falign(na_p))
            do ia = 1, na_p
               m%pops(ip)%falign(ia) = falign_powerlaw(m%pops(ip)%aeff(ia), &
                                                       f_max, a_align, alpha_align)
            end do
         end do
      end if

      m%align_fmax = f_max;  m%align_a = a_align;  m%align_alpha = alpha_align
      m%align_tabulated = .false.

   contains

      subroutine alignment_error(code, msg)
         integer,          intent(in) :: code
         character(len=*), intent(in) :: msg
         if (present(status)) then
            status = code
         else
            write(*,'(a,a)') 'dust_set_alignment: ', msg
            stop 1
         end if
      end subroutine alignment_error

   end subroutine dust_set_alignment


   ! Install an arbitrary tabulated alignment efficiency: the caller supplies
   ! (aeff_in [um], falign_in) pairs and each population's falign is
   ! interpolated from them onto its own radius grid. This is the route for
   ! prescriptions the HD23 power law cannot express -- a RAT-derived Rayleigh
   ! reduction factor R(a), or the GRADE-POL exponential
   ! f_max*[1 - exp(-(0.5a/a_align)**3)].
   !
   ! Interpolation is linear in log(a) and clamped at the ends of aeff_in,
   ! the same convention the optics tables are interpolated with, so radii
   ! outside the supplied range take the nearest tabulated efficiency.
   !
   ! RANGE POLICY. falign_in must lie in [-1, 1] and is REJECTED, not clamped,
   ! outside it. |f| <= 1 is the physical bound on a Rayleigh reduction
   ! factor, so a value beyond it is a caller error (wrong units, wrong
   ! column) that silent clamping would hide. Negative values are deliberately
   ! ALLOWED: a grain in the wrong internal alignment state has a negative
   ! reduction factor, which flips the polarization direction by 90 degrees
   ! and is the accepted origin of polarization parallel to B at millimeter
   ! wavelengths and perpendicular to it in the submillimeter. Clamping at
   ! zero would silently delete that effect. Interpolation of a bracket whose
   ! endpoints both lie in [-1, 1] stays in [-1, 1], so validating the input
   ! is enough.
   !
   ! status (optional): 0 = success; when absent a bad argument stops the run.
   !   1  aeff_in and falign_in differ in size, or fewer than 2 points
   !   2  aeff_in not positive and strictly increasing
   !   3  a falign_in value outside [-1, 1]
   subroutine dust_set_alignment_profile(m, aeff_in, falign_in, status)
      type(dust_model_t), intent(inout) :: m
      real(wp),           intent(in)    :: aeff_in(:)     ! (N) [um]
      real(wp),           intent(in)    :: falign_in(:)   ! (N)
      integer, optional,  intent(out)   :: status
      integer :: ip, ia, na_p, n, i

      if (present(status)) status = 0
      n = size(aeff_in)

      if (size(falign_in) /= n .or. n < 2) then
         call profile_error(1, 'need matching aeff_in/falign_in of length >= 2')
         return
      end if
      if (aeff_in(1) <= 0.0_wp) then
         call profile_error(2, 'aeff_in must be positive and strictly increasing')
         return
      end if
      do i = 2, n
         if (aeff_in(i) <= aeff_in(i-1)) then
            call profile_error(2, 'aeff_in must be positive and strictly increasing')
            return
         end if
      end do
      do i = 1, n
         if (falign_in(i) < -1.0_wp .or. falign_in(i) > 1.0_wp) then
            call profile_error(3, 'falign_in must lie in [-1, 1]')
            return
         end if
      end do

      if (allocated(m%pops)) then
         do ip = 1, size(m%pops)
            if (.not. polarizable(m%pops(ip))) cycle
            na_p = size(m%pops(ip)%aeff)
            if (.not. allocated(m%pops(ip)%falign)) allocate(m%pops(ip)%falign(na_p))
            do ia = 1, na_p
               m%pops(ip)%falign(ia) = falign_at_radius(m%pops(ip)%aeff(ia), &
                                                        aeff_in, falign_in)
            end do
         end do
      end if

      m%align_tabulated = .true.

   contains

      subroutine profile_error(code, msg)
         integer,          intent(in) :: code
         character(len=*), intent(in) :: msg
         if (present(status)) then
            status = code
         else
            write(*,'(a,a)') 'dust_set_alignment_profile: ', msg
            stop 1
         end if
      end subroutine profile_error

   end subroutine dust_set_alignment_profile


   ! Alignment efficiency at radius a [um] read off a tabulated profile,
   ! linear in log(a) and clamped at the ends of the table.
   pure function falign_at_radius(a, aeff_in, falign_in) result(f)
      real(wp), intent(in) :: a
      real(wp), intent(in) :: aeff_in(:), falign_in(:)
      real(wp)             :: f
      real(wp) :: loga, x_lo, x_hi, t
      integer  :: n, lo, hi, mid

      n = size(aeff_in)
      if (a <= 0.0_wp .or. a <= aeff_in(1)) then
         f = falign_in(1);  return
      end if
      if (a >= aeff_in(n)) then
         f = falign_in(n);  return
      end if
      loga = log(a)
      lo = 1;  hi = n
      do while (hi - lo > 1)
         mid = (lo + hi) / 2
         if (log(aeff_in(mid)) <= loga) then
            lo = mid
         else
            hi = mid
         end if
      end do
      x_lo = log(aeff_in(lo))
      x_hi = log(aeff_in(hi))
      t = (loga - x_lo) / (x_hi - x_lo)
      f = (1.0_wp - t) * falign_in(lo) + t * falign_in(hi)
   end function falign_at_radius


   ! .true. for a population that carries polarized optics on either the
   ! emission (Cpol) or the extinction (Cpol_ext) side, and so needs an
   ! alignment efficiency. Needs a radius grid to evaluate one on.
   pure logical function polarizable(p)
      type(grain_pop_t), intent(in) :: p
      polarizable = allocated(p%aeff) .and. &
                    (allocated(p%Cpol) .or. allocated(p%Cpol_ext))
   end function polarizable


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
      if (allocated(p%Cpol))      deallocate(p%Cpol)
      if (allocated(p%Cpol_ext))  deallocate(p%Cpol_ext)
      if (allocated(p%gsca))      deallocate(p%gsca)
      if (allocated(p%falign))    deallocate(p%falign)
      if (allocated(p%kappB))     deallocate(p%kappB)
      if (allocated(p%log_kappB)) deallocate(p%log_kappB)
      if (allocated(p%H))         deallocate(p%H)
      if (allocated(p%log_H))     deallocate(p%log_H)
      if (allocated(p%kappCMB))   deallocate(p%kappCMB)
   end subroutine free_pop

end module dust_model_mod
