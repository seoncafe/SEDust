program main_zubko
   ! Zubko, Dwek & Arendt (2004) ZDA BARE-GR-S driver: PAH + graphite +
   ! silicate with the ZDA size-distribution formula and the ZDA optics and
   ! calorimetry tables, heated by the Mathis ISRF at U = 1.
   !
   ! This is the emission counterpart of `./calc_kext.x zubko`, which only
   ! exercises the extinction size integral. It exists because the Zubko model
   ! is the one shipped model whose optics grid (1.0e-3 um, from the ZDA optics
   ! tables) reaches far past the band a transported radiation field occupies
   ! (0.0912 um for the Mathis ISRF). Anything in the stochastic-heating path
   ! that reads the grid where it should read the field shows up here and
   ! nowhere else, which is exactly how the regression recorded in
   ! docs/EUV_EXTENSION_HOST_REGRESSION.md escaped main_astrodust.x and
   ! main_dl07.x.
   !
   ! Usage:  ./main_zubko.x [heuristic | draine | qm | equil] [euv] [hardfield]
   !   solver    : stochastic solver; default 'heuristic' (the model default).
   !   euv       : build the model on the WHOLE ZDA optics range, 1.0e-3 um
   !               (1.24 keV) upward.  Without it the grid is cut at the Lyman
   !               limit, where the Mathis field stops, which is what the
   !               shipped kext_zubko_BARE_GR_S.dat is written on.
   !   hardfield : fill the band below the Lyman limit, where the Mathis field
   !               is identically zero, with a diluted 1e5 K blackbody, so that
   !               the field really does carry photons above 13.6 eV.  It
   !               implies `euv`, since a grid cut at the Lyman limit has no
   !               wavelength to put those photons on.  Without it the field
   !               stops at the Lyman limit and the single-photon bound must
   !               stay at 13.6 eV no matter how far the optics grid extends;
   !               with it the bound must follow the field up.  Running both is
   !               the regression test.
   !
   ! Output: output/zubko_sed_ours[_<solver>][_euv][_hardfield].dat
   !   columns: lambda[um]  lamI_total/NH  lamI_PAH/NH  lamI_GRA/NH  lamI_SIL/NH
   !            [erg s^-1 cm^-2 sr^-1 H^-1]

   use, intrinsic :: iso_fortran_env, only: real64
   use constants,         only: wp
   use radfield,          only: J_Mathis, use_mathis_corrected, bbody, &
                                hardest_photon_energy
   use sed_astrodust_mod, only: dust_model_t, build_zubko, dust_emission

   implicit none

   character(len=*), parameter :: F_QT_ZU    = '../data/zubko/sedust_zubko.h5'
   character(len=*), parameter :: F_ZDA_CFG = '../data/zubko/ZDA_BARE_GR_S_Config.dat'
   character(len=*), parameter :: D_ZUBKO   = '../data/zubko/'
   real(wp),         parameter :: U_ISRF = 1.0_wp        ! reference field: U = 1
   real(wp),         parameter :: T_LO = 2.7_wp, T_HI = 5.0e3_wp
   integer,          parameter :: NT_IN = 200
   character(len=*), parameter :: OUTDIR = 'output/zubko_sed_ours'

   ! Artificial hard component (only with the 'hardfield' argument): a diluted 1e5 K
   ! blackbody occupying the band below the Lyman limit, where the Mathis field
   ! is identically zero. 1e5 K puts the Wien cut-off of the added component at
   ! a few hundred eV, well above 13.6 eV and well inside the model's 1.24 keV
   ! grid, so the illuminated band ends strictly between the two -- which is
   ! what separates a field-based bound from a grid-based one and from a fixed
   ! 13.6 eV one in a single run. The dilution is chosen so the EUV band
   ! deposits roughly as much power as the whole Mathis field: enough for the
   ! hard photons to matter, little enough that the small grains still heat
   ! stochastically rather than settling into equilibrium.
   real(wp), parameter :: T_HARD = 1.0e5_wp
   real(wp), parameter :: W_HARD = 2.0e-19_wp
   real(wp), parameter :: LAM_LYMAN = 0.0912_wp

   type(dust_model_t)    :: m
   real(wp), allocatable :: J_lam(:), lamI_tot(:), lamI_chan(:,:)
   integer  :: k, ic, u, narg, iarg, status
   logical  :: add_hard_field, use_euv_grid
   character(len=32) :: solver, arg
   character(len=96) :: fname
   character(len=16) :: tag

   solver  = 'heuristic'
   add_hard_field = .false.;  use_euv_grid = .false.
   narg = command_argument_count()
   do iarg = 1, narg
      call get_command_argument(iarg, arg)
      select case (trim(arg))
      case ('euv')
         use_euv_grid = .true.
      case ('hardfield')
         add_hard_field = .true.
         use_euv_grid   = .true.    ! the photons need a grid to sit on
      case ('heuristic', 'draine', 'qm', 'equil')
         solver = trim(arg)
      case default
         write(*,'(a,a)') ' main_zubko: unknown argument ', trim(arg)
         write(*,'(a)')   ' usage: ./main_zubko.x [heuristic|draine|qm|equil] [euv] [hardfield]'
         stop 1
      end select
   end do

   write(*,'(a)') '=========================================================='
   write(*,'(a)') ' main_zubko: Zubko/ZDA BARE-GR-S dust emission'
   write(*,'(a)') '=========================================================='
   use_mathis_corrected = .true.
   write(*,'(a,f6.3)') ' U (Mathis)    : ', U_ISRF
   write(*,'(a,l1)')   ' mathis_corr   : ', use_mathis_corrected
   write(*,'(a,a)')    ' solver        : ', trim(solver)
   write(*,'(a,l1)')   ' EUV grid      : ', use_euv_grid
   write(*,'(a,l1)')   ' hard field    : ', add_hard_field

   call build_zubko(m, F_ZDA_CFG, D_ZUBKO, NT_IN, T_LO, T_HI, status, &
                    include_euv=use_euv_grid, qtable_path=F_QT_ZU)
   if (status /= 0) then
      write(*,'(a,i0)') ' main_zubko: build_zubko failed, status = ', status
      stop 1
   end if
   m%stoch_method = trim(solver)
   m%verbose      = .true.
   write(*,'(a,i0,a,es11.4,a,es11.4,a)') ' build_zubko done. NLAM=', m%NLAM, &
        ',  grid = ', m%lam(1), ' to ', m%lam(m%NLAM), ' um'

   allocate(J_lam(m%NLAM), lamI_tot(m%NLAM), lamI_chan(m%NLAM, m%n_channel))
   call J_Mathis(U_ISRF, m%lam, J_lam)
   if (add_hard_field) then
      do k = 1, m%NLAM
         if (m%lam(k) < LAM_LYMAN) J_lam(k) = W_HARD * bbody(T_HARD, m%lam(k))
      end do
   end if
   call report_field_extent()

   write(*,'(a)') ' solving Zubko SED (PAH + graphite + silicate) ...'
   call dust_emission(m, J_lam, lamI_tot, lamI_chan, status)
   if (status /= 0) then
      write(*,'(a,i0)') ' main_zubko: dust_emission failed, status = ', status
      stop 1
   end if

   tag = ''
   if (trim(solver) /= 'heuristic') tag = '_' // trim(solver)
   if (use_euv_grid)   tag = trim(tag) // '_euv'
   if (add_hard_field) tag = trim(tag) // '_hardfield'
   write(fname,'(a,a,a)') OUTDIR, trim(tag), '.dat'
   open(newunit=u, file=trim(fname), status='replace', action='write')
   write(u,'(a,f5.2)') '# Zubko/ZDA BARE-GR-S model SED (this work), Mathis ISRF U = ', U_ISRF
   write(u,'(a,a)')    '# solver = ', trim(solver)
   if (add_hard_field) write(u,'(a,es10.3,a,es10.3)') &
        '# artificial hard component below the Lyman limit: W*B_lambda(T), T = ', &
        T_HARD, ' K, W = ', W_HARD
   write(u,'(a)') '# columns: lambda[um]  lamI_total/NH  lamI_PAH/NH  lamI_GRA/NH  lamI_SIL/NH'
   write(u,'(a)') '#          [erg s^-1 cm^-2 sr^-1 H^-1]'
   ! e3 on the intensities: this model's optics grid already starts at 1e-3 um
   ! (1.24 keV), where lambda*I_lambda underflows to subnormals, and a two-digit
   ! exponent field drops the E there (4.29970510-319), which readers outside
   ! Fortran cannot parse.
   do k = 1, m%NLAM
      write(u,'(es14.6,4(1x,es16.8e3))') m%lam(k), lamI_tot(k), &
           (lamI_chan(k, ic), ic = 1, m%n_channel)
   end do
   close(u)
   write(*,'(a,a)') ' wrote ', trim(fname)
   write(*,'(a)') ' main_zubko: done.'

   deallocate(J_lam, lamI_tot, lamI_chan)

contains

   subroutine report_field_extent()
      ! Print the two competing definitions of "hardest photon" side by side:
      ! the short end of the model's optics grid, and the energy that
      ! hardest_photon_energy reads off the field. The stochastic solvers use
      ! the second; for this model the first is 91 times larger unless 'euv'
      ! is set, which is the whole content of the regression.
      real(wp), parameter :: ERG2EV = 1.0_wp / 1.602176565e-12_wp
      real(wp) :: u_field
      u_field = hardest_photon_energy(m%lam, J_lam)
      write(*,'(a,es11.4,a,es12.5,a)') ' optics grid short end : ', m%lam(1), &
           ' um  ->  ', 1.23984193_wp / m%lam(1), ' eV'
      write(*,'(a,es12.5,a)')          ' hardest photon in the field:   ', &
           u_field * ERG2EV, ' eV'
      write(*,'(a,es12.5,a)')          ' single-photon bound in use :   ', &
           max(13.6_wp / ERG2EV, u_field) * ERG2EV, ' eV'
   end subroutine report_field_extent

end program main_zubko
