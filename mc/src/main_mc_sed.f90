program main_mc_sed
   ! MC SED builder driver.  Mirrors calc_sed.x astrodust but uses the MC
   ! engine in place of the Guhathakurta-Draine matrix solver.
   !
   ! Inputs via namelist (positional arg 1):
   !
   !   &mc_sed_input
   !     qtable_path  = '../tmatrix/output/q_astrodust_P0.20_Fe0.00_1.400.dat'
   !     U_isrf       = 1.585           ! log U = 0.20
   !     ad_stage     = 'S2'            ! 'S1_C1' | 'S1_C2' | 'S2'
   !     N_events     = 20000
   !     base_seed    = 12345
   !     out_prefix   = 'output/sed_default'
   !     NT_init      = 250
   !     T_lo         = 1.0
   !     T_hi         = 3000.0
   !   /

   use omp_lib
   use iso_fortran_env, only: error_unit
   use constants,         only: wp
   use radfield,          only: J_Mathis, J_HD23, use_mathis_corrected
   use sed_astrodust_mod, only: sed_init, NLAM, lam
   use mc_sed,            only: mc_sed_solve_total
   implicit none

   character(len=512) :: qtable_path, out_prefix
   character(len=16)  :: ad_stage, mc_engine, radiation_field
   real(wp) :: U_isrf, T_lo, T_hi
   integer  :: N_events, base_seed, NT_init
   logical  :: use_mc_for_large_grain
   namelist /mc_sed_input/ qtable_path, U_isrf, ad_stage, &
                           N_events, base_seed, out_prefix, NT_init, T_lo, T_hi, &
                           use_mc_for_large_grain, mc_engine, radiation_field

   character(len=512) :: nml_path
   integer :: nargs, u, ios, k
   real(wp), allocatable :: J_lam(:), lamI_lam_ad(:), lamI_lam_pah(:), lamI_lam_tot(:)
   real(wp) :: t0, t1

   ! defaults
   qtable_path   = '../tmatrix/output/q_astrodust_P0.20_Fe0.00_1.400.dat'
   U_isrf        = 1.585_wp
   radiation_field = 'mathis'
   ad_stage      = 'S2'
   N_events      = 20000
   base_seed     = 12345
   out_prefix    = 'output/sed_default'
   NT_init       = 250
   T_lo          = 1.0_wp
   T_hi          = 3000.0_wp
   use_mc_for_large_grain        = .false.    ! default: calc_Teq for large grains (production)
   mc_engine                     = '2pass'    ! '2pass' (default) | 'buffered' | 'fixed'
                                              ! With all-grain MC, a fixed grid
                                              ! may resolve narrow distributions
                                              ! poorly. Such requests use '2pass'.

   nargs = command_argument_count()
   if (nargs < 1) then
      write(*,'(a)') 'Usage: main_mc_sed.x <input.nml>'
      stop 1
   end if
   call get_command_argument(1, nml_path)
   open(newunit=u, file=trim(nml_path), status='old', action='read', iostat=ios)
   if (ios /= 0) then
      write(error_unit,'(a,a)') 'Cannot open namelist: ', trim(nml_path)
      error stop 1
   end if
   read(u, nml=mc_sed_input, iostat=ios)
   if (ios /= 0) then
      write(error_unit,'(a,i0)') 'Failed to read &mc_sed_input, iostat=', ios
      error stop 1
   end if
   close(u)
   select case (trim(radiation_field))
   case ('mathis')
      use_mathis_corrected = .true.
   case ('mathis_orig')
      use_mathis_corrected = .false.
   case ('hd23')
      use_mathis_corrected = .true.
   case default
      write(*,'(a,a)') 'Invalid radiation_field: ', trim(radiation_field)
      stop 1
   end select

   ! Use an adaptive grid for all-grain MC as a conservative resolution
   ! choice. Fixed-grid adequacy depends on each recorded distribution.
   if (use_mc_for_large_grain .and. trim(mc_engine) == 'fixed') then
      write(*,'(a)') ' [warn] use_mc_for_large_grain=.true. with mc_engine="fixed"'
      write(*,'(a)') '        using adaptive grid to resolve narrow P(T); ' // &
                     'auto-promoting to "2pass".'
      mc_engine = '2pass'
   end if

   write(*,'(a)') '== mc_sed input =='
   write(*,'(a,a)') '  qtable    = ', trim(qtable_path)
   write(*,'(a,es12.4)') '  U_isrf    = ', U_isrf
   write(*,'(a,a)') '  radiation_field = ', trim(radiation_field)
   write(*,'(a,a)') '  ad_stage  = ', trim(ad_stage)
   write(*,'(a,l1,a)') '  use_mc_for_large = ', use_mc_for_large_grain, &
        merge('  (MC for large grains too)         ', &
              '  (calc_Teq for large grains, def)  ', use_mc_for_large_grain)
   write(*,'(a,a)') '  mc_engine = ', trim(mc_engine)
   write(*,'(a,i0)') '  N_events  = ', N_events
   write(*,'(a,i0)') '  base_seed = ', base_seed
   write(*,'(a,a)') '  out_prefix= ', trim(out_prefix)
   write(*,'(a,i0)') '  threads   = ', omp_get_max_threads()

   call sed_init(trim(qtable_path), NT_init, T_lo, T_hi)
   write(*,'(a,i0,a,i0)') '  loaded NLAM=', NLAM, '  NT=', NT_init

   allocate(J_lam(NLAM), lamI_lam_ad(NLAM), lamI_lam_pah(NLAM), lamI_lam_tot(NLAM))
   if (trim(radiation_field) == 'hd23') then
      call J_HD23(U_isrf, lam, J_lam)
   else
      call J_Mathis(U_isrf, lam, J_lam)
   end if
   open(newunit=u, file=trim(out_prefix)//'_field.dat', status='replace', action='write')
   write(u,'(a,a,a,es24.16e3)') '# radiation_field=', trim(radiation_field), ' U=', U_isrf
   write(u,'(a)') '# lambda[um] J_lambda[W m^-3 sr^-1]'
   do k = 1, NLAM
      write(u,'(2(es24.16e3,1x))') lam(k), J_lam(k)
   end do
   close(u)

   t0 = omp_get_wtime()
   call mc_sed_solve_total(J_lam, U_isrf, trim(ad_stage), N_events, base_seed, &
                           use_mc_for_large_grain, trim(mc_engine), &
                           lamI_lam_ad, lamI_lam_pah, lamI_lam_tot)
   t1 = omp_get_wtime()
   write(*,'(a,f8.2,a)') '  total wall: ', t1-t0, ' s'

   open(newunit=u, file=trim(out_prefix)//'_irem_mc.dat', status='replace', action='write')
   write(u,'(a,a,a,es24.16e3)') '# radiation_field=', trim(radiation_field), ' U=', U_isrf
   write(u,'(a)') '# lambda[um]  lamI_lam_total  lamI_lam_astrodust  lamI_lam_pah'
   do k = 1, NLAM
      write(u, '(es14.5e3,3(1x,es14.5e3))') lam(k), &
            lamI_lam_tot(k), lamI_lam_ad(k), lamI_lam_pah(k)
   end do
   close(u)
   write(*,'(a,a)') '  wrote ', trim(out_prefix)//'_irem_mc.dat'

   deallocate(J_lam, lamI_lam_ad, lamI_lam_pah, lamI_lam_tot)
end program main_mc_sed
