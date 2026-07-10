program main_mc_sed
   ! MC SED builder driver.  Mirrors main_astrodust.x but uses the MC
   ! engine in place of the Guhathakurta-Draine matrix solver.
   !
   ! Inputs via namelist (positional arg 1):
   !
   !   &mc_sed_input
   !     qtable_path  = '../tmatrix/output/q_astrodust_P0.20_Fe0.00_1.400.dat'
   !     sizedist_path = '../data/size_distribution.dat'
   !     U_isrf       = 1.585           ! log U = 0.20
   !     ad_stage     = 'S2'            ! 'S1_C1' | 'S1_C2' | 'S2'
   !     N_events     = 20000
   !     base_seed    = 12345
   !     out_prefix   = 'output/mc_sed_default'
   !     NT_init      = 250
   !     T_lo         = 1.0
   !     T_hi         = 3000.0
   !   /

   use, intrinsic :: omp_lib
   use constants,         only: wp
   use radfield,          only: J_Mathis
   use sed_astrodust_mod, only: sed_init, NLAM, lam
   use mc_sed,            only: mc_sed_solve_total
   implicit none

   character(len=512) :: qtable_path, sizedist_path, out_prefix
   character(len=16)  :: ad_stage, mc_engine
   real(wp) :: U_isrf, T_lo, T_hi
   integer  :: N_events, base_seed, NT_init
   logical  :: use_mc_for_large_grain
   namelist /mc_sed_input/ qtable_path, sizedist_path, U_isrf, ad_stage, &
                           N_events, base_seed, out_prefix, NT_init, T_lo, T_hi, &
                           use_mc_for_large_grain, mc_engine

   character(len=512) :: nml_path
   integer :: nargs, u, ios, k
   real(wp), allocatable :: J_lam(:), lamI_lam_ad(:), lamI_lam_pah(:), lamI_lam_tot(:)
   real(wp) :: t0, t1

   ! defaults
   qtable_path   = '../tmatrix/output/q_astrodust_P0.20_Fe0.00_1.400.dat'
   sizedist_path = '../data/size_distribution.dat'
   U_isrf        = 1.585_wp
   ad_stage      = 'S2'
   N_events      = 20000
   base_seed     = 12345
   out_prefix    = 'output/mc_sed_default'
   NT_init       = 250
   T_lo          = 1.0_wp
   T_hi          = 3000.0_wp
   use_mc_for_large_grain        = .false.    ! default: calc_Teq for large grains (production)
   mc_engine                     = '2pass'    ! '2pass' (default) | 'buffered' | 'fixed'
                                              ! Fixed grid is unsafe when
                                              ! use_mc_for_large_grain=.true.
                                              ! (P(T) collapses to one bin for
                                              ! a >~ 80 A); auto-promoted to
                                              ! '2pass' below if 'fixed' is
                                              ! requested with large-grain MC.

   nargs = command_argument_count()
   if (nargs < 1) then
      write(*,'(a)') 'Usage: main_mc_sed.x <input.nml>'
      stop 1
   end if
   call get_command_argument(1, nml_path)
   open(newunit=u, file=trim(nml_path), status='old', action='read', iostat=ios)
   if (ios /= 0) stop 'Cannot open namelist'
   read(u, nml=mc_sed_input, iostat=ios)
   if (ios /= 0) then
      write(*,'(a,i0)') 'Failed to read &mc_sed_input, iostat=', ios
      stop 1
   end if
   close(u)

   ! Safety: 'fixed' grid degenerates for equilibrium-regime grains
   ! (single-bin P(T) for a >~ 80 A), so it is unsafe when MC is forced on
   ! large grains.  Auto-promote to '2pass' with a warning.
   if (use_mc_for_large_grain .and. trim(mc_engine) == 'fixed') then
      write(*,'(a)') ' [warn] use_mc_for_large_grain=.true. with mc_engine="fixed"'
      write(*,'(a)') '        is unsafe (single-bin P(T) for large grains); ' // &
                     'auto-promoting to "2pass".'
      mc_engine = '2pass'
   end if

   write(*,'(a)') '== mc_sed input =='
   write(*,'(a,a)') '  qtable    = ', trim(qtable_path)
   write(*,'(a,a)') '  sizedist  = ', trim(sizedist_path)
   write(*,'(a,es12.4)') '  U_isrf    = ', U_isrf
   write(*,'(a,a)') '  ad_stage  = ', trim(ad_stage)
   write(*,'(a,l1,a)') '  use_mc_for_large = ', use_mc_for_large_grain, &
        merge('  (MC for large grains too)         ', &
              '  (calc_Teq for large grains, def)  ', use_mc_for_large_grain)
   write(*,'(a,a)') '  mc_engine = ', trim(mc_engine)
   write(*,'(a,i0)') '  N_events  = ', N_events
   write(*,'(a,i0)') '  base_seed = ', base_seed
   write(*,'(a,a)') '  out_prefix= ', trim(out_prefix)
   write(*,'(a,i0)') '  threads   = ', omp_get_max_threads()

   call sed_init(trim(qtable_path), trim(sizedist_path), NT_init, T_lo, T_hi)
   write(*,'(a,i0,a,i0)') '  loaded NLAM=', NLAM, '  NT=', NT_init

   allocate(J_lam(NLAM), lamI_lam_ad(NLAM), lamI_lam_pah(NLAM), lamI_lam_tot(NLAM))
   call J_Mathis(U_isrf, lam, J_lam)

   t0 = omp_get_wtime()
   call mc_sed_solve_total(J_lam, U_isrf, trim(ad_stage), N_events, base_seed, &
                           use_mc_for_large_grain, trim(mc_engine), &
                           lamI_lam_ad, lamI_lam_pah, lamI_lam_tot)
   t1 = omp_get_wtime()
   write(*,'(a,f8.2,a)') '  total wall: ', t1-t0, ' s'

   open(newunit=u, file=trim(out_prefix)//'_irem_mc.dat', status='replace', action='write')
   write(u,'(a)') '# lambda[um]  lamI_lam_total  lamI_lam_astrodust  lamI_lam_pah'
   do k = 1, NLAM
      write(u, '(es14.5e3,3(1x,es14.5e3))') lam(k), &
            lamI_lam_tot(k), lamI_lam_ad(k), lamI_lam_pah(k)
   end do
   close(u)
   write(*,'(a,a)') '  wrote ', trim(out_prefix)//'_irem_mc.dat'

   deallocate(J_lam, lamI_lam_ad, lamI_lam_pah, lamI_lam_tot)
end program main_mc_sed
