program main_pT_compare
   ! Single-grain P(T) cross-check: run both Monte Carlo and the
   ! Guhathakurta-Draine matrix solver (calc_P) on the same input
   ! (Cabs, J_lam, enthalpy) for a handful of representative grain sizes,
   ! and write each grain's two distributions to a file for plotting.
   !
   ! Namelist:
   !   &pT_compare_input
   !     qtable_path   = '...'
   !     sizedist_path = '...'
   !     U_isrf        = 1.585
   !     ad_stage      = 'S2'
   !     pop           = 'astrodust'   ! 'astrodust' or 'pah'
   !     a_um_list     = 0.001, 0.003, 0.01, 0.05  ! 4 grain sizes
   !     N_events      = 200000
   !     base_seed     = 12345
   !     out_prefix    = 'output/pT_compare'
   !     NT_init       = 250
   !     T_lo          = 1.0
   !     T_hi          = 3000.0
   !   /

   use, intrinsic :: omp_lib
   use constants,         only: wp
   use radfield,          only: J_Mathis
   use p_sub,             only: calc_P, calc_Teq
   use mc_heatcap,        only: C_of_T
   use mc_rng,            only: rng_t, rng_init
   use mc_grain_type,     only: mc_grain_t
   use mc_engine,         only: grain_setup_from_cabs, mc_run_engine, &
                                mc_run_engine_2pass, mc_run_engine_buffered
   use sed_astrodust_mod, only: sed_init, NLAM, NA, lam, aeff, T_first, &
                                Cabs, Cabs_pah, kappB_first, kappB_pah_first, &
                                H_first, H_pah_first, kappCMB, kappCMB_pah
   use enthalpy_astrodust_mod, only: s1_density_corrected
   implicit none

   integer, parameter :: NHIST_LOC = 600
   integer, parameter :: MAX_GRAINS = 8

   character(len=512) :: qtable_path, sizedist_path, out_prefix
   character(len=16)  :: ad_stage, pop
   real(wp) :: U_isrf, T_lo, T_hi
   real(wp) :: a_um_list(MAX_GRAINS)
   integer  :: N_events, base_seed, NT_init
   namelist /pT_compare_input/ qtable_path, sizedist_path, U_isrf, ad_stage, &
                               pop, a_um_list, N_events, base_seed, &
                               out_prefix, NT_init, T_lo, T_hi

   character(len=512) :: nml_path, fname
   character(len=16)  :: comp_for_mc
   integer :: nargs, u, ios, k, jg, ja_best, is_stage
   integer :: n_grains
   real(wp) :: aval, best_dist
   real(wp), allocatable :: J_lam(:)
   real(wp), allocatable :: Cabs_grain(:), kappB_grain(:), H_grain(:)
   real(wp), allocatable :: P_gd(:), lnP_gd(:)
   real(wp) :: T_edges(NHIST_LOC+1), dP_dT(NHIST_LOC), dP_dlnT(NHIST_LOC)
   real(wp) :: T_edges_A(NHIST_LOC+1), dP_dT_A(NHIST_LOC), dP_dlnT_A(NHIST_LOC)
   real(wp) :: T_edges_B(NHIST_LOC+1), dP_dT_B(NHIST_LOC), dP_dlnT_B(NHIST_LOC)
   real(wp) :: t_total, e_abs_tot, e_emit_tot, kappCMB_val, Teq_grain
   real(wp) :: t_total_A, e_abs_A, e_emit_A, Tmin_A, Tmax_A
   real(wp) :: t_total_B, e_abs_B, e_emit_B, Tmin_B, Tmax_B
   logical  :: is_log_A, is_log_B, buf_full_B
   integer  :: n_samp_B
   type(mc_grain_t) :: grain
   type(rng_t)      :: rng
   real(wp) :: t0, t1, t_A0, t_A1, t_B0, t_B1
   integer  :: N_burn, max_samp_B

   ! defaults
   qtable_path   = '../tmatrix/output/q_astrodust_P0.20_Fe0.00_1.400.dat'
   sizedist_path = '../data/release/size_distribution.dat'
   U_isrf        = 1.585_wp
   ad_stage      = 'S2'
   pop           = 'astrodust'
   a_um_list     = 0.0_wp
   a_um_list(1)  = 0.001_wp
   a_um_list(2)  = 0.003_wp
   a_um_list(3)  = 0.01_wp
   a_um_list(4)  = 0.05_wp
   N_events      = 200000
   base_seed     = 12345
   out_prefix    = 'output/pT_compare'
   NT_init       = 250
   T_lo          = 1.0_wp
   T_hi          = 3000.0_wp

   nargs = command_argument_count()
   if (nargs < 1) then
      write(*,'(a)') 'Usage: main_pT_compare.x <input.nml>'
      stop 1
   end if
   call get_command_argument(1, nml_path)
   open(newunit=u, file=trim(nml_path), status='old', action='read', iostat=ios)
   if (ios /= 0) stop 'Cannot open namelist'
   read(u, nml=pT_compare_input, iostat=ios)
   if (ios /= 0) then
      write(*,'(a,i0)') 'Failed to read &pT_compare_input, iostat=', ios
      stop 1
   end if
   close(u)

   ! Count non-zero entries in a_um_list
   n_grains = 0
   do k = 1, MAX_GRAINS
      if (a_um_list(k) > 0.0_wp) n_grains = n_grains + 1
   end do
   write(*,'(a,i0,a)') '== pT_compare input ==  (n_grains=', n_grains, ')'
   write(*,'(a,a)')        '  pop          = ', trim(pop)
   write(*,'(a,a)')        '  ad_stage     = ', trim(ad_stage)
   write(*,'(a,es12.4)')   '  U_isrf       = ', U_isrf
   write(*,'(a,i0)')       '  N_events     = ', N_events
   write(*,'(a,i0)')       '  base_seed    = ', base_seed
   write(*,'(a,a)')        '  out_prefix   = ', trim(out_prefix)

   ! Map ad_stage -> H_first stage index and mc comp tag.
   ! H_first now has only 2 stages: (1) Stage 1 silicate-only, (2) Stage 2
   ! silicate+carbonaceous.  The Stage-1 C1/C2 prefactor is selected by the
   ! s1_density_corrected toggle (set before sed_init builds H_first), not by
   ! a separate stage index.
   select case (trim(ad_stage))
   case ('S1_C1'); is_stage = 1; comp_for_mc = 'ad_s1_c1'; s1_density_corrected = .false.
   case ('S1_C2'); is_stage = 1; comp_for_mc = 'ad_s1_c2'; s1_density_corrected = .true.
   case ('S2');    is_stage = 2; comp_for_mc = 'ad_s2'
   case default; stop 'unknown ad_stage'
   end select
   if (trim(pop) == 'pah') comp_for_mc = 'pah'

   ! Load all pipeline data
   call sed_init(trim(qtable_path), trim(sizedist_path), NT_init, T_lo, T_hi)
   write(*,'(a,i0,a,i0,a,i0)') '  loaded NLAM=', NLAM, '  NA=', NA, '  NT=', NT_init

   allocate(J_lam(NLAM), Cabs_grain(NLAM), kappB_grain(NT_init), H_grain(NT_init))
   allocate(P_gd(NT_init), lnP_gd(NT_init))
   call J_Mathis(U_isrf, lam, J_lam)

   do jg = 1, n_grains
      aval = a_um_list(jg)

      ! Find closest aeff index (NA grid)
      ja_best = 1
      best_dist = abs(log(aeff(1)) - log(aval))
      do k = 2, NA
         if (abs(log(aeff(k)) - log(aval)) < best_dist) then
            ja_best = k
            best_dist = abs(log(aeff(k)) - log(aval))
         end if
      end do

      ! Pull this grain's Cabs, kappB, H, kappCMB depending on population
      if (trim(pop) == 'pah') then
         Cabs_grain  = Cabs_pah(:, ja_best)
         kappB_grain = kappB_pah_first(:, ja_best)
         H_grain     = H_pah_first(:, ja_best)
         kappCMB_val = kappCMB_pah(ja_best)
      else
         Cabs_grain  = Cabs(:, ja_best)
         kappB_grain = kappB_first(:, ja_best)
         H_grain     = H_first(:, ja_best, is_stage)
         kappCMB_val = kappCMB(ja_best)
      end if

      write(*,'(a,i2,a,f8.4,a,i4,a,f8.4,a)') &
         '  grain', jg, ': target a=', aval, &
         ' um  closest jaeff=', ja_best, ' (', aeff(ja_best), ' um)'

      ! --- Equilibrium temperature (for the side-by-side comparison) ---
      call calc_Teq(lam, Cabs_grain, J_lam, T_first, kappB_grain, Teq_grain)

      ! --- Guhathakurta-Draine calc_P ---
      call calc_P(lam, Cabs_grain, J_lam, T_first, kappB_grain, H_grain, &
                  P_gd, lnP_gd, kappCMB_val)

      ! --- Monte Carlo (3 engines, same seed) ---
      N_burn     = max(50, N_events / 20)
      max_samp_B = min(200000, 8 * N_events)

      ! Adaptive DA85 cutoff (same prescription as mc_sed_loop):
      ! lam_c = hc / (DT_THRESHOLD * Teq * C(Teq)), so photons whose
      ! jump dT = (hc/lam)/C(Teq) > DT_THRESHOLD * Teq per event remain
      ! stochastic, the rest fold into H_cont.  This keeps the trajectory
      ! near Teq for very large grains and prevents the CMB collapse seen
      ! with a static lam_c=1000 um.
      block
         real(wp), parameter :: DT_THRESHOLD = 0.01_wp
         real(wp), parameter :: HC_ERG_UM    = 1.98644582e-13_wp
         real(wp) :: C_Teq, lam_c_adapt
         C_Teq = C_of_T(Teq_grain, aeff(ja_best), comp_for_mc)
         if (C_Teq > 0.0_wp .and. Teq_grain > 0.0_wp) then
            lam_c_adapt = HC_ERG_UM / (DT_THRESHOLD * Teq_grain * C_Teq)
            if (lam_c_adapt < lam(2))      lam_c_adapt = lam(2)
            if (lam_c_adapt > lam(NLAM-1)) lam_c_adapt = lam(NLAM-1)
         else
            lam_c_adapt = 1000.0_wp
         end if
         write(*,'(a,es10.3,a)') '    adaptive lam_c = ', lam_c_adapt, ' um'

         ! engine 1: fixed-grid (original)
         call rng_init(rng, base_seed + 1009 * jg)
         call grain_setup_from_cabs(grain, aeff(ja_best), U_isrf, comp_for_mc, lam, &
              Cabs_grain / (3.141592653589793238_wp * (aeff(ja_best)*1.0e-4_wp)**2), &
              lam_c_in=lam_c_adapt)
      end block
      t0 = omp_get_wtime()
      call mc_run_engine(grain, rng, 2.725_wp, N_events, 0.0_wp, &
                         NHIST_LOC, T_edges, dP_dT, dP_dlnT, &
                         t_total, e_abs_tot, e_emit_tot)
      t1 = omp_get_wtime()

      ! engine 2: 2pass (Method A)
      call rng_init(rng, base_seed + 1009 * jg)
      t_A0 = omp_get_wtime()
      call mc_run_engine_2pass(grain, rng, 2.725_wp, N_events, -1, N_burn, 0.0_wp, &
                               NHIST_LOC, T_edges_A, dP_dT_A, dP_dlnT_A, &
                               t_total_A, e_abs_A, e_emit_A, &
                               is_log_out=is_log_A, &
                               T_min_obs_out=Tmin_A, T_max_obs_out=Tmax_A)
      t_A1 = omp_get_wtime()

      ! engine 3: buffered (Method B)
      call rng_init(rng, base_seed + 1009 * jg)
      t_B0 = omp_get_wtime()
      call mc_run_engine_buffered(grain, rng, 2.725_wp, N_events, N_burn, 0.0_wp, &
                                  NHIST_LOC, max_samp_B, &
                                  T_edges_B, dP_dT_B, dP_dlnT_B, &
                                  t_total_B, e_abs_B, e_emit_B, &
                                  is_log_out=is_log_B, &
                                  T_min_obs_out=Tmin_B, T_max_obs_out=Tmax_B, &
                                  n_samples_used=n_samp_B, &
                                  buffer_full_out=buf_full_B)
      t_B1 = omp_get_wtime()

      write(*,'(a,f7.2,a,f7.2,a,f7.2,a)') &
         '    walls (fixed / 2pass / buf) = ', t1-t0, ' / ', t_A1-t_A0, &
         ' / ', t_B1-t_B0, ' s'
      write(*,'(a,f7.4,a,f7.4,a,f7.4)') &
         '    emit/abs = ', e_emit_tot/max(e_abs_tot,1.0e-300_wp), &
         ' (fixed), ',     e_emit_A  /max(e_abs_A,  1.0e-300_wp), &
         ' (2pass), ',     e_emit_B  /max(e_abs_B,  1.0e-300_wp)
      write(*,'(a,l1,a,f7.2,a,f7.2,a,a,l1,a,i0,a,l1,a)') &
         '    2pass grid: log=', is_log_A, ' [', Tmin_A, ',', Tmax_A, ' K]', &
         '   buf grid: log=', is_log_B, ' (n_samp=', n_samp_B, ', full=', buf_full_B, ')'

      ! Write file: 4 blocks
      write(fname,'(a,a,i2.2,a,i6.6,a)') trim(out_prefix), '_', jg, '_a', &
            nint(aeff(ja_best)*1.0e4_wp), 'A.dat'
      open(newunit=u, file=trim(adjustl(fname)), status='replace', action='write')
      write(u,'(a,f10.5,a,a,a,f10.5,a)') '# aeff[um]=', aeff(ja_best), &
           '  comp=', trim(comp_for_mc), '  Teq[K]=', Teq_grain, ' (calc_Teq output)'
      write(u,'(a,i0,a,i0)') '# block 1 (GD): N=', NT_init, &
                             '   blocks 2-4 (MC fixed/2pass/buffered): N=', NHIST_LOC
      write(u,'(a)') '# === GD block: T[K]   P(T_i)   dP/dlnT'
      do k = 1, NT_init
         write(u,'(3(es14.5e3,1x))') T_first(k), P_gd(k), &
               P_gd(k) / max(log(T_first(min(k+1,NT_init))/T_first(max(k-1,1))) * 0.5_wp, 1.0e-30_wp)
      end do
      write(u,'(a)') '# === MC fixed-grid block: T_mid[K]   dP/dT [1/K]   dP/dlnT'
      do k = 1, NHIST_LOC
         write(u,'(3(es14.5e3,1x))') sqrt(T_edges(k)*T_edges(k+1)), &
                                      dP_dT(k), dP_dlnT(k)
      end do
      write(u,'(a,l1,a,es12.4,a,es12.4,a)') &
         '# === MC 2pass block (log=', is_log_A, ', range=[', Tmin_A, ',', Tmax_A, ' K])'
      write(u,'(a)') '#     T_mid[K]   dP/dT [1/K]   dP/dlnT'
      do k = 1, NHIST_LOC
         if (is_log_A) then
            write(u,'(3(es14.5e3,1x))') sqrt(T_edges_A(k)*T_edges_A(k+1)), &
                                         dP_dT_A(k), dP_dlnT_A(k)
         else
            write(u,'(3(es14.5e3,1x))') 0.5_wp*(T_edges_A(k)+T_edges_A(k+1)), &
                                         dP_dT_A(k), dP_dlnT_A(k)
         end if
      end do
      write(u,'(a,l1,a,es12.4,a,es12.4,a,i0,a,l1,a)') &
         '# === MC buffered block (log=', is_log_B, ', range=[', Tmin_B, ',', Tmax_B, &
         ' K], n_samp=', n_samp_B, ', buf_full=', buf_full_B, ')'
      write(u,'(a)') '#     T_mid[K]   dP/dT [1/K]   dP/dlnT'
      do k = 1, NHIST_LOC
         if (is_log_B) then
            write(u,'(3(es14.5e3,1x))') sqrt(T_edges_B(k)*T_edges_B(k+1)), &
                                         dP_dT_B(k), dP_dlnT_B(k)
         else
            write(u,'(3(es14.5e3,1x))') 0.5_wp*(T_edges_B(k)+T_edges_B(k+1)), &
                                         dP_dT_B(k), dP_dlnT_B(k)
         end if
      end do
      close(u)
      write(*,'(a,a)') '    wrote ', trim(adjustl(fname))

      ! Free grain
      if (allocated(grain%lam_grid))   deallocate(grain%lam_grid)
      if (allocated(grain%Q_grid))     deallocate(grain%Q_grid)
      if (allocated(grain%u_lam_grid)) deallocate(grain%u_lam_grid)
      if (allocated(grain%cdf_F))      deallocate(grain%cdf_F)
      if (allocated(grain%cdf_lam))    deallocate(grain%cdf_lam)
      if (allocated(grain%T_grid))     deallocate(grain%T_grid)
      if (allocated(grain%U_grid))     deallocate(grain%U_grid)
      if (allocated(grain%QPL_grid))   deallocate(grain%QPL_grid)
   end do

   deallocate(J_lam, Cabs_grain, kappB_grain, H_grain, P_gd, lnP_gd)
end program main_pT_compare
