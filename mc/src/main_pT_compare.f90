program main_pT_compare
   ! Single-grain P(T) cross-check: run both Monte Carlo and the
   ! Guhathakurta-Draine matrix solver (calc_P) on the same input
   ! (Cabs, J_lam, enthalpy) for a handful of representative grain sizes,
   ! and write each grain's two distributions to a file for plotting.
   !
   ! Namelist:
   !   &pT_compare_input
   !     qtable_path   = '...'
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

   use omp_lib
   use iso_fortran_env, only: error_unit
   use ieee_arithmetic, only: ieee_is_finite
   use constants,         only: wp
   use radfield,          only: J_Mathis, J_HD23, use_mathis_corrected
   use p_sub,             only: calc_P, calc_Teq
   use mc_heatcap,        only: C_of_T
   use mc_rng,            only: rng_t, rng_init
   use mc_grain_type,     only: mc_grain_t
   use mc_engine,         only: MC_OK, hc_erg_um, grain_setup_from_cabs, mc_run_engine, &
                                mc_run_engine_2pass, mc_run_engine_buffered, &
                                grain_radiative_equilibrium, mc_equilibrate
   use sed_astrodust_mod, only: sed_init, NLAM, NA, lam, aeff, T_first, &
                                Cabs, Cabs_pah, kappB_first, kappB_pah_first, &
                                H_first, H_pah_first
   use enthalpy_astrodust_mod, only: s1_density_corrected
   implicit none
   integer :: engine_status
   real(wp) :: omitted_time, floor_time, outside_time
   integer(kind=8) :: sample_count, samples_replay

   integer, parameter :: NHIST_LOC = 600
   integer, parameter :: MAX_GRAINS = 8

   character(len=512) :: qtable_path, out_prefix
   character(len=16)  :: ad_stage, pop, radiation_field
   real(wp) :: U_isrf, T_lo, T_hi
   real(wp) :: a_um_list(MAX_GRAINS)
   integer  :: N_events, base_seed, NT_init
   namelist /pT_compare_input/ qtable_path, U_isrf, ad_stage, &
                               pop, a_um_list, N_events, base_seed, &
                               out_prefix, NT_init, T_lo, T_hi, radiation_field

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
   real(wp) :: t_total, e_abs_tot, e_emit_tot, Teq_grain
   real(wp) :: T_eq_mc, T_record, p_abs, tau_relax, warm_time
   integer :: warm_events
   real(wp) :: t_total_A, e_abs_A, e_emit_A, Tmin_A, Tmax_A
   real(wp) :: t_total_B, e_abs_B, e_emit_B, Tmin_B, Tmax_B
   logical  :: is_log_A, is_log_B, buf_full_B
   integer  :: n_samp_B
   type(mc_grain_t) :: grain
   type(rng_t)      :: rng, rng_record
   real(wp) :: t0, t1, t_A0, t_A1, t_B0, t_B1
   integer  :: max_samp_B

   ! defaults
   qtable_path   = '../tmatrix/output/q_astrodust_P0.20_Fe0.00_1.400.dat'
   U_isrf        = 1.585_wp
   radiation_field = 'mathis'
   ad_stage      = 'S2'
   pop           = 'astrodust'
   ! Unspecified namelist elements retain this sentinel, so a short list
   ! cannot silently inherit the remaining default grain sizes.
   a_um_list     = -huge(1.0_wp)
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
   if (ios /= 0) then
      write(error_unit,'(a,a)') 'Cannot open namelist: ', trim(nml_path)
      error stop 1
   end if
   read(u, nml=pT_compare_input, iostat=ios)
   if (ios /= 0) then
      write(error_unit,'(a,i0)') 'Failed to read &pT_compare_input, iostat=', ios
      stop 1
   end if
   close(u)
   if (all(a_um_list == -huge(1.0_wp))) then
      a_um_list = 0.0_wp
      a_um_list(1:4) = [0.001_wp, 0.003_wp, 0.01_wp, 0.05_wp]
   end if
   if (.not. ieee_is_finite(U_isrf) .or. U_isrf < 0.0_wp) error stop 'U_isrf must be finite and nonnegative'
   if (N_events <= 0) error stop 'N_events must be positive'
   if (NT_init < 2) error stop 'NT_init must be at least 2'
   if (.not. ieee_is_finite(T_lo) .or. T_lo <= 0.0_wp) error stop 'T_lo must be finite and positive'
   if (.not. ieee_is_finite(T_hi) .or. T_hi <= T_lo) error stop 'T_hi must be finite and greater than T_lo'
   if (len_trim(qtable_path) == 0) error stop 'qtable_path must not be empty'
   if (len_trim(out_prefix) == 0) error stop 'out_prefix must not be empty'
   select case (trim(pop))
   case ('astrodust', 'pah')
   case default
      write(error_unit,'(a,a)') 'Invalid pop: ', trim(pop)
      error stop 1
   end select
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

   ! Compact the specified positive radii in their namelist order.
   n_grains = 0
   do k = 1, MAX_GRAINS
      if (a_um_list(k) == -huge(1.0_wp)) cycle
      if (.not. ieee_is_finite(a_um_list(k)) .or. a_um_list(k) < 0.0_wp) &
         error stop 'a_um_list must contain finite, nonnegative radii'
      if (a_um_list(k) == 0.0_wp) cycle
      n_grains = n_grains + 1
      a_um_list(n_grains) = a_um_list(k)
   end do
   if (n_grains == 0) error stop 'a_um_list must contain a positive radius'
   write(*,'(a,i0,a)') '== pT_compare input ==  (n_grains=', n_grains, ')'
   write(*,'(a,a)')        '  pop          = ', trim(pop)
   write(*,'(a,a)')        '  ad_stage     = ', trim(ad_stage)
   write(*,'(a,es12.4)')   '  U_isrf       = ', U_isrf
   write(*,'(a,a)')        '  radiation_field = ', trim(radiation_field)
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

   ! Load the optics, size grid and enthalpies (sed_init)
   call sed_init(trim(qtable_path), NT_init, T_lo, T_hi)
   write(*,'(a,i0,a,i0,a,i0)') '  loaded NLAM=', NLAM, '  NA=', NA, '  NT=', NT_init

   allocate(J_lam(NLAM), Cabs_grain(NLAM), kappB_grain(NT_init), H_grain(NT_init))
   allocate(P_gd(NT_init), lnP_gd(NT_init))
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

      ! Pull this grain's Cabs, kappB, and H depending on population.
      ! pop = 'pah' takes the f_ion-blended carbonaceous grain (Cabs_pah with
      ! kappB_pah_first built on that same Cabs_pah), so both solvers see
      ! one consistent grain. It is not a population
      ! sed_solve_pah solves -- that one runs the neutral and the cation
      ! separately -- but this program compares two P(T) solvers on the same
      ! input, which any self-consistent grain serves.
      if (trim(pop) == 'pah') then
         Cabs_grain  = Cabs_pah(:, ja_best)
         kappB_grain = kappB_pah_first(:, ja_best)
         H_grain     = H_pah_first(:, ja_best)
      else
         Cabs_grain  = Cabs(:, ja_best)
         kappB_grain = kappB_first(:, ja_best)
         H_grain     = H_first(:, ja_best, is_stage)
      end if

      write(*,'(a,i2,a,f8.4,a,i4,a,f8.4,a)') &
         '  grain', jg, ': target a=', aval, &
         ' um  closest jaeff=', ja_best, ' (', aeff(ja_best), ' um)'

      ! --- Equilibrium temperature (for the side-by-side comparison) ---
      call calc_Teq(lam, Cabs_grain, J_lam, T_first, kappB_grain, Teq_grain)

      ! --- Guhathakurta-Draine calc_P ---
      ! The constructed field includes the CMB, so no cooling offset is needed.
      call calc_P(lam, Cabs_grain, J_lam, T_first, kappB_grain, H_grain, &
                  P_gd, lnP_gd)

      ! --- Monte Carlo (3 engines, same seed) ---

      ! Adaptive DA85 cutoff (same prescription as mc_sed_loop):
      ! lam_c = hc / (DT_THRESHOLD * Teq * C(Teq)), so photons whose
      ! jump dT = (hc/lam)/C(Teq) > DT_THRESHOLD * Teq per event remain
      ! stochastic, the rest fold into H_cont.  This keeps the trajectory
      ! near Teq for very large grains and prevents the CMB collapse seen
      ! with a static lam_c=1000 um.
      block
         real(wp), parameter :: DT_THRESHOLD = 0.01_wp
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
              lam_c_in=lam_c_adapt, J_lam_in=J_lam)
      end block
      call grain_radiative_equilibrium(grain, T_eq_mc, p_abs, tau_relax, engine_status)
      if (engine_status /= MC_OK) error stop 'Radiative equilibrium failed'
      call mc_equilibrate(grain, rng, T_eq_mc, T_record, warm_time, warm_events, engine_status)
      if (engine_status /= MC_OK) error stop 'MC warmup failed'
      write(*,'(a,es24.16e3,a,es24.16e3,a,es24.16e3,a,es24.16e3,a,i0)') &
         '[mc_warmup] radius_um=', aeff(ja_best), ' initial_Teq_K=', T_eq_mc, &
         ' record_T_K=', T_record, ' warm_time_s=', warm_time, ' events=', warm_events
      rng_record = rng
      t0 = omp_get_wtime()
      call mc_run_engine(grain, rng, T_record, N_events, 0.0_wp, &
                         NHIST_LOC, T_edges, dP_dT, dP_dlnT, &
                         t_total, e_abs_tot, e_emit_tot, &
                         status_out=engine_status, omitted_time_out=omitted_time, &
                               floor_time_out=floor_time, out_of_grid_time_out=outside_time, &
                               sample_count_out=sample_count)
      if (engine_status /= MC_OK) then
         write(error_unit,'(a,i0,a,es24.16e3)') 'MC cooling failed: status=', engine_status, &
            ' incomplete_time_s=', omitted_time
         error stop 'MC interval incomplete; no grain comparison is written'
      end if
      write(*,'(a,3(es24.16e3,1x),i0)') '[mc_interval fixed] floor_s outside_grid_s omitted_s samples: ', &
         floor_time, outside_time, omitted_time, sample_count
      t1 = omp_get_wtime()

      ! engine 2: 2pass (Method A)
      rng = rng_record
      t_A0 = omp_get_wtime()
      call mc_run_engine_2pass(grain, rng, T_record, N_events, -1, 0.0_wp, &
                               NHIST_LOC, T_edges_A, dP_dT_A, dP_dlnT_A, &
                               t_total_A, e_abs_A, e_emit_A, &
                               is_log_out=is_log_A, &
                               T_min_obs_out=Tmin_A, T_max_obs_out=Tmax_A, &
                               status_out=engine_status, omitted_time_out=omitted_time, &
                               floor_time_out=floor_time, out_of_grid_time_out=outside_time, &
                               sample_count_out=sample_count)
      if (engine_status /= MC_OK) then
         write(error_unit,'(a,i0,a,es24.16e3)') 'MC cooling failed: status=', engine_status, &
            ' incomplete_time_s=', omitted_time
         error stop 'MC interval incomplete; no grain comparison is written'
      end if
      write(*,'(a,3(es24.16e3,1x),i0)') '[mc_interval 2pass] floor_s outside_grid_s omitted_s samples: ', &
         floor_time, outside_time, omitted_time, sample_count
      t_A1 = omp_get_wtime()

      ! The identical-seed replay measures the full sample capacity needed
      ! by the buffered calculation. At most 16 GB stores two real arrays.
      samples_replay = sample_count
      if (samples_replay <= 0_8 .or. samples_replay > 1000000000_8) &
         error stop 'Measured MC sample capacity exceeds the 16 GB comparison limit'
      max_samp_B = int(samples_replay)
      write(*,'(a,i0)') '    measured buffered capacity = ', max_samp_B

      ! engine 3: buffered (Method B)
      rng = rng_record
      t_B0 = omp_get_wtime()
      call mc_run_engine_buffered(grain, rng, T_record, N_events, 0.0_wp, &
                                  NHIST_LOC, max_samp_B, &
                                  T_edges_B, dP_dT_B, dP_dlnT_B, &
                                  t_total_B, e_abs_B, e_emit_B, &
                                  is_log_out=is_log_B, &
                                  T_min_obs_out=Tmin_B, T_max_obs_out=Tmax_B, &
                                  n_samples_used=n_samp_B, &
                                  buffer_full_out=buf_full_B, &
                                  status_out=engine_status, omitted_time_out=omitted_time, &
                               floor_time_out=floor_time, out_of_grid_time_out=outside_time, &
                               sample_count_out=sample_count)
      if (engine_status /= MC_OK) then
         write(error_unit,'(a,i0,a,es24.16e3)') 'MC cooling failed: status=', engine_status, &
            ' incomplete_time_s=', omitted_time
         error stop 'MC interval incomplete; no grain comparison is written'
      end if
      write(*,'(a,3(es24.16e3,1x),i0)') '[mc_interval buffered] floor_s outside_grid_s omitted_s samples: ', &
         floor_time, outside_time, omitted_time, sample_count
      if (sample_count /= samples_replay) error stop 'Buffered sample count differs from identical-seed replay'
      t_B1 = omp_get_wtime()

      write(*,'(a,f7.2,a,f7.2,a,f7.2,a)') &
         '    walls (fixed / 2pass / buf) = ', t1-t0, ' / ', t_A1-t_A0, &
         ' / ', t_B1-t_B0, ' s'
      write(*,'(a,es14.6e3,a,es14.6e3,a,es14.6e3)') &
         '    emit/event_abs = ', e_emit_tot/max(e_abs_tot,1.0e-300_wp), &
         ' (fixed), ',     e_emit_A  /max(e_abs_A,  1.0e-300_wp), &
         ' (2pass), ',     e_emit_B  /max(e_abs_B,  1.0e-300_wp)
      write(*,'(a,l1,a,f7.2,a,f7.2,a,a,l1,a,i0,a,l1,a)') &
         '    2pass grid: log=', is_log_A, ' [', Tmin_A, ',', Tmax_A, ' K]', &
         '   buf grid: log=', is_log_B, ' (n_samp=', n_samp_B, ', full=', buf_full_B, ')'

      ! Write file: 4 blocks
      write(fname,'(a,a,i2.2,a,i6.6,a)') trim(out_prefix), '_', jg, '_a', &
            nint(aeff(ja_best)*1.0e4_wp), 'A.dat'
      open(newunit=u, file=trim(adjustl(fname)), status='replace', action='write')
      write(u,'(a,a,a,es24.16e3)') '# radiation_field=', trim(radiation_field), ' U=', U_isrf
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
      if (allocated(grain%cdf_density)) deallocate(grain%cdf_density)
      if (allocated(grain%T_grid))     deallocate(grain%T_grid)
      if (allocated(grain%U_grid))     deallocate(grain%U_grid)
      if (allocated(grain%QPL_grid))   deallocate(grain%QPL_grid)
   end do

   deallocate(J_lam, Cabs_grain, kappB_grain, H_grain, P_gd, lnP_gd)
end program main_pT_compare
