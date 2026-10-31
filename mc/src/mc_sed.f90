module mc_sed
   ! MC-based SED builder.  Mirrors the structure of
   ! sed_astrodust_mod::sed_solve but replaces the calc_P matrix solver
   ! with the Monte Carlo engine.
   !
   ! Two populations are computed:
   !   - Astrodust grains    (Cabs, dn_ad, comp = 'ad_s1_c1'|'ad_s1_c2'|'ad_s2')
   !   - PAH grains          (Cabs_cneu/kappB_cneu, dn_cneu and
   !                          Cabs_cion/kappB_cion, dn_cion, comp = 'pah'),
   !                          the two charge states solved separately and
   !                          summed, as sed_solve_pah does
   !
   ! The size loop is OpenMP-parallel.  Each iteration creates its own
   ! grain_t and rng_t (both stack-local in the parallel region), runs
   ! MC, and atomically accumulates into the shared Jout(:).
   !
   ! Output convention matches sed_astrodust_mod::sed_solve:
   !   lamI_lam_out(:) = lambda * Jout * 1e-3
   ! so that the result can be directly compared with that subroutine.

   use iso_fortran_env,       only: error_unit
   use constants,             only: wp, pi, ev2erg, h_planck, k_boltz, c
   use radfield,              only: bbody
   use p_sub,                 only: calc_Teq
   use mc_heatcap,            only: U_of_T, C_of_T
   use mc_grain_type,         only: mc_grain_t
   use mc_rng,                only: rng_t, rng_init
   use mc_engine,             only: grain_setup_from_cabs, mc_run_engine, &
                                    mc_run_engine_2pass, mc_run_engine_buffered, &
                                    grain_radiative_equilibrium, mc_equilibrate, &
                                    grain_get_QPL, MC_OK, hc_erg_um
   use sed_astrodust_mod,     only: NLAM, NA, NT, lam, aeff, T_first, &
                                    dn_ad, Cabs, kappB_first, &
                                    Cabs_cneu, Cabs_cion, dn_cneu, dn_cion, &
                                    kappB_cneu, kappB_cion
   implicit none
   private
   public :: mc_sed_solve_astrodust, mc_sed_solve_pah, mc_sed_solve_total

   ! Equilibrium gate: if E_grain(Teq) >= EEQSS_ERG then treat grain as
   ! equilibrium (matches sed_astrodust_mod's threshold, 150 eV).
   real(wp), parameter :: EEQSS_ERG  = 150.0_wp * ev2erg
   real(wp), parameter :: sigma_sb = 2.0_wp*pi**5*k_boltz**4 / &
                                     (15.0_wp*h_planck**3*c**2)

contains

   subroutine mc_sed_solve_astrodust(J_lam, U_isrf, comp, N_events, base_seed, &
                                     use_mc_for_large_grain, mc_engine_kind, &
                                     lamI_lam_out)
      ! Stochastic grains always use MC.  The flag controls the equilibrium
      ! (large) grain branch:
      !   use_mc_for_large_grain = .false. (default): equilibrium grains use
      !       calc_Teq -- delta-function emission at Teq, the analytic
      !       path sed_solve takes.
      !   use_mc_for_large_grain = .true.: equilibrium grains also use MC.
      ! mc_engine_kind selects which engine handles the MC path:
      !   'fixed'    : mc_run_engine (log grid 0.5-5000 K, fast, low-resolution
      !                tails -- the original path)
      !   '2pass'    : mc_run_engine_2pass (auto log/linear adaptive grid)
      !   'buffered' : mc_run_engine_buffered (single-pass + sample buffer)
      real(wp),         intent(in)  :: J_lam(:), U_isrf
      character(len=*), intent(in)  :: comp, mc_engine_kind
      integer,          intent(in)  :: N_events, base_seed
      logical,          intent(in)  :: use_mc_for_large_grain
      real(wp),         intent(out) :: lamI_lam_out(:)
      call mc_sed_loop(Cabs, dn_ad, kappB_first, comp, 'astrodust', J_lam, U_isrf, &
                       N_events, base_seed, use_mc_for_large_grain, &
                       mc_engine_kind, lamI_lam_out)
   end subroutine mc_sed_solve_astrodust


   subroutine mc_sed_solve_pah(J_lam, U_isrf, N_events, base_seed, &
                               use_mc_for_large_grain, mc_engine_kind, lamI_lam_out)
      ! Carbonaceous (PAH) grains carry a size-dependent ionization
      ! fraction; the neutral and cation cross sections differ markedly
      ! in the near-IR (the cation has a Mattioda et al. 2005 NIR
      ! continuum that the neutral lacks).  The production matrix solver
      ! loops both charge states and sums them, so the MC builder must do
      ! the same -- using the single charge-averaged cross section biases
      ! the near-IR by ~20%.
      real(wp),         intent(in)  :: J_lam(:), U_isrf
      character(len=*), intent(in)  :: mc_engine_kind
      integer,          intent(in)  :: N_events, base_seed
      logical,          intent(in)  :: use_mc_for_large_grain
      real(wp),         intent(out) :: lamI_lam_out(:)
      real(wp), allocatable :: lamI_cion(:)
      allocate(lamI_cion(size(lamI_lam_out)))
      ! Each charge state is cooled by the kappB of its own Cabs, as in
      ! sed_solve_pah; the f_ion-blended kappB_pah_first would set Teq (the
      ! equilibrium gate, the equilibrium emission and the adaptive DA85
      ! split) of a grain that absorbs as one charge state and cools as the
      ! blend.
      call mc_sed_loop(Cabs_cneu, dn_cneu, kappB_cneu, 'pah', 'PAH_neutral', J_lam, U_isrf, &
                       N_events, base_seed, use_mc_for_large_grain, &
                       mc_engine_kind, lamI_lam_out)
      call mc_sed_loop(Cabs_cion, dn_cion, kappB_cion, 'pah', 'PAH_ionized', J_lam, U_isrf, &
                       N_events, base_seed + 5000, use_mc_for_large_grain, &
                       mc_engine_kind, lamI_cion)
      lamI_lam_out = lamI_lam_out + lamI_cion
      deallocate(lamI_cion)
   end subroutine mc_sed_solve_pah


   subroutine mc_sed_solve_total(J_lam, U_isrf, ad_stage, N_events, base_seed, &
                                 use_mc_for_large_grain, mc_engine_kind, &
                                 lamI_lam_ad, lamI_lam_pah, lamI_lam_tot)
      real(wp),         intent(in)  :: J_lam(:), U_isrf
      character(len=*), intent(in)  :: ad_stage, mc_engine_kind
      integer,          intent(in)  :: N_events, base_seed
      logical,          intent(in)  :: use_mc_for_large_grain
      real(wp),         intent(out) :: lamI_lam_ad(:), lamI_lam_pah(:), lamI_lam_tot(:)
      character(len=16) :: comp_for_mc
      integer :: u_diag
      select case (trim(ad_stage))
      case ('S1_C1');  comp_for_mc = 'ad_s1_c1'
      case ('S1_C2');  comp_for_mc = 'ad_s1_c2'
      case ('S2');     comp_for_mc = 'ad_s2'
      case default
         write(*,'(a,a)') 'mc_sed_solve_total: unknown ad_stage ', trim(ad_stage)
         stop 1
      end select
      open(newunit=u_diag, file='output/diag_Teq_vs_Trad.dat', &
           status='replace', action='write')
      write(u_diag,'(a)') '# a_um Teq_K Trad_K Trad_over_Teq histogram_mass population'
      close(u_diag)
      call mc_sed_solve_astrodust(J_lam, U_isrf, trim(comp_for_mc), &
                                  N_events, base_seed, use_mc_for_large_grain, &
                                  mc_engine_kind, lamI_lam_ad)
      call mc_sed_solve_pah(J_lam, U_isrf, N_events, base_seed + 10000, &
                            use_mc_for_large_grain, mc_engine_kind, lamI_lam_pah)
      lamI_lam_tot = lamI_lam_ad + lamI_lam_pah
   end subroutine mc_sed_solve_total


   ! ----- internal worker -----
   subroutine mc_sed_loop(Cabs_pop, dn_pop, kappB_pop, comp, population, J_lam, U_isrf, &
                          N_events, base_seed, use_mc_for_large_grain, &
                          mc_engine_kind, lamI_lam_out)
      ! Small (stochastic) grains always use MC.  When use_mc_for_large_grain
      ! is false, large (equilibrium) grains use calc_Teq -- the analytic
      ! production path that gives exact Planck emission at Teq.  When true,
      ! large grains also use MC.
      real(wp),         intent(in)  :: Cabs_pop(:,:)        ! (NLAM, NA)
      real(wp),         intent(in)  :: dn_pop(:)            ! (NA)
      real(wp),         intent(in)  :: kappB_pop(:,:)       ! (NT, NA)
      character(len=*), intent(in)  :: comp, population, mc_engine_kind
      real(wp),         intent(in)  :: J_lam(:), U_isrf
      integer,          intent(in)  :: N_events, base_seed
      logical,          intent(in)  :: use_mc_for_large_grain
      real(wp),         intent(out) :: lamI_lam_out(:)

      integer, parameter :: NHIST_LOC = 2000
      ! Buffered-engine sample capacity. Exhaustion is an explicit error;
      ! the default two-pass engine has no sample-storage limit.
      integer(kind=8), parameter :: MAX_SAMPLES_LOC = 1000000000_8
      integer  :: ekind   ! 1=fixed, 2=2pass, 3=buffered
      real(wp), allocatable :: Jout(:)
      real(wp), allocatable :: Jout_local(:)
      real(wp) :: sum_abs, sum_emit, sum_abs_cont
      real(wp) :: diagnostics(11,NA), prescribed_power, omitted_time, floor_time, outside_time, stored_energy
      real(wp) :: diagnostic_Teq(NA), diagnostic_Trad(NA), diagnostic_mass(NA)
      integer(kind=8) :: sample_count, samples_replay
      type(rng_t) :: rng_record
      integer :: max_samp
      real(wp) :: T_eq_mc, p_abs, tau_relax, warm_time
      integer :: warm_events
      integer :: engine_status
      type(mc_grain_t) :: grain
      type(rng_t)      :: rng
      real(wp) :: T_edges(NHIST_LOC+1), dP_dT(NHIST_LOC), dP_dlnT(NHIST_LOC)
      real(wp) :: t_total, e_abs_tot, e_emit_tot
      real(wp) :: Tmid, weight, BB, inv_pi_a2, Teq, EEQ, BB_eq
      real(wp), allocatable :: Q_ja(:)
      logical  :: is_log_grid
      integer  :: ja, ilam, ibin
      integer  :: n_equil, n_stoch
      real(wp) :: T_record
      ! Adaptive lam_c bookkeeping
      real(wp), parameter :: DT_THRESHOLD = 0.01_wp
      real(wp) :: C_Teq, lam_c_adapt
      ! Locals of the equilibrium-grain diagnostic and of the single-grain
      ! spectrum below, private in the parallel region.  They are declared
      ! here and listed in private(...) rather than declared in BLOCK
      ! constructs: a BLOCK-local variable inside a default(none) region is
      ! private only from OpenMP 5.1 on, and ifx rejects it without a clause.
      real(wp) :: w_sum, T4_sum, Trad_mc
      integer  :: u_diag
      real(wp) :: jgrain(NLAM), bol_grey, a_cm2, p_abs_grain, scale_ec

      select case (trim(mc_engine_kind))
      case ('fixed');    ekind = 1
      case ('2pass');    ekind = 2
      case ('buffered'); ekind = 3
      case default
         write(*,'(a,a)') 'mc_sed_loop: unknown mc_engine_kind=', trim(mc_engine_kind)
         stop 1
      end select
      write(*,'(a,a,a,i0,a)') '  mc_sed_loop: engine=', trim(mc_engine_kind), &
                              ' (ekind=', ekind, ')'

      allocate(Jout(NLAM))
      Jout         = 0.0_wp
      sum_abs      = 0.0_wp
      sum_abs_cont = 0.0_wp
      sum_emit     = 0.0_wp
      n_equil      = 0
      n_stoch      = 0
      diagnostics  = 0.0_wp
      diagnostic_Teq = 0.0_wp
      diagnostic_Trad = 0.0_wp
      diagnostic_mass = 0.0_wp

      !$omp parallel default(none) &
      !$omp   shared(Cabs_pop, dn_pop, kappB_pop, J_lam, U_isrf, comp, &
      !$omp          N_events, base_seed, lam, aeff, NLAM, NA, NT, T_first, &
      !$omp          Jout, sum_abs, sum_abs_cont, sum_emit, &
      !$omp          n_equil, n_stoch, diagnostics, diagnostic_Teq, diagnostic_Trad, diagnostic_mass, &
      !$omp          use_mc_for_large_grain, ekind) &
      !$omp   private(ja, ilam, ibin, grain, rng, T_edges, dP_dT, dP_dlnT, &
      !$omp          t_total, e_abs_tot, e_emit_tot, Tmid, weight, BB, BB_eq, &
      !$omp          Jout_local, Q_ja, inv_pi_a2, Teq, EEQ, is_log_grid, &
      !$omp          C_Teq, lam_c_adapt, w_sum, T4_sum, Trad_mc, &
      !$omp          jgrain, bol_grey, a_cm2, p_abs_grain, scale_ec, &
      !$omp          engine_status, omitted_time, prescribed_power, floor_time, outside_time, sample_count, T_record, stored_energy, &
      !$omp          T_eq_mc, p_abs, tau_relax, warm_time, warm_events, rng_record, samples_replay, max_samp)

      allocate(Jout_local(NLAM))
      Jout_local = 0.0_wp
      allocate(Q_ja(NLAM))

      !$omp do schedule(dynamic, 1)
      do ja = 1, NA
         ! ---- Equilibrium gate: compute Teq from absorbed = emitted, then ---
         ! ---- decide whether E_grain(Teq) is large enough that single   ---
         ! ---- photon absorptions barely perturb T. ----------------------
         call calc_Teq(lam, Cabs_pop(:, ja), J_lam, T_first, &
                       kappB_pop(:, ja), Teq)
         EEQ = U_of_T(Teq, aeff(ja), comp)

         if (EEQ >= EEQSS_ERG .and. .not. use_mc_for_large_grain) then
            ! Equilibrium path: delta-function emission at Teq, matching
            ! sed_solve's calc_Teq branch.  Skipped if user requested MC
            ! for large grains too.
            do ilam = 1, NLAM
               BB_eq = bbody(Teq, lam(ilam))
               Jout_local(ilam) = Jout_local(ilam) + &
                                  dn_pop(ja) * Cabs_pop(ilam, ja) * BB_eq
            end do
            !$omp atomic
            n_equil = n_equil + 1
         else
            ! Stochastic path: full MC.
            call rng_init(rng, base_seed + 1009 * ja)
            inv_pi_a2 = 1.0_wp / (PI * (aeff(ja)*1.0e-4_wp)**2)
            Q_ja = Cabs_pop(:, ja) * inv_pi_a2

            ! Adaptive DA85 cutoff: photons whose energy hc/lam produces a
            ! grain-temperature jump dT = (hc/lam)/C(Teq) below
            ! DT_THRESHOLD * Teq are too weak to perturb the grain and are
            ! folded into the continuous heating term H_cont; only
            ! shorter-wavelength photons remain stochastic events.  This
            ! approximates unresolved photon impulses near Teq. It does
            ! not impose a CMB lower boundary on a stochastic trajectory;
            ! continuous heating alone determines cooling between events.
            C_Teq = C_of_T(Teq, aeff(ja), comp)
            if (C_Teq > 0.0_wp .and. Teq > 0.0_wp) then
               lam_c_adapt = HC_ERG_UM / (DT_THRESHOLD * Teq * C_Teq)
               ! Clamp to the wavelength grid so the integration covers
               ! at least one bin on each side of the split.
               if (lam_c_adapt < lam(2))      lam_c_adapt = lam(2)
               if (lam_c_adapt > lam(NLAM-1)) lam_c_adapt = lam(NLAM-1)
            else
               lam_c_adapt = 1000.0_wp     ! the fixed DA85 split used before the adaptive one
            end if

            call grain_setup_from_cabs(grain, aeff(ja), U_isrf, comp, lam, Q_ja, &
                                       lam_c_in=lam_c_adapt, J_lam_in=J_lam)
            call grain_radiative_equilibrium(grain, T_eq_mc, p_abs, tau_relax, engine_status)
            if (engine_status /= MC_OK) error stop 'Radiative equilibrium failed'
            call mc_equilibrate(grain, rng, T_eq_mc, T_record, warm_time, warm_events, engine_status)
            if (engine_status /= MC_OK) error stop 'MC warmup failed'
            !$omp critical(mc_warmup_log)
            write(*,'(a,a,a,es24.16e3,a,es24.16e3,a,es24.16e3,a,es24.16e3,a,i0)') &
               '[mc_warmup] comp=', trim(comp), ' radius_um=', aeff(ja), &
               ' initial_Teq_K=', T_eq_mc, ' record_T_K=', T_record, &
               ' warm_time_s=', warm_time, ' events=', warm_events
            !$omp end critical(mc_warmup_log)
            rng_record = rng
            select case (ekind)
            case (1)   ! fixed-grid (original)
               call mc_run_engine(grain, rng, T_record, N_events, 0.0_wp, &
                                  NHIST_LOC, T_edges, dP_dT, dP_dlnT, &
                                  t_total, e_abs_tot, e_emit_tot, &
                                  status_out=engine_status, omitted_time_out=omitted_time, &
                                        floor_time_out=floor_time, out_of_grid_time_out=outside_time, &
                                        sample_count_out=sample_count, stored_energy_out=stored_energy)
               is_log_grid = .true.    ! fixed engine is always log 0.5-5000 K
            case (2)   ! 2pass adaptive
               call mc_run_engine_2pass(grain, rng, T_record, N_events, -1, &
                                        0.0_wp, NHIST_LOC, T_edges, dP_dT, dP_dlnT, &
                                        t_total, e_abs_tot, e_emit_tot, &
                                        is_log_out=is_log_grid, &
                                        status_out=engine_status, omitted_time_out=omitted_time, &
                                        floor_time_out=floor_time, out_of_grid_time_out=outside_time, &
                                        sample_count_out=sample_count, stored_energy_out=stored_energy)
            case (3)   ! buffered adaptive
               ! Serialize capacity measurement and allocation to bound process memory.
               !$omp critical(mc_buffered_run)
               call mc_run_engine_2pass(grain, rng, T_record, N_events, -1, &
                                        0.0_wp, NHIST_LOC, T_edges, dP_dT, dP_dlnT, &
                                        t_total, e_abs_tot, e_emit_tot, &
                                        status_out=engine_status, omitted_time_out=omitted_time, &
                                        sample_count_out=samples_replay)
               if (engine_status /= MC_OK) error stop 'Buffered capacity measurement failed'
               if (samples_replay <= 0_8 .or. samples_replay > MAX_SAMPLES_LOC) &
                  error stop 'Buffered capacity exceeds 1e9 samples (16 GB process limit)'
               max_samp = int(samples_replay)
               rng = rng_record
               call mc_run_engine_buffered(grain, rng, T_record, N_events, &
                                           0.0_wp, NHIST_LOC, max_samp, &
                                           T_edges, dP_dT, dP_dlnT, &
                                           t_total, e_abs_tot, e_emit_tot, &
                                           is_log_out=is_log_grid, &
                                        status_out=engine_status, omitted_time_out=omitted_time, &
                                        floor_time_out=floor_time, out_of_grid_time_out=outside_time, &
                                        sample_count_out=sample_count, stored_energy_out=stored_energy)
               if (sample_count /= samples_replay) error stop 'Buffered sample count differs from replay'
               !$omp end critical(mc_buffered_run)
            end select
            if (engine_status /= MC_OK) then
               !$omp critical(mc_failure)
               write(error_unit,'(a,a,a,es16.8,a,i0,a,es16.8)') &
                  'MC cooling failed: composition=', trim(comp), ' a_um=', aeff(ja), &
                  ' status=', engine_status, ' incomplete_time_s=', omitted_time
               !$omp end critical(mc_failure)
               error stop 'MC interval incomplete; no spectrum is accepted'
            end if
            ! Independent prescribed-field power, without a bolometric rescale.
            ! J is W m^-3 sr^-1, Cabs is cm^2, lambda is um: factor 1e-3
            ! converts the wavelength integral to erg s^-1.
            prescribed_power = 0.0_wp
            do ilam = 1, NLAM-1
               prescribed_power = prescribed_power + 0.5_wp * (lam(ilam+1)-lam(ilam)) * &
                  (Cabs_pop(ilam,ja)*J_lam(ilam) + Cabs_pop(ilam+1,ja)*J_lam(ilam+1))
            end do
            diagnostics(1,ja) = 4.0_wp * PI * 1.0e-3_wp * prescribed_power
            diagnostics(3,ja) = 4.0_wp * PI * (aeff(ja)*1.0e-4_wp)**2 * grain%H_cont
            if (t_total > 0.0_wp) then
               diagnostics(2,ja) = e_abs_tot / t_total
               diagnostics(4,ja) = e_emit_tot / t_total
            end if
            diagnostics(5,ja) = t_total
            diagnostics(6,ja) = omitted_time
            diagnostics(7,ja) = sum(dP_dT * (T_edges(2:)-T_edges(:NHIST_LOC)))
            diagnostics(8,ja) = floor_time
            diagnostics(9,ja) = outside_time
            diagnostics(10,ja) = real(sample_count,wp)
            diagnostics(11,ja) = stored_energy
            !$omp atomic
            sum_abs = sum_abs + e_abs_tot
            !$omp atomic
            sum_emit = sum_emit + e_emit_tot
            ! Continuous absorption rate is grain%H_cont; integrate
            ! over the grain surface (4 pi a^2) and the simulated time t_total
            ! to get the energy contribution from photons folded into the
            ! continuous heating term.  Add to the proper energy-balance
            ! denominator (events alone are misleading once lam_c is adaptive).
            !$omp atomic
            sum_abs_cont = sum_abs_cont &
                  + 4.0_wp * PI * (aeff(ja)*1.0e-4_wp)**2 * grain%H_cont * t_total
            ! Single-grain diagnostic: dump MC <T^4>^(1/4) vs calc_Teq to
            ! a side file when running MC-all mode.
            ! for the plot.
            if (use_mc_for_large_grain .and. EEQ >= EEQSS_ERG) then
               w_sum  = 0.0_wp
               T4_sum = 0.0_wp
               do ibin = 1, NHIST_LOC
                  if (dP_dT(ibin) <= 0.0_wp) cycle
                  if (is_log_grid) then
                     Tmid = sqrt(T_edges(ibin) * T_edges(ibin+1))
                  else
                     Tmid = 0.5_wp * (T_edges(ibin) + T_edges(ibin+1))
                  end if
                  weight = dP_dT(ibin) * (T_edges(ibin+1) - T_edges(ibin))
                  w_sum  = w_sum  + weight
                  T4_sum = T4_sum + weight * Tmid**4
               end do
               if (w_sum > 0.0_wp) then
                  Trad_mc = (T4_sum / w_sum) ** 0.25_wp
                  diagnostic_Teq(ja) = Teq
                  diagnostic_Trad(ja) = Trad_mc
                  diagnostic_mass(ja) = w_sum
               end if
            end if

            ! Build the emission spectrum from the temperature histogram and
            ! normalize its integral to the absorbed power. Evaluating Planck
            ! emission at bin midpoints introduces a quadrature error whose
            ! sign is not fixed by convexity alone. The normalization enforces
            ! stationary energy balance, not the accuracy of spectral shape
            ! or of the sampled temperature distribution. Event photons and
            ! continuous heating both enter the absorbed-power integral.
            a_cm2    = (aeff(ja) * 1.0e-4_wp)**2
            jgrain   = 0.0_wp
            bol_grey = 0.0_wp
            do ibin = 1, NHIST_LOC
               if (dP_dT(ibin) <= 0.0_wp) cycle
               if (is_log_grid) then
                  Tmid = sqrt(T_edges(ibin) * T_edges(ibin+1))
               else
                  Tmid = 0.5_wp * (T_edges(ibin) + T_edges(ibin+1))
               end if
               weight = dP_dT(ibin) * (T_edges(ibin+1) - T_edges(ibin))
               do ilam = 1, NLAM
                  jgrain(ilam) = jgrain(ilam) + &
                     dn_pop(ja) * Cabs_pop(ilam, ja) * bbody(Tmid, lam(ilam)) * weight
               end do
               bol_grey = bol_grey + weight * 4.0_wp * pi * a_cm2 * &
                          grain_get_QPL(grain, Tmid) * sigma_sb * Tmid**4
            end do
            ! grey absorbed power per unit time (events + continuous)
            if (t_total > 0.0_wp) then
               p_abs_grain = e_abs_tot / t_total + &
                             4.0_wp * pi * a_cm2 * grain%H_cont
            else
               p_abs_grain = 0.0_wp
            end if
            if (bol_grey > 0.0_wp) then
               scale_ec = p_abs_grain / bol_grey
            else
               scale_ec = 1.0_wp
            end if
            Jout_local = Jout_local + scale_ec * jgrain
            if (allocated(grain%lam_grid))  deallocate(grain%lam_grid)
            if (allocated(grain%Q_grid))    deallocate(grain%Q_grid)
            if (allocated(grain%u_lam_grid))deallocate(grain%u_lam_grid)
            if (allocated(grain%cdf_F))     deallocate(grain%cdf_F)
            if (allocated(grain%cdf_lam))   deallocate(grain%cdf_lam)
            if (allocated(grain%cdf_density)) deallocate(grain%cdf_density)
            if (allocated(grain%T_grid))    deallocate(grain%T_grid)
            if (allocated(grain%U_grid))    deallocate(grain%U_grid)
            if (allocated(grain%QPL_grid))  deallocate(grain%QPL_grid)
            !$omp atomic
            n_stoch = n_stoch + 1
         end if
      end do
      !$omp end do

      ! Combine thread-local accumulators
      !$omp critical
      Jout = Jout + Jout_local
      !$omp end critical

      deallocate(Jout_local, Q_ja)
      !$omp end parallel

      if (use_mc_for_large_grain) then
         open(newunit=u_diag, file='output/diag_Teq_vs_Trad.dat', &
              status='old', position='append', action='write')
         do ja = 1, NA
            if (diagnostic_mass(ja) <= 0.0_wp) cycle
            write(u_diag, '(5(es14.6,1x),a)') &
                 aeff(ja), diagnostic_Teq(ja), diagnostic_Trad(ja), &
                 diagnostic_Trad(ja)/diagnostic_Teq(ja), diagnostic_mass(ja), trim(population)
         end do
         close(u_diag)
      end if

      write(*,'(a,i4,a,i4,a,es11.3,a,es11.3,a,es11.3,a,es12.4)') &
         '  hybrid: n_equil=', n_equil, '  n_stoch=', n_stoch, &
         '  abs_evt=', sum_abs, &
         '  abs_cont=', sum_abs_cont, &
         '  emit=', sum_emit, &
         '  emit/abs_total=', sum_emit / max(sum_abs + sum_abs_cont, 1.0e-300_wp)

      ! Report each simulated grain separately: summing trajectories with
      ! different durations is not an independent radiation-field energy test.
      write(*,'(a)') '# mc_energy: comp a_um prescribed_erg_s event_erg_s continuous_erg_s ' // &
                    'emitted_erg_s elapsed_s omitted_s histogram_mass floor_s outside_grid_s samples stored_energy_erg'
      do ja = 1, NA
         if (diagnostics(5,ja) <= 0.0_wp) cycle
         write(*,'(a,1x,a,1x,12(es24.16e3,1x))') '[mc_energy]', trim(comp), aeff(ja), diagnostics(:,ja)
      end do

      ! Same unit convention as sed_astrodust_mod::sed_solve:
      !   lambda * I_lambda per H = lambda(cm) * Jout(CGS)
      !                           = (lambda_um * 1e-4) * (10 * Jout)
      !                           = lambda_um * Jout * 1e-3
      lamI_lam_out = lam * Jout * 1.0e-3_wp
      deallocate(Jout)
   end subroutine mc_sed_loop

end module mc_sed
