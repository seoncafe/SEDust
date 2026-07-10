module mc_sed
   ! MC-based SED builder.  Mirrors the structure of
   ! sed_astrodust_mod::sed_solve but replaces the calc_P matrix solver
   ! with the Monte Carlo engine.
   !
   ! Two populations are computed:
   !   - Astrodust grains    (Cabs, dn_ad, comp = 'ad_s1_c1'|'ad_s1_c2'|'ad_s2')
   !   - PAH grains          (Cabs_pah, dn_pah, comp = 'pah')
   !
   ! The size loop is OpenMP-parallel.  Each iteration creates its own
   ! grain_t and rng_t (both stack-local in the parallel region), runs
   ! MC, and atomically accumulates into the shared Jout(:).
   !
   ! Output convention matches sed_astrodust_mod::sed_solve:
   !   lamI_lam_out(:) = lambda * Jout * 1e-3
   ! so that the result can be directly compared with that subroutine.

   use, intrinsic :: omp_lib
   use constants,             only: wp
   use radfield,              only: bbody
   use p_sub,                 only: calc_Teq
   use mc_heatcap,            only: U_of_T, C_of_T
   use mc_grain_type,         only: mc_grain_t
   use mc_rng,                only: rng_t, rng_init
   use mc_engine,             only: grain_setup_from_cabs, mc_run_engine, &
                                    mc_run_engine_2pass, mc_run_engine_buffered, &
                                    grain_get_QPL
   use sed_astrodust_mod,     only: sed_init, NLAM, NA, NT, lam, aeff,    &
                                    dn_ad, dn_pah, Cabs, Cabs_pah, T_first, &
                                    kappB_first, kappB_pah_first, &
                                    Cabs_cneu, Cabs_cion, dn_cneu, dn_cion, &
                                    sed_solve, sed_solve_pah
   implicit none
   private
   public :: mc_sed_solve_astrodust, mc_sed_solve_pah, mc_sed_solve_total

   real(wp), parameter :: PI = 3.141592653589793238462643383279502884197_wp
   ! Equilibrium gate: if E_grain(Teq) >= EEQSS_ERG then treat grain as
   ! equilibrium (matches sed_astrodust_mod's threshold, 150 eV).
   real(wp), parameter :: EV_TO_ERG  = 1.60218e-12_wp
   real(wp), parameter :: EEQSS_ERG  = 150.0_wp * EV_TO_ERG

contains

   subroutine mc_sed_solve_astrodust(J_lam, U_isrf, comp, N_events, base_seed, &
                                     use_mc_for_large_grain, mc_engine_kind, &
                                     lamI_lam_out)
      ! Stochastic grains always use MC.  The flag controls the equilibrium
      ! (large) grain branch:
      !   use_mc_for_large_grain = .false. (default): equilibrium grains use
      !       calc_Teq -- delta-function emission at Teq, the analytic
      !       production-pipeline path.
      !   use_mc_for_large_grain = .true.: equilibrium grains also use MC.
      ! mc_engine_kind selects which engine handles the MC path:
      !   'fixed'    : mc_run_engine (log grid 2-5000 K, fast, low-resolution
      !                tails -- the original path)
      !   '2pass'    : mc_run_engine_2pass (auto log/linear adaptive grid)
      !   'buffered' : mc_run_engine_buffered (single-pass + sample buffer)
      real(wp),         intent(in)  :: J_lam(:), U_isrf
      character(len=*), intent(in)  :: comp, mc_engine_kind
      integer,          intent(in)  :: N_events, base_seed
      logical,          intent(in)  :: use_mc_for_large_grain
      real(wp),         intent(out) :: lamI_lam_out(:)
      call mc_sed_loop(Cabs, dn_ad, kappB_first, comp, J_lam, U_isrf, &
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
      call mc_sed_loop(Cabs_cneu, dn_cneu, kappB_pah_first, 'pah', J_lam, U_isrf, &
                       N_events, base_seed, use_mc_for_large_grain, &
                       mc_engine_kind, lamI_lam_out)
      call mc_sed_loop(Cabs_cion, dn_cion, kappB_pah_first, 'pah', J_lam, U_isrf, &
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
      select case (trim(ad_stage))
      case ('S1_C1');  comp_for_mc = 'ad_s1_c1'
      case ('S1_C2');  comp_for_mc = 'ad_s1_c2'
      case ('S2');     comp_for_mc = 'ad_s2'
      case default
         write(*,'(a,a)') 'mc_sed_solve_total: unknown ad_stage ', trim(ad_stage)
         stop 1
      end select
      call mc_sed_solve_astrodust(J_lam, U_isrf, trim(comp_for_mc), &
                                  N_events, base_seed, use_mc_for_large_grain, &
                                  mc_engine_kind, lamI_lam_ad)
      call mc_sed_solve_pah(J_lam, U_isrf, N_events, base_seed + 10000, &
                            use_mc_for_large_grain, mc_engine_kind, lamI_lam_pah)
      lamI_lam_tot = lamI_lam_ad + lamI_lam_pah
   end subroutine mc_sed_solve_total


   ! ----- internal worker -----
   subroutine mc_sed_loop(Cabs_pop, dn_pop, kappB_pop, comp, J_lam, U_isrf, &
                          N_events, base_seed, use_mc_for_large_grain, &
                          mc_engine_kind, lamI_lam_out)
      ! Small (stochastic) grains always use MC.  When use_mc_for_large_grain
      ! is false, large (equilibrium) grains use calc_Teq -- the analytic
      ! production path that gives exact Planck emission at Teq.  When true,
      ! large grains also use MC.
      real(wp),         intent(in)  :: Cabs_pop(:,:)        ! (NLAM, NA)
      real(wp),         intent(in)  :: dn_pop(:)            ! (NA)
      real(wp),         intent(in)  :: kappB_pop(:,:)       ! (NT, NA)
      character(len=*), intent(in)  :: comp, mc_engine_kind
      real(wp),         intent(in)  :: J_lam(:), U_isrf
      integer,          intent(in)  :: N_events, base_seed
      logical,          intent(in)  :: use_mc_for_large_grain
      real(wp),         intent(out) :: lamI_lam_out(:)

      integer, parameter :: NHIST_LOC = 2000
      ! Buffered-engine sub-step cap; trims memory while still over-sampling
      ! the stochastic excursion (8 sub-steps/event x N_events comfortably
      ! exceeds this for small grains, so the buffer fills mid-run).
      integer, parameter :: MAX_SAMPLES_LOC = 200000
      integer  :: ekind   ! 1=fixed, 2=2pass, 3=buffered
      real(wp), allocatable :: Jout(:)
      real(wp), allocatable :: Jout_local(:)
      real(wp) :: sum_abs, sum_emit, sum_abs_cont
      type(mc_grain_t) :: grain
      type(rng_t)      :: rng
      real(wp) :: T_edges(NHIST_LOC+1), dP_dT(NHIST_LOC), dP_dlnT(NHIST_LOC)
      real(wp) :: t_total, e_abs_tot, e_emit_tot
      real(wp) :: Tmid, weight, BB, inv_pi_a2, Teq, EEQ, BB_eq
      real(wp), allocatable :: Q_ja(:)
      logical  :: is_log_grid
      integer  :: ja, ilam, ibin, N_burn
      integer  :: n_equil, n_stoch
      real(wp) :: T_init
      ! Adaptive lam_c bookkeeping
      real(wp), parameter :: DT_THRESHOLD = 0.01_wp
      real(wp), parameter :: HC_ERG_UM    = 1.98644582e-13_wp  ! h*c in erg*um
      real(wp) :: C_Teq, lam_c_adapt

      T_init = 2.725_wp          ! CMB temperature (K)
      N_burn = max(50, N_events / 20)

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

      !$omp parallel default(none) &
      !$omp   shared(Cabs_pop, dn_pop, kappB_pop, J_lam, U_isrf, comp, &
      !$omp          N_events, base_seed, lam, aeff, NLAM, NA, NT, T_first, &
      !$omp          Jout, T_init, sum_abs, sum_abs_cont, sum_emit, &
      !$omp          n_equil, n_stoch, &
      !$omp          use_mc_for_large_grain, ekind, N_burn) &
      !$omp   private(ja, ilam, ibin, grain, rng, T_edges, dP_dT, dP_dlnT, &
      !$omp          t_total, e_abs_tot, e_emit_tot, Tmid, weight, BB, BB_eq, &
      !$omp          Jout_local, Q_ja, inv_pi_a2, Teq, EEQ, is_log_grid, &
      !$omp          C_Teq, lam_c_adapt)

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
            ! repairs the DA85 split for grains large enough that even
            ! UV photons cannot deliver an impulse comparable to a Teq
            ! fluctuation, preventing the trajectory from collapsing
            ! toward CMB between events.
            C_Teq = C_of_T(Teq, aeff(ja), comp)
            if (C_Teq > 0.0_wp .and. Teq > 0.0_wp) then
               lam_c_adapt = HC_ERG_UM / (DT_THRESHOLD * Teq * C_Teq)
               ! Clamp to the wavelength grid so the integration covers
               ! at least one bin on each side of the split.
               if (lam_c_adapt < lam(2))      lam_c_adapt = lam(2)
               if (lam_c_adapt > lam(NLAM-1)) lam_c_adapt = lam(NLAM-1)
            else
               lam_c_adapt = 1000.0_wp     ! fallback to original default
            end if

            call grain_setup_from_cabs(grain, aeff(ja), U_isrf, comp, lam, Q_ja, &
                                       lam_c_in=lam_c_adapt)
            select case (ekind)
            case (1)   ! fixed-grid (original)
               call mc_run_engine(grain, rng, T_init, N_events, 0.0_wp, &
                                  NHIST_LOC, T_edges, dP_dT, dP_dlnT, &
                                  t_total, e_abs_tot, e_emit_tot)
               is_log_grid = .true.    ! fixed engine is always log 2-5000 K
            case (2)   ! 2pass adaptive
               call mc_run_engine_2pass(grain, rng, T_init, N_events, -1, N_burn, &
                                        0.0_wp, NHIST_LOC, T_edges, dP_dT, dP_dlnT, &
                                        t_total, e_abs_tot, e_emit_tot, &
                                        is_log_out=is_log_grid)
            case (3)   ! buffered adaptive
               call mc_run_engine_buffered(grain, rng, T_init, N_events, N_burn, &
                                           0.0_wp, NHIST_LOC, MAX_SAMPLES_LOC, &
                                           T_edges, dP_dT, dP_dlnT, &
                                           t_total, e_abs_tot, e_emit_tot, &
                                           is_log_out=is_log_grid)
            end select
            !$omp atomic
            sum_abs = sum_abs + e_abs_tot
            !$omp atomic
            sum_emit = sum_emit + e_emit_tot
            ! Continuous absorption rate (per area) is grain%H_cont; integrate
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
               block
                  real(wp) :: w_sum, T4_sum, Trad_mc
                  integer  :: u_diag
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
                     !$omp critical
                     open(newunit=u_diag, &
                          file='output/diag_Teq_vs_Trad.dat', &
                          status='unknown', position='append', action='write')
                     write(u_diag, '(es14.6,1x,4(es14.6,1x))') &
                          aeff(ja), Teq, Trad_mc, Trad_mc/Teq, w_sum
                     close(u_diag)
                     !$omp end critical
                  end if
               end block
            end if

            ! Build the single-grain emission spectrum from the histogram and
            ! apply an energy-conservation rescale.  The histogram SED
            ! re-evaluates the Planck emission at each bin midpoint Tmid;
            ! because B(T) is convex (~T^4) this over-weights the hot tail
            ! relative to the trajectory-accurate emission, leaving the
            ! grain emitting more than it absorbs (most severe for the
            ! violently heated tiniest PAH grains).  A grain in steady
            ! state emits exactly what it absorbs (Kirchhoff), so we
            ! rescale the spectrum to the grey grain absorbed power
            ! (event photons + continuous heating), exactly as the matrix
            ! solver does by construction.
            block
               real(wp) :: jgrain(NLAM), bol_grey, a_cm2, p_abs_grain, scale_ec
               real(wp), parameter :: SIGMA_SB = 5.6703744e-5_wp
               real(wp), parameter :: PI_L = 3.141592653589793238_wp
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
                  bol_grey = bol_grey + weight * 4.0_wp * PI_L * a_cm2 * &
                             grain_get_QPL(grain, Tmid) * SIGMA_SB * Tmid**4
               end do
               ! grey absorbed power per unit time (events + continuous)
               if (t_total > 0.0_wp) then
                  p_abs_grain = e_abs_tot / t_total + &
                                4.0_wp * PI_L * a_cm2 * grain%H_cont
               else
                  p_abs_grain = 0.0_wp
               end if
               if (bol_grey > 0.0_wp) then
                  scale_ec = p_abs_grain / bol_grey
               else
                  scale_ec = 1.0_wp
               end if
               Jout_local = Jout_local + scale_ec * jgrain
            end block
            if (allocated(grain%lam_grid))  deallocate(grain%lam_grid)
            if (allocated(grain%Q_grid))    deallocate(grain%Q_grid)
            if (allocated(grain%u_lam_grid))deallocate(grain%u_lam_grid)
            if (allocated(grain%cdf_F))     deallocate(grain%cdf_F)
            if (allocated(grain%cdf_lam))   deallocate(grain%cdf_lam)
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

      write(*,'(a,i4,a,i4,a,es11.3,a,es11.3,a,es11.3,a,f7.4)') &
         '  hybrid: n_equil=', n_equil, '  n_stoch=', n_stoch, &
         '  abs_evt=', sum_abs, &
         '  abs_cont=', sum_abs_cont, &
         '  emit=', sum_emit, &
         '  emit/abs_total=', sum_emit / max(sum_abs + sum_abs_cont, 1.0e-300_wp)

      ! Same unit convention as sed_astrodust_mod::sed_solve:
      !   lambda * I_lambda per H = lambda(cm) * Jout(CGS)
      !                           = (lambda_um * 1e-4) * (10 * Jout)
      !                           = lambda_um * Jout * 1e-3
      lamI_lam_out = lam * Jout * 1.0e-3_wp
      deallocate(Jout)
   end subroutine mc_sed_loop

end module mc_sed
