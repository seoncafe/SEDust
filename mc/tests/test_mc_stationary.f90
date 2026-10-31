program test_mc_stationary
   use constants, only: wp
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
   use mc_grain_type, only: mc_grain_t
   use mc_heatcap, only: U_of_T
   use mc_rng, only: rng_t, rng_init
   use mc_engine, only: grain_setup_from_cabs, grain_radiative_equilibrium, mc_equilibrate, mc_run_engine, &
                        mc_run_engine_2pass, mc_run_engine_buffered, MC_OK, &
                        MC_STEP_BUDGET, MC_INVALID_STEP
   implicit none

   integer, parameter :: nh = 64, ne = 6
   real(wp), parameter :: lam(6) = [0.1_wp, 0.2_wp, 1.0_wp, 10.0_wp, 100.0_wp, 1000.0_wp]
   real(wp), parameter :: q(6) = 1.0e-3_wp, field(6) = 1.0e-8_wp
   type(mc_grain_t) :: g, bad
   type(rng_t) :: rng_a, rng_b, rng_record, rng_method, rng_zero
   real(wp) :: t_a, t_b, wt_a, wt_b, t_record, wt_record, teq, p_abs, tau, expected_time
   real(wp) :: edges(nh+1), density(nh), log_density(nh)
   real(wp) :: elapsed, absorbed, emitted, stored, pre(ne+1), post(ne+1), event_time(ne+1), dt(ne+1)
   real(wp) :: elapsed_ref, absorbed_ref, emitted_ref, stored_ref
   integer :: status, events_a, events_b, events_record, nrec, mode
   integer(kind=8) :: record_state

   call grain_setup_from_cabs(g, 3.55e-4_wp, 1.0_wp, 'ad_s2', lam, q, J_lam_in=field)
   g%rate_event = 1.0e-4_wp
   g%H_cont = 1.0e-3_wp
   call grain_radiative_equilibrium(g, teq, p_abs, tau, status)
   if (status /= MC_OK) error stop 'radiative equilibrium failed'
   expected_time = max(0.01_wp*tau, 4.0_wp/g%rate_event)

   call rng_init(rng_a, 20261030)
   rng_b = rng_a
   call mc_equilibrate(g, rng_a, 40.0_wp, t_a, wt_a, events_a, status, &
                       duration_factor=0.01_wp, expected_events=4, chunk_steps=1)
   if (status /= MC_OK) error stop 'warm start failed'
   call mc_equilibrate(g, rng_b, 40.0_wp, t_b, wt_b, events_b, status, &
                       duration_factor=0.01_wp, expected_events=4, chunk_steps=1000)
   if (status /= MC_OK) error stop 'warm start with large chunks failed'
   if (t_a /= t_b .or. wt_a /= wt_b .or. events_a /= events_b .or. &
       rng_a%state /= rng_b%state) error stop 'warm start depends on chunk size'
   if (wt_a /= expected_time) error stop 'warm start missed fixed duration'
   if (events_a < 0) error stop 'negative warm event count'
   if (.not. ieee_is_finite(t_a) .or. t_a < 0.5_wp) error stop 'invalid warm temperature'

   ! Rewind the same seed for each recording method. Its accounting begins at
   ! the warmed temperature and RNG state, independently of the warm interval.
   call rng_init(rng_record, 20261030)
   call mc_equilibrate(g, rng_record, 40.0_wp, t_record, wt_record, events_record, status, &
                       duration_factor=0.01_wp, expected_events=4, chunk_steps=13)
   if (status /= MC_OK) error stop 'warm replay failed'
   if (t_record /= t_a .or. wt_record /= wt_a .or. &
       events_record /= events_a .or. rng_record%state /= rng_a%state) &
      error stop 'warm replay differs'
   record_state = rng_record%state

   do mode = 1, 3
      rng_method = rng_record
      select case (mode)
      case (1)
         call mc_run_engine(g, rng_method, t_record, ne, 0.0_wp, nh, edges, density, log_density, &
                            elapsed, absorbed, emitted, event_time, pre, post, dt, nrec, &
                            status_out=status, stored_energy_out=stored)
         if (nrec /= ne+1) error stop 'fixed recording count includes warm events'
         if (event_time(1) /= 0.0_wp .or. event_time(2) <= 0.0_wp .or. event_time(2) >= elapsed) &
            error stop 'event time is not relative to recording start'
         elapsed_ref = elapsed
         absorbed_ref = absorbed
         emitted_ref = emitted
         stored_ref = stored
      case (2)
         call mc_run_engine_2pass(g, rng_method, t_record, ne, ne, 0.0_wp, nh, &
                                  edges, density, log_density, elapsed, absorbed, emitted, &
                                  status_out=status, stored_energy_out=stored)
      case (3)
         call mc_run_engine_buffered(g, rng_method, t_record, ne, 0.0_wp, nh, 100000, &
                                     edges, density, log_density, elapsed, absorbed, emitted, &
                                     status_out=status, stored_energy_out=stored)
      end select
      if (status /= MC_OK) then
         write(*,*) 'recording method and status:', mode, status
         error stop 'recording after warm start failed'
      end if
      if (rng_method%state == record_state) error stop 'recording did not advance RNG'
      if (elapsed <= 0.0_wp .or. absorbed <= 0.0_wp) error stop 'empty recording interval'
      if (abs(elapsed-elapsed_ref) > 1.0e-12_wp*elapsed_ref .or. &
          abs(absorbed-absorbed_ref) > 1.0e-12_wp*absorbed_ref .or. &
          abs(emitted-emitted_ref) > 1.0e-12_wp*max(tiny(1.0_wp),abs(emitted_ref)) .or. &
          abs(stored-stored_ref) > 1.0e-12_wp*max(tiny(1.0_wp),abs(stored_ref))) &
          error stop 'recording totals differ between methods'
      if (abs(stored-(U_of_T(post(ne+1),g%a_um,g%comp)-U_of_T(t_record,g%a_um,g%comp))) > &
          1.0e-12_wp*max(tiny(1.0_wp),abs(stored))) &
          error stop 'stored heat includes warm interval'
      if (abs(sum(density*(edges(2:)-edges(:nh)))-1.0_wp) > 1.0e-9_wp) &
         error stop 'recording histogram is not normalized'
   end do

   call rng_init(rng_a, 1234)
   call mc_equilibrate(g, rng_a, 40.0_wp, t_a, wt_a, events_a, status, &
                       duration_factor=1.0_wp, expected_events=4, max_events=1)
   if (status /= MC_STEP_BUDGET) error stop 'warm event limit was ignored'
   call rng_init(rng_a, 1234)
   call mc_equilibrate(g, rng_a, 40.0_wp, t_a, wt_a, events_a, status, &
                       duration_factor=1.0_wp, expected_events=4, max_steps_total=1)
   if (status /= MC_STEP_BUDGET) error stop 'warm step limit was ignored'
   bad = g
   bad%H_cont = ieee_value(0.0_wp, ieee_quiet_nan)
   call rng_init(rng_a, 1234)
   call mc_equilibrate(bad, rng_a, 40.0_wp, t_a, wt_a, events_a, status, &
                       duration_factor=0.01_wp, expected_events=4)
   if (status /= MC_INVALID_STEP) error stop 'nonfinite warm heating was accepted'
   bad = g
   bad%rate_event = 0.0_wp
   call rng_init(rng_a, 1234)
   rng_zero = rng_a
   call mc_equilibrate(bad, rng_a, 40.0_wp, t_a, wt_a, events_a, status, &
                       duration_factor=0.01_wp, expected_events=4)
   if (status /= MC_OK .or. events_a /= 0 .or. wt_a /= 0.01_wp*tau .or. &
       rng_a%state /= rng_zero%state) error stop 'zero-rate warmup sampled an event'
   bad = g
   bad%rate_event = 1.0e-12_wp
   call rng_init(rng_a, 1234)
   rng_zero = rng_a
   call mc_equilibrate(bad, rng_a, 40.0_wp, t_a, wt_a, events_a, status, &
                       duration_factor=0.01_wp, expected_events=0)
   if (status /= MC_OK .or. events_a /= 0 .or. wt_a /= 0.01_wp*tau .or. &
       rng_a%state == rng_zero%state) error stop 'fixed boundary includes crossing photon'
   write(*,'(a)') 'MC stationary warm-start tests passed'
end program test_mc_stationary
