program test_mc_cooling_interval
   use constants, only: wp
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_finite
   use mc_grain_type, only: mc_grain_t
   use mc_heatcap, only: U_of_T
   use mc_rng, only: rng_t, rng_init
   use mc_engine, only: grain_setup_from_cabs, grain_T_from_U, mc_run_engine, &
                        mc_run_engine_2pass, mc_run_engine_buffered, &
                        hc_erg_um, MC_OK, MC_STEP_BUDGET, MC_INVALID_STEP, MC_INVALID_WAIT, &
                        MC_HISTOGRAM_RANGE, MC_SAMPLE_CAPACITY, &
                        mc_traj_enable, mc_traj_disable, mc_traj_get
   implicit none
   integer, parameter :: nh = 80, nr = 8
   real(wp), parameter :: lam(6) = [0.1_wp, 0.2_wp, 1.0_wp, 10.0_wp, 100.0_wp, 1000.0_wp]
   real(wp), parameter :: q(6) = 1.0e-3_wp, field(6) = 1.0e-8_wp
   type(mc_grain_t) :: g

   ! Synthetic optical data exercise the actual 3.55e-4 um radius.
   ! A 1 um photon has energy 1.239841984332 eV. This checks the
   ! erg*um constant shared by event energies and the wavelength split.
   if (abs(hc_erg_um / 1.602176634e-12_wp - 1.239841984332_wp) > 1.0e-10_wp) &
      error stop 'Photon energy and wavelength units are inconsistent'
   call grain_setup_from_cabs(g, 3.55e-4_wp, 1.0_wp, 'ad_s2', lam, q, J_lam_in=field)
   call compare_case('cooling', g, 100.0_wp, 1.0e-4_wp, 0.0_wp, 0.0_wp, 3)
   call compare_case('continuous heating', g, 10.0_wp, 1.0e-4_wp, 1.0e-3_wp, 0.0_wp, 3)
   call compare_case('floor residence', g, 0.5_wp, 1.0e-8_wp, 0.0_wp, 0.0_wp, 1)
   call compare_case('long wait', g, 0.5_wp, 1.0e-12_wp, 0.0_wp, 0.0_wp, 1)
   call compare_case('finite time boundary', g, 100.0_wp, 1.0e-12_wp, 0.0_wp, 100.0_wp, 3)
   call compare_case('zero event rate', g, 100.0_wp, 0.0_wp, 0.0_wp, 100.0_wp, 3)
   call check_enthalpy_inverse()
   call check_discrete_enthalpy_target()
   call check_budget(g)
   call check_invalid(g)
   call check_histogram_failures(g)
   call check_trajectory_capacity(g)
   write(*,'(a)') 'MC cooling interval tests passed'

contains

   subroutine compare_case(label, grain, initial, rate, heating, limit, events)
      character(len=*), intent(in) :: label
      type(mc_grain_t), intent(in) :: grain
      real(wp), intent(in) :: initial, rate, heating, limit
      integer, intent(in) :: events
      type(mc_grain_t) :: local
      type(rng_t) :: rng_a, rng_b
      real(wp) :: edge_a(nh+1), edge_b(nh+1), density_a(nh), density_b(nh)
      real(wp) :: log_a(nh), log_b(nh), elapsed_a, elapsed_b
      real(wp) :: absorbed_a, absorbed_b, emitted_a, emitted_b, omitted_a, omitted_b
      real(wp) :: time_a(nr), time_b(nr), pre_a(nr), pre_b(nr)
      real(wp) :: temp_a(nr), temp_b(nr), dtime_a(nr), dtime_b(nr)
      real(wp) :: mass_a, mass_b, floor_a, floor_b, outside_a, outside_b
      real(wp) :: elapsed_fixed, absorbed_fixed, emitted_fixed, floor_fixed, stored_fixed
      real(wp) :: stored_a, stored_b
      integer(kind=8) :: samples_a, samples_b, samples_fixed, state_fixed
      integer :: mode, status_a, status_b, count_a, count_b

      samples_fixed = -1_8
      state_fixed = -1_8
      local = grain
      local%rate_event = rate
      local%H_cont = heating
      do mode = 1, 3
         call rng_init(rng_a, 271828)
         call rng_init(rng_b, 271828)
         select case (mode)
         case (1)
            call mc_run_engine(local, rng_a, initial, events, limit, nh, edge_a, density_a, log_a, &
                 elapsed_a, absorbed_a, emitted_a, time_a, pre_a, temp_a, dtime_a, count_a, &
                 status_out=status_a, chunk_steps=1, omitted_time_out=omitted_a, &
                 floor_time_out=floor_a, out_of_grid_time_out=outside_a, sample_count_out=samples_a, &
                 stored_energy_out=stored_a)
            call mc_run_engine(local, rng_b, initial, events, limit, nh, edge_b, density_b, log_b, &
                 elapsed_b, absorbed_b, emitted_b, time_b, pre_b, temp_b, dtime_b, count_b, &
                 status_out=status_b, chunk_steps=4000, omitted_time_out=omitted_b, &
                 floor_time_out=floor_b, out_of_grid_time_out=outside_b, sample_count_out=samples_b, &
                 stored_energy_out=stored_b)
            if (count_a /= count_b) error stop 'fixed record count differs'
            if (any(temp_a(1:count_a) < pre_a(1:count_a))) &
               error stop 'photon event lowered grain temperature'
            call same_vector('event time', time_a(1:count_a), time_b(1:count_b))
            call same_vector('pre-event T', pre_a(1:count_a), pre_b(1:count_b))
            call same_vector('post-event T', temp_a(1:count_a), temp_b(1:count_b))
            call same_vector('event interval', dtime_a(1:count_a), dtime_b(1:count_b))
         case (2)
            call mc_run_engine_2pass(local, rng_a, initial, events, events, limit, nh, &
                 edge_a, density_a, log_a, elapsed_a, absorbed_a, emitted_a, &
                 status_out=status_a, chunk_steps=2, omitted_time_out=omitted_a, &
                 floor_time_out=floor_a, out_of_grid_time_out=outside_a, sample_count_out=samples_a, &
                 stored_energy_out=stored_a)
            call mc_run_engine_2pass(local, rng_b, initial, events, events, limit, nh, &
                 edge_b, density_b, log_b, elapsed_b, absorbed_b, emitted_b, &
                 status_out=status_b, chunk_steps=300, omitted_time_out=omitted_b, &
                 floor_time_out=floor_b, out_of_grid_time_out=outside_b, sample_count_out=samples_b, &
                 stored_energy_out=stored_b)
         case (3)
            call mc_run_engine_buffered(local, rng_a, initial, events, limit, nh, 100000, &
                 edge_a, density_a, log_a, elapsed_a, absorbed_a, emitted_a, &
                 status_out=status_a, chunk_steps=2, omitted_time_out=omitted_a, &
                 floor_time_out=floor_a, out_of_grid_time_out=outside_a, sample_count_out=samples_a, &
                 stored_energy_out=stored_a)
            call mc_run_engine_buffered(local, rng_b, initial, events, limit, nh, 100000, &
                 edge_b, density_b, log_b, elapsed_b, absorbed_b, emitted_b, &
                 status_out=status_b, chunk_steps=300, omitted_time_out=omitted_b, &
                 floor_time_out=floor_b, out_of_grid_time_out=outside_b, sample_count_out=samples_b, &
                 stored_energy_out=stored_b)
         end select
         if (status_a /= MC_OK .or. status_b /= MC_OK) then
            write(*,*) label, mode, status_a, status_b
            error stop 'interval integration failed'
         end if
         if (rng_a%state /= rng_b%state) error stop 'RNG state differs'
         if (mode == 1) then
            elapsed_fixed = elapsed_a
            absorbed_fixed = absorbed_a
            emitted_fixed = emitted_a
            floor_fixed = floor_a
            stored_fixed = stored_a
            samples_fixed = samples_a
            state_fixed = rng_a%state
         else
            call same_scalar('method elapsed time', elapsed_a, elapsed_fixed)
            call same_scalar('method absorbed energy', absorbed_a, absorbed_fixed)
            call same_scalar('method emitted energy', emitted_a, emitted_fixed)
            call same_scalar('method floor time', floor_a, floor_fixed)
            call same_scalar('method stored energy', stored_a, stored_fixed)
            if (samples_a /= samples_fixed) error stop 'method sample counts differ'
            if (rng_a%state /= state_fixed) error stop 'method RNG states differ'
         end if
         call same_scalar('elapsed time', elapsed_a, elapsed_b)
         call same_scalar('absorbed energy', absorbed_a, absorbed_b)
         call same_scalar('emitted energy', emitted_a, emitted_b)
         call same_scalar('omitted time', omitted_a, omitted_b)
         call same_scalar('floor time', floor_a, floor_b)
         call same_scalar('stored energy', stored_a, stored_b)
         call same_scalar('out-of-grid time', outside_a, outside_b)
         if (samples_a /= samples_b) error stop 'sample counts differ across chunks'
         if (samples_a <= 0_8) error stop 'missing cooling samples'
         call same_vector('bin edges', edge_a, edge_b)
         call same_vector('temperature density', density_a, density_b)
         call same_vector('log temperature density', log_a, log_b)
         mass_a = sum(density_a * (edge_a(2:) - edge_a(:nh)))
         mass_b = sum(density_b * (edge_b(2:) - edge_b(:nh)))
         call same_scalar('histogram mass', mass_a, mass_b)
         if (abs(mass_a - 1.0_wp) > 1.0e-9_wp) then
            write(*,*) label, mode, mass_a
            error stop 'histogram lost interval residence'
         end if
         if (omitted_a /= 0.0_wp) error stop 'successful interval has omitted time'
         if (outside_a /= 0.0_wp) error stop 'successful interval has out-of-grid time'
         if (mode == 1) then
            call same_scalar('fixed lower edge', edge_a(1), 0.5_wp)
            call same_scalar('fixed upper edge', edge_a(nh+1), 5000.0_wp)
         end if
         if (initial == 0.5_wp .and. heating == 0.0_wp) &
            call same_scalar('floor residence', floor_a, elapsed_a)
         if (limit > 0.0_wp) then
            call same_scalar('time boundary', elapsed_a, limit)
            if (absorbed_a /= 0.0_wp) error stop 'absorption after time boundary'
         end if
      end do
   end subroutine compare_case

   subroutine check_budget(grain)
      type(mc_grain_t), intent(in) :: grain
      type(mc_grain_t) :: local
      type(rng_t) :: rng
      real(wp) :: edges(nh+1), density(nh), log_density(nh)
      real(wp) :: elapsed, absorbed, emitted, omitted
      integer :: mode, status
      local = grain
      local%rate_event = 1.0e-12_wp
      local%H_cont = 0.0_wp
      do mode = 1, 3
         call rng_init(rng, 314159)
         select case (mode)
         case (1)
            call mc_run_engine(local, rng, 100.0_wp, 1, 0.0_wp, nh, edges, density, &
                 log_density, elapsed, absorbed, emitted, status_out=status, &
                 max_steps_total=1, chunk_steps=1, omitted_time_out=omitted)
         case (2)
            call mc_run_engine_2pass(local, rng, 100.0_wp, 1, 1, 0.0_wp, nh, &
                 edges, density, log_density, elapsed, absorbed, emitted, &
                 status_out=status, max_steps_total=1, chunk_steps=1, omitted_time_out=omitted)
         case (3)
            call mc_run_engine_buffered(local, rng, 100.0_wp, 1, 0.0_wp, nh, 10000, &
                 edges, density, log_density, elapsed, absorbed, emitted, &
                 status_out=status, max_steps_total=1, chunk_steps=1, omitted_time_out=omitted)
         end select
         if (status /= MC_STEP_BUDGET) then
            write(*,*) mode, status
            error stop 'step budget was not reported'
         end if
         if (absorbed /= 0.0_wp) error stop 'absorption after budget failure'
         if (omitted <= 0.0_wp) then
            write(*,*) 'mode, omitted, elapsed:', mode, omitted, elapsed
            error stop 'missing omitted interval remainder'
         end if
      end do

   end subroutine check_budget

   subroutine check_invalid(grain)
      type(mc_grain_t), intent(in) :: grain
      type(mc_grain_t) :: local
      type(rng_t) :: rng
      real(wp) :: edges(nh+1), density(nh), log_density(nh)
      real(wp) :: elapsed, absorbed, emitted, omitted
      integer :: mode, defect, status, expected
      do defect = 1, 3
         local = grain
         local%rate_event = 1.0e-4_wp
         local%H_cont = 0.0_wp
         expected = MC_INVALID_STEP
         select case (defect)
         case (1)
            local%H_cont = ieee_value(0.0_wp, ieee_quiet_nan)
         case (2)
            local%QPL_grid = ieee_value(0.0_wp, ieee_quiet_nan)
         case (3)
            local%rate_event = ieee_value(0.0_wp, ieee_quiet_nan)
            expected = MC_INVALID_WAIT
         end select
         do mode = 1, 3
            call rng_init(rng, 161803)
            select case (mode)
            case (1)
               call mc_run_engine(local, rng, 100.0_wp, 1, 0.0_wp, nh, edges, density, &
                    log_density, elapsed, absorbed, emitted, status_out=status, &
                    chunk_steps=1, omitted_time_out=omitted)
            case (2)
               call mc_run_engine_2pass(local, rng, 100.0_wp, 1, 1, 0.0_wp, nh, &
                    edges, density, log_density, elapsed, absorbed, emitted, &
                    status_out=status, chunk_steps=1, omitted_time_out=omitted)
            case (3)
               call mc_run_engine_buffered(local, rng, 100.0_wp, 1, 0.0_wp, nh, 10000, &
                    edges, density, log_density, elapsed, absorbed, emitted, &
                    status_out=status, chunk_steps=1, omitted_time_out=omitted)
            end select
            if (status /= expected) then
               write(*,*) 'defect, mode, status:', defect, mode, status
               error stop 'invalid input was not reported'
            end if
            if (absorbed /= 0.0_wp) error stop 'absorption after invalid interval'
            if (elapsed /= 0.0_wp) then
               write(*,*) 'defect, mode, elapsed:', defect, mode, elapsed
               error stop 'invalid interval advanced time'
            end if
         end do
      end do
   end subroutine check_invalid

   subroutine check_histogram_failures(grain)
      type(mc_grain_t), intent(in) :: grain
      type(mc_grain_t) :: local
      type(rng_t) :: rng
      real(wp) :: edges(nh+1), density(nh), log_density(nh)
      real(wp) :: elapsed, absorbed, emitted, outside_time
      integer(kind=8) :: samples
      integer :: status

      local = grain
      local%rate_event = 1.0e-4_wp
      local%H_cont = 0.0_wp
      call rng_init(rng, 271828)
      call mc_run_engine(local, rng, 6000.0_wp, 1, 0.0_wp, nh, edges, density, &
           log_density, elapsed, absorbed, emitted, status_out=status, &
           out_of_grid_time_out=outside_time, sample_count_out=samples)
      if (status /= MC_HISTOGRAM_RANGE .or. outside_time <= 0.0_wp) &
         error stop 'fixed grid accepted residence above 5000 K'
      if (any(density /= 0.0_wp) .or. any(log_density /= 0.0_wp)) &
         error stop 'range failure returned an accepted distribution'
      if (absorbed /= 0.0_wp) error stop 'absorption after range failure'

      ! One cool interval probes the floor; the next follows a photon jump.
      call rng_init(rng, 271828)
      call mc_run_engine_2pass(local, rng, 0.5_wp, 2, 1, 0.0_wp, nh, &
           edges, density, log_density, elapsed, absorbed, emitted, &
           status_out=status, out_of_grid_time_out=outside_time, sample_count_out=samples)
      if (status /= MC_HISTOGRAM_RANGE .or. outside_time <= 0.0_wp) &
         error stop 'truncated probe accepted out-of-grid residence'
      if (any(density /= 0.0_wp) .or. any(log_density /= 0.0_wp)) &
         error stop 'replay range failure returned an accepted distribution'

      call rng_init(rng, 271828)
      call mc_run_engine_buffered(local, rng, 100.0_wp, 1, 0.0_wp, nh, 1, &
           edges, density, log_density, elapsed, absorbed, emitted, &
           status_out=status, sample_count_out=samples)
      if (status /= MC_SAMPLE_CAPACITY .or. samples <= 1_8) &
         error stop 'buffer capacity failure was not reported'
      if (any(density /= 0.0_wp) .or. any(log_density /= 0.0_wp)) &
         error stop 'buffer capacity failure returned an accepted distribution'
   end subroutine check_histogram_failures

   subroutine check_trajectory_capacity(grain)
      ! A full visualization buffer must fail explicitly without changing
      ! the cooling trajectory or the result when recording is disabled.
      type(mc_grain_t), intent(in) :: grain
      type(mc_grain_t) :: local
      type(rng_t) :: rng_base, rng_over, rng_reset
      real(wp) :: edges(nh+1), density(nh), log_density(nh)
      real(wp) :: elapsed_base, absorbed_base, emitted_base
      real(wp) :: elapsed_over, absorbed_over, emitted_over
      real(wp) :: elapsed_reset, absorbed_reset, emitted_reset
      real(wp) :: trace_time(2), trace_temp(2)
      integer :: status_base, status_over, status_reset, n_trace

      local = grain
      local%rate_event = 1.0e-12_wp
      local%H_cont = 0.0_wp
      call rng_init(rng_base, 12345)
      call rng_init(rng_over, 12345)
      call rng_init(rng_reset, 12345)
      call mc_traj_disable()
      call mc_run_engine(local, rng_base, 0.5_wp, 1, 100.0_wp, nh, edges, density, &
                         log_density, elapsed_base, absorbed_base, emitted_base, &
                         status_out=status_base)
      if (status_base /= MC_OK) error stop 'disabled trajectory failed'
      call same_scalar('disabled trajectory histogram mass', &
                       sum(density*(edges(2:)-edges(:nh))), 1.0_wp)

      call mc_traj_enable(0.0_wp, 1)
      call mc_run_engine(local, rng_over, 0.5_wp, 1, 100.0_wp, nh, edges, density, &
                         log_density, elapsed_over, absorbed_over, emitted_over, &
                         status_out=status_over)
      if (status_over /= MC_SAMPLE_CAPACITY) error stop 'trajectory overflow was silent'
      if (any(density /= 0.0_wp) .or. any(log_density /= 0.0_wp)) &
         error stop 'overflow returned an accepted histogram'
      call same_scalar('overflow elapsed time', elapsed_over, elapsed_base)
      call same_scalar('overflow absorbed energy', absorbed_over, absorbed_base)
      call same_scalar('overflow emitted energy', emitted_over, emitted_base)
      if (rng_over%state /= rng_base%state) error stop 'trajectory overflow changed RNG sequence'

      call mc_traj_enable(1.0e30_wp, 2)
      call mc_run_engine(local, rng_reset, 0.5_wp, 1, 100.0_wp, nh, edges, density, &
                         log_density, elapsed_reset, absorbed_reset, emitted_reset, &
                         status_out=status_reset)
      if (status_reset /= MC_OK) error stop 'trajectory overflow flag was not reset'
      call mc_traj_get(n_trace, trace_time, trace_temp)
      if (n_trace /= 2) error stop 'trajectory capacity counted skipped samples'
      call same_scalar('recorded trajectory endpoint', trace_time(2), elapsed_reset)
      call same_scalar('reset elapsed time', elapsed_reset, elapsed_base)
      call same_scalar('reset emitted energy', emitted_reset, emitted_base)
      call mc_traj_disable()
   end subroutine check_trajectory_capacity

   subroutine check_enthalpy_inverse()
      type(mc_grain_t) :: thermal
      real(wp), parameter :: temperatures(10) = [0.5_wp, 0.75_wp, 1.0_wp, 2.725_wp, &
           10.0_wp, 100.0_wp, 300.0_wp, 1000.0_wp, 5000.0_wp, 12000.0_wp]
      real(wp) :: temperature, enthalpy, photon, target, returned_temperature
      real(wp) :: returned_enthalpy, enthalpy_tolerance, temperature_tolerance
      integer :: material, i, j, status

      do material = 1, 4
         select case (material)
         case (1)
            thermal%a_um = 3.55e-4_wp
            thermal%comp = 'ad_s2'
         case (2)
            thermal%a_um = 5.0_wp
            thermal%comp = 'ad_s2'
         case (3)
            thermal%a_um = 5.0e-4_wp
            thermal%comp = 'pah'
         case (4)
            thermal%a_um = 5.0e-3_wp
            thermal%comp = 'pah'
         end select
         thermal%NT = 2000
         if (allocated(thermal%T_grid)) deallocate(thermal%T_grid, thermal%U_grid)
         allocate(thermal%T_grid(thermal%NT), thermal%U_grid(thermal%NT))
         do i = 1, thermal%NT
            thermal%T_grid(i) = (1.0e4_wp)**(real(i-1,wp)/real(thermal%NT-1,wp))
            thermal%U_grid(i) = U_of_T(thermal%T_grid(i), thermal%a_um, thermal%comp)
         end do
         do i = 1, size(temperatures)
            temperature = temperatures(i)
            enthalpy = U_of_T(temperature, thermal%a_um, thermal%comp)
            if (.not. ieee_is_finite(enthalpy) .or. enthalpy <= 0.0_wp) &
               error stop 'invalid analytic grain enthalpy'
            do j = 1, 4
               select case (j)
               case (1)
                  photon = 0.0_wp
               case (2)
                  photon = 16.0_wp * spacing(enthalpy)
               case (3)
                  photon = hc_erg_um / 1000.0_wp
               case (4)
                  photon = hc_erg_um
               end select
               target = enthalpy + photon
               returned_temperature = grain_T_from_U(thermal, target, status_out=status)
               if (status /= MC_OK .or. .not. ieee_is_finite(returned_temperature)) &
                  error stop 'finite enthalpy inversion failed'
               temperature_tolerance = 64.0_wp * epsilon(temperature) * max(1.0_wp, temperature)
               if (returned_temperature < temperature - temperature_tolerance) then
                  write(*,*) 'inverse cooled grain:', material, temperature, photon, returned_temperature
                  error stop 'nonnegative photon energy lowered temperature'
               end if
               returned_enthalpy = U_of_T(returned_temperature, thermal%a_um, thermal%comp)
               enthalpy_tolerance = max(2.0e-12_wp*abs(target), &
                                        256.0_wp*epsilon(target)*abs(target), 1.0e-8_wp*photon)
               if (.not. ieee_is_finite(returned_enthalpy) .or. &
                   abs(returned_enthalpy-target) > enthalpy_tolerance) then
                  write(*,*) 'inverse enthalpy residual:', material, temperature, photon, &
                             returned_enthalpy-target, enthalpy_tolerance
                  error stop 'photon enthalpy jump was not conserved'
               end if
            end do
         end do
         target = 0.99_wp * U_of_T(0.5_wp, thermal%a_um, thermal%comp)
         returned_temperature = grain_T_from_U(thermal, target, status_out=status)
         if (status /= MC_INVALID_STEP .or. returned_temperature /= 0.0_wp) &
            error stop 'enthalpy below temperature floor was accepted'
         target = ieee_value(0.0_wp, ieee_quiet_nan)
         returned_temperature = grain_T_from_U(thermal, target, status_out=status)
         if (status /= MC_INVALID_STEP .or. returned_temperature /= 0.0_wp) &
            error stop 'nonfinite enthalpy was accepted'
      end do
   end subroutine check_enthalpy_inverse

   subroutine check_discrete_enthalpy_target()
      type(mc_grain_t) :: thermal
      real(wp), parameter :: target = 6.1689123306119163e-13_wp
      real(wp), parameter :: lower_temperature = 356.34773826465414_wp
      real(wp) :: upper_temperature, returned_temperature, returned_enthalpy
      real(wp) :: half_enthalpy_step
      integer :: i, status

      ! This target lies between adjacent floating-point temperatures for
      ! the first HD23 astrodust radius; neither endpoint gives exact U.
      thermal%a_um = 10.0_wp**(-3.45_wp)
      thermal%comp = 'ad_s2'
      thermal%NT = 2000
      allocate(thermal%T_grid(thermal%NT), thermal%U_grid(thermal%NT))
      do i = 1, thermal%NT
         thermal%T_grid(i) = (1.0e4_wp)**(real(i-1,wp)/real(thermal%NT-1,wp))
         thermal%U_grid(i) = U_of_T(thermal%T_grid(i), thermal%a_um, thermal%comp)
      end do
      upper_temperature = nearest(lower_temperature, 1.0_wp)
      half_enthalpy_step = 0.5_wp * abs(U_of_T(upper_temperature, thermal%a_um, thermal%comp) - &
                                         U_of_T(lower_temperature, thermal%a_um, thermal%comp))
      returned_temperature = grain_T_from_U(thermal, target, status_out=status)
      if (status /= MC_OK .or. .not. ieee_is_finite(returned_temperature)) &
         error stop 'finite discrete enthalpy target was rejected'
      returned_enthalpy = U_of_T(returned_temperature, thermal%a_um, thermal%comp)
      if (abs(returned_enthalpy-target) > half_enthalpy_step + &
          8.0_wp*epsilon(target)*abs(target)) &
         error stop 'discrete enthalpy inverse missed closest temperature'
   end subroutine check_discrete_enthalpy_target

   subroutine same_scalar(label, a, b)
      character(len=*), intent(in) :: label
      real(wp), intent(in) :: a, b
      if (.not. ieee_is_finite(a) .or. .not. ieee_is_finite(b)) then
         write(*,*) label, a, b
         error stop 'nonfinite scalar'
      end if
      if (abs(a-b) > 2.0e-12_wp * max(1.0e-100_wp, abs(a), abs(b))) then
         write(*,*) label, a, b
         error stop 'scalar mismatch'
      end if
   end subroutine same_scalar

   subroutine same_vector(label, a, b)
      character(len=*), intent(in) :: label
      real(wp), intent(in) :: a(:), b(:)
      integer :: i
      if (size(a) /= size(b)) error stop 'vector length differs'
      do i = 1, size(a)
         call same_scalar(label, a(i), b(i))
      end do
   end subroutine same_vector
end program test_mc_cooling_interval
