program main_mc
   ! Standalone single-grain MC.  Uses the thread-safe mc_engine API with
   ! a single mc_grain_t and rng_t.  Inputs are a Fortran namelist file.
   !
   ! Example input.nml:
   !   &mc_input
   !     a_um = 0.01, U_isrf = 1.0, comp = 'gra_dl01'
   !     N_events = 50000, t_max = 0.0
   !     lam_c_um = 1000.0, seed = 1, out_prefix = 'output/run01'
   !   /
   ! The default records a stationary distribution after equilibration.
   ! An explicit positive T_init requests a transient experiment.

   use iso_fortran_env, only: error_unit
   use ieee_arithmetic, only: ieee_is_finite
   use constants,      only: wp, pi
   use mc_rng,         only: rng_t, rng_init
   use mc_grain_type,  only: mc_grain_t
   use mc_engine,      only: MC_OK, grain_setup_graphite, mc_run_engine, &
                             mc_run_engine_2pass, mc_run_engine_buffered, &
                             mc_traj_enable, mc_traj_disable, mc_traj_get, &
                             grain_radiative_equilibrium, mc_equilibrate
   implicit none
   integer :: engine_status
   real(wp) :: omitted_time, floor_time, outside_time, stored_energy, continuous_energy, absorbed_energy
   integer(kind=8) :: sample_count, samples_replay
   type(rng_t) :: rng_record
   real(wp) :: T_record, T_eq, p_abs, tau_relax, warm_time
   integer :: warm_events, max_samp

   integer, parameter :: NHIST = 600

   real(wp)            :: a_um, U_isrf, T_init, t_max, lam_c_um
   integer             :: N_events, seed
   character(len=16)   :: comp, mc_engine_kind
   character(len=128)  :: out_prefix
   logical             :: record_trajectory
   namelist /mc_input/ a_um, U_isrf, comp, T_init, N_events, t_max, &
                       lam_c_um, seed, out_prefix, mc_engine_kind, &
                       record_trajectory

   character(len=512) :: nml_path
   integer            :: nargs, u, ios, k
   type(mc_grain_t)   :: grain
   type(rng_t)        :: rng
   real(wp) :: T_edges(NHIST+1), dP_dT(NHIST), dP_dlnT(NHIST)
   real(wp) :: t_total, e_abs_tot, e_emit_tot
   real(wp), allocatable :: t_rec(:), temp_pre_rec(:), temp_rec(:), dt_rec(:)
   integer  :: n_rec
   ! Sub-step trajectory buffer (filled by cool_segment via mc_engine globals)
   integer,  parameter :: TRAJ_CAP = 200000
   real(wp), allocatable :: traj_time(:), traj_temp(:)
   integer  :: traj_n, traj_capacity
   real(wp) :: traj_stride

   a_um            = 0.01_wp
   U_isrf          = 1.0_wp
   comp            = 'gra_dl01'
   T_init          = -1.0_wp         ! negative: stationary; positive: transient
   N_events        = 50000
   t_max           = 0.0_wp
   lam_c_um        = 1000.0_wp
   seed            = 1
   out_prefix      = 'output/run01'
   mc_engine_kind  = '2pass'        ! '2pass' | 'buffered' | 'fixed'
   record_trajectory = .false.      ! .true. -> additionally write _Tt.dat
                                     ! (forces engine='fixed' for the run)

   nargs = command_argument_count()
   if (nargs < 1) then
      write(*,'(a)') 'Usage: main_mc.x <input.nml>'
      stop 1
   end if
   call get_command_argument(1, nml_path)
   open(newunit=u, file=trim(nml_path), status='old', action='read', iostat=ios)
   if (ios /= 0) then
      write(error_unit,'(a,a)') 'Cannot open namelist: ', trim(nml_path)
      error stop 1
   end if
   read(u, nml=mc_input, iostat=ios)
   close(u)
   if (ios /= 0) then
      write(error_unit,'(a,i0)') 'Invalid &mc_input namelist, iostat=', ios
      error stop 1
   end if
   if (.not. ieee_is_finite(a_um) .or. a_um <= 0.0_wp) error stop 'a_um must be finite and positive'
   if (.not. ieee_is_finite(U_isrf) .or. U_isrf < 0.0_wp) error stop 'U_isrf must be finite and nonnegative'
   if (.not. ieee_is_finite(lam_c_um) .or. lam_c_um <= 0.0_wp) &
      error stop 'lam_c_um must be finite and positive'
   if (.not. ieee_is_finite(t_max) .or. t_max < 0.0_wp) error stop 't_max must be finite and nonnegative'
   if (N_events <= 0) error stop 'N_events must be positive'
   if (len_trim(out_prefix) == 0) error stop 'out_prefix must not be empty'
   if (.not. ieee_is_finite(T_init)) error stop 'T_init must be finite'
   if (T_init == 0.0_wp) then
      write(error_unit,'(a)') 'T_init must be positive for a transient or negative for stationary sampling'
      error stop 1
   end if
   select case (trim(mc_engine_kind))
   case ('fixed', '2pass', 'buffered')
   case default
      write(error_unit,'(a,a)') 'Unknown mc_engine_kind: ', trim(mc_engine_kind)
      error stop 1
   end select

   write(*,'(a)') '== mc_pT input =='
   write(*,'(a,es12.4,a)') '  a_um     = ', a_um, ' um'
   write(*,'(a,es12.4)')   '  U_isrf   = ', U_isrf
   write(*,'(a,a)')        '  comp     = ', trim(comp)
   write(*,'(a,i12)')      '  N_events = ', N_events
   write(*,'(a,i0)')       '  seed     = ', seed

   call rng_init(rng, seed)
   call grain_setup_graphite(grain, a_um, U_isrf, trim(comp), lam_c_in=lam_c_um)

   if (T_init < 0.0_wp) then
      call grain_radiative_equilibrium(grain, T_eq, p_abs, tau_relax, engine_status)
      if (engine_status /= MC_OK) error stop 'Radiative equilibrium failed'
      call mc_equilibrate(grain, rng, T_eq, T_record, warm_time, warm_events, engine_status)
      if (engine_status /= MC_OK) error stop 'MC warmup failed'
      write(*,'(a,es24.16e3,a,es24.16e3,a,es24.16e3,a,es24.16e3,a,i0)') &
         '[mc_warmup] radius_um=', a_um, ' initial_Teq_K=', T_eq, &
         ' record_T_K=', T_record, ' warm_time_s=', warm_time, ' events=', warm_events
      write(*,'(a)') '  sampling_mode = stationary; elapsed time starts at recording'
   else
      T_record = T_init
      write(*,'(a,es24.16e3,a)') '  sampling_mode = transient; T_init=', T_init, ' K'
   end if
   rng_record = rng

   write(*,'(a)') '== setup =='
   if (record_trajectory) then
      mc_engine_kind = 'fixed'
      write(*,'(a)') '  record_trajectory=.true. -> engine forced to fixed for _Tt.dat output'
      ! Reserve the time-sampled records plus three forced records for
      ! each event (segment start, pre-absorption, and post-absorption).
      ! An exhausted recording capacity is an explicit engine failure.
      if (N_events > (huge(traj_capacity)-TRAJ_CAP-4)/3) error stop 'Trajectory capacity exceeds integer range'
      traj_capacity = TRAJ_CAP + 3*N_events + 4
      if (t_max > 0.0_wp) then
         traj_stride = t_max / real(TRAJ_CAP, wp)
      else
         traj_stride = 0.0_wp     ! record every sub-step (capacity-limited)
      end if
      call mc_traj_enable(traj_stride, traj_capacity)
   end if
   write(*,'(a,es12.4)') '  rate_event = ', grain%rate_event
   write(*,'(a,es12.4)') '  H_cont     = ', grain%H_cont
   write(*,'(a,a)')      '  engine     = ', trim(mc_engine_kind)

   select case (trim(mc_engine_kind))
   case ('2pass')
      call mc_run_engine_2pass(grain, rng, T_record, N_events, -1, t_max, &
                               NHIST, T_edges, dP_dT, dP_dlnT, &
                               t_total, e_abs_tot, e_emit_tot, &
                               status_out=engine_status, omitted_time_out=omitted_time, &
                               floor_time_out=floor_time, out_of_grid_time_out=outside_time, &
                               sample_count_out=sample_count, stored_energy_out=stored_energy)
   case ('buffered')
      call mc_run_engine_2pass(grain, rng, T_record, N_events, -1, t_max, &
                               NHIST, T_edges, dP_dT, dP_dlnT, &
                               t_total, e_abs_tot, e_emit_tot, &
                               status_out=engine_status, omitted_time_out=omitted_time, &
                               sample_count_out=samples_replay)
      if (engine_status /= MC_OK) error stop 'Buffered capacity measurement failed'
      if (samples_replay <= 0_8 .or. samples_replay > 1000000000_8) &
         error stop 'Buffered sample capacity exceeds 1e9 samples (16 GB)'
      max_samp = int(samples_replay)
      rng = rng_record
      call mc_run_engine_buffered(grain, rng, T_record, N_events, t_max, &
                                  NHIST, max_samp, &
                                  T_edges, dP_dT, dP_dlnT, &
                                  t_total, e_abs_tot, e_emit_tot, &
                               status_out=engine_status, omitted_time_out=omitted_time, &
                               floor_time_out=floor_time, out_of_grid_time_out=outside_time, &
                               sample_count_out=sample_count, stored_energy_out=stored_energy)
   case ('fixed')
      if (record_trajectory) then
         allocate(t_rec(N_events+1), temp_pre_rec(N_events+1), &
                  temp_rec(N_events+1), dt_rec(N_events+1))
         call mc_run_engine(grain, rng, T_record, N_events, t_max, &
                            NHIST, T_edges, dP_dT, dP_dlnT, &
                            t_total, e_abs_tot, e_emit_tot, &
                            time_rec_out=t_rec, temp_pre_rec_out=temp_pre_rec, &
                            temp_rec_out=temp_rec, &
                            dtime_rec_out=dt_rec, n_rec_out=n_rec, &
                            status_out=engine_status, omitted_time_out=omitted_time, &
                               floor_time_out=floor_time, out_of_grid_time_out=outside_time, &
                               sample_count_out=sample_count, stored_energy_out=stored_energy)
      else
         call mc_run_engine(grain, rng, T_record, N_events, t_max, &
                            NHIST, T_edges, dP_dT, dP_dlnT, &
                            t_total, e_abs_tot, e_emit_tot, &
                               status_out=engine_status, omitted_time_out=omitted_time, &
                               floor_time_out=floor_time, out_of_grid_time_out=outside_time, &
                               sample_count_out=sample_count, stored_energy_out=stored_energy)
      end if
   end select

   if (engine_status /= MC_OK) then
      write(error_unit,'(a,i0,a,es24.16e3)') 'MC cooling failed: status=', engine_status, &
         ' incomplete_time_s=', omitted_time
      error stop 'MC interval incomplete; no distribution is written'
   end if

   if (trim(mc_engine_kind) == 'buffered') then
      if (sample_count /= samples_replay) error stop 'Buffered sample count differs from replay'
   end if
   write(*,'(a,3(es24.16e3,1x),i0)') '[mc_interval] floor_s outside_grid_s omitted_s samples: ', &
      floor_time, outside_time, omitted_time, sample_count
   write(*,'(a)') '== run =='
   write(*,'(a,es12.4)') '  t_total    = ', t_total
   write(*,'(a,es12.4)') '  e_abs_tot  = ', e_abs_tot
   write(*,'(a,es12.4)') '  e_emit_tot = ', e_emit_tot
   if (e_abs_tot > 0.0_wp) then
      write(*,'(a,es14.6e3)') '  emit/event_abs = ', e_emit_tot/e_abs_tot
   end if

   continuous_energy = 4.0_wp*pi*(a_um*1.0e-4_wp)**2 * grain%H_cont * t_total
   absorbed_energy = e_abs_tot + continuous_energy
   write(*,'(a)') '# mc_energy_single: event_abs_erg continuous_abs_erg emitted_erg stored_erg ' // &
                 'elapsed_s closure_residual_erg histogram_mass'
   write(*,'(a,7(es24.16e3,1x))') '[mc_energy_single] ', e_abs_tot, continuous_energy, &
      e_emit_tot, stored_energy, t_total, absorbed_energy-e_emit_tot-stored_energy, &
      sum(dP_dT*(T_edges(2:)-T_edges(:NHIST)))

   open(newunit=u, file=trim(out_prefix)//'_PT.dat', status='replace', action='write')
   write(u, '(a)') '# columns: T_lo[K]  T_hi[K]  T_mid[K]  dP/dT[1/K]  dP/dlnT'
   do k = 1, NHIST
      write(u, '(5(es14.5e3,1x))') T_edges(k), T_edges(k+1), &
            sqrt(T_edges(k)*T_edges(k+1)), dP_dT(k), dP_dlnT(k)
   end do
   close(u)

   if (record_trajectory) then
      ! ---- write the boundary record of each event (T_pre, T_post) ----
      open(newunit=u, file=trim(out_prefix)//'_evt.dat', status='replace', action='write')
      write(u, '(a)') '# columns: t[s]  T_pre[K]  T_post[K]  dt_event[s]'
      write(u, '(a)') '#   T_pre  : cooling-segment endpoint (just BEFORE the photon jump)'
      write(u, '(a)') '#   T_post : grain T immediately AFTER the photon jump'
      do k = 1, n_rec
         write(u, '(4(es14.7e3,1x))') t_rec(k), temp_pre_rec(k), temp_rec(k), dt_rec(k)
      end do
      close(u)
      write(*,'(a,i0,a)') '  wrote ', n_rec, ' event-boundary records to _evt.dat'

      ! ---- write the actual sub-step trajectory (stride-sampled) ----
      allocate(traj_time(traj_capacity), traj_temp(traj_capacity))
      call mc_traj_get(traj_n, traj_time, traj_temp)
      open(newunit=u, file=trim(out_prefix)//'_Tt.dat', status='replace', action='write')
      write(u, '(a)')   '# Sub-step trajectory dumped from cool_segment.'
      write(u, '(a,es10.3,a)') '# stride = ', traj_stride, ' s (min time between samples)'
      write(u, '(a)')   '# columns: t[s]   T[K]'
      do k = 1, traj_n
         write(u, '(2(es14.7e3,1x))') traj_time(k), traj_temp(k)
      end do
      close(u)
      write(*,'(a,i0,a)') '  wrote ', traj_n, ' sub-step samples to _Tt.dat'
      call mc_traj_disable()
      deallocate(t_rec, temp_pre_rec, temp_rec, dt_rec, traj_time, traj_temp)
   end if

   write(*,'(a)') '== done =='
end program main_mc
