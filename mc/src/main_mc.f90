program main_mc
   ! Standalone single-grain MC.  Uses the thread-safe mc_engine API with
   ! a single mc_grain_t and rng_t.  Inputs are a Fortran namelist file.
   !
   ! Example input.nml:
   !   &mc_input
   !     a_um = 0.01, U_isrf = 1.0, comp = 'gra_dl01'
   !     T_init = 2.725, N_events = 50000, t_max = 0.0
   !     lam_c_um = 1000.0, seed = 1, out_prefix = 'output/run01'
   !   /
   ! Defaults: T_init = T_CMB = 2.725 K (grain starts in equilibrium with
   ! the cosmic background; faster convergence than starting too hot).

   use constants,      only: wp
   use mc_rng,         only: rng_t, rng_init
   use mc_grain_type,  only: mc_grain_t
   use mc_engine,      only: grain_setup_graphite, mc_run_engine, &
                             mc_run_engine_2pass, mc_run_engine_buffered, &
                             mc_traj_enable, mc_traj_disable, mc_traj_get
   implicit none

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
   integer  :: traj_n
   real(wp) :: traj_stride

   a_um            = 0.01_wp
   U_isrf          = 1.0_wp
   comp            = 'gra_dl01'
   T_init          = 2.725_wp        ! CMB temperature (K)
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
   if (ios /= 0) stop 'Cannot open namelist'
   read(u, nml=mc_input, iostat=ios)
   close(u)

   write(*,'(a)') '== mc_pT input =='
   write(*,'(a,es12.4,a)') '  a_um     = ', a_um, ' um'
   write(*,'(a,es12.4)')   '  U_isrf   = ', U_isrf
   write(*,'(a,a)')        '  comp     = ', trim(comp)
   write(*,'(a,i12)')      '  N_events = ', N_events
   write(*,'(a,i0)')       '  seed     = ', seed

   call rng_init(rng, seed)
   call grain_setup_graphite(grain, a_um, U_isrf, trim(comp), lam_c_in=lam_c_um)

   write(*,'(a)') '== setup =='
   if (record_trajectory) then
      mc_engine_kind = 'fixed'
      write(*,'(a)') '  record_trajectory=.true. -> engine forced to fixed for _Tt.dat output'
      ! Enable the sub-step trajectory buffer.  Stride = t_max / TRAJ_CAP
      ! so the buffer covers the full window at uniform time resolution.
      if (t_max > 0.0_wp) then
         traj_stride = t_max / real(TRAJ_CAP, wp)
      else
         traj_stride = 0.0_wp     ! record every sub-step (capacity-limited)
      end if
      call mc_traj_enable(traj_stride, TRAJ_CAP)
   end if
   write(*,'(a,es12.4)') '  rate_event = ', grain%rate_event
   write(*,'(a,es12.4)') '  H_cont     = ', grain%H_cont
   write(*,'(a,a)')      '  engine     = ', trim(mc_engine_kind)

   select case (trim(mc_engine_kind))
   case ('2pass')
      call mc_run_engine_2pass(grain, rng, T_init, N_events, -1, -1, t_max, &
                               NHIST, T_edges, dP_dT, dP_dlnT, &
                               t_total, e_abs_tot, e_emit_tot)
   case ('buffered')
      call mc_run_engine_buffered(grain, rng, T_init, N_events, -1, t_max, &
                                  NHIST, 200000, &
                                  T_edges, dP_dT, dP_dlnT, &
                                  t_total, e_abs_tot, e_emit_tot)
   case default       ! 'fixed'
      if (record_trajectory) then
         allocate(t_rec(N_events+1), temp_pre_rec(N_events+1), &
                  temp_rec(N_events+1), dt_rec(N_events+1))
         call mc_run_engine(grain, rng, T_init, N_events, t_max, &
                            NHIST, T_edges, dP_dT, dP_dlnT, &
                            t_total, e_abs_tot, e_emit_tot, &
                            time_rec_out=t_rec, temp_pre_rec_out=temp_pre_rec, &
                            temp_rec_out=temp_rec, &
                            dtime_rec_out=dt_rec, n_rec_out=n_rec)
      else
         call mc_run_engine(grain, rng, T_init, N_events, t_max, &
                            NHIST, T_edges, dP_dT, dP_dlnT, &
                            t_total, e_abs_tot, e_emit_tot)
      end if
   end select

   write(*,'(a)') '== run =='
   write(*,'(a,es12.4)') '  t_total    = ', t_total
   write(*,'(a,es12.4)') '  e_abs_tot  = ', e_abs_tot
   write(*,'(a,es12.4)') '  e_emit_tot = ', e_emit_tot
   if (e_abs_tot > 0.0_wp) then
      write(*,'(a,f8.4)') '  emit/abs   = ', e_emit_tot/e_abs_tot
   end if

   open(newunit=u, file=trim(out_prefix)//'_PT.dat', status='replace', action='write')
   write(u, '(a)') '# columns: T_lo[K]  T_hi[K]  T_mid[K]  dP/dT[1/K]  dP/dlnT'
   do k = 1, NHIST
      write(u, '(5(es14.5e3,1x))') T_edges(k), T_edges(k+1), &
            sqrt(T_edges(k)*T_edges(k+1)), dP_dT(k), dP_dlnT(k)
   end do
   close(u)

   if (record_trajectory) then
      ! ---- write the per-event boundary record (T_pre, T_post) ----
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
      allocate(traj_time(TRAJ_CAP), traj_temp(TRAJ_CAP))
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
