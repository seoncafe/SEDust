program test_mc_stationary_start
   use constants, only: wp
   use mc_grain_type, only: mc_grain_t
   use mc_rng, only: rng_t, rng_init
   use mc_engine, only: grain_setup_from_cabs, grain_radiative_equilibrium, &
                        mc_equilibrate, MC_OK
   implicit none
   integer, parameter :: nseed=1000
   real(wp), parameter :: lam(6)=[0.1_wp,0.2_wp,1.0_wp,10.0_wp,100.0_wp,1000.0_wp]
   real(wp), parameter :: q(6)=1.0e-3_wp, field(6)=1.0e-8_wp
   type(mc_grain_t) :: grain
   type(rng_t) :: rng
   real(wp) :: temperature, elapsed, teq, p_abs, tau, expected_time
   integer :: seed, events, status, event_count(nseed)

   call grain_setup_from_cabs(grain,3.55e-4_wp,1.0_wp,'ad_s2',lam,q,J_lam_in=field)
   grain%rate_event=1.0e-4_wp
   grain%H_cont=1.0e-3_wp
   call grain_radiative_equilibrium(grain,teq,p_abs,tau,status)
   if(status/=MC_OK) error stop 'radiative equilibrium failed'
   expected_time=max(0.001_wp*tau,4.0_wp/grain%rate_event)
   do seed=1,nseed
      call rng_init(rng,seed)
      call mc_equilibrate(grain,rng,40.0_wp,temperature,elapsed,events,status, &
                          duration_factor=0.001_wp,expected_events=4)
      if(status/=MC_OK .or. elapsed/=expected_time .or. events<0) &
         error stop 'fixed-time warm start failed'
      event_count(seed)=events
   end do
   ! The requested count determines the observation duration, not a minimum
   ! realized event count. Both sides of the Poisson count must occur.
   if(count(event_count<4)==0 .or. count(event_count>4)==0) &
      error stop 'fixed-time warm start waits for an absorption event'
   if(abs(real(sum(event_count),wp)/nseed-4.0_wp)>0.5_wp) &
      error stop 'warm event count inconsistent with a Poisson mean of four'
   write(*,'(a,f9.5,a,i0,a,i0)') 'Fixed-time mean events=', &
      real(sum(event_count),wp)/nseed,' below expected=',count(event_count<4), &
      ' above expected=',count(event_count>4)
end program test_mc_stationary_start
