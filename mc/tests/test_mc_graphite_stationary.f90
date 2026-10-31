program test_mc_graphite_stationary
 use constants, only: wp
 use mc_grain_type, only: mc_grain_t
 use mc_rng, only: rng_t, rng_init
 use mc_engine, only: grain_setup_graphite, grain_radiative_equilibrium, grain_get_QPL, &
                      mc_equilibrate, mc_run_engine, grain_T_from_U, MC_OK
 use mc_heatcap, only: U_of_T
 implicit none
 type(mc_grain_t) :: g
 type(rng_t) :: rng, warm_rng(2,2), trial_rng
 integer, parameter :: nh=600, nb=8, ne=250
 real(wp), parameter :: sigma=5.6703744e-5_wp
 real(wp) :: a,teq,pabs,tau,trec,wt,lo,hi,mid,root,edges(nh+1),dens(nh),ldens(nh)
 real(wp) :: time,ea,ee,stored,meanT(nb),massbelow(nb),coolmass(nb),mass,mean1,mean2,se
 real(wp) :: cdf(nb,nh),ks,maxz
 real(wp) :: warm_T(2,2),warm_elapsed(2,2),edge_compare(nh+1)
 real(wp) :: density_compare(2,nh),log_compare(nh),elapsed_compare,abs_compare,emit_compare
 integer :: warm_count(2,2),init_case
 real(wp) :: tstart
 integer :: isize,iinit,ifactor,nev,status,b,j
 do isize=1,3
   select case(isize)
   case(1); a=0.020_wp
   case(2); a=0.010_wp
   case(3); a=0.005_wp
   end select
   call grain_setup_graphite(g,a,1.0_wp,'gra_dl01',lam_c_in=1000.0_wp)
   call grain_radiative_equilibrium(g,teq,pabs,tau,status)
   if(status/=MC_OK) error stop 'radiative equilibrium failed'
   lo=0.5_wp; hi=teq
   do j=1,100
     mid=0.5_wp*(lo+hi)
     if(sigma*grain_get_QPL(g,mid)*mid**4 < g%H_cont) then
       lo=mid
     else
       hi=mid
     end if
   end do
   root=0.5_wp*(lo+hi)
   print '(a,f6.0,a,es12.4,a,f10.5,a,f10.5,a,es12.4)', &
      'radius_A=',a*1.e4_wp,' rate=',g%rate_event,' Teq_K=',teq,' Tcont_K=',root,' tau_s=',tau
   do iinit=1,2
     do ifactor=1,2
       call rng_init(rng,82)
       tstart=teq
       if(iinit==1) tstart=2.725_wp
       call mc_equilibrate(g,rng,tstart,trec,wt,nev,status, &
            duration_factor=real(ifactor*10,wp))
       if(status/=MC_OK) then
          print *, 'warm failed',isize,iinit,ifactor,status
          error stop
       end if
       warm_T(iinit,ifactor)=trec
       warm_elapsed(iinit,ifactor)=wt
       if(wt/=max(real(ifactor*10,wp)*tau,100.0_wp/g%rate_event)) &
          error stop 'graphite warmup missed fixed endpoint'
       warm_count(iinit,ifactor)=nev
       warm_rng(iinit,ifactor)=rng
       print '(a,f7.3,a,i0,a,f10.5,a,es12.4,a,i0)', &
          ' initial_K=',tstart,' factor=',ifactor*10,' record_K=',trec, &
          ' warm_s=',wt,' events=',nev
       if(isize/=1 .or. iinit/=2 .or. ifactor/=1) cycle
       do b=1,nb
         call mc_run_engine(g,rng,trec,ne,0.0_wp,nh,edges,dens,ldens, &
              time,ea,ee,status_out=status,stored_energy_out=stored)
         if(status/=MC_OK) error stop 'recording failed'
         mass=sum(dens*(edges(2:)-edges(:nh)))
         if(abs(mass-1.0_wp)>1.e-8_wp) error stop 'histogram mass failed'
         meanT(b)=sum(0.5_wp*(edges(2:)+edges(:nh))*dens*(edges(2:)-edges(:nh)))
         massbelow(b)=0.0_wp
         coolmass(b)=0.0_wp
         do j=1,nh
           cdf(b,j)=sum(dens(:j)*(edges(2:j+1)-edges(:j)))
           if(edges(j+1)<=root) massbelow(b)=massbelow(b)+dens(j)*(edges(j+1)-edges(j))
           if(edges(j+1)<=teq) coolmass(b)=coolmass(b)+dens(j)*(edges(j+1)-edges(j))
         end do
         trec=grain_T_from_U(g,U_of_T(trec,g%a_um,g%comp)+stored,status)
         if(status/=MC_OK) error stop 'end state inverse failed'
       end do
       mean1=sum(meanT(1:nb/2))/real(nb/2,wp)
       mean2=sum(meanT(nb/2+1:nb))/real(nb/2,wp)
       se=sqrt(sum((meanT(1:nb/2)-mean1)**2)/real((nb/2)*(nb/2-1),wp) + &
               sum((meanT(nb/2+1:nb)-mean2)**2)/real((nb/2)*(nb/2-1),wp))
       print '(a,8f10.5)', ' block_meanT_K=',meanT
       print '(a,f10.5,a,f10.5,a,f10.5,a,f10.5)', &
          ' first_half_K=',mean1,' second_half_K=',mean2,' SE_K=',se, &
          ' z=',abs(mean1-mean2)/se
       if(se <= 0.0_wp .or. abs(mean1-mean2)>3.0_wp*se) &
          error stop 'temperature histogram moment is not stationary within sampled noise'
       print '(a,es12.4)', ' mass_below_cont_root=',sum(massbelow)/nb
       if(sum(massbelow)>1.0e-12_wp) error stop 'histogram extends below continuous-heating root'
       mean1=sum(coolmass(1:nb/2))/real(nb/2,wp)
       mean2=sum(coolmass(nb/2+1:nb))/real(nb/2,wp)
       se=sqrt(sum((coolmass(1:nb/2)-mean1)**2)/real((nb/2)*(nb/2-1),wp) + &
               sum((coolmass(nb/2+1:nb)-mean2)**2)/real((nb/2)*(nb/2-1),wp))
       print '(a,8f10.5)', ' block_P_below_Teq=',coolmass
       print '(a,f10.5,a,f10.5,a,f10.5,a,f10.5)', &
          ' first_half_P=',mean1,' second_half_P=',mean2,' SE_P=',se, &
          ' z=',abs(mean1-mean2)/se
       if(se <= 0.0_wp .or. abs(mean1-mean2)>3.0_wp*se) &
          error stop 'temperature histogram occupancy is not stationary within sampled noise'
       ks=0.0_wp
       maxz=0.0_wp
       do j=1,nh
         mean1=sum(cdf(1:nb/2,j))/real(nb/2,wp)
         mean2=sum(cdf(nb/2+1:nb,j))/real(nb/2,wp)
         ks=max(ks,abs(mean1-mean2))
         se=sqrt(sum((cdf(1:nb/2,j)-mean1)**2)/real((nb/2)*(nb/2-1),wp) + &
                 sum((cdf(nb/2+1:nb,j)-mean2)**2)/real((nb/2)*(nb/2-1),wp))
         if(se>0.0_wp) then
           maxz=max(maxz,abs(mean1-mean2)/se)
         else if(abs(mean1-mean2)>1.0e-12_wp) then
           error stop 'CDF halves differ without block variation'
         end if
       end do
       print '(a,f10.5,a,f10.5)', ' CDF_KS=',ks,' max_block_SE_units=',maxz
       if(maxz>4.0_wp) &
          error stop 'full CDF halves differ beyond block-estimated sampling noise'
     end do
   end do
   do ifactor=1,2
     if(warm_elapsed(1,ifactor)/=warm_elapsed(2,ifactor) .or. &
        warm_count(1,ifactor)/=warm_count(2,ifactor) .or. &
        warm_rng(1,ifactor)%state/=warm_rng(2,ifactor)%state) &
        error stop 'warm trajectory schedule depends on initial temperature'
     print '(a,f6.0,a,i0,a,es12.4)', ' radius_A=',a*1.e4_wp, &
          ' factor=',ifactor*10,' warm_T_difference_K=',abs(warm_T(1,ifactor)-warm_T(2,ifactor))
     if(abs(warm_T(1,ifactor)-warm_T(2,ifactor))>1.e-6_wp*teq) &
        error stop 'CMB and radiative equilibrium warm starts have different final temperatures'
     if(warm_T(1,ifactor)/=warm_T(2,ifactor)) then
       do init_case=1,2
         trial_rng=warm_rng(init_case,ifactor)
         call mc_run_engine(g,trial_rng,warm_T(init_case,ifactor),ne,0.0_wp,nh, &
              edge_compare,density_compare(init_case,:),log_compare,elapsed_compare, &
              abs_compare,emit_compare,status_out=status)
         if(status/=MC_OK) error stop 'warm initial-state recording failed'
       end do
       cdf(1,:)=0.0_wp
       cdf(2,:)=0.0_wp
       do j=1,nh
         cdf(1,j)=sum(density_compare(1,:j)*(edge_compare(2:j+1)-edge_compare(:j)))
         cdf(2,j)=sum(density_compare(2,:j)*(edge_compare(2:j+1)-edge_compare(:j)))
       end do
       if(maxval(abs(cdf(1,:)-cdf(2,:)))>1.e-4_wp) &
          error stop 'recorded CDF depends on warmed initial temperature'
     end if
   end do
 end do
 print '(a)', 'MC graphite stationary diagnostics passed'
end program
