program test_photon_support
   use constants, only: wp, pi
   use mc_engine, only: grain_setup_from_cabs, grain_setup_graphite, sample_photon_lam, hc_erg_um
   use mc_grain_type, only: mc_grain_t
   use mc_rng, only: rng_t, rng_init, rng_uniform
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_invalid, ieee_get_flag, ieee_set_flag
   use radfield, only: J_Mathis
   implicit none
   integer, parameter :: n=400000
   real(wp), parameter :: x(9)=[0.1_wp,0.2_wp,0.4_wp,0.8_wp,1.6_wp,3.2_wp,6.4_wp,12.8_wp,25.6_wp]
   real(wp), parameter :: p(9)=[0.0_wp,0.0_wp,1.0_wp,0.0_wp,0.0_wp,0.0_wp,2.0_wp,0.0_wp,0.0_wp]
   real(wp) :: field(9),q(9),z,e,mean,var,expect,second,area,a,b,p0,p1,slope,offset,cut
   real(wp) :: lam(30),jfield(30),actual(1),one(1), uv(2),uvj(2),expected_rate
   integer :: i,k
   logical :: invalid
   type(mc_grain_t) :: g,legacy,explicit_field,boundary
   type(rng_t) :: rng, boundary_rng
   q=1.0_wp
   field=p/x
   do k=1,2
      cut=1000.0_wp
      if(k==2)cut=5.0_wp
      call grain_setup_from_cabs(g,0.01_wp,1.0_wp,'ad_s2',x,q,lam_c_in=cut,J_lam_in=field)
      area=0.0_wp;expect=0.0_wp;second=0.0_wp
      do i=1,8
         a=x(i);b=min(x(i+1),cut)
         if(b<=a)cycle
         slope=(p(i+1)-p(i))/(x(i+1)-x(i));offset=p(i)-slope*a
         p0=p(i);p1=p(i)+slope*(b-a)
         area=area+0.5_wp*(p0+p1)*(b-a)
         expect=expect+offset*log(b/a)+slope*(b-a)
         second=second+offset*(1/a-1/b)+slope*log(b/a)
      end do
      expect=hc_erg_um*expect/area;second=hc_erg_um**2*second/area
      call rng_init(rng,137)
      mean=0.0_wp;var=0.0_wp
      do i=1,n
         z=sample_photon_lam(g,rng)
         if(.not. ((z>0.2_wp.and.z<0.8_wp).or.(z>3.2_wp.and.z<12.8_wp))) &
            error stop 'sample outside positive field support'
         if(z>cut)error stop 'sample outside stochastic cutoff'
         e=hc_erg_um/z;mean=mean+e;var=var+e**2
      end do
      mean=mean/n;var=max(0.0_wp,second-expect**2)/n
      write(*,'(a,i0,a,3es22.13)')'support case ',k,' mean/expected/standard_error ',mean,expect,sqrt(var)
      if(abs(mean-expect)>6*sqrt(var))error stop 'photon mean energy differs from analytic expectation'
      deallocate(g%lam_grid,g%Q_grid,g%u_lam_grid,g%cdf_F,g%cdf_lam,g%cdf_density, &
                 g%T_grid,g%U_grid,g%QPL_grid)
   end do
   ! Force an exact interior-plateau probability, which ordinary random
   ! draws almost never encounter. Check both support and IEEE validity.
   call rng_init(rng,137)
   boundary_rng=rng
   z=rng_uniform(rng)
   boundary%NCDF=7
   allocate(boundary%cdf_F(7),boundary%cdf_lam(7),boundary%cdf_density(7))
   boundary%cdf_F=[0.0_wp,z,z,z,0.5_wp*(1.0_wp+z),1.0_wp,1.0_wp]
   boundary%cdf_lam=[1.0_wp,2.0_wp,3.0_wp,4.0_wp,5.0_wp,6.0_wp,7.0_wp]
   boundary%cdf_density=[2*z,0.0_wp,0.0_wp,0.0_wp,1-z,0.0_wp,0.0_wp]
   call ieee_set_flag(ieee_invalid,.false.)
   e=sample_photon_lam(boundary,boundary_rng)
   call ieee_get_flag(ieee_invalid,invalid)
   if(invalid.or..not.ieee_is_finite(e))error stop 'invalid exact-plateau inverse'
   if(e<=4.0_wp.or.e>=5.0_wp)error stop 'exact-plateau sample outside positive support'
   write(*,'(a)')'Exact plateau probability: PASS'
   ! Both setup routes must respect the actual hard UV support of Mathis.
   call grain_setup_graphite(g,0.001_wp,1.0_wp,'gra_dl01')
   call rng_init(rng,982)
   do i=1,n
      one(1)=sample_photon_lam(g,rng)
      call J_Mathis(1.0_wp,one,actual)
      if(actual(1)<=0.0_wp.or.one(1)<0.0912_wp)error stop 'graphite photon in zero UV field'
   end do
   uv=[0.0912_wp,0.1_wp]
   call J_Mathis(1.0_wp,uv,uvj)
   expected_rate=0.5_wp*(uv(2)-uv(1))*sum(uv*uvj) &
      *4.0_wp*pi**2*(0.01_wp*1.0e-4_wp)**2*1.0e-3_wp/hc_erg_um
   call grain_setup_from_cabs(legacy,0.01_wp,1.0_wp,'ad_s2',[0.08_wp,0.2_wp], &
      [1.0_wp,1.0_wp],lam_c_in=0.1_wp)
   if(abs(legacy%rate_event/expected_rate-1.0_wp)>1.0e-13_wp) &
      error stop 'Lyman-edge cutoff photon rate is inconsistent'
   write(*,'(a,es22.13)') 'Lyman-edge cutoff rate/reference ',legacy%rate_event/expected_rate
   do i=1,30
      lam(i)=0.0912_wp*(1000.0_wp/0.0912_wp)**(real(i-1,wp)/29)
   end do
   call J_Mathis(1.0_wp,lam,jfield)
   call grain_setup_from_cabs(legacy,0.01_wp,1.0_wp,'ad_s2',lam,lam*0+1)
   call grain_setup_from_cabs(explicit_field,0.01_wp,1.0_wp,'ad_s2',lam,lam*0+1,J_lam_in=jfield)
   if(any(legacy%cdf_lam/=explicit_field%cdf_lam).or.any(legacy%cdf_F/=explicit_field%cdf_F)) &
      error stop 'explicit and default fields disagree above Lyman edge'
   write(*,'(a)')'Photon support and analytic mean energy: PASS'
end program test_photon_support
