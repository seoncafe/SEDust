program test_radiation_support
   use constants, only: wp
   use radfield, only: hardest_photon_energy, radiation_field_support, &
      isolated_radiation_support_count, RADIATION_SUPPORT_ISOLATED, &
      RADIATION_SUPPORT_EMPTY, RADIATION_SUPPORT_OK
   use p_sub, only: radiation_absorption, p_sub_setup, calc_P, calc_Teq
   implicit none
   real(wp), parameter :: hc=6.62606957e-27_wp*2.99792458e10_wp*1e4_wp
   real(wp) :: lam(5), cross(5), field(5), t(2), cooling(2), h(2)
   real(wp) :: p(2), lnp(2), reference(2), reference_log(2), teq, expected, integral
   integer :: first, last, support_status
   integer(kind(isolated_radiation_support_count)) :: warnings_before
   lam=[0.5_wp,1.0_wp,2.0_wp,4.0_wp,8.0_wp]
   cross=1.0_wp
   field=[0.0_wp,1.0_wp,1.0_wp,0.0_wp,0.0_wp]
   expected=1.5_wp*log(2.0_wp)
   integral=radiation_absorption(lam,cross,field)
   if(abs(integral/expected-1.0_wp)>1e-14_wp) error stop 'absorption support'
   call p_sub_setup(lam)
   if(radiation_absorption(lam,cross,field)/=integral) error stop 'cache changed integral'
   ! A different caller grid must not use the cached wavelengths.
   if(abs(radiation_absorption(2*lam,cross,field)/(2*expected)-1.0_wp)>1e-14_wp) &
      error stop 'different wavelength grid'
   t=[1.0_wp,100.0_wp]
   cooling=[0.1_wp,10.0_wp]
   call calc_Teq(lam,cross,field,t,cooling,teq)
   if(abs(teq/(10*expected)-1.0_wp)>1e-14_wp) error stop 'equilibrium support'
   h=[hc,hc+hc/1.5_wp]
   call calc_P(lam,cross,field,t,cooling,h,p,lnp)
   call calc_P(lam(2:3),cross(2:3),field(2:3),t,cooling,h,reference,reference_log)
   if(any(abs(p-reference)>1e-14_wp)) error stop 'exterior zero extension'
   ! With two bins every absorbed photon raises the grain into the top bin
   ! (Draine & Li 2001 eqs. 16 and 28), so 4 pi P(2) cooling(2) = P(1) times
   ! the absorbed power (hence the check on P(2)/P(1)): only the support
   ! [1,2] um counts, whatever the gap.
   h=[hc,hc+hc/10.0_wp]
   call calc_P(lam,cross,field,t,cooling,h,p,lnp)
   expected=1.5_wp*log(2.0_wp)
   integral=p(2)/p(1)*cooling(2)
   if(abs(integral/expected-1.0_wp)>1e-14_wp) error stop 'top-bin support'
   warnings_before=isolated_radiation_support_count
   field=0.0_wp
   call radiation_field_support(lam,field,first,last,support_status)
   if(support_status/=RADIATION_SUPPORT_EMPTY) error stop 'dark support status'
   if(hardest_photon_energy(lam,field)/=0.0_wp) error stop 'dark photon bound'
   call calc_P(lam,cross,field,t,cooling,h,p,lnp)
   if(p(1)/=1.0_wp.or.p(2)/=0.0_wp) error stop 'dark-field probability'
   if(radiation_absorption(lam,cross,field)/=0.0_wp) error stop 'dark absorption'
   if(isolated_radiation_support_count/=warnings_before) error stop 'dark field warning'
   field(3)=1.0_wp
   call radiation_field_support(lam,field,first,last,support_status)
   if(support_status/=RADIATION_SUPPORT_ISOLATED.or.first/=0.or.last/=0) &
      error stop 'isolated support status'
   if(isolated_radiation_support_count/=warnings_before) error stop 'pure query changed counter'
   if(hardest_photon_energy(lam,field)/=0.0_wp) error stop 'isolated photon bound' 
   call calc_P(lam,cross,field,t,cooling,h,p,lnp)
   if(p(1)/=1.0_wp.or.p(2)/=0.0_wp) error stop 'isolated node bandwidth'
   if(radiation_absorption(lam,cross,field)/=0.0_wp) error stop 'isolated absorption'
   if(isolated_radiation_support_count/=warnings_before+2) error stop 'shared warning count'
   ! The photon bound retains its original support for unordered arrays.
   call radiation_field_support([2.0_wp,1.0_wp],[1.0_wp,1.0_wp],first,last,support_status)
   if(support_status/=RADIATION_SUPPORT_OK.or.first/=2.or.last/=1) &
      error stop 'unordered support endpoints'
   if(hardest_photon_energy([2.0_wp,1.0_wp],[1.0_wp,1.0_wp])/=hc) &
      error stop 'unordered photon bound'
   if(isolated_radiation_support_count/=warnings_before+2) error stop 'ordinary support warning'
   write(*,'(a)') 'PASS: absorption support, grid cache, equilibrium, transitions, dark field'
end program test_radiation_support
