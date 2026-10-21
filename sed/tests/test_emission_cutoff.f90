program test_emission_cutoff
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
   use constants, only: wp
   use dust_lib, only: dust_model_t, dust_emission, emission_bin_power_cutoff, &
      radiation_absorption, isolated_radiation_support_count
   implicit none
   type(dust_model_t) :: m
   real(wp) :: j(1), out(1)
   integer :: status
   integer(kind(isolated_radiation_support_count)) :: warnings_before
   real(wp) :: absorption
   j=0.0_wp
   emission_bin_power_cutoff=-1.0_wp
   call dust_emission(m,j,out,status=status)
   if (status/=4) stop 1
   emission_bin_power_cutoff=2.0_wp
   call dust_emission(m,j,out,status=status)
   if (status/=4) stop 2
   emission_bin_power_cutoff=ieee_value(0.0_wp,ieee_quiet_nan)
   call dust_emission(m,j,out,status=status)
   if (status/=4) stop 3
   emission_bin_power_cutoff=0.0_wp
   call dust_emission(m,j,out,status=status)
   if (status/=3) stop 4
   warnings_before=isolated_radiation_support_count
   absorption=radiation_absorption([1.0_wp,2.0_wp,3.0_wp], &
      [1.0_wp,1.0_wp,1.0_wp],[0.0_wp,1.0_wp,0.0_wp])
   if(absorption/=0.0_wp) stop 5
   if(isolated_radiation_support_count/=warnings_before+1) stop 6
   print *, 'CUTOFF_STATUS_OK' 
end program test_emission_cutoff
