program test_mc_supplied_field
   use constants, only: wp
   use radfield, only: J_Mathis, J_HD23, use_mathis_corrected
   use mc_engine, only: grain_setup_from_cabs
   use mc_grain_type, only: mc_grain_t
   implicit none
   real(wp), parameter :: pi = 3.141592653589793238462643383279502884197_wp
   real(wp), parameter :: c = 2.99792458e8_wp
   real(wp), parameter :: lam(6) = [0.1_wp, 0.2_wp, 1.0_wp, 10.0_wp, 100.0_wp, 1000.0_wp]
   real(wp), parameter :: q(6) = 1.0_wp
   real(wp) :: field(6), legacy(6), factor
   type(mc_grain_t) :: supplied_a, supplied_b, default_grain, explicit_grain

   factor = (4.0_wp*pi/c)*1.0e-5_wp
   call J_HD23(1.6_wp, lam, field)
   call grain_setup_from_cabs(supplied_a, 0.01_wp, 1.0_wp, 'ad_s2', lam, q, J_lam_in=field)
   call grain_setup_from_cabs(supplied_b, 0.01_wp, 12.0_wp, 'ad_s2', lam, q, J_lam_in=field)
   if (any(abs(supplied_a%u_lam_grid-field*factor) > 1.0e-13_wp*maxval(field*factor))) &
      error stop 'supplied field was not used'
   if (any(supplied_a%u_lam_grid /= supplied_b%u_lam_grid)) &
      error stop 'supplied field depends on U_isrf'

   use_mathis_corrected = .true.
   call J_Mathis(1.6_wp, lam, legacy)
   call grain_setup_from_cabs(default_grain, 0.01_wp, 1.6_wp, 'ad_s2', lam, q)
   call grain_setup_from_cabs(explicit_grain, 0.01_wp, 1.6_wp, 'ad_s2', lam, q, J_lam_in=legacy)
   if (any(default_grain%u_lam_grid /= explicit_grain%u_lam_grid)) &
      error stop 'default Mathis differs from explicit Mathis'
   write(*,'(a)') 'MC supplied field and legacy Mathis paths passed'
end program test_mc_supplied_field
