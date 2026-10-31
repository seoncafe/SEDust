program test_qm_invalid_bins
   use constants, only: wp
   use stoch_qm_mod, only: qm_solve_grain, qm_nstate_default, qm_nisrf_max, qm_verbose
   implicit none
   real(wp), parameter :: lam(3) = [0.1_wp, 0.2_wp, 0.3_wp]
   real(wp), parameter :: cabs(3) = 1.0e-21_wp
   real(wp), parameter :: field(3) = 1.0e-12_wp
   real(wp) :: emission(3), short_emission(1)
   logical :: solved

   qm_verbose = .false.
   qm_nstate_default = 4
   qm_nisrf_max = 200
   emission = -1.0_wp
   call qm_solve_grain(3,lam,cabs,field,1.0e-15_wp, &
        2.4e-10_wp,1.0e-7_wp,'pah',emission,solved)
   if (solved .or. any(emission /= 0.0_wp)) error stop 'nstate<5 accepted'

   qm_nstate_default = 200
   qm_nisrf_max = 1
   emission = -1.0_wp
   call qm_solve_grain(3,lam,cabs,field,1.0e-15_wp, &
        2.4e-10_wp,1.0e-7_wp,'pah',emission,solved)
   if (solved .or. any(emission /= 0.0_wp)) error stop 'nisrf<2 accepted'

   qm_nisrf_max = 200
   short_emission = -1.0_wp
   call qm_solve_grain(1,lam(:1),cabs(:1),field(:1), &
        1.0e-15_wp,2.4e-10_wp,1.0e-7_wp,'pah',short_emission,solved)
   if (solved .or. any(short_emission /= 0.0_wp)) error stop 'nlam<2 accepted'
   write(*,'(a)') 'invalid QM bin and wavelength counts rejected'
end program test_qm_invalid_bins
