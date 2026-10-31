program test_cmb_population_invariance
   ! The population solver receives the complete radiation field. Changing
   ! the stored partial CMB integral must not change its emitted spectrum.
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use constants, only: wp
   use radfield, only: J_Mathis
   use dust_lib, only: dust_model_t, build_astrodust, dust_emission, sed_set_data_root
   use omp_lib, only: omp_set_num_threads
   implicit none
   type(dust_model_t) :: model
   character(len=512) :: data_dir, qtable_path
   character(len=9), parameter :: methods(2) = [character(len=9) :: 'heuristic', 'draine']
   real(wp), allocatable :: field(:), reference(:), changed(:), background(:,:)
   integer :: status, imethod, ipop

   call get_command_argument(1, data_dir)
   if (len_trim(data_dir) == 0) error stop 'data directory argument required'
   call sed_set_data_root(trim(data_dir))
   qtable_path = trim(data_dir)//'/astrodust/sedust_astrodust.h5'
   call omp_set_num_threads(1)
   call build_astrodust(model, trim(qtable_path), 120, 2.725_wp, 5000.0_wp, &
                        status=status, include_euv=.false., &
                        load_polarized_optics=.false.)
   if (status /= 0) error stop 'build_astrodust failed'
   allocate(field(model%NLAM), reference(model%NLAM), changed(model%NLAM))
   allocate(background(model%NA, size(model%pops)))
   do ipop = 1, size(model%pops)
      background(:, ipop) = model%pops(ipop)%kappCMB
   end do
   call J_Mathis(1.0_wp, model%lam, field)

   do imethod = 1, size(methods)
      model%stoch_method = methods(imethod)
      call dust_emission(model, field, reference, status=status)
      if (status /= 0) error stop 'baseline dust_emission failed'
      if (.not. all(ieee_is_finite(reference)) .or. .not. any(reference > 0.0_wp)) &
         error stop 'invalid baseline emission'

      ! A negative sentinel keeps the old cooling rates finite while making
      ! any accidental use of this stored value plainly visible.
      do ipop = 1, size(model%pops)
         model%pops(ipop)%kappCMB = -model%pops(ipop)%kappB(model%NT/2, :)
      end do
      call dust_emission(model, field, changed, status=status)
      if (status /= 0) error stop 'perturbed dust_emission failed'
      if (.not. all(ieee_is_finite(changed))) error stop 'nonfinite emission'
      if (any(changed /= reference)) error stop 'emission depends on stored CMB integral'
      write(*,'(a,a)') trim(methods(imethod)), ': unchanged after CMB-integral perturbation'
      do ipop = 1, size(model%pops)
         model%pops(ipop)%kappCMB = background(:, ipop)
      end do
   end do
end program test_cmb_population_invariance
