program test_qm_rate_closure
   ! Exercise the public QM solver on the Zubko U=1 hard-field PAH optics.
   ! Pass the absolute path to data/zubko, a PAH size index, the cooling
   ! method, and the transition-grid count (150 or 0 for the full grid).
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use constants, only: wp
   use sed_astrodust_mod, only: dust_model_t, build_zubko
   use sed_apply_options, only: add_hard_euv_component
   use sed_mathlib, only: interp
   use p_sub, only: calc_Teq
   use radfield, only: J_Mathis
   use stoch_qm_mod, only: qm_solve_grain, qm_method, qm_nstate_default, &
                           qm_nisrf_max, qm_verbose, QM_SOLVED, QM_INVALID_INPUT
   implicit none
   type(dust_model_t) :: model
   real(wp), allocatable :: field(:), emission(:), invalid_cabs(:)
   real(wp) :: teq, eeq, absorbed, emitted, dlnl, wavelength
   integer :: status, grain_index, grid_count, raw_negative_count, i
   logical :: solved
   character(len=512) :: data_dir, arg
   character(len=5) :: method, actual_cooling

   if (command_argument_count() /= 4) error stop &
      'usage: test_qm_rate_closure.x DATA_ZUBKO SIZE_INDEX METHOD GRID_COUNT'
   call get_command_argument(1, data_dir)
   call get_command_argument(2, arg)
   read(arg, *) grain_index
   call get_command_argument(3, arg)
   method = trim(arg)
   call get_command_argument(4, arg)
   read(arg, *) grid_count
   if (grain_index < 1 .or. grid_count < 0) error stop 'invalid test arguments'
   select case (method)
   case ('dbcon', 'dbdis', 'stati')
   case default
      error stop 'unsupported cooling method'
   end select

   call build_zubko(model, trim(data_dir)//'/ZDA_BARE_GR_S_Config.dat', &
                    trim(data_dir)//'/', 200, 2.7_wp, 5000.0_wp, status, &
                    include_euv=.true., &
                    qtable_path=trim(data_dir)//'/sedust_zubko.h5')
   if (status /= 0) error stop 'Zubko model construction failed'
   if (grain_index > size(model%pops(1)%dn)) error stop 'missing PAH size'
   if (model%pops(1)%dn(grain_index) <= 0.0_wp) error stop 'unpopulated PAH size'

   allocate(field(model%NLAM), emission(model%NLAM))
   call J_Mathis(1.0_wp, model%lam, field)
   call add_hard_euv_component(model%lam, field)
   call calc_Teq(model%lam, model%pops(1)%Cabs(:, grain_index), field, &
                 model%T_first, model%pops(1)%kappB(:, grain_index), teq)
   call interp(model%T_first, model%pops(1)%H(:, grain_index), teq, eeq)

   qm_method = method
   qm_nstate_default = 200
   qm_nisrf_max = grid_count
   qm_verbose = .false.
   call qm_solve_grain(model%NLAM, model%lam, &
                       model%pops(1)%Cabs(:, grain_index), field, &
                       eeq, 150.0_wp * 1.602176634e-12_wp, &
                       model%pops(1)%aeff(grain_index) * 1.0e-4_wp, 'pah', &
                       emission, solved, status, actual_cooling, &
                       raw_negative_count)
   if (.not. solved .or. status /= QM_SOLVED) error stop 'QM grain not solved'
   if (raw_negative_count /= 0) error stop 'negative raw transition rate'
   if (.not. all(ieee_is_finite(emission))) error stop 'nonfinite QM emission'
   if (any(emission < 0.0_wp)) error stop 'negative QM emission'
   if (method /= 'stati' .and. actual_cooling /= method) &
      error stop 'unexpected cooling method'
   if (method == 'stati' .and. actual_cooling /= 'stati' .and. &
       actual_cooling /= 'dbdis') error stop 'unexpected actual cooling method'

   absorbed = 0.0_wp
   emitted = 0.0_wp
   do i = 1, model%NLAM
      if (i == 1) then
         dlnl = 0.5_wp * log(model%lam(2) / model%lam(1))
      else if (i == model%NLAM) then
         dlnl = 0.5_wp * log(model%lam(i) / model%lam(i-1))
      else
         dlnl = 0.5_wp * log(model%lam(i+1) / model%lam(i-1))
      end if
      wavelength = model%lam(i) * 1.0e-4_wp
      absorbed = absorbed + model%pops(1)%Cabs(i, grain_index) * &
                 (4.0_wp * acos(-1.0_wp) * 10.0_wp * field(i)) * &
                 wavelength * dlnl
      emitted = emitted + emission(i) * dlnl
   end do
   if (absorbed <= 0.0_wp) error stop 'zero absorbed power'
   write(*,'(a,i0,1x,a,a,1x,a,i0,1x,a,a,1x,a,es16.8,1x,a,i0)') &
      'size=', grain_index, 'method=', trim(method), 'grid=', grid_count, &
      'cooling=', trim(actual_cooling), 'Pem/Pabs=', emitted / absorbed, &
      'negative_rates=', raw_negative_count
   if (abs(emitted / absorbed - 1.0_wp) > 0.01_wp) &
      error stop 'QM emission differs from absorbed power by more than 1%'

   allocate(invalid_cabs(model%NLAM))
   invalid_cabs = model%pops(1)%Cabs(:, grain_index)
   invalid_cabs(1) = -invalid_cabs(1)
   raw_negative_count = -1
   emission = -1.0_wp
   call qm_solve_grain(model%NLAM, model%lam, invalid_cabs, field, &
                       eeq, 150.0_wp * 1.602176634e-12_wp, &
                       model%pops(1)%aeff(grain_index) * 1.0e-4_wp, 'pah', &
                       emission, solved, status, actual_cooling, &
                       raw_negative_count)
   if (solved .or. status /= QM_INVALID_INPUT .or. raw_negative_count /= 0 &
       .or. any(emission /= 0.0_wp)) &
      error stop 'negative absorption accepted as a physical input'
end program test_qm_rate_closure
