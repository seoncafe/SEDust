program main_astrodust
   ! Production driver.
   !
   ! Computes the astrodust thermal-emission SED for a single Mathis-ISRF
   ! cell for the two enthalpy stages (S1, S2) and the PAH population. Writes:
   !   output/astrodust_irem_ours_S1.dat
   !   output/astrodust_irem_ours_S2.dat
   !   output/astrodust_irem_ours_PAH.dat
   ! Each: header + (lambda[um], lambda*I_lambda / N_H) per line.
   !
   ! Heating: U_mathis = 1.585 (log U = 0.20, HD23 best fit) by default.
   !
   ! Optional CLI arguments (any position):
   !   euv       solve on the _euv Q table, whose grid runs to 1.0e-4 um
   !             (12398 eV) instead of stopping at the Lyman limit; tags
   !             filenames (euv_). The Mathis field is zero below that limit,
   !             so the added wavelengths carry no photon and change the
   !             emergent SED only through the quadrature.
   !   logU=X    set U_mathis = 10^X; tags output filenames (logUX_).
   !   qm        use the energy-space transition-matrix stochastic solver
   !             (thermal-discrete) instead of the default GD narrowing
   !             solver; tags filenames (qm_).
   !   qm_dbcon  energy-space solver in thermal-continuous cooling mode;
   !             tags filenames (qm_dbcon_).
   !   induced   multiply the emergent SED by the stimulated-emission factor
   !             (1 + J_lam/B_env); diagnostic only, must stay off for
   !             production. Tags filenames (induced_).

   use constants,         only: wp
   use radfield,          only: J_Mathis, use_mathis_corrected
   use sed_astrodust_mod, only: sed_init, sed_solve, sed_solve_pah, &
                                sed_solve_qm_batch, stoch_method, NLAM, lam, &
                                use_induced_emission, gd_photon_cutoff
   use stoch_qm_mod,      only: qm_method, qm_nstate_default, qm_nisrf_max
   use enthalpy_astrodust_mod, only: s1_density_corrected
   ! The astrodust EUV band, when one is asked for, is solved on the oblate
   ! spheroid of the Q table rather than on a volume-equivalent sphere. That
   ! calculation is a separate module so that the library links without the
   ! T-matrix; naming it here makes this driver able to carry a grid below
   ! 0.0912 um the moment sed_init is given a lam_min.
   use euv_astrodust_tmatrix, only: use_tmatrix_euv_band_optics
   implicit none

   ! The two shipped scalar Q tables.  The model's wavelength grid is whichever
   ! one this driver is given, so the `euv` argument decides whether the SED is
   ! solved over the ionizing band as well.  Default is the non-ionizing table,
   ! matching the field this driver illuminates with: J_Mathis is exactly zero
   ! below the Lyman limit, so the 633 extra wavelengths carry no photon and
   ! only cost solver time.  `euv` is there for a run that wants the band
   ! resolved anyway, and tags its output files with it.
   character(len=*), parameter :: F_QTAB =  &
      '../data/astrodust/q_astrodust_P0.20_Fe0.00_1.400.dat'
   character(len=*), parameter :: F_QTAB_EUV =  &
      '../data/astrodust/q_astrodust_P0.20_Fe0.00_1.400_euv.dat'
   character(len=*), parameter :: F_SIZE = '../data/release/size_distribution.dat'
   character(len=8), parameter :: STAGES(2) = ['S1      ', 'S2      ']
   character(len=*), parameter :: OUTDIR = 'output/astrodust_irem_ours_'
   real(wp),         parameter :: T_LO  = 2.7_wp
   real(wp),         parameter :: T_HI  = 5.0e3_wp
   integer,          parameter :: NT_IN = 200

   real(wp)              :: U_MATHIS = 1.585_wp   ! log U = 0.20
   real(wp)              :: logU_val
   real(wp), allocatable :: J_lam(:), lamI_lam(:)
   integer               :: is, narg, iarg, n_solver
   character(len=64)     :: arg, suffix, solver_tag
   character(len=24)     :: nstag, nitag, logutag
   logical               :: set_logu, set_nstate, set_nisrf
   ! Solve over the ionizing band as well, on the _euv Q table.
   logical               :: use_euv_grid = .false.

   ! Parse the optional CLI arguments, then build the filename tag.
   !
   ! The two steps are kept apart on purpose. Parsing only records WHAT was
   ! asked for; the tag is assembled afterwards in ONE fixed order, so that
   ! options given in any order produce the same filename -- which is what
   ! lets a non-default run coexist with the production one instead of
   ! overwriting it. (Tagging as the arguments are read cannot do that: it
   ! either lets a later option erase an earlier one's tag, or gives the same
   ! combination two different names depending on how it was typed.)
   !
   ! The order runs from the coarsest distinction to the finest -- solver,
   ! then enthalpy stage, then radiation field, then the emission term, then
   ! numerical settings -- so that a directory listing groups runs by solver.
   solver_tag = '';  logutag = '';  nstag = '';  nitag = ''
   n_solver   = 0
   set_logu   = .false.;  set_nstate = .false.;  set_nisrf = .false.

   narg = command_argument_count()
   do iarg = 1, narg
      call get_command_argument(iarg, arg)
      if (trim(arg) == 'qm') then
         stoch_method = 'qm'
         solver_tag = 'qm_';        n_solver = n_solver + 1
      else if (trim(arg) == 'qm_dbcon') then
         stoch_method = 'qm'
         qm_method    = 'dbcon'       ! thermal-continuous (GD-collapse) QM matrix
         solver_tag = 'qm_dbcon_';   n_solver = n_solver + 1
      else if (trim(arg) == 'draine') then
         stoch_method = 'draine'      ! Draine's original GD (default is 'heuristic')
         solver_tag = 'draine_';     n_solver = n_solver + 1
      else if (index(arg, 'logU=') > 0) then
         read(arg(index(arg,'=')+1:), *) logU_val
         U_MATHIS = 10.0_wp ** logU_val
         set_logu = .true.
      else if (index(arg, 'nstate=') > 0) then
         read(arg(index(arg,'=')+1:), *) qm_nstate_default
         set_nstate = .true.
      else if (index(arg, 'nisrf=') > 0) then
         read(arg(index(arg,'=')+1:), *) qm_nisrf_max
         set_nisrf = .true.
      else if (trim(arg) == 'induced') then
         use_induced_emission = .true.
      else if (trim(arg) == 'photcut') then
         gd_photon_cutoff = .true.    ! enable dbdis photon cutoff in GD emission
      else if (trim(arg) == 'c2') then
         s1_density_corrected = .true.   ! Stage-1 density-corrected prefactor
      else if (trim(arg) == 'mathis_orig') then
         use_mathis_corrected = .false.
      else if (trim(arg) == 'euv') then
         use_euv_grid = .true.        ! solve on the _euv Q table's 1762 grid
      end if
   end do

   ! qm, qm_dbcon and draine each select the solver, so naming two of them
   ! is a contradiction. Silently keeping the last one would write a file
   ! whose name claims a solver the run did not use.
   if (n_solver > 1) then
      write(*,'(a)') ' main_astrodust: give at most one of qm / qm_dbcon / draine'
      stop 1
   end if

   if (set_logu)   write(logutag, '(a,f0.2,a)') 'logU', logU_val, '_'
   if (set_nstate) write(nstag,   '(a,i0,a)')   'ns',   qm_nstate_default, '_'
   if (set_nisrf)  write(nitag,   '(a,i0,a)')   'nisrf', qm_nisrf_max, '_'

   ! The wavelength grid comes first: it is the coarsest distinction of all,
   ! since it decides how many points every other setting is evaluated on.
   suffix = ''
   if (use_euv_grid)            suffix = 'euv_'
   suffix = trim(suffix)//trim(solver_tag)
   if (s1_density_corrected)    suffix = trim(suffix)//'c2_'
   if (.not. use_mathis_corrected) suffix = trim(suffix)//'morig_'
   suffix = trim(suffix)//trim(logutag)
   if (use_induced_emission)    suffix = trim(suffix)//'induced_'
   suffix = trim(suffix)//trim(nstag)//trim(nitag)
   if (gd_photon_cutoff)        suffix = trim(suffix)//'photcut_'

   write(*,'(a)') '=========================================================='
   write(*,'(a)') ' main_astrodust: production driver'
   write(*,'(a)') '=========================================================='
   if (use_euv_grid) then
      write(*,'(a,a)') ' Q table     : ', F_QTAB_EUV
   else
      write(*,'(a,a)') ' Q table     : ', F_QTAB
   end if
   write(*,'(a,a)')    ' size_dist   : ', F_SIZE
   write(*,'(a,f8.3)') ' U_mathis    : ', U_MATHIS
   write(*,'(a,i0)')   ' NT (T grid) : ', NT_IN
   write(*,'(a,a)')    ' stoch_method: ', trim(stoch_method)
   if (trim(stoch_method) == 'qm') &
      write(*,'(a,a)') ' qm_method   : ', trim(qm_method)
   ! Say where the run will land before it starts, so a mistyped option is
   ! caught now rather than after the solve.
   write(*,'(a,a,a)')  ' output files: ', OUTDIR, trim(suffix)//'<stage>.dat'

   call use_tmatrix_euv_band_optics()
   if (use_euv_grid) then
      call sed_init(F_QTAB_EUV, F_SIZE, NT_IN, T_LO, T_HI)
   else
      call sed_init(F_QTAB, F_SIZE, NT_IN, T_LO, T_HI)
   end if
   write(*,'(a,i0,a)') ' sed_init done. NLAM=', NLAM, ' wavelengths cached.'
   write(*,'(a)') ''

   allocate(J_lam(NLAM), lamI_lam(NLAM))
   call J_Mathis(U_MATHIS, lam, J_lam)

   if (trim(stoch_method) == 'qm') then
      ! QM batch mode: all 4 grain types (S1, S2, PAH-neutral, PAH-cation)
      ! are processed in one parallel region for maximum thread utilisation.
      block
         real(wp), allocatable :: lamI_stages(:,:), lamI_pah_b(:)
         allocate(lamI_stages(NLAM, 2), lamI_pah_b(NLAM))

         write(*,'(a)') ' solving all stages (QM batch) ...'
         call sed_solve_qm_batch(J_lam, lamI_stages, lamI_pah_b)

         do is = 1, 2
            call write_sed(trim(STAGES(is)), lamI_stages(:, is), .false.)
         end do
         call write_sed('PAH', lamI_pah_b, .true.)

         deallocate(lamI_stages, lamI_pah_b)
      end block
   else
      ! Default GD (Draine-narrowing) solver, sequential per grain type.
      do is = 1, 2
         write(*,'(a,a,a)', advance='no') ' solving stage ', trim(STAGES(is)), ' ... '
         call sed_solve(J_lam, trim(STAGES(is)), lamI_lam)
         call write_sed(trim(STAGES(is)), lamI_lam, .false.)
      end do

      write(*,'(a)', advance='no') ' solving stage PAH ... '
      call sed_solve_pah(J_lam, lamI_lam)
      call write_sed('PAH', lamI_lam, .true.)
   end if

   write(*,'(a)') ''
   write(*,'(a)') ' main_astrodust: done.'

contains

   subroutine write_sed(stage, lamI, is_pah)
      ! Write one SED file: output/astrodust_irem_ours_<suffix><stage>.dat
      character(len=*), intent(in) :: stage
      real(wp),         intent(in) :: lamI(:)
      logical,          intent(in) :: is_pah
      integer :: uu, kk
      character(len=96) :: fn
      write(fn,'(a,a,a,a)') OUTDIR, trim(suffix), trim(stage), '.dat'
      open(newunit=uu, file=trim(fn), status='replace', action='write')
      if (is_pah) then
         write(uu,'(a)') '# DH21 PAH SED for Mathis ISRF, U = 1.585'
         write(uu,'(a)') '# DL07 PAH cross sections, neutral + cation mixed by f_ion(a)'
      else
         write(uu,'(a)') '# DH21 astrodust SED for Mathis ISRF, U = 1.585'
         write(uu,'(a,a)') '# Enthalpy stage: ', trim(stage)
      end if
      write(uu,'(a)') '# columns: lambda[um]    lambda*I_lambda / N_H [erg s^-1 cm^-2 sr^-1 H^-1]'
      ! e3 on the intensity: carried into the EUV and X-ray, lambda*I_lambda
      ! underflows to subnormals (~1e-320) on the Wien side of a 20 K grain, and
      ! a two-digit exponent field drops the E there -- Fortran writes
      ! 4.29970510-319, which its own list-directed read accepts but numpy, IDL
      ! and awk do not.  The wavelength column needs no widening: the grid spans
      ! 1e-4 to 4e4 um, two exponent digits throughout.
      do kk = 1, NLAM
         write(uu,'(es14.6,1x,es16.8e3)') lam(kk), lamI(kk)
      end do
      close(uu)
      write(*,'(a,a)') 'wrote ', trim(fn)
   end subroutine write_sed

end program main_astrodust
