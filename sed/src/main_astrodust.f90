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
   !   logU=X    set U_mathis = 10^X; tags output filenames (logUX_).
   !   qm        use the energy-space transition-matrix stochastic solver
   !             (thermal-discrete) instead of the default GD narrowing
   !             solver; tags filenames (qm_).
   !   qm_dbcon  energy-space solver in thermal-continuous cooling mode;
   !             tags filenames (qm_dbcon_).

   use constants,         only: wp
   use radfield,          only: J_Mathis, use_mathis_corrected
   use sed_astrodust_mod, only: sed_init, sed_solve, sed_solve_pah, &
                                sed_solve_qm_batch, stoch_method, NLAM, lam, &
                                use_induced_emission, gd_photon_cutoff
   use stoch_qm_mod,      only: qm_method, qm_nstate_default, qm_nisrf_max
   use enthalpy_astrodust_mod, only: s1_density_corrected
   implicit none

   character(len=*), parameter :: F_QTAB =  &
      '../tmatrix/output/q_astrodust_P0.20_Fe0.00_1.400.dat'
   character(len=*), parameter :: F_SIZE = '../data/release/size_distribution.dat'
   character(len=8), parameter :: STAGES(2) = ['S1      ', 'S2      ']
   character(len=*), parameter :: OUTDIR = 'output/astrodust_irem_ours_'
   real(wp),         parameter :: T_LO  = 2.7_wp
   real(wp),         parameter :: T_HI  = 5.0e3_wp
   integer,          parameter :: NT_IN = 200

   real(wp)              :: U_MATHIS = 1.585_wp   ! log U = 0.20
   real(wp)              :: logU_val
   real(wp), allocatable :: J_lam(:), lamI_lam(:)
   integer               :: is, narg, iarg
   character(len=64)     :: arg, suffix, logutag

   ! Parse optional CLI arguments (any position). The solver toggle (qm /
   ! qm_dbcon) and the logU= override each contribute a filename tag so a
   ! non-default run does not clobber the production output.
   suffix  = ''
   logutag = ''
   narg = command_argument_count()
   do iarg = 1, narg
      call get_command_argument(iarg, arg)
      if (trim(arg) == 'qm') then
         stoch_method = 'qm'
         suffix = 'qm_'
      else if (trim(arg) == 'qm_dbcon') then
         stoch_method = 'qm'
         qm_method    = 'dbcon'       ! thermal-continuous (GD-collapse) QM matrix
         suffix = 'qm_dbcon_'
      else if (trim(arg) == 'draine') then
         stoch_method = 'draine'      ! Draine's original GD (default is 'heuristic')
         suffix = 'draine_'
      else if (index(arg, 'logU=') > 0) then
         read(arg(index(arg,'=')+1:), *) logU_val
         U_MATHIS = 10.0_wp ** logU_val
         write(logutag, '(a,f0.2,a)') 'logU', logU_val, '_'
      else if (index(arg, 'nstate=') > 0) then
         read(arg(index(arg,'=')+1:), *) qm_nstate_default
         block
            character(len=12) :: nstag
            write(nstag, '(a,i0,a)') 'ns', qm_nstate_default, '_'
            suffix = trim(suffix)//trim(nstag)
         end block
      else if (index(arg, 'nisrf=') > 0) then
         read(arg(index(arg,'=')+1:), *) qm_nisrf_max
      else if (trim(arg) == 'induced') then
         use_induced_emission = .true.
         suffix = trim(suffix)//'ind_'
      else if (trim(arg) == 'photcut') then
         gd_photon_cutoff = .true.    ! enable dbdis photon cutoff in GD emission
         suffix = trim(suffix)//'photcut_'
      else if (trim(arg) == 'c2') then
         s1_density_corrected = .true.   ! Stage-1 density-corrected prefactor
         suffix = trim(suffix)//'c2_'
      else if (trim(arg) == 'mathis_orig') then
         use_mathis_corrected = .false.
         suffix = trim(suffix)//'morig_'
      end if
   end do
   suffix = trim(suffix)//trim(logutag)

   write(*,'(a)') '=========================================================='
   write(*,'(a)') ' main_astrodust: production driver'
   write(*,'(a)') '=========================================================='
   write(*,'(a,a)')    ' Q table     : ', F_QTAB
   write(*,'(a,a)')    ' size_dist   : ', F_SIZE
   write(*,'(a,f8.3)') ' U_mathis    : ', U_MATHIS
   write(*,'(a,i0)')   ' NT (T grid) : ', NT_IN
   write(*,'(a,a)')    ' stoch_method: ', trim(stoch_method)
   if (trim(stoch_method) == 'qm') &
      write(*,'(a,a)') ' qm_method   : ', trim(qm_method)

   call sed_init(F_QTAB, F_SIZE, NT_IN, T_LO, T_HI)
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
      do kk = 1, NLAM
         write(uu,'(es14.6,1x,es16.8)') lam(kk), lamI(kk)
      end do
      close(uu)
      write(*,'(a,a)') 'wrote ', trim(fn)
   end subroutine write_sed

end program main_astrodust
