program calc_sed
   !====================================================================
   ! Thermal-emission SED of any dust model this library can build, for a
   ! single illuminated cell.
   !
   !   ./calc_sed.x astrodust [settings ...]
   !   ./calc_sed.x dl07 [submodel] [settings ...]
   !   ./calc_sed.x mrn [settings ...]
   !   ./calc_sed.x zubko [settings ...]
   !   ./calc_sed.x themis [settings ...]
   !   ./calc_sed.x g18d [settings ...]
   !
   ! The model name is required; without one, or with one this program does
   ! not know, it prints the usage and stops.  Everything after it is the
   ! shared run-setting vocabulary of sed_run_options -- solver, wavelength
   ! grid, radiation field, emission term, transition-matrix sizes -- plus the
   ! settings only some models have a referent for: the graphite of the PAH
   ! xi blend (astrodust, dl07), the Stage-1 enthalpy prefactor c2 (astrodust),
   ! and the size-distribution submodel (dl07).  A setting the named model has
   ! no referent for is refused by name rather than accepted and ignored.
   !
   ! Output, under output/, with every non-default setting in the tag:
   !   sed_astrodust[_<tag>]_{S1,S2,PAH}.dat
   !     lambda[um]  lambda*I_lambda / N_H            (per enthalpy stage)
   !   sed_dl07_<submodel>[_<tag>].dat
   !     lambda[um]  total  silicate  carbonaceous
   !   sed_mrn[_<tag>].dat
   !     lambda[um]  total  graphite  silicate
   !   sed_zubko[_<tag>].dat
   !     lambda[um]  total  PAH  graphite  silicate
   !   sed_themis[_<tag>].dat, sed_g18d[_<tag>].dat
   !     lambda[um]  total  <one column per DustEM grain population>
   ! all in erg s^-1 cm^-2 sr^-1 H^-1 (the HD23 / DL07 convention).
   !
   ! This is the emission counterpart of calc_kext.x, which does the size
   ! integral for the transport optics of the same models.
   !====================================================================
   use constants,         only: wp
   use radfield,          only: J_Mathis, hardest_photon_energy
   use sed_astrodust_mod, only: sed_init, sed_solve, sed_solve_pah, &
                                sed_solve_qm_batch, stoch_method, NLAM, lam, &
                                sed_init_dl07, sed_solve_dl07, &
                                dust_model_t, build_zubko, build_mrn, build_dustem, &
                                dust_emission
   ! Which graphite optics the carbonaceous xi blend takes; each model sets its
   ! own before the arguments are applied on top.
   use qpah,              only: qpah_graphite_source, nc_coeff, nc_integer
   use grain_dist_mod,    only: gd_apply_d03_reduction
   use sed_run_options,   only: run_options_t, declare_run_options, &
                                widen_run_options, read_run_subject, &
                                read_run_option, check_run_options, &
                                run_options_tag, report_run_options, &
                                write_run_option_usage, T_HARD_EUV, W_HARD_EUV
   use sed_apply_options, only: apply_run_options, add_hard_euv_component
   ! The astrodust EUV band, when one is asked for, is solved on the oblate
   ! spheroid of the Q table rather than on a volume-equivalent sphere.  That
   ! calculation is a separate module so that the library links without the
   ! T-matrix; naming it here makes this driver able to carry a grid below
   ! 0.0912 um the moment sed_init is given a lam_min.
   use euv_astrodust_tmatrix, only: use_tmatrix_euv_band_optics
   implicit none

   ! ---- model inputs, the same products calc_kext.x builds from ----------
   ! Each model's own product carries ONE wavelength axis and the index where
   ! the non-ionizing part of it begins, so `euv` selects a view of the same
   ! file rather than a second one.
   character(len=*), parameter :: F_QT_AD = '../data/astrodust/sedust_astrodust.h5'
   character(len=*), parameter :: F_QT_DL = '../data/dl07/sedust_dl07.h5'
   character(len=*), parameter :: F_QT_MRN = '../data/mrn/sedust_mrn.h5'
   character(len=*), parameter :: F_QT_ZU = '../data/zubko/sedust_zubko.h5'
   character(len=*), parameter :: F_ZDA_CFG = '../data/zubko/ZDA_BARE_GR_S_Config.dat'
   character(len=*), parameter :: D_ZUBKO   = '../data/zubko/'
   ! The two models defined by DustEM input files: THEMIS (Jones et al. 2013,
   ! Koehler et al. 2014) and Guillet et al. (2018) Model D, which Hensley &
   ! Draine (2023) Sec. 6.2.2 compare the astrodust model against.
   character(len=*), parameter :: F_GRAIN_THEMIS = '../data/themis/GRAIN_J13.DAT'
   character(len=*), parameter :: D_THEMIS       = '../data/themis/'
   character(len=*), parameter :: F_GRAIN_G18D   = '../data/g18d/GRAIN_G17_ModelD.DAT'
   character(len=*), parameter :: D_G18D         = '../data/g18d/'

   character(len=8), parameter :: STAGES(2) = ['S1      ', 'S2      ']
   real(wp),         parameter :: T_LO  = 2.7_wp
   real(wp),         parameter :: T_HI  = 5.0e3_wp
   integer,          parameter :: NT_IN = 200

   type(run_options_t) :: opt
   real(wp)            :: U_field
   integer             :: narg, iarg, n_submodel
   logical             :: taken
   character(len=64)   :: arg
   character(len=32)   :: model, submodel
   character(len=160)  :: tag

   ! ---- the axes this program has a referent for -------------------------
   ! Declared before a single argument is read.  Every model here is heated by
   ! a radiation field and solved by a stochastic-heating solver, so those axes
   ! belong to all three; the graphite of the PAH xi blend and the astrodust
   ! Stage-1 enthalpy exist only for the models that compute them, and are
   ! added below once the model is known.
   call declare_run_options(opt, program='calc_sed', &
        subjects=[character(len=16):: 'astrodust', 'dl07', 'mrn', 'zubko', &
                                      'themis', 'g18d'], &
        solver=.true., grid=.true., field=.true., emission=.true., qm_size=.true.)

   narg = command_argument_count()
   if (narg < 1) then
      call print_usage();  stop 1
   end if
   call get_command_argument(1, model)
   call read_run_subject(opt, trim(model), taken)
   if (.not. taken) then
      write(*,'(a,a)') ' calc_sed: unknown model ', trim(model)
      call print_usage();  stop 1
   end if
   select case (trim(model))
   case ('astrodust')
      call widen_run_options(opt, graphite=.true., stage1=.true.)
      U_field = 1.585_wp                ! log U = 0.20, the HD23 best fit
   case ('dl07')
      call widen_run_options(opt, graphite=.true., pah_xsec=.true.)
      U_field = 1.0_wp                  ! the DL07 reference intensity
   case ('mrn')
      U_field = 1.0_wp                  ! the Mathis ISRF at its nominal strength
   case ('zubko')
      U_field = 1.0_wp                  ! the reference intensity of its published SED
   case ('themis', 'g18d')
      ! DustEM's own reference run for both is the Mathis ISRF at G0 = 1, the
      ! scaling their GRAIN files state.
      U_field = 1.0_wp
   end select

   ! ---- settings, in any order ------------------------------------------
   ! Reading only records WHAT was asked for; the filename tag is assembled
   ! from the settings afterwards, in one fixed order, so that the order they
   ! are typed in does not change where the run lands.
   submodel = 'mw31_60';  n_submodel = 0
   do iarg = 2, narg
      call get_command_argument(iarg, arg)
      call read_run_option(trim(arg), opt, taken)
      if (.not. taken) then
         ! Not a word of any axis.  For dl07 that is the size-distribution
         ! submodel; the other two models have no positional argument left.
         if (trim(model) == 'dl07') then
            submodel = trim(arg);  n_submodel = n_submodel + 1
         else
            write(*,'(a,a,a,a)') ' calc_sed: ', trim(model), &
               ' has no setting called ', trim(arg)
            call print_usage();  stop 1
         end if
      end if
   end do
   ! The size distribution is an axis like any other: two names of it is a
   ! contradiction, not a request for the second.
   if (n_submodel > 1) then
      write(*,'(a)') ' calc_sed: give at most one dl07 submodel name'
      stop 1
   end if
   call check_run_options(opt)

   ! `logU=X` scales whichever field the model is illuminated with.
   if (opt%set_logU) U_field = 10.0_wp ** opt%logU

   tag = run_options_tag(opt)

   select case (trim(model))
   case ('astrodust');  call solve_astrodust()
   case ('dl07');       call solve_dl07()
   case ('mrn');        call solve_mrn()
   case ('zubko');      call solve_zubko()
   case ('themis');     call solve_dustem(F_GRAIN_THEMIS, D_THEMIS, 'themis', &
                             'THEMIS (Jones et al. 2013 / Koehler et al. 2014)')
   case ('g18d');       call solve_dustem(F_GRAIN_G18D, D_G18D, 'g18d', &
                             'Guillet et al. (2018) Model D')
   end select

contains

   ! ===================================================================
   subroutine print_usage()
      write(*,'(a)') ' usage:'
      write(*,'(a)') '   ./calc_sed.x astrodust [settings ...]'
      write(*,'(a)') '   ./calc_sed.x dl07 [submodel] [settings ...]'
      write(*,'(a)') '   ./calc_sed.x mrn [settings ...]'
      write(*,'(a)') '   ./calc_sed.x zubko [settings ...]'
      write(*,'(a)') '   ./calc_sed.x themis [settings ...]'
      write(*,'(a)') '   ./calc_sed.x g18d [settings ...]'
      write(*,'(a)') ''
      write(*,'(a)') '   dl07 submodel : mw31_00..mw31_60 (default mw31_60),'// &
                     ' lmc2_00, lmc2_05, lmc2_10, smc'
      write(*,'(a)') '   graphite  : gra_d03_sphere | gra_d16_sphere |'// &
                     ' gra_d16_spheroid'
      write(*,'(a)') '               the graphite of the PAH xi blend'// &
                     ' (astrodust, dl07)'
      write(*,'(a)') '   enthalpy  : c2          astrodust Stage-1'// &
                     ' density-corrected prefactor'
      call write_run_option_usage(opt)
      write(*,'(a)') ''
      write(*,'(a)') ' Writes output/sed_<model>[_<submodel>][_<tag>]'// &
                     '[_<stage>].dat, every'
      write(*,'(a)') ' non-default setting appearing in the tag so a variant'// &
                     ' run cannot'
      write(*,'(a)') ' overwrite the production one.'
   end subroutine print_usage


   subroutine report_header(title)
      character(len=*), intent(in) :: title
      write(*,'(a)') '=========================================================='
      write(*,'(a,a)') ' calc_sed: ', title
      write(*,'(a)') '=========================================================='
   end subroutine report_header


   function tagged(stem) result(fn)
      ! output/<stem>[_<tag>] -- the settings tag, appended only when the run
      ! is not the default one, so the production filenames carry no tag.
      character(len=*), intent(in) :: stem
      character(len=192) :: fn
      if (len_trim(tag) > 0) then
         fn = 'output/'//trim(stem)//'_'//trim(tag)
      else
         fn = 'output/'//trim(stem)
      end if
   end function tagged


   ! ===================================================================
   subroutine solve_astrodust()
      ! Hensley & Draine (2023) astrodust + PAH, for a single Mathis-ISRF cell,
      ! solved for the two enthalpy stages (S1, S2) and the PAH population.
      real(wp), allocatable :: J_lam(:), lamI_lam(:)
      character(len=192)    :: stem
      integer :: is

      stem = trim(tagged('sed_astrodust'))

      call report_header('Hensley & Draine (2023) astrodust + PAH SED')
      write(*,'(a,a)')    ' Q table     : ', F_QT_AD
      write(*,'(a)')      ' size_dist   : HD23 eqs. (17), (24), analytic'
      write(*,'(a,f8.3)') ' U_mathis    : ', U_field
      write(*,'(a,i0)')   ' NT (T grid) : ', NT_IN
      call apply_run_options(opt)
      write(*,'(a,a)')    ' PAH graphite: ', trim(qpah_graphite_source)
      write(*,'(a,l1)')   ' Stage-1 c2  : ', opt%stage1_density_corrected
      call report_run_options(opt)
      ! Say where the run will land before it starts, so a mistyped setting is
      ! caught now rather than after the solve.
      write(*,'(a,a,a)')  ' output files: ', trim(stem), '_<stage>.dat'

      call use_tmatrix_euv_band_optics()
      if (opt%euv) then
         call sed_init(F_QT_AD, NT_IN, T_LO, T_HI, include_euv=.true.)
      else
         call sed_init(F_QT_AD, NT_IN, T_LO, T_HI)
      end if
      write(*,'(a,i0,a)') ' sed_init done. NLAM=', NLAM, ' wavelengths cached.'
      write(*,'(a)') ''

      allocate(J_lam(NLAM), lamI_lam(NLAM))
      call J_Mathis(U_field, lam, J_lam)
      if (opt%hard_euv_field) call add_hard_euv_component(lam, J_lam)

      if (trim(stoch_method) == 'qm') then
         ! QM batch mode: all 4 grain types (S1, S2, PAH-neutral, PAH-cation)
         ! are processed in one parallel region for maximum thread utilisation.
         block
            real(wp), allocatable :: lamI_stages(:,:), lamI_pah_b(:)
            allocate(lamI_stages(NLAM, 2), lamI_pah_b(NLAM))
            write(*,'(a)') ' solving all stages (QM batch) ...'
            call sed_solve_qm_batch(J_lam, lamI_stages, lamI_pah_b)
            do is = 1, 2
               call write_astrodust_stage(trim(stem), trim(STAGES(is)), &
                                          lamI_stages(:, is), .false.)
            end do
            call write_astrodust_stage(trim(stem), 'PAH', lamI_pah_b, .true.)
            deallocate(lamI_stages, lamI_pah_b)
         end block
      else
         ! Temperature-window solvers (and the equilibrium branch), sequential
         ! per grain type.
         do is = 1, 2
            write(*,'(a,a,a)', advance='no') ' solving stage ', trim(STAGES(is)), ' ... '
            call sed_solve(J_lam, trim(STAGES(is)), lamI_lam)
            call write_astrodust_stage(trim(stem), trim(STAGES(is)), lamI_lam, .false.)
         end do
         write(*,'(a)', advance='no') ' solving stage PAH ... '
         call sed_solve_pah(J_lam, lamI_lam)
         call write_astrodust_stage(trim(stem), 'PAH', lamI_lam, .true.)
      end if

      write(*,'(a)') ''
      write(*,'(a)') ' calc_sed: done.'
      deallocate(J_lam, lamI_lam)
   end subroutine solve_astrodust


   subroutine write_astrodust_stage(stem, stage, lamI, is_pah)
      character(len=*), intent(in) :: stem, stage
      real(wp),         intent(in) :: lamI(:)
      logical,          intent(in) :: is_pah
      integer :: uu, kk
      character(len=224) :: fn
      fn = trim(stem)//'_'//trim(stage)//'.dat'
      open(newunit=uu, file=trim(fn), status='replace', action='write')
      ! The intensity is written out rather than quoted from the default, so
      ! that a logU= run does not claim the production value it did not use.
      if (is_pah) then
         write(uu,'(a,f0.3)') '# DH21 PAH SED for Mathis ISRF, U = ', U_field
         write(uu,'(a)') '# DL07 PAH cross sections, neutral + cation mixed by f_ion(a)'
      else
         write(uu,'(a,f0.3)') '# DH21 astrodust SED for Mathis ISRF, U = ', U_field
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
   end subroutine write_astrodust_stage


   ! ===================================================================
   subroutine solve_dl07()
      ! Draine & Li (2007) silicate + carbonaceous (PAH + graphite) with the
      ! WD01 size distributions.  Reproduces the DL07spec reference spectra
      ! at U = 1.
      real(wp), allocatable :: J_lam(:), lamI_tot(:), lamI_sil(:), lamI_carb(:)
      integer  :: sd_index, k, u
      logical  :: own_optics
      character(len=192) :: fname

      select case (trim(submodel))
      case ('mw31_00'); sd_index = 1
      case ('mw31_10'); sd_index = 2
      case ('mw31_20'); sd_index = 3
      case ('mw31_30'); sd_index = 4
      case ('mw31_40'); sd_index = 5
      case ('mw31_50'); sd_index = 6
      case ('mw31_60'); sd_index = 7
      case ('lmc2_00'); sd_index = 29
      case ('lmc2_05'); sd_index = 30
      case ('lmc2_10'); sd_index = 31
      case ('smc');     sd_index = 32
      case default
         write(*,'(a,a)') ' calc_sed: unknown dl07 submodel or setting ', trim(submodel)
         call print_usage();  stop 1
      end select

      ! DL07 / WD MW R_V=3.1_60 = the 2003 model: WD01 abundances reduced by
      ! the Draine-2003a factor 0.93, heated by the MMP83 field.  The
      ! "mathis_corrected" field (w_4000 = 1.65e-13) IS the canonical MMP83
      ! (it matches u_star = 8.64e-13 to 0.4%); the literal 1e-13 of the paper
      ! under-normalizes the optical band by ~7%, and `mathis_orig` asks for
      ! that one.
      qpah_graphite_source   = 'd03_sphere'  ! DL07 carbonaceous uses D03 graphite
      gd_apply_d03_reduction = .true.        ! 0.93 abundance reduction = the 2003 model
      nc_coeff   = 470.0d0                   ! DL07 Nc coefficient (rho~2.2)
      nc_integer = .true.                    ! Nc is rounded to an integer
      call apply_run_options(opt)

      ! Naming a graphite source is a sensitivity test on the xi blend, and the
      ! stored cross sections cannot answer it: they were computed once, with
      ! this model's own 'd03_sphere' graphite, and reading them back would
      ! leave the request with no effect at all.  So a named graphite solves
      ! every optic from the dielectric functions instead -- which is what a
      ! comparison of two graphite choices needs anyway, both variants then
      ! taking the same route.
      ! The same holds for the LD01 carbonaceous vintage: the stored tables
      ! are the DL07 one, so asking for the other has to solve the optics.
      own_optics = (len_trim(opt%graphite) > 0 .or. trim(opt%pah_xsec) /= 'dl07')

      fname = trim(tagged('sed_dl07_'//trim(submodel)))//'.dat'

      call report_header('Draine & Li (2007) silicate + carbonaceous SED')
      write(*,'(a,a)')   ' submodel      : ', trim(submodel)
      write(*,'(a,i0)')  ' WD01 index    : ', sd_index
      write(*,'(a,f8.3)')' U (Mathis)    : ', U_field
      write(*,'(a,a)')   ' graphite      : ', trim(qpah_graphite_source)
      write(*,'(a,a)')   ' PAH xsec      : ', trim(opt%pah_xsec)
      if (own_optics) write(*,'(a)') ' optics        : solved from the'// &
         ' dielectric functions (no stored table)'
      write(*,'(a,l1)')  ' d03_reduction : ', gd_apply_d03_reduction
      write(*,'(a,f6.1,a,l1)') ' Nc_coeff      : ', nc_coeff, '   Nc_integer: ', nc_integer
      call report_run_options(opt)
      write(*,'(a,a)')   ' output file   : ', trim(fname)

      ! One product, two views: include_euv picks which part of its axis.
      if (own_optics) then
         call sed_init_dl07(F_QT_DL, sd_index, U_field, NT_IN, T_LO, T_HI, &
                            include_euv=opt%euv, stored_q_dir='')
      else
         call sed_init_dl07(F_QT_DL, sd_index, U_field, NT_IN, T_LO, T_HI, &
                            include_euv=opt%euv)
      end if
      write(*,'(a,i0,a)') ' sed_init_dl07 done. NLAM=', NLAM, '.'

      allocate(J_lam(NLAM), lamI_tot(NLAM), lamI_sil(NLAM), lamI_carb(NLAM))
      call J_Mathis(U_field, lam, J_lam)
      if (opt%hard_euv_field) call add_hard_euv_component(lam, J_lam)

      write(*,'(a)') ' solving DL07 SED (silicate + carbonaceous) ...'
      call sed_solve_dl07(J_lam, lamI_tot, lamI_sil, lamI_carb)

      open(newunit=u, file=trim(fname), status='replace', action='write')
      ! f5.2 holds every intensity the reference runs use; a field scaled far
      ! up by logU= would fill it with asterisks, so that one takes an exponent.
      if (U_field < 1000.0_wp) then
         write(u,'(a,a,a,f5.2)')  '# DL07 model SED (this work), ', trim(submodel), &
              ', Mathis ISRF U = ', U_field
      else
         write(u,'(a,a,a,es10.3)')'# DL07 model SED (this work), ', trim(submodel), &
              ', Mathis ISRF U = ', U_field
      end if
      ! Name the graphite the blend actually took, so a sensitivity run does
      ! not carry the model's own choice in its header.
      if (trim(qpah_graphite_source) == 'd03_sphere') then
         write(u,'(a)') '# silicate (D03) + carbonaceous (DL07 PAH Nc=470 + D03 graphite blend)'
      else
         write(u,'(a,a,a)') '# silicate (D03) + carbonaceous (DL07 PAH Nc=470 + ', &
              trim(qpah_graphite_source), ' graphite blend)'
      end if
      if (trim(opt%pah_xsec) /= 'dl07') &
         write(u,'(a,a)') '# carbonaceous absorption vintage: ', trim(opt%pah_xsec)
      write(u,'(a)') '# PAH ionization computed via WD01b grain charging (pah_ionfrac)'
      write(u,'(a)') '# columns: lambda[um]  lamI_total/NH  lamI_sil/NH  lamI_carb/NH'
      write(u,'(a)') '#          [erg s^-1 cm^-2 sr^-1 H^-1]'
      ! e3 on the intensities: on a grid carried into the EUV, lambda*I_lambda
      ! underflows to subnormals, and a two-digit exponent field drops the E
      ! there (4.29970510-319), which readers outside Fortran cannot parse.
      do k = 1, NLAM
         write(u,'(es14.6,3(1x,es16.8e3))') lam(k), lamI_tot(k), lamI_sil(k), lamI_carb(k)
      end do
      close(u)
      write(*,'(a,a)') ' wrote ', trim(fname)
      write(*,'(a)') ' calc_sed: done.'

      deallocate(J_lam, lamI_tot, lamI_sil, lamI_carb)
   end subroutine solve_dl07



   ! ===================================================================
   subroutine solve_mrn()
      ! Mathis, Rumpl & Nordsieck (1977) graphite + silicate spheres, one
      ! a^-3.5 power law per material over 0.005 - 0.25 um.  No PAHs: the
      ! model has none, so the SED carries no aromatic features and the
      ! shortest-wavelength emission is that of a stochastically heated 50 A
      ! graphite sphere.
      type(dust_model_t)    :: m
      real(wp), allocatable :: J_lam(:), lamI_tot(:), lamI_chan(:,:)
      integer  :: k, ic, u, status
      character(len=192) :: fname

      fname = trim(tagged('sed_mrn'))//'.dat'

      call report_header('Mathis, Rumpl & Nordsieck (1977) graphite + silicate emission')
      write(*,'(a,f8.3)') ' U (Mathis)    : ', U_field
      ! Before the build: the field convention fixes the CMB temperature that
      ! goes into the cooling term, and the solver choice is stamped on the
      ! model.
      call apply_run_options(opt)
      call report_run_options(opt)
      write(*,'(a,a)')    ' output file   : ', trim(fname)

      ! One product, two views: include_euv picks which part of its axis.
      call build_mrn(m, F_QT_MRN, NT_IN, T_LO, T_HI, status=status, &
                     include_euv=opt%euv)
      if (status /= 0) then
         write(*,'(a,i0)') ' calc_sed: build_mrn failed, status = ', status
         stop 1
      end if
      m%verbose = .true.
      write(*,'(a,i0,a,es11.4,a,es11.4,a)') ' build_mrn done. NLAM=', m%NLAM, &
           ',  grid = ', m%lam(1), ' to ', m%lam(m%NLAM), ' um'

      allocate(J_lam(m%NLAM), lamI_tot(m%NLAM), lamI_chan(m%NLAM, m%n_channel))
      call J_Mathis(U_field, m%lam, J_lam)
      if (opt%hard_euv_field) call add_hard_euv_component(m%lam, J_lam)

      write(*,'(a)') ' solving MRN SED (graphite + silicate) ...'
      call dust_emission(m, J_lam, lamI_tot, lamI_chan, status)
      if (status /= 0) then
         write(*,'(a,i0)') ' calc_sed: dust_emission failed, status = ', status
         stop 1
      end if

      open(newunit=u, file=trim(fname), status='replace', action='write')
      ! f5.2 holds every intensity the reference runs use; a field scaled far
      ! up by logU= would fill it with asterisks, so that one takes an exponent.
      if (U_field < 1000.0_wp) then
         write(u,'(a,f5.2)')  '# MRN (1977) model SED (this work),'// &
              ' Mathis ISRF U = ', U_field
      else
         write(u,'(a,es10.3)')'# MRN (1977) model SED (this work),'// &
              ' Mathis ISRF U = ', U_field
      end if
      write(u,'(a)') '# graphite + silicate spheres, dn/da = A_i a^-3.5,'// &
           ' 0.005 - 0.25 um, Mie on the D03 dielectric functions'
      write(u,'(a)') '# normalization: Draine & Lee (1984), log10 A ='// &
           ' -25.16 (graphite), -25.11 (silicate) [cm^2.5/H]'
      write(u,'(a,a)') '# solver = ', trim(m%stoch_method)
      if (opt%hard_euv_field) write(u,'(a,es10.3,a,es10.3)') &
           '# artificial hard component below the Lyman limit: W*B_lambda(T), T = ', &
           T_HARD_EUV, ' K, W = ', W_HARD_EUV
      write(u,'(a)') '# columns: lambda[um]  lamI_total/NH  lamI_GRA/NH  lamI_SIL/NH'
      write(u,'(a)') '#          [erg s^-1 cm^-2 sr^-1 H^-1]'
      ! e3 on the intensities: carried into the EUV, lambda*I_lambda underflows
      ! to subnormals and a two-digit exponent field drops the E there
      ! (4.29970510-319), which readers outside Fortran cannot parse.
      do k = 1, m%NLAM
         write(u,'(es14.6,3(1x,es16.8e3))') m%lam(k), lamI_tot(k), &
              (lamI_chan(k, ic), ic = 1, m%n_channel)
      end do
      close(u)
      write(*,'(a,a)') ' wrote ', trim(fname)
      write(*,'(a)') ' calc_sed: done.'

      deallocate(J_lam, lamI_tot, lamI_chan)
   end subroutine solve_mrn


   ! ===================================================================
   subroutine solve_zubko()
      ! Zubko, Dwek & Arendt (2004) BARE-GR-S: PAH + graphite + silicate with
      ! the ZDA size-distribution formula and the ZDA optics and calorimetry
      ! tables.
      !
      ! This is the one shipped model whose optics grid (1.0e-3 um, 1.24 keV)
      ! reaches far past the band a transported radiation field occupies
      ! (0.0912 um for the Mathis ISRF), so anything in the stochastic-heating
      ! path that reads the grid where it should read the field shows up here
      ! and nowhere else -- the regression recorded in
      ! docs/EUV_EXTENSION_HOST_REGRESSION.md.  `euv` builds it on the whole
      ! ZDA range and `hardfield` puts photons above 13.6 eV into the field;
      ! running both is that regression.
      type(dust_model_t)    :: m
      real(wp), allocatable :: J_lam(:), lamI_tot(:), lamI_chan(:,:)
      integer  :: k, ic, u, status
      character(len=192) :: fname

      fname = trim(tagged('sed_zubko'))//'.dat'

      call report_header('Zubko/ZDA BARE-GR-S dust emission')
      write(*,'(a,f8.3)') ' U (Mathis)    : ', U_field
      ! Before the build: the field convention fixes the CMB temperature that
      ! goes into the cooling term, and the solver choice is stamped on the
      ! model.
      call apply_run_options(opt)
      call report_run_options(opt)
      write(*,'(a,a)')    ' output file   : ', trim(fname)

      call build_zubko(m, F_ZDA_CFG, D_ZUBKO, NT_IN, T_LO, T_HI, status, &
                       include_euv=opt%euv, qtable_path=F_QT_ZU)
      if (status /= 0) then
         write(*,'(a,i0)') ' calc_sed: build_zubko failed, status = ', status
         stop 1
      end if
      m%verbose = .true.
      write(*,'(a,i0,a,es11.4,a,es11.4,a)') ' build_zubko done. NLAM=', m%NLAM, &
           ',  grid = ', m%lam(1), ' to ', m%lam(m%NLAM), ' um'

      allocate(J_lam(m%NLAM), lamI_tot(m%NLAM), lamI_chan(m%NLAM, m%n_channel))
      call J_Mathis(U_field, m%lam, J_lam)
      if (opt%hard_euv_field) call add_hard_euv_component(m%lam, J_lam)
      call report_field_extent(m%lam, J_lam)

      write(*,'(a)') ' solving Zubko SED (PAH + graphite + silicate) ...'
      call dust_emission(m, J_lam, lamI_tot, lamI_chan, status)
      if (status /= 0) then
         write(*,'(a,i0)') ' calc_sed: dust_emission failed, status = ', status
         stop 1
      end if

      open(newunit=u, file=trim(fname), status='replace', action='write')
      ! f5.2 holds every intensity the reference runs use; a field scaled far
      ! up by logU= would fill it with asterisks, so that one takes an exponent.
      if (U_field < 1000.0_wp) then
         write(u,'(a,f5.2)')  '# Zubko/ZDA BARE-GR-S model SED (this work),'// &
              ' Mathis ISRF U = ', U_field
      else
         write(u,'(a,es10.3)')'# Zubko/ZDA BARE-GR-S model SED (this work),'// &
              ' Mathis ISRF U = ', U_field
      end if
      write(u,'(a,a)')    '# solver = ', trim(m%stoch_method)
      if (opt%hard_euv_field) write(u,'(a,es10.3,a,es10.3)') &
           '# artificial hard component below the Lyman limit: W*B_lambda(T), T = ', &
           T_HARD_EUV, ' K, W = ', W_HARD_EUV
      write(u,'(a)') '# columns: lambda[um]  lamI_total/NH  lamI_PAH/NH  lamI_GRA/NH  lamI_SIL/NH'
      write(u,'(a)') '#          [erg s^-1 cm^-2 sr^-1 H^-1]'
      ! e3 on the intensities: this model's optics grid already starts at
      ! 1e-3 um (1.24 keV), where lambda*I_lambda underflows to subnormals, and
      ! a two-digit exponent field drops the E there (4.29970510-319), which
      ! readers outside Fortran cannot parse.
      do k = 1, m%NLAM
         write(u,'(es14.6,4(1x,es16.8e3))') m%lam(k), lamI_tot(k), &
              (lamI_chan(k, ic), ic = 1, m%n_channel)
      end do
      close(u)
      write(*,'(a,a)') ' wrote ', trim(fname)
      write(*,'(a)') ' calc_sed: done.'

      deallocate(J_lam, lamI_tot, lamI_chan)
   end subroutine solve_zubko


   ! ===================================================================
   subroutine solve_dustem(grain_path, data_dir, stem, title)
      ! A model defined by DustEM input files: THEMIS, or Guillet et al. (2018)
      ! Model D.  Both are built from the GRAIN_*.DAT that defines them plus the
      ! optics and calorimetry tables under their own oprop/ and hcap/, all read
      ! unchanged, so the emission solved here rests on the same numbers the
      ! transport curve calc_kext.x writes for them does.
      !
      ! Like the Zubko model, their optics grid starts far shortward of the band
      ! a transported radiation field occupies (0.04 um against 0.0912 um), so
      ! `euv` builds them on the whole grid and its absence cuts it at the Lyman
      ! limit.
      !
      ! One output column per grain population, in the order of the GRAIN file,
      ! which is the layout of DustEM's own SED_*.RES.
      character(len=*), intent(in) :: grain_path, data_dir, stem, title
      type(dust_model_t)    :: m
      real(wp), allocatable :: J_lam(:), lamI_tot(:), lamI_chan(:,:)
      integer  :: k, ic, u, status
      character(len=192) :: fname
      character(len=256) :: gmiss, cols

      fname = trim(tagged('sed_'//trim(stem)))//'.dat'

      call report_header(trim(title)//' dust emission')
      write(*,'(a,f8.3)') ' U (Mathis)    : ', U_field
      call apply_run_options(opt)
      call report_run_options(opt)
      write(*,'(a,a)')    ' output file   : ', trim(fname)

      ! The optics come from the model's own product, as they do for zubko, so
      ! that this SED and the transport curve in the same file are the same
      ! numbers.  A tree that has not run ./calc_qtable.x for this model yet
      ! has no product, and build_dustem reads the DustEM text tables instead;
      ! the two routes agree bit for bit.
      call build_dustem(m, grain_path, data_dir, NT_IN, T_LO, T_HI, status, &
                        include_euv=opt%euv, gsca_missing=gmiss, &
                        qtable_path=trim(data_dir)//'sedust_'//trim(stem)//'.h5')
      if (status /= 0) then
         write(*,'(a,i0)') ' calc_sed: build_dustem failed, status = ', status
         stop 1
      end if
      m%verbose = .true.
      if (.not. m%gsca_complete) then
         write(*,'(a,a)') ' *** no scattering asymmetry for: ', trim(gmiss)
         write(*,'(a)')   ' *** (the emission below does not use it; the transport' // &
                          ' curve does)'
      end if
      write(*,'(a,i0,a,es11.4,a,es11.4,a)') ' build_dustem done. NLAM=', m%NLAM, &
           ',  grid = ', m%lam(1), ' to ', m%lam(m%NLAM), ' um'

      allocate(J_lam(m%NLAM), lamI_tot(m%NLAM), lamI_chan(m%NLAM, m%n_channel))
      call J_Mathis(U_field, m%lam, J_lam)
      if (opt%hard_euv_field) call add_hard_euv_component(m%lam, J_lam)
      call report_field_extent(m%lam, J_lam)

      write(*,'(a,i0,a)') ' solving the SED over ', m%n_channel, ' grain populations ...'
      call dust_emission(m, J_lam, lamI_tot, lamI_chan, status)
      if (status /= 0) then
         write(*,'(a,i0)') ' calc_sed: dust_emission failed, status = ', status
         stop 1
      end if

      open(newunit=u, file=trim(fname), status='replace', action='write')
      if (U_field < 1000.0_wp) then
         write(u,'(a,a,f5.2)') '# ', trim(title)//' model SED (this work),'// &
              ' Mathis ISRF U = ', U_field
      else
         write(u,'(a,a,es10.3)') '# ', trim(title)//' model SED (this work),'// &
              ' Mathis ISRF U = ', U_field
      end if
      write(u,'(a,a)')    '# solver = ', trim(m%stoch_method)
      if (opt%hard_euv_field) write(u,'(a,es10.3,a,es10.3)') &
           '# artificial hard component below the Lyman limit: W*B_lambda(T), T = ', &
           T_HARD_EUV, ' K, W = ', W_HARD_EUV
      cols = '# columns: lambda[um]  lamI_total/NH'
      do ic = 1, m%n_channel
         cols = trim(cols)//'  lamI_'//trim(m%channel_name(ic))//'/NH'
      end do
      write(u,'(a)') trim(cols)
      write(u,'(a)') '#          [erg s^-1 cm^-2 sr^-1 H^-1]'
      ! e3 on the intensities: this grid starts at 0.04 um, where lambda*I_lambda
      ! underflows to subnormals and a two-digit exponent field drops the E.
      do k = 1, m%NLAM
         write(u,'(es14.6,100(1x,es16.8e3))') m%lam(k), lamI_tot(k), &
              (lamI_chan(k, ic), ic = 1, m%n_channel)
      end do
      close(u)
      write(*,'(a,a)') ' wrote ', trim(fname)
      write(*,'(a)') ' calc_sed: done.'

      deallocate(J_lam, lamI_tot, lamI_chan)
   end subroutine solve_dustem


   subroutine report_field_extent(lam_m, J_lam)
      ! Print the two competing definitions of "hardest photon" side by side:
      ! the short end of the model's optics grid, and the energy that
      ! hardest_photon_energy reads off the field.  The stochastic solvers use
      ! the second; for the Zubko model the first is 91 times larger unless
      ! 'euv' is set, which is the whole content of the regression.
      real(wp), intent(in) :: lam_m(:), J_lam(:)
      real(wp), parameter  :: ERG2EV = 1.0_wp / 1.602176565e-12_wp
      real(wp) :: u_field
      u_field = hardest_photon_energy(lam_m, J_lam)
      write(*,'(a,es11.4,a,es12.5,a)') ' optics grid short end : ', lam_m(1), &
           ' um  ->  ', 1.23984193_wp / lam_m(1), ' eV'
      write(*,'(a,es12.5,a)')          ' hardest photon in the field:   ', &
           u_field * ERG2EV, ' eV'
      write(*,'(a,es12.5,a)')          ' single-photon bound in use :   ', &
           max(13.6_wp / ERG2EV, u_field) * ERG2EV, ' eV'
   end subroutine report_field_extent

end program calc_sed
