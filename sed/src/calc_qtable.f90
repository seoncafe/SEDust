program calc_qtable
   ! Writes the (lambda, a_eff) cross-section tables of every shipped model, in
   ! two forms: the text products under data/qtable/, one on the non-ionizing
   ! grid and one carried into the ionizing band, and one HDF5 file per model,
   !
   !   data/astrodust/sedust_astrodust.h5   data/dl07/sedust_dl07.h5
   !   data/mrn/sedust_mrn.h5               data/zubko/sedust_zubko.h5
   !   data/themis/sedust_themis.h5         data/g18d/sedust_g18d.h5
   !
   ! which holds the WIDER of the two wavelength axes once, together with the
   ! index at which the non-ionizing part of it begins (see sedust_h5.f90).
   ! The narrow product is then a slice rather than a second copy, so the pair
   ! cannot disagree.  A build without HDF5 writes the text alone.
   !
   !   ./calc_qtable.x [astrodust | dl07 | mrn | zubko | themis | g18d | all]
   !                                                            (default: all)
   !
   ! The optics are computed HERE, from the dielectric functions and the PAH
   ! cross sections, by the same routines the SED solver calls -- these are our
   ! own tables, not a copy of anyone else's.  For DL07 that is what the model
   ! already does at build time; the tables record it instead of repeating it.
   ! For Zubko it is a genuine recomputation: the ZDA optics tables
   ! distributed with the Camps et al. (2015) benchmark are a
   ! homogeneous-sphere Mie calculation on Draine's optical constants, and
   ! these files are what the same calculation gives here, on the 2003
   ! revision of those constants (the one this tree carries and DL07 uses).
   ! The two are NOT identical -- docs/SEDust_user_manual, "The ZDA optics are
   ! a homogeneous-sphere Mie calculation", measures the difference.
   !
   ! THEMIS and G18 Model D are the exception, and are so by definition: those
   ! models ARE their published DustEM Q_/G_ tables, and there is no dielectric
   ! function behind them to solve anything from.  What their products hold is
   ! those tables interpolated onto each population's own model radii, which is
   ! the one step between the distribution's files and the numbers the size
   ! integral uses.  No text form is written beside them: the distribution's
   ! own tables are text already.
   !
   ! Grids.  Each model keeps its own, because they are the grids its optics
   ! are defined on:
   !
   !   astrodust lambda = the T-matrix Q table's, 1762 (EUV) or 1129
   !             a_eff  = that table's own 169 radii
   !   DL07    lambda = the astrodust Q table's, 1129 (non-EUV) or 1823 (EUV,
   !                    carried to 6.205e-5 um by the D03 optical constants)
   !           a_eff  = Draine's analytic 84-point grid, 3.548e-4 - 5.012 um
   !   MRN     lambda = the DL07 axis, these optics being the same Mie on the
   !                    same D03 dielectric functions
   !           a_eff  = the model's own 70-point grid, 0.005 - 0.25 um, the
   !                    sharp cutoffs of the power law ON its ends
   !   Zubko   lambda = the ZDA optics tables', 1201 (EUV, from 1e-3 um) or
   !                    865 (non-EUV, cut at the Lyman limit)
   !           a_eff  = the ZDA tables' own, 121 radii for silicate and
   !                    graphite, 28 for the PAHs
   !   THEMIS  lambda = oprop/LAMBDA.DAT, the 800-point DustEM grid, from
   !   G18D             0.04 um (non-EUV: 750, cut at the Lyman limit)
   !           a_eff  = the MODEL's own radii, the nsize points the GRAIN line
   !                    spaces evenly in ln a between its amin and amax:
   !                    4 x 25 for THEMIS, 10 + 25 + 25 for G18D
   !
   ! Layout follows the astrodust Q table: lambda-major, a_eff inner, one row
   ! per cell.  There is no regime flag column -- that belongs to the T-matrix,
   ! which selects between Rayleigh, T-matrix and geometric optics; everything
   ! here is one method over the whole table.
   use, intrinsic :: iso_fortran_env, only: real64
   use constants,         only: wp, PI
   use q_silicate_mod,    only: q_silicate_full
   use q_graphite_mod,    only: q_graphite_full
   use qpah,              only: qpah_dl07, qpah_graphite_source, nc_coeff, nc_integer
   use grain_dist_mod,    only: gd_apply_d03_reduction
   use zubko_io,          only: read_zubko_optics
   use q_table_mod,       only: load_q_table, n_lam_qt => n_lam, n_aeff_qt => n_aeff, &
                                lam_qt => lam_t, aeff_qt => aeff_t, &
                                qext_qt => qext, qabs_qt => qabs, qsca_qt => qsca, &
                                gpar_qt => gpar, flag_qt => flag
   use sedust_h5
   use sed_astrodust_mod, only: sed_init_dl07, build_mrn, dust_model_t, &
                                NLAM, NA, lam, aeff, &
                                RHO_ASTROSIL, RHO_GRAPHITE, UM2CM, &
                                d03_euv_lambda_floor
   use dustem_io,         only: dustem_pop_t, DUSTEM_MAXPOP, read_dustem_grain, &
                                dustem_size_distribution, read_dustem_wavelengths, &
                                read_dustem_qtable, read_dustem_gtable, &
                                optics_at_radii, dustem_population_name
   use sedust_product_mod, only: lyman_index
   use enthalpy_astrodust_mod, only: RHO_AD, RHO_PAH
   use sed_run_options, only: run_options_t, declare_run_options, &
                              read_run_subject, read_run_option
   implicit none

   ! This program WRITES data/astrodust/sedust_astrodust.h5, so it cannot read
   ! the astrodust Q from it -- it would be reading back its own output.  It
   ! reads the text table beside it, which the T-matrix sweep produces
   ! (cd tmatrix && make && ./run_tmatrix.x, then make lyman_cut, and install
   ! the two files here).  That table is not shipped, for the reason
   ! .gitignore states; a checkout that only RUNS the models never needs it,
   ! because every other program reads the HDF5 product.
   character(len=*), parameter :: F_QT  = '../data/astrodust/q_astrodust_P0.20_Fe0.00_1.400.dat'
   character(len=*), parameter :: F_QT_EUV = '../data/astrodust/q_astrodust_P0.20_Fe0.00_1.400_euv.dat'
   character(len=*), parameter :: D_ZUBKO = '../data/zubko/'
   ! Every product a model owns goes into that model's own directory, so the
   ! text tables and the HDF5 file land beside each other and beside the
   ! model's extinction curves.  There is no shared output directory: a host
   ! ships data/<model>/ for what it runs.
   character(len=*), parameter :: D_ASTRODUST = '../data/astrodust/'
   character(len=*), parameter :: D_DL07      = '../data/dl07/'
   character(len=*), parameter :: D_MRN       = '../data/mrn/'
   ! The two models defined by DustEM input files.  Each is its GRAIN_*.DAT
   ! plus the optics and calorimetry tables under its own oprop/ and hcap/,
   ! and the GRAIN file is the only one of them that belongs to one model
   ! alone -- the oprop tables are shared by every model naming the same gtype.
   character(len=*), parameter :: D_THEMIS       = '../data/themis/'
   character(len=*), parameter :: F_GRAIN_THEMIS = D_THEMIS//'GRAIN_J13.DAT'
   character(len=*), parameter :: D_G18D         = '../data/g18d/'
   character(len=*), parameter :: F_GRAIN_G18D   = D_G18D//'GRAIN_G17_ModelD.DAT'
   ! DL07 reference model, as calc_sed.x dl07 runs it.
   integer,  parameter :: SD_INDEX = 7
   real(wp), parameter :: U_ISRF = 1.0_wp
   integer,  parameter :: NT_IN = 100
   real(wp), parameter :: T_LO = 1.0_wp, T_HI = 3000.0_wp
   real(wp), parameter :: LAM_LYMAN = 0.0912_wp
   ! Where the optics come from.  Named here so that the provenance an HDF5
   ! group carries is the path the routine below actually opens, and cannot
   ! drift from it the way a header written by hand does.
   character(len=*), parameter :: F_SIL_D03 = '../data/dielectric/index_silD03'
   character(len=*), parameter :: F_GRA_D03 = &
        '../data/dielectric/index_CpaD03 + index_CpeD03'

   character(len=32)  :: which
   integer :: narg, iarg
   logical :: taken
   type(run_options_t) :: o
   character(len=64)   :: arg

   ! Which model's tables to compute is the only axis here.  These tables ARE
   ! the model definition -- the graphite of the xi blend, the PAH cross-section
   ! vintage, the size distribution are all fixed by which model this is -- so
   ! no other axis has a referent, and a word naming one is refused by name.
   call declare_run_options(o, program='calc_qtable', &
        subjects=[character(len=16):: 'astrodust', 'dl07', 'mrn', 'zubko', &
                                      'themis', 'g18d', 'all'])
   narg = command_argument_count()
   which = 'all'
   if (narg >= 1) then
      call get_command_argument(1, which)
      call read_run_subject(o, trim(which), taken)
      if (.not. taken) then
         write(*,'(a,a)') ' calc_qtable: unknown model ', trim(which)
         write(*,'(a)')   ' usage: ./calc_qtable.x [astrodust | dl07 | mrn |' // &
                          ' zubko | themis | g18d | all]'
         stop 1
      end if
   end if
   do iarg = 2, narg
      call get_command_argument(iarg, arg)
      call read_run_option(trim(arg), o, taken)
      if (.not. taken) then
         write(*,'(a,a)') ' calc_qtable: unknown argument ', trim(arg)
         write(*,'(a)')   ' usage: ./calc_qtable.x [astrodust | dl07 | mrn |' // &
                          ' zubko | themis | g18d | all]'
         stop 1
      end if
   end do

   select case (trim(which))
   case ('astrodust');  call write_astrodust_tables()
   case ('dl07');       call write_dl07_tables()
   case ('mrn');        call write_mrn_tables()
   case ('zubko');      call write_zubko_tables()
   case ('themis');     call write_themis_tables()
   case ('g18d');       call write_g18d_tables()
   case ('all');        call write_astrodust_tables()
                        call write_dl07_tables();  call write_mrn_tables()
                        call write_zubko_tables()
                        call write_themis_tables();  call write_g18d_tables()
   end select

contains


   ! ===================================================================
   subroutine write_astrodust_tables()
      ! The astrodust model is TWO materials, not three: one composite grain
      ! (silicate + carbon + iron in a single porous spheroid, which is the
      ! premise of HD23) and the PAHs.  It does not split into silicate and
      ! graphite the way DL07 and ZDA do.
      !
      ! The composite grain's Q comes from the T-matrix sweep, already stored
      ! as text under tmatrix/output/; this reads it rather than re-solving
      ! hours of T-matrix.  The PAHs are the DL07 absorption cross sections on
      ! the same 169 radii, which until now were the one part of any model
      ! recomputed at every build with nothing written down.
      real(wp), allocatable :: Qn(:,:), Qi(:,:)
      real(wp) :: Q1
      integer  :: jw, ja, u1, u2

      call load_q_table(F_QT_EUV)
      allocate(Qn(n_lam_qt, n_aeff_qt), Qi(n_lam_qt, n_aeff_qt))
      qpah_graphite_source = 'd16_sphere'    ! HD23: Nc = 417, D16 graphite
      nc_coeff   = 417.0_wp
      nc_integer = .false.
      do ja = 1, n_aeff_qt
         do jw = 1, n_lam_qt
            call qpah_dl07(0, aeff_qt(ja), lam_qt(jw), Q1);  Qn(jw, ja) = Q1
            call qpah_dl07(1, aeff_qt(ja), lam_qt(jw), Q1);  Qi(jw, ja) = Q1
         end do
      end do

      call open_table(u1, D_ASTRODUST//'q_astrodust_pah_neu.dat', &
           'HD23 astrodust PAH, neutral: the DL07 absorption cross sections (not Mie)', &
           n_lam_qt, n_aeff_qt, .false.)
      call open_table(u2, D_ASTRODUST//'q_astrodust_pah_ion.dat', &
           'HD23 astrodust PAH, cation: the DL07 absorption cross sections (not Mie)', &
           n_lam_qt, n_aeff_qt, .false.)
      do jw = 1, n_lam_qt
         do ja = 1, n_aeff_qt
            call put_abs_row(u1, lam_qt(jw), aeff_qt(ja), Qn(jw, ja))
            call put_abs_row(u2, lam_qt(jw), aeff_qt(ja), Qi(jw, ja))
         end do
      end do
      close(u1);  close(u2)
      write(*,'(a,i0,a,i0,a)') ' astrodust PAH: ', n_lam_qt, ' x ', n_aeff_qt, ' cells written'

      call h5_model_file('astrodust', lam_qt, lyman_index(lam_qt))
      call h5_add_component('astrodust', 'astrodust', aeff_qt, qext_qt, qabs_qt, qsca_qt, &
                            gpar_qt, RHO_AD, &
                            'random-orientation T-matrix on the DH21 astrodust spheroid', &
                            F_QT_EUV, flag_in = flag_qt)
      call h5_add_component('astrodust', 'pah_neu', aeff_qt, absonly = Qn, &
                            rho = RHO_PAH, &
                            method = 'DL07 PAH absorption cross sections, neutral', &
                            source = 'qpah.f90 (Nc = 417, HD23)')
      call h5_add_component('astrodust', 'pah_ion', aeff_qt, absonly = Qi, &
                            rho = RHO_PAH, &
                            method = 'DL07 PAH absorption cross sections, cation', &
                            source = 'qpah.f90 (Nc = 417, HD23)')
      deallocate(Qn, Qi)
   end subroutine write_astrodust_tables



   ! ===================================================================
   subroutine write_dl07_tables()
      ! The DL07 optics are Mie on the D03 dielectric functions plus the DL07
      ! PAH cross sections, evaluated by the very routines sed_init_dl07 calls,
      ! so these tables are that model's own optics written down.
      call dl07_settings()
      ! stored_q_dir = '' : this program WRITES the stored tables, so it must
      ! not read them -- a model built from them would dump back what it just
      ! read, and they would freeze at whatever they already held instead of
      ! following the optics routines.  F_QT is a text table, so no HDF5
      ! product is read either.  Solve from the dielectric functions, always.
      call sed_init_dl07(F_QT, SD_INDEX, U_ISRF, NT_IN, T_LO, T_HI, &
                         stored_q_dir='')
      call dl07_dump('', .false.)
      call dl07_settings()
      ! The floor comes from the model, not from a number written here: the
      ! extinction curve of the same model is built on it too, and the two must
      ! land on the same wavelength grid to share one axis in the HDF5 file.
      call sed_init_dl07(F_QT_EUV, SD_INDEX, U_ISRF, NT_IN, T_LO, T_HI, &
                         lam_min=d03_euv_lambda_floor(), stored_q_dir='')
      ! The HDF5 file carries the WIDER of the two axes and the index where the
      ! non-ionizing part of it begins, so the narrow product is a slice of it
      ! rather than a second copy to keep in step.  Hence only this pass writes.
      call dl07_dump('_euv', .true.)
   end subroutine write_dl07_tables


   subroutine dl07_settings()
      ! The 2003 DL07 model, exactly as calc_sed.x dl07 sets it up.
      qpah_graphite_source   = 'd03_sphere'
      gd_apply_d03_reduction = .true.
      nc_coeff   = 470.0_wp
      nc_integer = .true.
   end subroutine dl07_settings


   subroutine dl07_dump(tag, to_h5)
      character(len=*), intent(in) :: tag
      ! Write the HDF5 file as well as the text tables.  Set on the wide pass
      ! only; see write_dl07_tables.
      logical, intent(in) :: to_h5
      real(wp) :: Qe, Qs, Qa, g, Qn, Qi
      integer  :: u1, u2, u3, u4, jw, ja
      ! The same cells the text rows carry, kept so that the HDF5 datasets are
      ! written from one calculation rather than a second sweep of the optics.
      real(wp), allocatable :: se(:,:), sa(:,:), ss(:,:), sg(:,:)
      real(wp), allocatable :: ge(:,:), ga(:,:), gs(:,:), gg(:,:)
      real(wp), allocatable :: pn(:,:), pc(:,:)

      allocate(se(NLAM,NA), sa(NLAM,NA), ss(NLAM,NA), sg(NLAM,NA))
      allocate(ge(NLAM,NA), ga(NLAM,NA), gs(NLAM,NA), gg(NLAM,NA))
      allocate(pn(NLAM,NA), pc(NLAM,NA))

      call open_table(u1, D_DL07//'q_dl07_sil'//trim(tag)//'.dat', &
           'DL07 silicate: Mie on the D03 astrosilicate dielectric function', &
           NLAM, NA, .true.)
      call open_table(u2, D_DL07//'q_dl07_gra'//trim(tag)//'.dat', &
           'DL07 graphite: Mie on the D03 graphite dielectric functions, ' // &
           '1/3 E||c + 2/3 E-perp-c', NLAM, NA, .true.)
      call open_table(u3, D_DL07//'q_dl07_pah_neu'//trim(tag)//'.dat', &
           'DL07 PAH, neutral: the DL07 absorption cross sections (not Mie)', &
           NLAM, NA, .false.)
      call open_table(u4, D_DL07//'q_dl07_pah_ion'//trim(tag)//'.dat', &
           'DL07 PAH, cation: the DL07 absorption cross sections (not Mie)', &
           NLAM, NA, .false.)

      do jw = 1, NLAM
         do ja = 1, NA
            call q_silicate_full(aeff(ja), lam(jw), Qe, Qs, Qa, g)
            call put_row(u1, lam(jw), aeff(ja), Qe, Qa, Qs, g)
            se(jw,ja) = Qe;  sa(jw,ja) = Qa;  ss(jw,ja) = Qs;  sg(jw,ja) = g
            call q_graphite_full(aeff(ja), lam(jw), Qe, Qs, Qa, g)
            call put_row(u2, lam(jw), aeff(ja), Qe, Qa, Qs, g)
            ge(jw,ja) = Qe;  ga(jw,ja) = Qa;  gs(jw,ja) = Qs;  gg(jw,ja) = g
            call qpah_dl07(0, aeff(ja), lam(jw), Qn)
            call qpah_dl07(1, aeff(ja), lam(jw), Qi)
            call put_abs_row(u3, lam(jw), aeff(ja), Qn)
            call put_abs_row(u4, lam(jw), aeff(ja), Qi)
            pn(jw,ja) = Qn;  pc(jw,ja) = Qi
         end do
      end do
      close(u1);  close(u2);  close(u3);  close(u4)
      write(*,'(a,i0,a,i0,a)') ' dl07'//trim(tag)//': ', NLAM, ' x ', NA, ' cells written'

      if (to_h5) then
         call h5_model_file('dl07', lam, lyman_index(lam))
         call h5_add_component('dl07', 'sil', aeff, se, sa, ss, sg, RHO_ASTROSIL, &
              'Mie on the D03 astrosilicate dielectric function', F_SIL_D03)
         call h5_add_component('dl07', 'gra', aeff, ge, ga, gs, gg, RHO_GRAPHITE, &
              'Mie on the D03 graphite dielectric functions, 1/3 E||c + 2/3 E-perp-c', &
              F_GRA_D03)
         ! The carbonaceous charge states differ in absorption alone; their
         ! scattering IS the graphite component above, which is why they carry
         ! Q_abs only and the graphite density rather than one of their own.
         call h5_add_component('dl07', 'pah_neu', aeff, absonly = pn, &
              rho = RHO_GRAPHITE, &
              method = 'DL07 PAH absorption cross sections, neutral', &
              source = 'qpah.f90 (Nc = 470, DL07)')
         call h5_add_component('dl07', 'pah_ion', aeff, absonly = pc, &
              rho = RHO_GRAPHITE, &
              method = 'DL07 PAH absorption cross sections, cation', &
              source = 'qpah.f90 (Nc = 470, DL07)')
      end if
      deallocate(se, sa, ss, sg, ge, ga, gs, gg, pn, pc)
   end subroutine dl07_dump



   ! ===================================================================
   subroutine write_mrn_tables()
      ! The MRN optics are Mie on the D03 dielectric functions -- the SAME
      ! calculation the DL07 silicate and graphite come from -- evaluated on
      ! this model's own 70-point radius grid.  Nothing is copied from the DL07
      ! tables: those are on Draine's 84-point grid, which runs from 3.5 A to
      ! 5 um and shares no node with 0.005 - 0.25 um.
      !
      ! The grids come from build_mrn itself, with stored_q_dir = '' so that it
      ! reads none of the tables this program is about to write; that is what
      ! makes these files the model's own optics and not a copy of whatever
      ! they already held.
      call mrn_dump('', .false.)
      call mrn_dump('_euv', .true.)
   end subroutine write_mrn_tables


   subroutine mrn_dump(tag, wide)
      character(len=*), intent(in) :: tag
      ! The ionizing pass.  It alone writes the HDF5 file: that file carries
      ! the WIDER axis and the index where the non-ionizing part of it begins,
      ! so the narrow product is a slice rather than a second copy.
      logical, intent(in) :: wide
      type(dust_model_t) :: m
      real(wp) :: Qe, Qs, Qa, g
      integer  :: u1, u2, jw, ja
      real(wp), allocatable :: se(:,:), sa(:,:), ss(:,:), sg(:,:)
      real(wp), allocatable :: ge(:,:), ga(:,:), gs(:,:), gg(:,:)

      ! The axis is the astrodust Q table's, as the DL07 one is: these optics
      ! are the same Mie on the same dielectric functions, so the two models
      ! that share a calculation share a wavelength grid.  The floor of the
      ! wide pass comes from those functions rather than from a number written
      ! here, so this table and the extinction curve of the same model land on
      ! one grid and can share an axis on disk.
      if (wide) then
         call build_mrn(m, F_QT_EUV, NT_IN, T_LO, T_HI, &
                        lam_min=d03_euv_lambda_floor(), stored_q_dir='')
      else
         call build_mrn(m, F_QT, NT_IN, T_LO, T_HI, stored_q_dir='')
      end if

      allocate(se(NLAM,NA), sa(NLAM,NA), ss(NLAM,NA), sg(NLAM,NA))
      allocate(ge(NLAM,NA), ga(NLAM,NA), gs(NLAM,NA), gg(NLAM,NA))

      call open_table(u1, D_MRN//'q_mrn_sil'//trim(tag)//'.dat', &
           'MRN silicate: Mie on the D03 astrosilicate dielectric function', &
           NLAM, NA, .true.)
      call open_table(u2, D_MRN//'q_mrn_gra'//trim(tag)//'.dat', &
           'MRN graphite: Mie on the D03 graphite dielectric functions, ' // &
           '1/3 E||c + 2/3 E-perp-c', NLAM, NA, .true.)

      do jw = 1, NLAM
         do ja = 1, NA
            call q_silicate_full(aeff(ja), lam(jw), Qe, Qs, Qa, g)
            call put_row(u1, lam(jw), aeff(ja), Qe, Qa, Qs, g)
            se(jw,ja) = Qe;  sa(jw,ja) = Qa;  ss(jw,ja) = Qs;  sg(jw,ja) = g
            call q_graphite_full(aeff(ja), lam(jw), Qe, Qs, Qa, g)
            call put_row(u2, lam(jw), aeff(ja), Qe, Qa, Qs, g)
            ge(jw,ja) = Qe;  ga(jw,ja) = Qa;  gs(jw,ja) = Qs;  gg(jw,ja) = g
         end do
      end do
      close(u1);  close(u2)
      write(*,'(a,i0,a,i0,a)') ' mrn'//trim(tag)//': ', NLAM, ' x ', NA, ' cells written'

      if (wide) then
         call h5_model_file('mrn', lam, lyman_index(lam))
         call h5_add_component('mrn', 'gra', aeff, ge, ga, gs, gg, RHO_GRAPHITE, &
              'Mie on the D03 graphite dielectric functions, 1/3 E||c + 2/3 E-perp-c', &
              F_GRA_D03)
         call h5_add_component('mrn', 'sil', aeff, se, sa, ss, sg, RHO_ASTROSIL, &
              'Mie on the D03 astrosilicate dielectric function', F_SIL_D03)
      end if
      deallocate(se, sa, ss, sg, ge, ga, gs, gg)
   end subroutine mrn_dump


   ! ===================================================================
   subroutine write_zubko_tables()
      ! TWO sets of optics go into this model's product, and the file says which
      ! is which.
      !
      !   /qtable/{sil,gra,pah}           the tables the Camps et al. (2015)
      !                                   benchmark distributes, verbatim.  This
      !                                   is the DEFAULT the model is built on:
      !                                   the seven codes that benchmark compares
      !                                   read these files, so a run against it
      !                                   measures the solver and not an optics
      !                                   difference.
      !   /qtable/{sil,gra,pah}_mie_d03   this tree's own recomputation, Mie on
      !                                   the Draine (2003) optical constants,
      !                                   selectable through build_zubko's
      !                                   optics argument.
      !
      ! The text q_zubko_*.dat products are the recomputation, as before.
      !
      ! WHAT THE DISTRIBUTED FILES ARE.  Their own headers: the silicate and
      ! graphite are Zubko's multilayer-sphere code (1997-2002), converted by
      ! Misselt in 2002; the PAH file is Misselt's 2009 implementation of Li &
      ! Draine (2001) / Draine & Li (2007).  Measured against Gra_121_1201.dat
      ! cell by cell, that PAH file holds graphite Q_sca and <cos> at every size
      ! and all 1201 wavelengths, and graphite Q_abs shortward of 0.0585 um
      ! (21.2 eV, the DL07 PAH-to-graphite transition, 1/17.25 um in qpah's
      ! form) -- so the ZDA PAH population scatters as the graphite sphere of
      ! its own radius.  The recomputation below does the same, by calling
      ! q_graphite_full at the PAH radii; it used to store no scattering for
      ! that population at all.
      !
      ! Draine (2003) is what the recomputation is on, and it is evidently what
      ! the shipped tables were computed on too: this route reproduces their
      ! silicate to a mean 2.4e-6, whatever the label in their header says.
      real(wp), allocatable :: a_sil(:), a_gra(:), a_pah(:), lz(:)
      ! The distributed Q columns.
      real(wp), allocatable :: za_sil(:,:), zs_sil(:,:), zg_sil(:,:)
      real(wp), allocatable :: za_gra(:,:), zs_gra(:,:), zg_gra(:,:)
      real(wp), allocatable :: za_pah(:,:), zs_pah(:,:), zg_pah(:,:)
      ! This tree's recomputation, on the same axes.
      real(wp), allocatable :: se(:,:), sa(:,:), ss(:,:), sg(:,:)
      real(wp), allocatable :: ge(:,:), ga(:,:), gs(:,:), gg(:,:)
      real(wp), allocatable :: pe(:,:), pa(:,:), ps(:,:), pg(:,:)
      real(wp) :: r_sil, r_gra, r_pah
      integer  :: n_sil, n_gra, n_pah, nw, nwc, k0

      ! The ZDA PAHs are the Li & Draine (2001) / Draine & Li (2007) cross
      ! sections -- their own file says so, and states the density 2.24 g/cm^3
      ! that N_C = 470 (a/nm)^3 implies.  So this component is solved on the
      ! DL07 carbon-atom count and the D03 graphite continuum, the same setting
      ! the DL07 tables above use, and NOT on the HD23 astrodust one (N_C =
      ! 417, rho = 2.0, D16 turbostratic graphite), which would contradict the
      ! density the table itself carries.
      !
      ! Stated here rather than inherited: these are saved module variables, so
      ! without this line the tables this pass writes would depend on whether
      ! the DL07 pass had run first in the same process -- ./calc_qtable.x all
      ! and ./calc_qtable.x zubko wrote PAH tables differing by a mean 17%.
      call dl07_settings()

      ! Read them whole: the axes and densities the recomputation needs, and
      ! the Q columns that are now a stored set in their own right.
      call read_zubko_optics(D_ZUBKO//'suvSil_121_1201.dat', n_sil, nw, a_sil, lz, &
                             za_sil, zs_sil, r_sil, gpar=zg_sil)
      call read_zubko_optics(D_ZUBKO//'Gra_121_1201.dat', n_gra, nw, a_gra, lz, &
                             za_gra, zs_gra, r_gra, gpar=zg_gra)
      call read_zubko_optics(D_ZUBKO//'PAH_28_1201_neu.dat', n_pah, nw, a_pah, lz, &
                             za_pah, zs_pah, r_pah, gpar=zg_pah)

      ! Recompute once, on the WIDE grid; the narrow product is a row slice of
      ! it, so the two cannot drift apart.
      call zubko_recompute(nw, lz, a_sil, a_gra, a_pah, se, sa, ss, sg, &
                           ge, ga, gs, gg, pe, pa, ps, pg)

      ! The non-ionizing set is the same grid cut where lyman_index cuts it,
      ! which is also where build_zubko's include_euv cuts it, so the stored
      ! pair matches the pair a run can build.
      k0  = lyman_index(lz(1:nw))
      nwc = nw - k0 + 1

      call zubko_text('_euv', nw,  lz(1:nw), a_sil, a_gra, a_pah, r_sil, r_gra, r_pah, &
                      se, sa, ss, sg, ge, ga, gs, gg, pe, pa, ps, pg)
      call zubko_text('',    nwc, lz(k0:nw), a_sil, a_gra, a_pah, r_sil, r_gra, r_pah, &
                      se(k0:,:), sa(k0:,:), ss(k0:,:), sg(k0:,:), &
                      ge(k0:,:), ga(k0:,:), gs(k0:,:), gg(k0:,:), &
                      pe(k0:,:), pa(k0:,:), ps(k0:,:), pg(k0:,:))

      ! Only the wide pass goes into HDF5; the file records where the narrow
      ! view starts.  See write_dl07_tables.
      call h5_model_file('zubko', lz(1:nw), k0)

      ! The model as the benchmark distributes it -- the default.
      call h5_add_component('zubko', 'sil', a_sil, za_sil+zs_sil, za_sil, zs_sil, zg_sil, &
           r_sil, 'ZDA BARE-GR-S as distributed: homogeneous-sphere Mie', &
           D_ZUBKO//'suvSil_121_1201.dat (Zubko multilayer-sphere code, 1997-2002)')
      call h5_add_component('zubko', 'gra', a_gra, za_gra+zs_gra, za_gra, zs_gra, zg_gra, &
           r_gra, 'ZDA BARE-GR-S as distributed: homogeneous-sphere Mie', &
           D_ZUBKO//'Gra_121_1201.dat (Zubko multilayer-sphere code, 1997-2002)')
      call h5_add_component('zubko', 'pah', a_pah, za_pah+zs_pah, za_pah, zs_pah, zg_pah, &
           r_pah, 'ZDA BARE-GR-S as distributed: LD01/DL07 absorption above ' // &
           '0.0585 um, graphite Mie below it, graphite Q_sca and <cos> throughout', &
           D_ZUBKO//'PAH_28_1201_neu.dat (Misselt 2009, from Li & Draine 2001 / ' // &
           'Draine & Li 2007)')

      ! This tree's recomputation, on the same axes -- the optics argument.
      call h5_add_component('zubko', 'sil_mie_d03', a_sil, se, sa, ss, sg, r_sil, &
           'Mie on the D03 astrosilicate dielectric function', &
           F_SIL_D03//'; grids from '//D_ZUBKO//'suvSil_121_1201.dat')
      call h5_add_component('zubko', 'gra_mie_d03', a_gra, ge, ga, gs, gg, r_gra, &
           'Mie on the D03 graphite dielectric functions, 1/3 E||c + 2/3 ' // &
           'E-perp-c, free-electron term at each radius', &
           F_GRA_D03//'; grids from '//D_ZUBKO//'Gra_121_1201.dat')
      call h5_add_component('zubko', 'pah_mie_d03', a_pah, pe, pa, ps, pg, r_pah, &
           'DL07 PAH absorption cross sections (neutral) with the scattering ' // &
           'of the graphite sphere of the same radius, as the distributed ' // &
           'table carries it', &
           'qpah.f90 + '//F_GRA_D03//'; grids from '//D_ZUBKO//'PAH_28_1201_neu.dat')

      deallocate(za_sil, zs_sil, zg_sil, za_gra, zs_gra, zg_gra, za_pah, zs_pah, zg_pah)
      deallocate(se, sa, ss, sg, ge, ga, gs, gg, pe, pa, ps, pg)
   end subroutine write_zubko_tables


   subroutine zubko_recompute(nw, lz, a_sil, a_gra, a_pah, se, sa, ss, sg, &
                              ge, ga, gs, gg, pe, pa, ps, pg)
      ! This tree's own optics on the ZDA axes: Mie on the D03 constants for
      ! the two solid components, and for the PAHs the DL07 absorption with the
      ! scattering of the graphite sphere of the same radius -- which is what
      ! the distributed PAH table holds, measured cell by cell.
      integer,  intent(in) :: nw
      real(wp), intent(in) :: lz(:), a_sil(:), a_gra(:), a_pah(:)
      real(wp), allocatable, intent(out) :: se(:,:), sa(:,:), ss(:,:), sg(:,:)
      real(wp), allocatable, intent(out) :: ge(:,:), ga(:,:), gs(:,:), gg(:,:)
      real(wp), allocatable, intent(out) :: pe(:,:), pa(:,:), ps(:,:), pg(:,:)
      real(wp) :: Qe, Qs, Qa, g, Qn
      integer  :: jw, ja

      allocate(se(nw,size(a_sil)), sa(nw,size(a_sil)), ss(nw,size(a_sil)), sg(nw,size(a_sil)))
      allocate(ge(nw,size(a_gra)), ga(nw,size(a_gra)), gs(nw,size(a_gra)), gg(nw,size(a_gra)))
      allocate(pe(nw,size(a_pah)), pa(nw,size(a_pah)), ps(nw,size(a_pah)), pg(nw,size(a_pah)))

      do jw = 1, nw
         do ja = 1, size(a_sil)
            call q_silicate_full(a_sil(ja), lz(jw), Qe, Qs, Qa, g)
            se(jw,ja) = Qe;  sa(jw,ja) = Qa;  ss(jw,ja) = Qs;  sg(jw,ja) = g
         end do
         do ja = 1, size(a_gra)
            call q_graphite_full(a_gra(ja), lz(jw), Qe, Qs, Qa, g)
            ge(jw,ja) = Qe;  ga(jw,ja) = Qa;  gs(jw,ja) = Qs;  gg(jw,ja) = g
         end do
         do ja = 1, size(a_pah)
            call qpah_dl07(0, a_pah(ja), lz(jw), Qn)
            ! Scattering from the graphite sphere at the SAME radius the PAH
            ! table is tabulated on.  Its Q_abs is discarded: the absorption is
            ! qpah's, which already hands over to graphite shortward of
            ! 1/17.25 um itself.
            call q_graphite_full(a_pah(ja), lz(jw), Qe, Qs, Qa, g)
            pa(jw,ja) = Qn;  ps(jw,ja) = Qs;  pg(jw,ja) = g
            pe(jw,ja) = Qn + Qs
         end do
      end do
   end subroutine zubko_recompute


   subroutine zubko_text(tag, nw, lz, a_sil, a_gra, a_pah, r_sil, r_gra, r_pah, &
                         se, sa, ss, sg, ge, ga, gs, gg, pe, pa, ps, pg)
      ! The text products, for reading with an eye and for diff.  They hold the
      ! RECOMPUTATION; the distributed tables are already text, beside them.
      character(len=*), intent(in) :: tag
      integer,  intent(in) :: nw
      real(wp), intent(in) :: lz(:), a_sil(:), a_gra(:), a_pah(:)
      real(wp), intent(in) :: r_sil, r_gra, r_pah
      real(wp), intent(in) :: se(:,:), sa(:,:), ss(:,:), sg(:,:)
      real(wp), intent(in) :: ge(:,:), ga(:,:), gs(:,:), gg(:,:)
      real(wp), intent(in) :: pe(:,:), pa(:,:), ps(:,:), pg(:,:)
      integer :: u1, u2, u3, jw, ja

      call open_table(u1, D_ZUBKO//'q_zubko_sil'//trim(tag)//'.dat', &
           'ZDA BARE-GR-S silicate: Mie on the D03 astrosilicate dielectric function', &
           nw, size(a_sil), .true., rho_in = r_sil)
      call open_table(u2, D_ZUBKO//'q_zubko_gra'//trim(tag)//'.dat', &
           'ZDA BARE-GR-S graphite: Mie on the D03 graphite dielectric functions, ' // &
           '1/3 E||c + 2/3 E-perp-c, free-electron term per radius', &
           nw, size(a_gra), .true., rho_in = r_gra)
      call open_table(u3, D_ZUBKO//'q_zubko_pah'//trim(tag)//'.dat', &
           'ZDA BARE-GR-S PAH, neutral: DL07 absorption with the scattering of ' // &
           'the graphite sphere of the same radius, as the distributed table has it', &
           nw, size(a_pah), .true., rho_in = r_pah)

      do jw = 1, nw
         do ja = 1, size(a_sil)
            call put_row(u1, lz(jw), a_sil(ja), se(jw,ja), sa(jw,ja), ss(jw,ja), sg(jw,ja))
         end do
         do ja = 1, size(a_gra)
            call put_row(u2, lz(jw), a_gra(ja), ge(jw,ja), ga(jw,ja), gs(jw,ja), gg(jw,ja))
         end do
         do ja = 1, size(a_pah)
            call put_row(u3, lz(jw), a_pah(ja), pe(jw,ja), pa(jw,ja), ps(jw,ja), pg(jw,ja))
         end do
      end do
      close(u1);  close(u2);  close(u3)
      write(*,'(a,i0,a)') ' zubko'//trim(tag)//': ', nw, ' wavelengths written'
   end subroutine zubko_text



   ! ===================================================================
   subroutine write_themis_tables()
      ! THEMIS -- Jones et al. (2013) for the carbon grains, Koehler et al.
      ! (2014) for the silicates, Ysard et al. (2015) for the densities and
      ! mass ratios -- as the DustEM distribution carries it.
      !
      ! WHAT GOES INTO THE PRODUCT.  /grid is oprop/LAMBDA.DAT, the 800-point
      ! DustEM wavelength grid every population shares.  /qtable/<population>
      ! is Q_abs and Q_sca from oprop/Q_<gtype>.DAT and <cos theta> from
      ! oprop/G_<gtype>.DAT, put on the MODEL's own radii -- the nsize points
      ! GRAIN_J13.DAT spaces evenly in ln a between its amin and amax --
      ! linearly in radius, which is what DustEM does with the same tables and
      ! what build_dustem does when it reads them itself.  Q_ext is the sum of
      ! the two, the tables carrying absorption and scattering and not
      ! extinction.  Nothing is recomputed from a dielectric function: this
      ! model IS those tables.
      !
      ! GRAIN_J13.DAT names CM20 TWICE -- the large a-C:H/a-C grains as a power
      ! law with exponential decay, the small a-C grains as a log-normal -- and
      ! the two have different radius ranges, so their stored optics are
      ! different arrays and are named apart; see dustem_population_name.
      call write_dustem_product('themis', D_THEMIS, F_GRAIN_THEMIS, &
           'THEMIS (Jones et al. 2013, Koehler et al. 2014, Ysard et al. 2015)')
   end subroutine write_themis_tables


   ! ===================================================================
   subroutine write_g18d_tables()
      ! Guillet et al. (2018) Model D, as the DustEM distribution carries it.
      ! Three populations: PAHs, and amorphous-carbon and silicate PROLATE
      ! spheroids whose Q tables are already averaged over orientation, so
      ! they enter this scalar product exactly as a sphere table would.
      !
      ! NO g FOR THE TWO LARGE POPULATIONS.  The distribution ships no
      ! G_amCBE_0.3333x.DAT and no G_aSil2001BE6pctG_0.4x.DAT, so those two
      ! groups get Q_ext, Q_abs and Q_sca and NO g dataset.  The absence is the
      ! record: a zero written there would be indistinguishable from a measured
      ! zero, and a model built from this product would report an asymmetry
      ! parameter the model does not have.  build_dustem reads the absence back
      ! and comes to the same gsca_complete = .false. it reaches from the text
      ! tables.
      call write_dustem_product('g18d', D_G18D, F_GRAIN_G18D, &
           'Guillet et al. (2018) Model D')
   end subroutine write_g18d_tables


   subroutine write_dustem_product(model, ddir, grain_file, what)
      ! The two routines above, once.  What differs between the models is the
      ! GRAIN_*.DAT and where its files hang, both of which are arguments; the
      ! steps from those files to the stored arrays are the model's definition
      ! and are the same for any DustEM model.
      character(len=*), intent(in) :: model, ddir, grain_file, what
      type(dustem_pop_t) :: pops(DUSTEM_MAXPOP)
      real(wp) :: G0
      integer  :: npop, ip, nwave, ntab, nsize, k0
      logical  :: rok, has_g
      real(wp), allocatable :: lamd(:), a_tab(:), qa_tab(:,:), qs_tab(:,:), gg_tab(:,:)
      real(wp), allocatable :: ac_tab(:)
      real(wp), allocatable :: a_cm(:), lna(:), dnda(:), a_um(:)
      real(wp), allocatable :: qa(:,:), qs(:,:), qe(:,:), gfac(:,:)
      character(len=512) :: qpath, gpath
      character(len=140) :: pname

      call read_dustem_grain(grain_file, G0, npop, pops, ok=rok)
      if (.not. rok) then
         write(*,'(a,a)') ' calc_qtable: cannot use ', grain_file;  stop 1
      end if
      call read_dustem_wavelengths(ddir//'oprop/LAMBDA.DAT', nwave, lamd, ok=rok)
      if (.not. rok) then
         write(*,'(a,a)') ' calc_qtable: cannot read ', ddir//'oprop/LAMBDA.DAT'
         stop 1
      end if

      ! The whole DustEM grid goes in, with the index where its non-ionizing
      ! part begins; the narrow view is a slice of it, as it is for every other
      ! model.
      k0 = lyman_index(lamd)
      call h5_model_file(model, lamd, k0)

      do ip = 1, npop
         nsize = pops(ip)%nsize
         call dustem_size_distribution(pops(ip), a_cm, lna, dnda, ok=rok)
         if (.not. rok) then
            write(*,'(a,a)') ' calc_qtable: bad size distribution for ', &
                 trim(pops(ip)%gtype)
            stop 1
         end if
         allocate(a_um(nsize))
         a_um = a_cm / UM2CM

         qpath = ddir//'oprop/Q_'//trim(pops(ip)%gtype)//'.DAT'
         call read_dustem_qtable(trim(qpath), nwave, ntab, a_tab, qa_tab, qs_tab, ok=rok)
         if (.not. rok) then
            write(*,'(a,a)') ' calc_qtable: cannot read ', trim(qpath);  stop 1
         end if
         allocate(qa(nwave, nsize), qs(nwave, nsize), qe(nwave, nsize))
         call optics_at_radii(a_tab, qa_tab, a_um, qa, ok=rok, what='Q_abs')
         if (rok) call optics_at_radii(a_tab, qs_tab, a_um, qs, ok=rok, what='Q_sca')
         if (.not. rok) then
            write(*,'(a,a)') ' calc_qtable: optics do not cover the sizes of ', &
                 trim(pops(ip)%gtype)
            stop 1
         end if
         qe = qa + qs

         gpath = ddir//'oprop/G_'//trim(pops(ip)%gtype)//'.DAT'
         inquire(file=trim(gpath), exist=has_g)
         if (has_g) then
            call read_dustem_gtable(trim(gpath), nwave, ntab, ac_tab, gg_tab, ok=rok)
            if (.not. rok) then
               write(*,'(a,a)') ' calc_qtable: cannot read ', trim(gpath);  stop 1
            end if
            allocate(gfac(nwave, nsize))
            call optics_at_radii(ac_tab, gg_tab, a_um, gfac, ok=rok, what='<cos theta>')
            if (.not. rok) then
               write(*,'(a,a)') ' calc_qtable: the g table does not cover the sizes of ', &
                    trim(pops(ip)%gtype)
               stop 1
            end if
            deallocate(ac_tab, gg_tab)
         end if

         pname = dustem_population_name(pops, npop, ip)
         if (has_g) then
            call h5_add_component(model, trim(pname), a_um, qe, qa, qs, gfac, &
                 pops(ip)%rho, &
                 trim(what)//': the DustEM Q_ and G_ tables as distributed, ' // &
                 'interpolated onto this population''s model radii linearly in radius', &
                 trim(qpath)//' + '//trim(gpath)//'; radii from '//grain_file)
         else
            call h5_add_component(model, trim(pname), a_um, qe, qa, qs, &
                 rho = pops(ip)%rho, &
                 method = trim(what)//': the DustEM Q_ table as distributed, ' // &
                 'interpolated onto this population''s model radii linearly in radius; ' // &
                 'the distribution ships no G_ file, so this population has no <cos theta>', &
                 source = trim(qpath)//'; radii from '//grain_file)
         end if
         write(*,'(a,a,a,i0,a,i0,a)') ' ', trim(model)//' '//trim(pname), ': ', &
              nwave, ' x ', nsize, ' cells written'

         deallocate(a_cm, lna, dnda, a_um, a_tab, qa_tab, qs_tab, qa, qs, qe)
         if (allocated(gfac)) deallocate(gfac)
      end do
      deallocate(lamd)
   end subroutine write_dustem_product



   ! ===================================================================
   ! HDF5 side.  One file per model, opened once by h5_model_file (which lays
   ! down the wavelength axis and the Lyman index) and then extended one
   ! component at a time.  The text products above are written unchanged; this
   ! is the primary form, that one is for reading with an eye.
   ! ===================================================================
   function h5_path(model) result(p)
      ! The model's product, in the model's own directory.
      character(len=*), intent(in) :: model
      character(len=128) :: p
      p = '../data/'//model//'/sedust_'//model//'.h5'
   end function h5_path


   subroutine h5_model_file(model, lam, i_lyman)
      character(len=*), intent(in) :: model
      real(real64),     intent(in) :: lam(:)
      integer,          intent(in) :: i_lyman
      integer(h5id_k) :: fid, gid
      logical :: ok
      if (.not. sedust_has_hdf5) then
         write(*,'(a)') ' calc_qtable: built without HDF5, writing text only'
         return
      end if
      call h5_begin(ok);  if (.not. ok) return
      call h5_create_file(h5_path(model), fid, ok)
      if (.not. ok) then
         write(*,'(a,a)') ' calc_qtable: cannot create ', h5_path(model)
         call h5_end();  return
      end if
      call h5_put_attr_s(fid, 'model', model)
      call h5_put_attr_s(fid, 'generator', 'SEDust sed/calc_qtable.x')
      call h5_put_attr_s(fid, 'layout', &
           'cross sections are (n_a, n_lam) in h5py; lambda is the fastest Fortran axis')
      call h5_group(fid, 'grid', gid, ok)
      call h5_write_1d(gid, 'lambda', lam, units='um', &
                       long_name='wavelength grid of the model, increasing')
      call h5_put_attr_i(gid, 'i_lyman', i_lyman)
      call h5_put_attr_d(gid, 'lambda_lyman_um', LAM_LYMAN)
      call h5_put_attr_s(gid, 'i_lyman_meaning', &
           'last index at or shortward of the Lyman limit; include_euv=.false. starts here')
      call h5_group_close(gid)
      call h5_close_file(fid)
      call h5_end()
      write(*,'(a,a)') ' wrote ', h5_path(model)
      ! The file is replaced, not extended: this pass lays down the wavelength
      ! axis, and an extinction curve integrated on the previous one would be
      ! filed against wavelengths it was not computed at.  So /kext goes with
      ! it, and calc_kext.x has to run again -- which is also the only order in
      ! which the two can be consistent.
      write(*,'(a,a,a)') '   /kext is gone with it; run  ./calc_kext.x ', model, ' euv'
   end subroutine h5_model_file


   subroutine h5_add_component(model, comp, aeff, qe, qa, qs, gg, rho, method, source, &
                               flag_in, absonly)
      ! One grain population.  A population that only absorbs -- the PAHs, whose
      ! cross sections are a prescription and not a Mie solution -- is given
      ! through absonly and gets Q_abs alone, so a reader cannot mistake an
      ! absent scattering term for a computed zero.  A population that
      ! scatters but has no PUBLISHED asymmetry parameter omits gg and gets no
      ! g dataset, for the same reason.
      character(len=*), intent(in) :: model, comp
      real(real64),     intent(in) :: aeff(:)
      real(real64),     intent(in), optional :: qe(:,:), qa(:,:), qs(:,:), gg(:,:)
      real(real64),     intent(in) :: rho
      character(len=*), intent(in) :: method, source
      integer,          intent(in), optional :: flag_in(:,:)
      real(real64),     intent(in), optional :: absonly(:,:)
      integer(h5id_k) :: fid, gid
      logical :: ok
      real(real64), allocatable :: rflag(:,:)
      if (.not. sedust_has_hdf5) return
      call h5_begin(ok);  if (.not. ok) return
      call h5fopen_or_skip(h5_path(model), fid, ok)
      if (.not. ok) then;  call h5_end();  return;  end if
      call h5_group(fid, 'qtable/'//comp, gid, ok)
      call h5_write_1d(gid, 'a_eff', aeff, units='um', long_name='grain effective radius')
      if (present(absonly)) then
         call h5_write_2d(gid, 'Q_abs', absonly, units='1', &
                          long_name='absorption efficiency')
         call h5_put_attr_s(gid, 'scatters', 'no -- absorption prescription only')
      else
         call h5_write_2d(gid, 'Q_ext', qe, units='1', long_name='extinction efficiency')
         call h5_write_2d(gid, 'Q_abs', qa, units='1', long_name='absorption efficiency')
         call h5_write_2d(gid, 'Q_sca', qs, units='1', long_name='scattering efficiency')
         if (present(gg)) then
            call h5_write_2d(gid, 'g',  gg, units='1', &
                             long_name='scattering asymmetry <cos>')
            call h5_put_attr_s(gid, 'scatters', 'yes')
         else
            ! Scattering with no asymmetry parameter: the two large G18D
            ! populations, for which the DustEM distribution ships no G_ file.
            ! No g dataset is written, because a zero here would be a
            ! measurement and this is the absence of one.
            call h5_put_attr_s(gid, 'scatters', &
                 'yes -- but no asymmetry parameter is published for it, so no g')
         end if
      end if
      if (present(flag_in)) then
         allocate(rflag(size(flag_in,1), size(flag_in,2)))
         rflag = real(flag_in, real64)
         call h5_write_2d(gid, 'regime', rflag, units='1', &
              long_name='0 T-matrix, 10 Rayleigh dipole (x<0.1), 20 geometric optics (x>50)')
         deallocate(rflag)
      end if
      call h5_put_attr_d(gid, 'rho_bulk_g_cm3', rho)
      call h5_put_attr_s(gid, 'method', method)
      call h5_put_attr_s(gid, 'source', source)
      call h5_group_close(gid)
      call h5_close_file(fid)
      call h5_end()
   end subroutine h5_add_component


   subroutine h5fopen_or_skip(path, fid, ok)
      ! Open for append.  h5_open_file is read-only, so use the library call.
      character(len=*), intent(in)  :: path
      integer(h5id_k),  intent(out) :: fid
      logical,          intent(out) :: ok
      call h5_open_rw(path, fid, ok)
   end subroutine h5fopen_or_skip

   ! ===================================================================
   subroutine open_table(u, path, what, nw, na_, with_sca, note_in, rho_in)
      integer,          intent(out) :: u
      character(len=*), intent(in)  :: path, what
      integer,          intent(in)  :: nw, na_
      logical,          intent(in)  :: with_sca
      ! Solid density of the material [g/cm^3].  A model built from these
      ! tables needs it for the grain mass and gets it here rather than from
      ! the file the optics used to come from.
      real(wp), intent(in), optional :: rho_in
      ! Anything the reader has to know about this particular table beyond
      ! what and on which grid; blank when there is nothing.
      character(len=*), intent(in), optional :: note_in
      character(len=512) :: note
      note = ''
      if (present(note_in)) note = note_in
      open(newunit=u, file=path, status='replace', action='write')
      write(u,'(a,a)')     '# ', what
      write(u,'(a,i0)')    '# a_eff stride: every 1 of ', na_
      write(u,'(a,i0)')    '# lambda stride: every 1 of ', nw
      write(u,'(a)')       '# Written by sed/calc_qtable.x from this tree''s own optics'
      write(u,'(a)')       '#   routines; see src/calc_qtable.f90 for the grids and inputs.'
      if (present(rho_in)) write(u,'(a,es13.6)') '# density [g/cm^3]: ', rho_in
      if (len_trim(note) > 0) then
         write(u,'(a,a)') '# ', note
      end if
      if (with_sca) then
         write(u,'(a)') '#   lambda[um]     a_eff[um]      Q_ext         Q_abs' // &
                        '          Q_sca         albedo         g'
      else
         write(u,'(a)') '#   lambda[um]     a_eff[um]      Q_abs'
      end if
   end subroutine open_table


   subroutine put_row(u, l, a, Qe, Qa, Qs, g)
      integer,  intent(in) :: u
      real(wp), intent(in) :: l, a, Qe, Qa, Qs, g
      real(wp) :: alb
      alb = 0.0_wp
      if (Qe > 0.0_wp) alb = Qs / Qe
      write(u,'(7(1x,es14.6e3))') l, a, Qe, Qa, Qs, alb, g
   end subroutine put_row


   subroutine put_abs_row(u, l, a, Qa)
      integer,  intent(in) :: u
      real(wp), intent(in) :: l, a, Qa
      write(u,'(3(1x,es14.6e3))') l, a, Qa
   end subroutine put_abs_row


   ! ===================================================================
end program calc_qtable
