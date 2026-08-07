program make_qtable
   ! Writes the (lambda, a_eff) cross-section tables of every shipped model, in
   ! two forms: the text products under data/qtable/, one on the non-ionizing
   ! grid and one carried into the ionizing band, and one HDF5 file per model,
   !
   !   data/sedust_astrodust.h5   data/sedust_dl07.h5   data/sedust_zubko.h5
   !
   ! which holds the WIDER of the two wavelength axes once, together with the
   ! index at which the non-ionizing part of it begins (see sedust_h5.f90).
   ! The narrow product is then a slice rather than a second copy, so the pair
   ! cannot disagree.  A build without HDF5 writes the text alone.
   !
   !   ./make_qtable.x [astrodust | dl07 | zubko | all]     (default: all)
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
   ! Grids.  Each model keeps its own, because they are the grids its optics
   ! are defined on:
   !
   !   astrodust lambda = the T-matrix Q table's, 1762 (EUV) or 1129
   !             a_eff  = that table's own 169 radii
   !   DL07    lambda = the astrodust Q table's, 1129 (non-EUV) or 1823 (EUV,
   !                    carried to 6.205e-5 um by the D03 optical constants)
   !           a_eff  = Draine's analytic 84-point grid, 3.548e-4 - 5.012 um
   !   Zubko   lambda = the ZDA optics tables', 1201 (EUV, from 1e-3 um) or
   !                    865 (non-EUV, cut at the Lyman limit)
   !           a_eff  = the ZDA tables' own, 121 radii for silicate and
   !                    graphite, 28 for the PAHs
   !
   ! Layout follows the astrodust Q table: lambda-major, a_eff inner, one row
   ! per cell.  There is no regime flag column -- that belongs to the T-matrix,
   ! which selects between Rayleigh, T-matrix and geometric optics; everything
   ! here is one method over the whole table.
   use, intrinsic :: iso_fortran_env, only: real64
   use constants,         only: wp, PI
   use q_silicate_mod,    only: q_silicate_full
   use q_graphite_mod,    only: q_graphite_full
   use qpah,              only: qpah_dl07, qpah_use_d03_graphite, nc_coeff, nc_integer
   use grain_dist_mod,    only: gd_apply_d03_reduction
   use zubko_io,          only: read_zubko_optics
   use q_table_mod,       only: load_q_table, n_lam_qt => n_lam, n_aeff_qt => n_aeff, &
                                lam_qt => lam_t, aeff_qt => aeff_t, &
                                qext_qt => qext, qabs_qt => qabs, qsca_qt => qsca, &
                                gpar_qt => gpar, flag_qt => flag
   use sedust_h5
   use sed_astrodust_mod, only: sed_init_dl07, NLAM, NA, lam, aeff, &
                                use_stored_q_tables, RHO_ASTROSIL, RHO_GRAPHITE, &
                                dl07_euv_lambda_floor
   use enthalpy_astrodust_mod, only: RHO_AD, RHO_PAH
   implicit none

   character(len=*), parameter :: F_QT  = '../data/astrodust/q_astrodust_P0.20_Fe0.00_1.400.dat'
   character(len=*), parameter :: F_QT_EUV = '../data/astrodust/q_astrodust_P0.20_Fe0.00_1.400_euv.dat'
   character(len=*), parameter :: F_SD  = '../data/release/size_distribution.dat'
   character(len=*), parameter :: D_ZUBKO = '../data/zubko/'
   ! Every product a model owns goes into that model's own directory, so the
   ! text tables and the HDF5 file land beside each other and beside the
   ! model's extinction curves.  There is no shared output directory: a host
   ! ships data/<model>/ for what it runs.
   character(len=*), parameter :: D_ASTRODUST = '../data/astrodust/'
   character(len=*), parameter :: D_DL07      = '../data/dl07/'
   ! DL07 reference model, as main_dl07.x runs it.
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
   integer :: narg

   ! This program WRITES the stored tables, so it must not read them: a model
   ! built from them would dump back what it just read, and the tables would
   ! freeze at whatever they already held instead of following the optics
   ! routines.  Solve from the dielectric functions here, always.
   use_stored_q_tables = .false.

   narg = command_argument_count()
   which = 'all'
   if (narg >= 1) call get_command_argument(1, which)

   select case (trim(which))
   case ('astrodust');  call write_astrodust_tables()
   case ('dl07');       call write_dl07_tables()
   case ('zubko');      call write_zubko_tables()
   case ('all');        call write_astrodust_tables()
                        call write_dl07_tables();  call write_zubko_tables()
   case default
      write(*,'(a)') ' usage: ./make_qtable.x [astrodust | dl07 | zubko | all]'
      stop 1
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
      qpah_use_d03_graphite = .false.        ! HD23: Nc = 417, D16 graphite
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

      call h5_model_file('astrodust', lam_qt, i_lyman_of(lam_qt))
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


   pure integer function i_lyman_of(lam) result(i)
      ! First index at or longward of the Lyman limit.  A reader given
      ! include_euv = .false. starts here; everything shortward is the ionizing
      ! band, where an interstellar radiation field carries no photon.
      real(wp), intent(in) :: lam(:)
      integer :: k
      i = 1
      do k = 1, size(lam)
         if (lam(k) >= LAM_LYMAN) then;  i = k;  exit;  end if
      end do
   end function i_lyman_of

   ! ===================================================================
   subroutine write_dl07_tables()
      ! The DL07 optics are Mie on the D03 dielectric functions plus the DL07
      ! PAH cross sections, evaluated by the very routines sed_init_dl07 calls,
      ! so these tables are that model's own optics written down.
      call dl07_settings()
      call sed_init_dl07(F_QT, F_SD, SD_INDEX, U_ISRF, NT_IN, T_LO, T_HI)
      call dl07_dump('', .false.)
      call dl07_settings()
      ! The floor comes from the model, not from a number written here: the
      ! extinction curve of the same model is built on it too, and the two must
      ! land on the same wavelength grid to share one axis in the HDF5 file.
      call sed_init_dl07(F_QT_EUV, F_SD, SD_INDEX, U_ISRF, NT_IN, T_LO, T_HI, &
                         lam_min=dl07_euv_lambda_floor())
      ! The HDF5 file carries the WIDER of the two axes and the index where the
      ! non-ionizing part of it begins, so the narrow product is a slice of it
      ! rather than a second copy to keep in step.  Hence only this pass writes.
      call dl07_dump('_euv', .true.)
   end subroutine write_dl07_tables


   subroutine dl07_settings()
      ! The 2003 DL07 model, exactly as main_dl07.x sets it up.
      qpah_use_d03_graphite  = .true.
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
         call h5_model_file('dl07', lam, i_lyman_of(lam))
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
   subroutine write_zubko_tables()
      ! Recomputed on the ZDA tables' own grids, by the SAME routines the DL07
      ! model uses -- Mie on the Draine (2003) optical constants this tree
      ! already carries.  Two things follow from that choice.  The D03 files
      ! span 6.2e-5 to 1.24e5 um for the silicate and 1.24e6 for the graphite,
      ! so the ZDA grid is covered outright and nothing is extrapolated; and
      ! q_graphite_full rebuilds the graphite index for every radius, so the
      ! size-dependent free-electron term is carried rather than frozen at one
      ! grain size.  D03 is the 2003 revision of the same astronomical silicate
      ! and graphite ZDA cited (Draine & Lee 1984, Laor & Draine 1993,
      ! Weingartner & Draine 2001), and it is evidently what the shipped tables
      ! were computed on: this route reproduces their silicate to a mean
      ! 2.4e-6, whatever the label in their header says.
      real(wp), allocatable :: a_sil(:), a_gra(:), a_pah(:), lz(:), qa2(:,:), qs2(:,:)
      real(wp) :: r_sil, r_gra, r_pah
      integer  :: n_sil, n_gra, n_pah, nw, nwc, i

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
      ! the DL07 pass had run first in the same process -- ./make_qtable.x all
      ! and ./make_qtable.x zubko wrote PAH tables differing by a mean 17%.
      call dl07_settings()

      ! The ZDA tables are read for their AXES and densities only -- the radii
      ! each component is tabulated on and the shared wavelength grid.  Their Q
      ! columns are what this program recomputes.
      call read_zubko_optics(D_ZUBKO//'suvSil_121_1201.dat', n_sil, nw, a_sil, lz, &
                             qa2, qs2, r_sil)
      deallocate(qa2, qs2)
      call read_zubko_optics(D_ZUBKO//'Gra_121_1201.dat', n_gra, nw, a_gra, lz, &
                             qa2, qs2, r_gra)
      deallocate(qa2, qs2)
      call read_zubko_optics(D_ZUBKO//'PAH_28_1201_neu.dat', n_pah, nw, a_pah, lz, &
                             qa2, qs2, r_pah)
      deallocate(qa2, qs2)

      ! The non-ionizing set is the same grid cut at the Lyman limit, which is
      ! also where build_zubko's lam_min cuts it, so the stored pair matches
      ! the pair a run can build.
      nwc = 0
      do i = 1, nw
         if (lz(i) >= LAM_LYMAN) nwc = nwc + 1
      end do

      ! Only the wide pass writes HDF5; the narrow text product is the same
      ! rows and the file records where they start.  See write_dl07_tables.
      call zubko_dump('_euv', nw,  lz(nw-nw+1:nw),  a_sil, a_gra, a_pah, &
                      r_sil, r_gra, r_pah, .true.)
      call zubko_dump('',    nwc, lz(nw-nwc+1:nw), a_sil, a_gra, a_pah, &
                      r_sil, r_gra, r_pah, .false.)
   end subroutine write_zubko_tables


   subroutine zubko_dump(tag, nw, lz, a_sil, a_gra, a_pah, r_sil, r_gra, r_pah, to_h5)
      character(len=*), intent(in) :: tag
      integer,  intent(in) :: nw
      real(wp), intent(in) :: lz(:), a_sil(:), a_gra(:), a_pah(:)
      real(wp), intent(in) :: r_sil, r_gra, r_pah
      logical,  intent(in) :: to_h5
      real(wp) :: Qe, Qs, Qa, g, Qn
      integer  :: u1, u2, u3, jw, ja
      real(wp), allocatable :: se(:,:), sa(:,:), ss(:,:), sg(:,:)
      real(wp), allocatable :: ge(:,:), ga(:,:), gs(:,:), gg(:,:)
      real(wp), allocatable :: pa(:,:)

      allocate(se(nw,size(a_sil)), sa(nw,size(a_sil)), ss(nw,size(a_sil)), sg(nw,size(a_sil)))
      allocate(ge(nw,size(a_gra)), ga(nw,size(a_gra)), gs(nw,size(a_gra)), gg(nw,size(a_gra)))
      allocate(pa(nw,size(a_pah)))

      call open_table(u1, D_ZUBKO//'q_zubko_sil'//trim(tag)//'.dat', &
           'ZDA BARE-GR-S silicate: Mie on the D03 astrosilicate dielectric function', &
           nw, size(a_sil), .true., rho_in = r_sil)
      call open_table(u2, D_ZUBKO//'q_zubko_gra'//trim(tag)//'.dat', &
           'ZDA BARE-GR-S graphite: Mie on the D03 graphite dielectric functions, ' // &
           '1/3 E||c + 2/3 E-perp-c, free-electron term per radius', &
           nw, size(a_gra), .true., rho_in = r_gra)
      call open_table(u3, D_ZUBKO//'q_zubko_pah'//trim(tag)//'.dat', &
           'ZDA BARE-GR-S PAH, neutral: the DL07 absorption cross sections (not Mie)', &
           nw, size(a_pah), .false., rho_in = r_pah)

      do jw = 1, nw
         do ja = 1, size(a_sil)
            call q_silicate_full(a_sil(ja), lz(jw), Qe, Qs, Qa, g)
            call put_row(u1, lz(jw), a_sil(ja), Qe, Qa, Qs, g)
            se(jw,ja) = Qe;  sa(jw,ja) = Qa;  ss(jw,ja) = Qs;  sg(jw,ja) = g
         end do
         do ja = 1, size(a_gra)
            call q_graphite_full(a_gra(ja), lz(jw), Qe, Qs, Qa, g)
            call put_row(u2, lz(jw), a_gra(ja), Qe, Qa, Qs, g)
            ge(jw,ja) = Qe;  ga(jw,ja) = Qa;  gs(jw,ja) = Qs;  gg(jw,ja) = g
         end do
         do ja = 1, size(a_pah)
            call qpah_dl07(0, a_pah(ja), lz(jw), Qn)
            call put_abs_row(u3, lz(jw), a_pah(ja), Qn)
            pa(jw,ja) = Qn
         end do
      end do
      close(u1);  close(u2);  close(u3)
      write(*,'(a,i0,a)') ' zubko'//trim(tag)//': ', nw, ' wavelengths written'

      if (to_h5) then
         ! The densities are the ZDA tables' own, read from their headers, not
         ! the DL07 numbers -- these are that model's materials.
         call h5_model_file('zubko', lz(1:nw), i_lyman_of(lz(1:nw)))
         call h5_add_component('zubko', 'sil', a_sil, se, sa, ss, sg, r_sil, &
              'Mie on the D03 astrosilicate dielectric function', &
              F_SIL_D03//'; grids from '//D_ZUBKO//'suvSil_121_1201.dat')
         call h5_add_component('zubko', 'gra', a_gra, ge, ga, gs, gg, r_gra, &
              'Mie on the D03 graphite dielectric functions, 1/3 E||c + 2/3 ' // &
              'E-perp-c, free-electron term at each radius', &
              F_GRA_D03//'; grids from '//D_ZUBKO//'Gra_121_1201.dat')
         call h5_add_component('zubko', 'pah', a_pah, absonly = pa, rho = r_pah, &
              method = 'DL07 PAH absorption cross sections, neutral', &
              source = 'qpah.f90; grids from '//D_ZUBKO//'PAH_28_1201_neu.dat')
      end if
      deallocate(se, sa, ss, sg, ge, ga, gs, gg, pa)
   end subroutine zubko_dump



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
         write(*,'(a)') ' make_qtable: built without HDF5, writing text only'
         return
      end if
      call h5_begin(ok);  if (.not. ok) return
      call h5_create_file(h5_path(model), fid, ok)
      if (.not. ok) then
         write(*,'(a,a)') ' make_qtable: cannot create ', h5_path(model)
         call h5_end();  return
      end if
      call h5_put_attr_s(fid, 'model', model)
      call h5_put_attr_s(fid, 'generator', 'SEDust sed/make_qtable.x')
      call h5_put_attr_s(fid, 'layout', &
           'cross sections are (n_a, n_lam) in h5py; lambda is the fastest Fortran axis')
      call h5_group(fid, 'grid', gid, ok)
      call h5_write_1d(gid, 'lambda', lam, units='um', &
                       long_name='wavelength grid of the model, increasing')
      call h5_put_attr_i(gid, 'i_lyman', i_lyman)
      call h5_put_attr_d(gid, 'lambda_lyman_um', LAM_LYMAN)
      call h5_put_attr_s(gid, 'i_lyman_meaning', &
           'first index at or longward of the Lyman limit; include_euv=.false. starts here')
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
      ! absent scattering term for a computed zero.
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
         call h5_write_2d(gid, 'g',     gg, units='1', long_name='scattering asymmetry <cos>')
         call h5_put_attr_s(gid, 'scatters', 'yes')
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
      write(u,'(a)')       '# Written by sed/make_qtable.x from this tree''s own optics'
      write(u,'(a)')       '#   routines; see src/make_qtable.f90 for the grids and inputs.'
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
end program make_qtable
