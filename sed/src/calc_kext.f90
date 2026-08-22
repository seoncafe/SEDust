program calc_kext
   !====================================================================
   ! Size-integrated transport optics -- extinction, absorption and
   ! scattering cross section per H, albedo and scattering asymmetry <cos> --
   ! for any dust model this library can build.
   !
   !   ./calc_kext.x astrodust [euv | <lam_min in um>]
   !   ./calc_kext.x dl07 [ld01] [euv | <lam_min in um>]
   !   ./calc_kext.x zubko     [formula | table] [euv]
   !   ./calc_kext.x from_files <descriptor> [data_dir]
   !
   ! `ld01` asks for the DL07 model with the Li & Draine (2001) carbonaceous
   ! absorption in place of the Draine & Li (2007) one -- everything else the
   ! same.  Both are curves of the same model, so both land in ../data/dl07/,
   ! as kext_dl07_MW[_euv].dat and kext_ld01_MW[_euv].dat, and inside that
   ! model's own product as /kext and /kext_ld01.
   !
   ! Every model goes through the same two calls -- one build_* to set the
   ! model up, then size_integrated_extinction to do the size integral over its
   ! optics from first principles.  The tables written here are what an RT host
   ! then gets back from the library: dust_extinction reads these files rather
   ! than repeating the integral during a transport run.  What differs between
   ! models is only (a) which builder is called, (b) what sets the
   ! short-wavelength floor, and (c) the model-specific text in the file header.
   !
   ! The extreme-ultraviolet extension (`euv`, or an explicit floor) carries
   ! the wavelength grid below the T-matrix Q table's own short-wavelength end,
   ! for a host that transports ionizing radiation.  The floor is never
   ! hard-coded: it is asked of the dielectric function that will supply the
   ! optics there, because past its own tabulation the refractive index would
   ! freeze at the boundary value.  The Zubko model needs no such extension and
   ! could not have one -- its optics table IS the model definition and already
   ! reaches 0.001 um (1.24 keV) -- so for it `euv` means the opposite: it keeps
   ! the whole tabulated range, while the default cuts the grid at the Lyman
   ! limit, where an interstellar radiation field stops carrying photons.
   !
   ! MASS ABSORPTION COEFFICIENT.  Every product carries K_abs [cm^2/g] =
   ! C_abs/H divided by the model's dust mass per H, which the model itself
   ! reports (dust_mass_per_H) from its own size distribution and the solid
   ! density of each grain population.  It is one wavelength-independent
   ! constant, so it is as well defined in the EUV as in the optical, and the
   ! header of each file states its value and where the densities came from.
   !
   ! POLARIZED EXTINCTION.  ../data/astrodust/kext_astrodust_MW.dat carries the dichroic
   ! extinction as an eighth column, because that file is tracked and hosts
   ! read it; it comes from the same size_integrated_extinction call, as
   ! Cpol_ext.  (No column of it is read back by dust_extinction, which keeps
   ! computing the polarized optics from the model so that they follow the
   ! host's runtime alignment state.)  Both astrodust products carry it, on the
   ! 1129-wavelength grid and on the 1762-wavelength one; no other model has
   ! polarized optics, so no other product has the column.
   !
   ! The orientation-resolved table the column comes from stops at 0.0912 um,
   ! and the companion that would carry the band below it (run_q_jori.x euv)
   ! is not shipped, so on the _euv grid build_Cpol leaves the 633 rows below
   ! that wavelength at zero and announces the deficit on stderr at every run.
   ! That zero is an omission, not a physical value.  Letting a reader mistake
   ! it for a measurement would be exactly the laundering these products
   ! refuse, so header_polarized_column states the wavelengths and the row
   ! count of that run, taken at run time from the polarized table's own
   ! coverage.  On the non-EUV grid there is no such run: the model grid and
   ! the polarized table are the same 1129 wavelengths.
   !====================================================================
   use constants,         only: wp
   use sed_astrodust_mod, only: build_astrodust, build_dl07, build_mrn, build_zubko, &
                                build_from_files, size_integrated_extinction, &
                                dust_mass_per_H, dust_model_t, &
                                d03_euv_lambda_floor, astrodust_euv_lambda_floor
   use q_astrodust_mod,   only: astrodust_index_lambda_range, get_astrodust_index_path
   use q_silicate_mod,    only: silicate_index_lambda_range
   use q_graphite_mod,    only: graphite_index_lambda_range
   ! Coverage of the orientation-resolved polarized table, read after the model
   ! is built so the header can state which rows of the eighth column are zero
   ! for want of it.  Same basis build_Cpol works from: the model grid is that
   ! table's own wavelength axis, possibly extended shortward, and nj_lam_euv
   ! is non-zero only when an EUV companion table filled the extension.
   use q_table_jori_mod,  only: nj_lam, lam_j, nj_lam_euv
   ! Should a requested floor fall below the astrodust Q table, that band is
   ! solved on the oblate spheroid of the table rather than on a
   ! volume-equivalent sphere.  The calculation is a separate module so that
   ! the library links without the T-matrix; this driver registers it because
   ! `./calc_kext.x astrodust euv` is what could ask for such a band.  The Q
   ! table shipped with this release already reaches into the ionizing band,
   ! so with it nothing is left for this route to solve.
   use euv_astrodust_tmatrix, only: use_tmatrix_euv_band_optics
   use sedust_h5
   ! One command-line system for every program here: this one declares the
   ! wavelength-grid axis, the DL07 model's carbonaceous vintage, and (for
   ! zubko) that model's size-distribution and optics-set axes, and nothing
   ! else.  A word of any other axis -- a solver, a radiation field, an
   ! emission term -- is refused by name, because a size-integrated transport
   ! curve has no referent for one.  Nor is the graphite of the PAH xi blend
   ! offered: these products ARE the model definition.
   use sed_run_options,   only: run_options_t, declare_run_options, &
                                widen_run_options, read_run_subject, &
                                read_run_option, check_run_options, &
                                write_run_option_usage
   implicit none

   ! ---- model inputs shared with the SED drivers -------------------------
   ! The plain product preserves the historical non-ionizing wavelength grid.
   ! The _euv product prepends DH21 nodes below the Lyman limit and is used
   ! only when the caller explicitly asks for the EUV table.
   ! Each model's own product carries ONE wavelength axis and the index where
   ! the non-ionizing part of it begins, so `euv` selects a view of the same
   ! file rather than a second one.
   character(len=*), parameter :: F_QT   = '../data/astrodust/sedust_astrodust.h5'
   character(len=*), parameter :: F_QT_DL = '../data/dl07/sedust_dl07.h5'
   character(len=*), parameter :: F_SD  = '../data/release/size_distribution.dat'
   character(len=*), parameter :: F_EXT = '../data/release/extinction.dat'
   character(len=*), parameter :: F_SCA = '../data/release/scattering.dat'
   character(len=*), parameter :: F_REF_DL07 = &
      '../data/release/kext_albedo_WD_MW_3.1_60_D03.all_2003'
   ! Draine's published curve for the MRN model, the size integral of the very
   ! power law build_mrn evaluates: his parameter file for it states the grain
   ! volumes per H, 2.49e-27 and 2.79e-27 cm^3/H, which are
   ! log10 A = -25.16 and -25.11 over 0.005 - 0.25 um.
   character(len=*), parameter :: F_REF_MRN = '../data/release/kext_albedo_MRN'
   character(len=*), parameter :: F_QT_MRN   = '../data/mrn/sedust_mrn.h5'
   character(len=*), parameter :: F_QT_ZU    = '../data/zubko/sedust_zubko.h5'
   character(len=*), parameter :: F_ZDA_CFG  = '../data/zubko/ZDA_BARE_GR_S_Config.dat'
   character(len=*), parameter :: F_ZDA_DESC = '../data/zubko/zubko_descriptor.txt'
   character(len=*), parameter :: D_ZUBKO    = '../data/zubko/'

   ! DL07 model: WD01 MW R_V = 3.1, b_C = 6e-5 (Draine's "60" model), MMP83 field.
   integer,  parameter :: SD_INDEX_DL07 = 7
   real(wp), parameter :: U_ISRF_DL07   = 1.0_wp
   ! Thermal-table grid.  It sets H(T,a) / kappB, which the extinction never
   ! touches; the builders need it all the same.
   integer,  parameter :: NT_IN = 100
   real(wp), parameter :: T_LO = 1.0_wp, T_HI = 3000.0_wp
   ! The Lyman limit, 13.6 eV.  It is where the non-ionizing products stop,
   ! because an interstellar radiation field is exactly zero below it.
   real(wp), parameter :: LAM_LYMAN = 0.0912_wp

   integer, parameter :: MAXNOTE = 120
   character(len=200) :: note(MAXNOTE)
   integer            :: nnote

   type(dust_model_t)    :: m
   real(wp), allocatable :: Cext(:), Cabs(:), Csca(:), alb(:), gbar(:)
   ! Dichroic (polarized) extinction per H.  Only the astrodust model fills
   ! it; every other builder leaves it zero, which is why only the astrodust
   ! product writes it out.
   real(wp), allocatable :: Cpolext(:)
   real(wp) :: lam_min, lam_lo, lam_hi, lam_lo2, Mdust_H
   integer  :: narg, ios, nlam_out, iarg
   logical  :: euv, zubko_formula
   ! Which published PAH absorption cross section the DL07-model curve is the
   ! size integral of, and which model's data directory and HDF5 product the
   ! curve belongs in.  The two are not the same thing: 'dl07' and 'ld01' are
   ! two vintages of ONE model, so both land in ../data/dl07/ and in that
   ! model's own product, distinguished by the file name and the group.
   character(len=8)   :: pah_xsec
   character(len=16)  :: product
   character(len=32)  :: model
   type(run_options_t) :: o
   logical            :: taken
   integer            :: n_free
   character(len=16)  :: zubko_optics
   character(len=256) :: opt, fout, desc, ddir, arg

   call use_tmatrix_euv_band_optics()

   ! ---- command line -----------------------------------------------------
   ! The axes this program has a referent for, declared before any argument is
   ! read; the DL07 model and zubko add their own once the subject is known.
   call declare_run_options(o, program='calc_kext', &
        subjects=[character(len=16):: 'astrodust', 'dl07', 'mrn', 'zubko', &
                                      'from_files'], &
        grid=.true.)

   narg = command_argument_count()
   if (narg < 1) then
      call print_usage();  stop
   end if
   call get_command_argument(1, model)

   call read_run_subject(o, trim(model), taken)
   if (.not. taken) then
      write(*,'(a,a)') ' calc_kext: unknown model ', trim(model)
      call print_usage();  stop 1
   end if
   if (trim(model) == 'zubko') &
      call widen_run_options(o, zubko_sizedist=.true., zubko_optics=.true.)
   ! The DL07 model is defined with the DL07 (2007) carbonaceous absorption and
   ! can be built with the earlier LD01 (2001) one instead; both are curves of
   ! the same model, filed beside each other.
   if (trim(model) == 'dl07') call widen_run_options(o, pah_xsec=.true.)

   euv = .false.;  lam_min = 0.0_wp
   pah_xsec = 'dl07';  product = trim(model)
   opt = '';  desc = F_ZDA_DESC;  ddir = D_ZUBKO
   ! Every word after the model, through the shared parser.  What it does not
   ! recognize is this program's own positional argument: a wavelength floor
   ! for astrodust and the DL07 model, a descriptor path and a data directory
   ! for from_files.
   n_free = 0
   do iarg = 2, narg
      call get_command_argument(iarg, arg)
      call read_run_option(trim(arg), o, taken)
      if (.not. taken) then
         n_free = n_free + 1
         select case (n_free)
         case (1);  opt  = trim(arg)
         case (2);  ddir = trim(arg)
         case default
            write(*,'(a,a)') ' calc_kext: unexpected argument ', trim(arg)
            call print_usage();  stop 1
         end select
      end if
   end do
   call check_run_options(o)
   euv = o%euv
   ! Which of the two stored zubko optics sets this curve is the size integral
   ! of.  The default matches build_dust's, so the shipped curve is the one a
   ! host that names nothing is served.
   zubko_formula = o%zubko_formula
   zubko_optics  = o%zubko_optics
   pah_xsec      = o%pah_xsec

   nnote = 0
   Mdust_H = 0.0_wp

   select case (trim(model))

   ! ===================================================================
   case ('astrodust')
      ! A bare number is a wavelength floor in um, and asking for one is asking
      ! for the band below the Q table's own short end.
      if (len_trim(opt) > 0) then
         euv = .true.
         read(opt, *, iostat=ios) lam_min
         if (ios /= 0 .or. lam_min <= 0.0_wp) then;  call print_usage();  stop 1;  end if
      end if
      if (euv) then
         ! The requested floor is tested against the DH21 astrodust dielectric
         ! function's own short-wavelength end: whatever the Q table does not
         ! already cover would have to be computed from that function, and past
         ! its tabulation (n, k) would freeze at the boundary value.
         call astrodust_index_lambda_range(lam_lo, lam_hi)
         if (lam_min <= 0.0_wp) lam_min = astrodust_euv_lambda_floor()
         write(*,'(a,es12.5,a,es12.5,a)') ' DH21 astrodust dielectric function covers ', &
            lam_lo, ' - ', lam_hi, ' um'
         write(*,'(a,es12.5,a)') ' requested lam_min = ', lam_min, ' um'
         call build_astrodust(m, F_QT, F_SD, NT_IN, T_LO, T_HI, lam_min=lam_min, &
                              include_euv=.true.)
         fout = '../data/astrodust/kext_astrodust_MW_euv.dat'
      else
         call build_astrodust(m, F_QT, F_SD, NT_IN, T_LO, T_HI)
         fout = '../data/astrodust/kext_astrodust_MW.dat'
      end if

   ! ===================================================================
   case ('dl07')
      ! One model, two published carbonaceous absorption vintages.  Everything
      ! else -- the WD01 size distribution, the D03 silicate Mie optics, the D03
      ! graphite scattering, the charge mixing -- is shared, so the pair
      ! isolates that one cross section.  `ld01` is the vintage Draine's 2003
      ! kext_albedo table was computed with, which is why the comparison this
      ! program prints at the end is the sharper test for it.
      product = 'dl07'
      if (len_trim(opt) > 0) then
         euv = .true.
         read(opt, *, iostat=ios) lam_min
         if (ios /= 0 .or. lam_min <= 0.0_wp) then;  call print_usage();  stop 1;  end if
      end if
      if (euv) then
         ! This model's optics are Mie on the D03 dielectric functions at every
         ! wavelength, and it needs BOTH the silicate and the graphite table, so
         ! the floor is the longer of the two tables' short-wavelength ends.
         call silicate_index_lambda_range(lam_lo, lam_hi)
         call graphite_index_lambda_range(lam_lo2, lam_hi)
         write(*,'(a,es12.5,a)') ' D03 silicate dielectric function starts at ', lam_lo,  ' um'
         write(*,'(a,es12.5,a)') ' D03 graphite dielectric function starts at ', lam_lo2, ' um'
         lam_lo = max(lam_lo, lam_lo2)
         if (lam_min <= 0.0_wp) lam_min = d03_euv_lambda_floor()
         write(*,'(a,es12.5,a)') ' requested lam_min = ', lam_min, ' um'
         call build_dl07(m, F_QT_DL, F_SD, SD_INDEX_DL07, U_ISRF_DL07, NT_IN, T_LO, T_HI, &
                         lam_min=lam_min, include_euv=.true., pah_xsec=pah_xsec)
         fout = '../data/dl07/kext_'//trim(pah_xsec)//'_MW_euv.dat'
      else
         call build_dl07(m, F_QT_DL, F_SD, SD_INDEX_DL07, U_ISRF_DL07, NT_IN, T_LO, T_HI, &
                         pah_xsec=pah_xsec)
         fout = '../data/dl07/kext_'//trim(pah_xsec)//'_MW.dat'
      end if

   ! ===================================================================
   case ('mrn')
      ! Graphite + silicate spheres on one a^-3.5 power law each.  Its optics
      ! are Mie on the D03 dielectric functions at every wavelength, exactly as
      ! the DL07 model's are, so `euv` extends the grid down to the longer of
      ! the two functions' short-wavelength ends and nothing is extrapolated.
      product = 'mrn'
      if (len_trim(opt) > 0) then
         euv = .true.
         read(opt, *, iostat=ios) lam_min
         if (ios /= 0 .or. lam_min <= 0.0_wp) then;  call print_usage();  stop 1;  end if
      end if
      if (euv) then
         call silicate_index_lambda_range(lam_lo, lam_hi)
         call graphite_index_lambda_range(lam_lo2, lam_hi)
         write(*,'(a,es12.5,a)') ' D03 silicate dielectric function starts at ', lam_lo,  ' um'
         write(*,'(a,es12.5,a)') ' D03 graphite dielectric function starts at ', lam_lo2, ' um'
         lam_lo = max(lam_lo, lam_lo2)
         if (lam_min <= 0.0_wp) lam_min = d03_euv_lambda_floor()
         write(*,'(a,es12.5,a)') ' requested lam_min = ', lam_min, ' um'
         call build_mrn(m, F_QT_MRN, NT_IN, T_LO, T_HI, lam_min=lam_min, &
                        include_euv=.true.)
         fout = '../data/mrn/kext_mrn_euv.dat'
      else
         call build_mrn(m, F_QT_MRN, NT_IN, T_LO, T_HI, include_euv=.false.)
         fout = '../data/mrn/kext_mrn.dat'
      end if

   ! ===================================================================
   case ('zubko')
      ! Two routes to the same ZDA BARE-GR-S model, differing only in where the
      ! size distribution comes from; the optics tables are the same files.
      ! Orthogonal to that, `euv` keeps the whole ZDA range and its absence cuts
      ! the grid at the Lyman limit, exactly as for the other two models -- with
      ! the difference that here the wide product is the tables' own range and
      ! the narrow one is that range cut, rather than the other way round.
      ! The optics come from the model's own product, so that /kext and
      ! /qtable of one file are the same numbers.  Left to the text tables this
      ! took the axis from the product and the optics from their seven written
      ! digits, and filed the two together as though they matched.
      if (euv) then
         if (zubko_formula) then
            call build_zubko(m, F_ZDA_CFG, D_ZUBKO, NT_IN, T_LO, T_HI, &
                             include_euv=.true., qtable_path=F_QT_ZU, &
                             optics=trim(zubko_optics))
         else
            call build_from_files(m, F_ZDA_DESC, D_ZUBKO, NT_IN, T_LO, T_HI, &
                                  include_euv=.true.)
         end if
         fout = '../data/zubko/kext_zubko_BARE_GR_S'//trim(kext_tag())//'_euv.dat'
      else
         ! The narrow product is an INDEX cut at the Lyman limit, as it is for
         ! the other two models -- not a wavelength floor.  The ZDA grid has no
         ! node at 0.0912 um, so a floor there used to start this curve 1.2e-5
         ! of itself inside the limit and a host on its own Lyman-limit grid
         ! was refused it.
         if (zubko_formula) then
            call build_zubko(m, F_ZDA_CFG, D_ZUBKO, NT_IN, T_LO, T_HI, &
                             include_euv=.false., qtable_path=F_QT_ZU, &
                             optics=trim(zubko_optics))
         else
            call build_from_files(m, F_ZDA_DESC, D_ZUBKO, NT_IN, T_LO, T_HI, &
                                  include_euv=.false.)
         end if
         fout = '../data/zubko/kext_zubko_BARE_GR_S'//trim(kext_tag())//'.dat'
      end if

   ! ===================================================================
   case ('from_files')
      if (len_trim(opt) == 0) then;  call print_usage();  stop 1;  end if
      desc = opt
      if (n_free < 2) ddir = dirname_of(desc)
      call build_from_files(m, trim(desc), trim(ddir), NT_IN, T_LO, T_HI)
      fout = '../data/' // trim(m%name) // '/kext_' // trim(m%name) // '.dat'

   case default
      call print_usage();  stop 1
   end select

   ! ---- size integral, identical call for every model --------------------
   nlam_out = m%NLAM
   allocate(Cext(nlam_out), Cabs(nlam_out), Csca(nlam_out), alb(nlam_out), gbar(nlam_out), &
            Cpolext(nlam_out))
   call size_integrated_extinction(m, Cext, Cabs, Csca, gbar=gbar, Cpol_ext=Cpolext, &
                                   albedo=alb)
   write(*,'(a,i0,a,es12.5,a,es12.5,a)') ' size-integrated ', nlam_out, &
      ' wavelengths, ', m%lam(1), ' - ', m%lam(nlam_out), ' um'

   ! ---- dust mass per H, the K_abs normalization, also model-independent --
   ! It comes from the model object: every builder states the solid density of
   ! each population's material, so this is the same size distribution the size
   ! integral above just used, weighted by those densities.
   Mdust_H = dust_mass_per_H(m)
   write(*,'(a,es12.5,a)') ' M_dust/H = ', Mdust_H, ' g/H'

   ! ---- header text, the one part that is model-specific -----------------
   select case (trim(model))
   case ('astrodust')
      if (euv) then
         call header_astrodust_euv()
      else
         call header_astrodust_qtable_grid()
      end if
   case ('dl07');       call header_dl07()
   case ('mrn');        call header_mrn()
   case ('zubko');      call header_zubko()
   case ('from_files'); call header_from_files()
   end select

   ! ---- write ------------------------------------------------------------
   ! Every product carries K_abs, because every model states its grain densities
   ! and so has a dust mass per H, and the header of each records the
   ! M_dust/N_H the column is normalized by. The normalization is one
   ! wavelength-independent constant, as well defined in the EUV as in the
   ! optical. The tracked ../data/astrodust/kext_astrodust_MW.dat is a regression
   ! reference, so it alone keeps the original narrow-field format; every other
   ! product uses the wider default fields. The dichroic C_polext column goes on
   ! BOTH astrodust products, because astrodust is the only model with polarized
   ! optics and the quantity is defined over its whole grid -- where the
   ! orientation-resolved table does not reach, the header states the deficit
   ! rather than leaving the reader to infer it from a missing column.
   if (Mdust_H > 0.0_wp) then
      call write_kext_curve(trim(fout), kabs_norm=Mdust_H, &
                            legacy_format=(trim(model) == 'astrodust' .and. .not. euv), &
                            with_polext_column=(trim(model) == 'astrodust'))
   else
      call write_kext_curve(trim(fout))
   end if
   write(*,'(a,a)') ' wrote ', trim(fout)

   ! ---- the same curve into the model's HDF5 file ------------------------
   ! Only the WIDE run of a coded model writes: data/<model>/sedust_<model>.h5 holds
   ! one wavelength axis, and a reader asked for the non-ionizing part slices
   ! it at /grid/i_lyman.  The narrow text product is that slice, so writing it
   ! here as well would be a second copy of the same numbers on a shorter axis.
   if (euv) call write_kext_h5(trim(product))

   ! ---- checks -----------------------------------------------------------
   call check_internal_consistency()
   select case (trim(model))
   case ('astrodust');  call compare_hd23_release()
   case ('dl07');       call compare_draine_kext(F_REF_DL07, &
                             'Draine kext_albedo_WD_MW_3.1_60_D03.all (2003)')
   case ('mrn');        call compare_draine_kext(F_REF_MRN, &
                             'Draine kext_albedo_MRN')
   end select

contains

   function kext_tag() result(t)
      ! The suffix a curve carries when it is NOT the model's default optics:
      ! the second stored optics set for zubko, the earlier cross-section
      ! vintage for the DL07 model.  Each model's own choice is unmarked, so
      ! the shipped file names and the /kext group do not move.
      character(len=16) :: t
      t = ''
      if (trim(model) == 'zubko' .and. trim(zubko_optics) /= 'zda') &
         t = '_'//trim(zubko_optics)
      if (trim(product) == 'dl07' .and. trim(pah_xsec) /= 'dl07') &
         t = '_'//trim(pah_xsec)
   end function kext_tag


   function h5_kext_group() result(g)
      ! Where that curve goes inside the product: /kext for the default set,
      ! /kext_<optics> for another.
      character(len=32) :: g
      g = 'kext'//trim(kext_tag())
   end function h5_kext_group


   ! ===================================================================
   subroutine print_usage()
      write(*,'(a)') ' usage:'
      write(*,'(a)') '   ./calc_kext.x astrodust [euv | <lam_min in um>]'
      write(*,'(a)') '   ./calc_kext.x dl07 [ld01] [euv | <lam_min in um>]'
      write(*,'(a)') '   ./calc_kext.x mrn       [euv | <lam_min in um>]'
      write(*,'(a)') '   ./calc_kext.x zubko     [formula | table] [euv]'
      write(*,'(a)') '   ./calc_kext.x from_files <descriptor> [data_dir]'
      write(*,'(a)') ''
      write(*,'(a)') ' Writes lambda / albedo / <cos> / C_ext per H (+ C_abs, C_sca, and'
      write(*,'(a)') ' K_abs = C_abs/H normalized by the model dust mass per H) under'
      write(*,'(a)') ' ../data/.  `euv` asks for the ionizing band: for astrodust and'
      write(*,'(a)') ' dl07 it extends the grid down to whatever the model dielectric'
      write(*,'(a)') ' function itself covers; for zubko, whose own tables already reach'
      write(*,'(a)') ' 1.24 keV, it keeps the band that the default cuts at the Lyman limit.'
      write(*,'(a)') ' `dl07 ld01` is the DL07 model with the Li & Draine (2001) carbonaceous'
      write(*,'(a)') ' absorption in place of the Draine & Li (2007) one -- everything else'
      write(*,'(a)') ' the same -- and lands beside it as ../data/dl07/kext_ld01_MW[_euv].dat.'
      write(*,'(a)') ''
      call write_run_option_usage(o)
      write(*,'(a)') ' kext_astrodust_MW.dat alone also carries the dichroic C_polext/H.'
   end subroutine print_usage


   function dirname_of(p) result(d)
      ! Directory part of a path, with its trailing separator; '.' if there is none.
      character(len=*), intent(in) :: p
      character(len=256) :: d
      integer :: k
      k = index(trim(p), '/', back=.true.)
      if (k > 0) then
         d = p(1:k)
      else
         d = './'
      end if
   end function dirname_of


   subroutine add_note(s)
      character(len=*), intent(in) :: s
      if (nnote >= MAXNOTE) then
         write(*,'(a)') ' calc_kext: header note array too small';  stop 1
      end if
      nnote = nnote + 1
      note(nnote) = s
   end subroutine add_note


   ! ===================================================================
   ! Shared writer.  Columns are
   !   lambda[um]  albedo  <cos>  C_ext/H[cm^2/H]  C_abs/H[cm^2/H]  C_sca/H[cm^2/H]
   ! -- the first four in the order of Draine's kext_albedo tables and of
   ! MoCHII's par%ion_dust_kext reader, which takes four reals from each row
   ! and uses lambda, albedo and C_ext.  A trailing K_abs [cm^2/g] column is
   ! added when the model supplies its dust mass per H, and the frozen format
   ! adds C_polext/H after it.
   subroutine write_kext_curve(path, kabs_norm, legacy_format, with_polext_column)
      character(len=*),   intent(in) :: path
      ! Dust mass per H [g/H]; when present a K_abs = C_abs/H / this column is added.
      real(wp), optional, intent(in) :: kabs_norm
      ! The frozen on-disk field widths of ../data/astrodust/kext_astrodust_MW.dat, kept
      ! so that tracked file is reproduced byte for byte.  It is also the only
      ! product carrying the dichroic C_polext/H column.  Everything else uses
      ! the wider default fields.
      logical,  optional, intent(in) :: legacy_format
      ! Write the dichroic C_polext/H as an eighth column.  Every astrodust
      ! product carries it: the quantity is defined over the whole grid, and
      ! where the orientation-resolved table does not reach, the header says so
      ! in as many words rather than leaving the reader to infer it from a
      ! missing column.
      logical,  optional, intent(in) :: with_polext_column
      integer :: u, jw
      logical :: legacy, with_kabs, with_polext
      legacy    = .false.;  if (present(legacy_format)) legacy = legacy_format
      with_polext = legacy
      if (present(with_polext_column)) with_polext = with_polext_column
      with_kabs = present(kabs_norm)
      if (legacy .and. .not. with_kabs) then
         write(*,'(a)') ' calc_kext: the frozen format carries K_abs and so needs' // &
                        ' the dust mass per H'
         stop 1
      end if

      open(newunit=u, file=path, status='replace', action='write')
      do jw = 1, nnote
         write(u,'(a)') trim(note(jw))
      end do
      if (legacy) then
         ! The wavelength column carries 8 significant digits, not the 6 this
         ! format historically used.  The astrodust grid resolves each X-ray
         ! absorption edge with a pair of points ~1.2e-6 apart in relative
         ! wavelength, which 6 digits round to the same string: the column then
         ! repeats a value, and load_kext_table rejects the file for not being
         ! strictly ascending -- the model builds but has no extinction to
         ! serve.  8 digits separate the closest pair by 18 units of the last
         ! digit (7 would leave under 2).  The header is written through the
         ! same field widths so the two cannot drift apart.
         write(u,'(a1,a14,2(1x,a10),5(1x,a15))') '#', 'lambda', 'albedo', '<cos>', &
              'C_ext/H', 'C_abs/H', 'C_sca/H', 'K_abs', 'C_polext/H'
         write(u,'(a1,a14,2(1x,a10),5(1x,a15))') '#', '(micron)', ' ', ' ', &
              '(cm^2/H)', '(cm^2/H)', '(cm^2/H)', '(cm^2/g)', '(cm^2/H)'
         do jw = 1, nlam_out
            write(u,'(es15.7e3,2(1x,f10.6),5(1x,es15.7e3))') &
               m%lam(jw), alb(jw), gbar(jw), Cext(jw), Cabs(jw), Csca(jw), &
               Cabs(jw)/kabs_norm, Cpolext(jw)
         end do
      else if (with_kabs .and. with_polext) then
         write(u,'(a1,a19,7(1x,a20))') '#', 'lambda', 'albedo', '<cos>', &
              'C_ext/H', 'C_abs/H', 'C_sca/H', 'K_abs', 'C_polext/H'
         write(u,'(a1,a19,7(1x,a20))') '#', '(micron)', ' ', ' ', &
              '(cm^2/H)', '(cm^2/H)', '(cm^2/H)', '(cm^2/g)', '(cm^2/H)'
         do jw = 1, nlam_out
            write(u,'(es20.12e3,7(1x,es20.12e3))') &
               m%lam(jw), alb(jw), gbar(jw), Cext(jw), Cabs(jw), Csca(jw), &
               Cabs(jw)/kabs_norm, Cpolext(jw)
         end do
      else if (with_kabs) then
         write(u,'(a1,a19,6(1x,a20))') '#', 'lambda', 'albedo', '<cos>', &
              'C_ext/H', 'C_abs/H', 'C_sca/H', 'K_abs'
         write(u,'(a1,a19,6(1x,a20))') '#', '(micron)', ' ', ' ', &
              '(cm^2/H)', '(cm^2/H)', '(cm^2/H)', '(cm^2/g)'
         do jw = 1, nlam_out
            write(u,'(es20.12e3,6(1x,es20.12e3))') &
               m%lam(jw), alb(jw), gbar(jw), Cext(jw), Cabs(jw), Csca(jw), Cabs(jw)/kabs_norm
         end do
      else
         write(u,'(a1,a19,5(1x,a20))') '#', 'lambda', 'albedo', '<cos>', &
              'C_ext/H', 'C_abs/H', 'C_sca/H'
         write(u,'(a1,a19,5(1x,a20))') '#', '(micron)', ' ', ' ', &
              '(cm^2/H)', '(cm^2/H)', '(cm^2/H)'
         do jw = 1, nlam_out
            write(u,'(es20.12e3,5(1x,es20.12e3))') &
               m%lam(jw), alb(jw), gbar(jw), Cext(jw), Cabs(jw), Csca(jw)
         end do
      end if
      close(u)
   end subroutine write_kext_curve


   ! ===================================================================
   subroutine write_kext_h5(model)
      ! The MODEL whose product this is -- 'dl07' for both vintages of the DL07
      ! model -- not the word the command line used.  Which curve inside that
      ! product is kext_tag()'s business.
      !
      ! Append the size-integrated curve to ../data/<model>/sedust_<model>.h5, beside
      ! the (lambda, a_eff) tables calc_qtable.x wrote into the same file.
      !
      ! The file holds ONE wavelength axis and every wavelength-indexed array
      ! is filed against it, so this refuses to write unless the grid the size
      ! integral just ran on IS that axis.  A curve stored under a wavelength
      ! it was not computed at is precisely the stale label this format carries
      ! provenance attributes to prevent, and the two grids are built by two
      ! programs, so the check is worth its cost.
      character(len=*), intent(in) :: model
      character(len=256) :: path, sdfile
      integer(h5id_k)    :: fid, gid
      logical  :: ok
      real(wp), allocatable :: gl(:), kabs(:)
      real(wp) :: dmax
      integer  :: jw

      if (.not. sedust_has_hdf5) then
         write(*,'(a)') ' calc_kext: built without HDF5, wrote the text product only'
         return
      end if
      select case (model)
      case ('astrodust', 'dl07')
         sdfile = F_SD
      case ('mrn')
         ! This model's size distribution is a formula with no file behind it,
         ! so the provenance attribute records the power law rather than a path.
         sdfile = 'dn/da = A a^-3.5, 0.005-0.25 um, log10 A = -25.16 (gra),'// &
                  ' -25.11 (sil) [Draine & Lee 1984]'
      case ('zubko')
         ! Two routes to the same model: the ZDA size-distribution formula from
         ! the config, or the tabulated dn/da the descriptor names.
         if (zubko_formula) then
            sdfile = F_ZDA_CFG
         else
            sdfile = desc
         end if
      case default
         return                      ! from_files: no shipped file to extend
      end select

      path = '../data/'//model//'/sedust_'//model//'.h5'
      call h5_begin(ok);  if (.not. ok) return
      call h5_open_rw(trim(path), fid, ok)
      if (.not. ok) then
         write(*,'(a,a)') ' calc_kext: no HDF5 file to extend at ', trim(path)
         write(*,'(a,a)') '            write it first with  ./calc_qtable.x ', model
         call h5_end();  return
      end if

      call h5_read_1d(fid, 'grid/lambda', gl, ok)
      if (.not. ok) then
         write(*,'(a,a)') ' calc_kext: no /grid/lambda in ', trim(path)
         call h5_close_file(fid);  call h5_end();  return
      end if
      if (size(gl) /= nlam_out) then
         write(*,'(a,i0,a,i0,a)') ' calc_kext: /grid/lambda has ', size(gl), &
            ' points and this curve has ', nlam_out, ' -- not written'
         deallocate(gl);  call h5_close_file(fid);  call h5_end();  return
      end if
      dmax = 0.0_wp
      do jw = 1, nlam_out
         dmax = max(dmax, abs(m%lam(jw)/gl(jw) - 1.0_wp))
      end do
      deallocate(gl)
      ! Both grids come from the same construction in the same library, so they
      ! agree to the last bit or they are not the same grid at all.
      if (dmax > 1.0e-12_wp) then
         write(*,'(a,es9.2,a)') ' calc_kext: this curve''s grid differs from ' // &
            '/grid/lambda by ', dmax, ' -- not written'
         call h5_close_file(fid);  call h5_end();  return
      end if

      ! /kext for the model as it is built by default; a second set of optics
      ! gets its own group, so the curve a host is served is the size integral
      ! of the very optics its model was built on and not of another set.
      if (h5_has(fid, trim(h5_kext_group()))) call h5_unlink(fid, trim(h5_kext_group()))
      call h5_group(fid, trim(h5_kext_group()), gid, ok)
      if (.not. ok) then
         write(*,'(a,a)') ' calc_kext: cannot create /kext in ', trim(path)
         call h5_close_file(fid);  call h5_end();  return
      end if
      call h5_write_1d(gid, 'albedo', alb,  units='1', &
                       long_name='scattering albedo C_sca/C_ext')
      call h5_write_1d(gid, 'g',      gbar, units='1', &
                       long_name='scattering asymmetry <cos>')
      call h5_write_1d(gid, 'C_ext',  Cext, units='cm^2/H', &
                       long_name='extinction cross section per H nucleon')
      call h5_write_1d(gid, 'C_abs',  Cabs, units='cm^2/H', &
                       long_name='absorption cross section per H nucleon')
      call h5_write_1d(gid, 'C_sca',  Csca, units='cm^2/H', &
                       long_name='scattering cross section per H nucleon')
      ! The dichroic term, which only the astrodust model has: the extinction
      ! difference an aligned grain presents to the two linear polarizations.
      ! It is written for that model alone rather than as a column of zeros,
      ! so that a reader can tell a model that has no polarized optics from one
      ! whose polarization happens to vanish.
      if (model == 'astrodust') then
         call h5_write_1d(gid, 'C_polext', Cpolext, units='cm^2/H', &
              long_name='polarized (dichroic) extinction cross section per H nucleon')
         call h5_put_attr_d(gid, 'pol_valid_from', pol_valid_from())
      end if
      if (Mdust_H > 0.0_wp) then
         allocate(kabs(nlam_out))
         kabs = Cabs / Mdust_H
         call h5_write_1d(gid, 'K_abs', kabs, units='cm^2/g', &
                          long_name='absorption mass opacity, C_abs/H / (M_dust/N_H)')
         deallocate(kabs)
         call h5_put_attr_d(gid, 'M_dust_per_H', Mdust_H)
      end if
      call h5_put_attr_s(gid, 'size_dist_file', trim(sdfile))
      call h5_put_attr_s(gid, 'text_product', trim(fout))
      call h5_put_attr_s(gid, 'generator', 'SEDust sed/calc_kext.x')
      call h5_put_attr_s(gid, 'method', &
           'size integral over the built model''s populations, ' // &
           'the same call dust_extinction serves from')
      call h5_group_close(gid)
      call h5_close_file(fid)
      call h5_end()
      write(*,'(a,a)') ' wrote /kext into ', trim(path)
   end subroutine write_kext_h5


   real(wp) function pol_valid_from() result(l0)
      ! Shortest wavelength at which C_polext is a computed value rather than
      ! an omitted zero.  Below the orientation-resolved table's own first node
      ! there is nothing to read, and build_Cpol leaves the rows at zero -- a
      ! b/a = 1.4 spheroid has a nonzero dichroic extinction there, so that zero
      ! is an omission and a reader must be able to tell.  An EUV companion
      ! table fills the band, and then the whole grid is valid.
      if (nj_lam > 0 .and. nj_lam_euv == 0) then
         l0 = lam_j(1)
      else
         l0 = m%lam(1)
      end if
   end function pol_valid_from



   subroutine header_common_format_and_scope(density_source)
      ! Format and scope statements every new product carries, ending with the
      ! K_abs normalization and where that model's grain densities come from.
      ! (The tracked ../data/astrodust/kext_astrodust_MW.dat predates these lines and is
      ! frozen, so header_astrodust_qtable_grid does not call this routine and
      ! that file keeps its own wording.)
      character(len=*), intent(in) :: density_source
      character(len=200) :: s
      call add_note('#')
      call add_note('# Column order follows Draine''s kext_albedo tables and MoCHII''s')
      call add_note('#   par%ion_dust_kext reader: the first four columns are')
      call add_note('#   lambda, albedo, <cos>, C_ext/H, and a reader may stop there.')
      call add_note('#   Draine''s 5th column is K_abs [cm^2/g]; ours is C_abs/H [cm^2/H],')
      call add_note('#   so no dust-mass normalization is folded into the first six')
      call add_note('#   columns; K_abs is a separate trailing column instead.')
      call add_note('#   A host running the Zubko model through an ionizing band needs')
      call add_note('#   exactly such a file: that model has no EUV extension of its own')
      call add_note('#   (its optics table IS the model definition), so the host reads the')
      call add_note('#   ionizing-band optics from a curve file instead.')
      call add_note('#')
      call add_note('# SCOPE: these are size-integrated TRANSPORT optics, not a dust model.')
      call add_note('#   There is no size-resolved Cabs(lambda, a), no grain enthalpy and no')
      call add_note('#   Planck integral here, so stochastic heating and thermal emission')
      call add_note('#   cannot be computed from this file.  A host that needs emission must')
      call add_note('#   build the model (build_astrodust / build_dl07 / build_zubko /')
      call add_note('#   build_from_files).  The builder reads this very curve back in, and')
      call add_note('#   dust_extinction serves it on the model wavelength grid while')
      call add_note('#   dust_emission supplies the re-emitted spectrum, so the absorbed')
      call add_note('#   power and the emission refer to one and the same grain.')
      call add_note('#')
      call add_note('# Units: lambda micron; cross sections cm^2 per H nucleon;')
      call add_note('#   albedo and <cos> dimensionless.')
      if (Mdust_H > 0.0_wp) then
         call add_note('#')
         write(s,'(a,es13.6,a)') '# K_abs [cm^2/g] = C_abs/H / (M_dust/N_H), ' // &
              'M_dust/N_H = ', Mdust_H, ' g/H.'
         call add_note(trim(s))
         call add_note('#   The dust mass per H is this model''s own size distribution')
         call add_note('#   weighted by the solid density of each grain material:')
         call add_note('#   M_dust/N_H = sum_pop rho * sum_a (4/3) pi a^3 dn(a).')
         call add_note('#   ' // density_source)
      end if
   end subroutine header_common_format_and_scope


   subroutine header_wavelength_range(what_sets_the_floor)
      character(len=*), intent(in) :: what_sets_the_floor
      character(len=200) :: s
      call add_note('#')
      write(s,'(a,es13.6,a,es13.6,a,i0,a)') '# Wavelength range: ', m%lam(1), ' - ', &
           m%lam(nlam_out), ' um, ', nlam_out, ' points.'
      call add_note(trim(s))
      call add_note('#   ' // what_sets_the_floor)
   end subroutine header_wavelength_range


   ! ===================================================================
   subroutine header_astrodust_qtable_grid()
      ! FROZEN down to the dust-mass line.  Those lines reproduce
      ! ../data/astrodust/kext_astrodust_MW.dat as it has been tracked since before this
      ! file carried a polarized column; do not reword them.  The C_polext/H
      ! block below them describes the eighth column and is not frozen.
      character(len=200) :: s
      call add_note('# Extinction, albedo, and scattering asymmetry for the')
      call add_note('# Hensley & Draine (2023) astrodust+PAH Milky-Way model.')
      call add_note('#')
      call add_note('# Astrodust grains: random-orientation T-matrix optics')
      call add_note('#   (Draine & Hensley 2021 dielectric, P=0.20, fFe=0, b/a=1.4);')
      call add_note('#   Cabs, Csca and <cos> from the T-matrix Q table.')
      call add_note('# PAH: charge-resolved DL07 absorption (neutral+cation by f_ion,')
      call add_note('#   100 A -> 100% ionized); scattering from the graphite fraction')
      call add_note('#   xi_gra(a) of HD23 eq. 15, which is all of the mixture that')
      call add_note('#   scatters (3.5% of tau_sca at 0.1 um, nothing past 0.3 um).')
      call add_note('# Size integral over the HD23 release size distribution (per H).')
      call add_note('# Computed from the SED pipeline optics (sed_init), all components.')
      write(s,'(a,es13.6,a)') '# Dust mass per H, M_dust/N_H = ', Mdust_H, &
           ' g/H  (rho_Ad=2.74, rho_PAH=2.0 g/cm^3; K_abs = C_abs/H / this).'
      call add_note(trim(s))
      call add_note('#')
      call header_polarized_column()
      call add_note('#')
   end subroutine header_astrodust_qtable_grid


   subroutine header_polarized_column()
      ! The C_polext/H column, described the same way in both astrodust
      ! products.  Whether it ends in a DEFICIT paragraph depends on the grid:
      ! the model grid is the polarized table's own wavelength axis, optionally
      ! carried shortward, and the rows below that table's first node are the
      ! ones build_Cpol cannot fill.  On the non-EUV grid there are none; on the
      ! EUV grid there are 633, unless an EUV companion table (nj_lam_euv > 0)
      ! filled them.
      character(len=200) :: s
      integer :: npol0
      call add_note('# Column 8, C_polext/H, is the dichroic (polarized) extinction')
      call add_note('#   sum_a dn_Ad(a) * C_pol,ext(lambda,a) * f_align(a), with')
      call add_note('#   C_pol,ext = 0.5*(Q_ext(E perp a) - Q_ext(E || a)) * pi a^2 at')
      call add_note('#   k perp a, from the orientation-resolved DH21 spheroid table.')
      call add_note('#   PAHs are unaligned (f_align = 0) and do not contribute.')
      call add_note('#   This is the MAXIMUM polarized extinction: it assumes the')
      call add_note('#   symmetry axis of every aligned grain lies in the plane of the')
      call add_note('#   sky, i.e. it is quoted BEFORE the sin^2(gamma) geometry factor')
      call add_note('#   and before any turbulent depolarization. A radiative transfer')
      call add_note('#   host must apply both. Codes that read only the first seven')
      call add_note('#   columns are unaffected.')
      npol0 = 0
      if (nj_lam > 0 .and. nj_lam_euv == 0) npol0 = max(nlam_out - nj_lam, 0)
      if (npol0 > 0) then
         call add_note('#')
         write(s,'(a,i0,a)') '#   DEFICIT: C_polext/H is exactly zero in the first ', &
              npol0, ' rows of'
         call add_note(trim(s))
         write(s,'(a,es13.6,a,es13.6,a)') '#   this file, lambda ', m%lam(1), ' - ', &
              m%lam(npol0), ' um.  The orientation-resolved'
         call add_note(trim(s))
         write(s,'(a,es13.6,a)') '#   table starts at ', lam_j(1), &
              ' um and there is nothing below it to read,'
         call add_note(trim(s))
         call add_note('#   so that zero is an OMISSION, not physics: a b/a = 1.4 spheroid')
         call add_note('#   has a nonzero dichroic extinction there -- only a sphere gives')
         call add_note('#   exactly zero.  build_Cpol says the same on stderr at every run.')
         call add_note('#   To fill the band, compute the companion polarized table with')
         call add_note('#   tmatrix/driver/run_q_jori.x (`euv` for the whole band, or')
         call add_note('#   `lam L1 L2 ...` for selected wavelengths) and hand it to')
         call add_note('#   sed_init through qpol_euv_path / qpol_euv_wave_path.')
      end if
   end subroutine header_polarized_column


   subroutine header_astrodust_euv()
      character(len=200) :: s
      call add_note('# model = astrodust_MW_euv')
      call add_note('# Size-integrated extinction, albedo and scattering asymmetry per H')
      call add_note('# for the Hensley & Draine (2023) astrodust+PAH Milky-Way model,')
      call add_note('# carried into the extreme ultraviolet.')
      call add_note('#')
      call add_note('# Model definition: Hensley & Draine (2023), ApJ 948, 55')
      call add_note('#   (astrodust + PAH, Milky Way R_V = 3.1 fit).')
      call add_note('# Size distribution: ' // F_SD)
      call add_note('# Optics:')
      call add_note('#   astrodust : the random-orientation T-matrix Q table, read at EVERY')
      call add_note('#     wavelength of the range stated below, the ionizing band included:')
      call add_note('#     ' // F_QT)
      call add_note('#     (oblate spheroid b/a = 1.4, porosity P = 0.20, f_Fe = 0, on the')
      call add_note('#     Draine & Hensley (2021) dielectric function,')
      call add_note('#     ' // trim(get_astrodust_index_path()) // ').')
      call add_note('#     Where the size parameter leaves the range the T-matrix is solved')
      call add_note('#     in, the table itself carries the Rayleigh dipole limit (x < 0.1)')
      call add_note('#     or the geometric optics limit (x > 50).  No astrodust cross')
      call add_note('#     section is solved on the fly by this program, and no sphere is')
      call add_note('#     substituted for the spheroid at any wavelength.')
      call add_note('#   PAH : charge-resolved DL07 absorption (neutral + cation mixed by')
      call add_note('#     f_ion, a > 100 A fully ionized).  There is no PAH scattering')
      call add_note('#     cross section, so what scatters is the graphite fraction')
      call add_note('#     xi_gra(a) of HD23 eq. 15, on the D03 dielectric functions by')
      call add_note('#     Mie: 3.5% of tau_sca at 0.1 um and nothing past 0.3 um.')
      call header_wavelength_range('The grid starts at whichever is shorter, ' // &
           'the Q table''s own short end')
      call add_note('#   or the requested floor.  A floor below the DH21 astrodust')
      write(s,'(a,es13.6,a)') '#   dielectric table''s own end, ', lam_lo, &
           ' um, is refused, because (n, k)'
      call add_note(trim(s))
      call add_note('#   would freeze at the boundary value there.')
      call add_note('#')
      call add_note('# How far to trust the EUV band of THIS model:')
      call add_note('#   The astrodust optics are the b/a = 1.4 oblate spheroid at every')
      call add_note('#   wavelength, read from one table, so the grain does not change')
      call add_note('#   shape and the solution does not change character across the band.')
      call add_note('#   There is NO external reference for this band: the HD23 author')
      call add_note('#   release (extinction.dat, scattering.dat) stops at 12.4 eV.')
      call add_note('#   Above 21.4 eV (1/lambda = 17.25 um^-1) the DL07 PAH cross section')
      call add_note('#   is zero by construction, so the carbonaceous absorption there is')
      call add_note('#   entirely D03 graphite -- Mie on the Draine 2003 dielectric')
      call add_note('#   functions, random-orientation 1/3 parallel + 2/3 perpendicular.')
      call header_common_format_and_scope('Densities: rho_Ad = 2.74, ' // &
           'rho_PAH = 2.0 g/cm^3 (HD23 conventions).')
      call add_note('#')
      call header_polarized_column()
      call add_note('#')
   end subroutine header_astrodust_euv


   subroutine header_dl07()
      character(len=200) :: s
      if (euv) then
         call add_note('# model = '//trim(pah_xsec)//'_MW_euv')
      else
         call add_note('# model = '//trim(pah_xsec)//'_MW')
      end if
      call add_note('# Size-integrated extinction, albedo and scattering asymmetry per H')
      call add_note('# for the Draine & Li (2007) / Weingartner & Draine (2001) Milky-Way')
      call add_note('# dust model.')
      call add_note('#')
      call add_note('# Model definition: WD01 R_V = 3.1, b_C = 6e-5 size distribution')
      call add_note('#   (Draine''s "MW 3.1_60"), with the Draine (2003a) x0.93 abundance')
      call add_note('#   reduction applied, heated by the MMP83 field at U = 1 for the')
      call add_note('#   PAH ionization fraction (WD01b charging).')
      call add_note('# Size grid: Draine''s 84-point log grid, 3.548 A - 5.012 um.')
      call add_note('# Optics:')
      call add_note('#   silicate     : Mie on the D03 astrosilicate dielectric function')
      call add_note('#     ' // '../data/dielectric/index_silD03')
      if (trim(pah_xsec) == 'ld01') then
         call add_note('#   carbonaceous : absorption from the PAH <-> graphite xi-blend in')
         call add_note('#     the Li & Draine (2001) cross-section vintage (qpah_ld01),')
         call add_note('#     N_C = 468 (a/10A)^3, WITHOUT the PAH cation near-infrared')
         call add_note('#     additions Draine & Li (2007) eq. 2 makes after Mattioda et al.')
         call add_note('#     (2005) -- the 1.05 and 1.26 um resonances, the 1.905 um')
         call add_note('#     negative term and the NIR continuum.  This is the vintage the')
         call add_note('#     2003 reference table below was computed with.')
      else
         call add_note('#   carbonaceous : absorption from the DL07 PAH <-> graphite xi-blend')
      end if
      call add_note('#     (neutral and cation mixed by the ionization fraction, a > 100 A')
      call add_note('#     fully ionized); scattering and <cos> from random-orientation D03')
      call add_note('#     graphite (1/3 parallel + 2/3 perpendicular), which is the')
      call add_note('#     standard DL07 treatment -- PAH scattering is Rayleigh-negligible.')
      call add_note('#     ' // '../data/dielectric/index_CpaD03, index_CpeD03')
      ! Said on the LD01 side only, so that the DL07 products keep the header
      ! they were tracked with.
      if (trim(pah_xsec) == 'ld01') then
         call add_note('#   The two vintages differ in the carbonaceous ABSORPTION only, so')
         call add_note('#   this file and kext_dl07_MW*.dat share their C_sca and <cos>')
         call add_note('#   exactly and a ratio of the two isolates that one cross section.')
      end if
      if (euv) then
         call header_wavelength_range('The floor is set by the D03 dielectric ' // &
              'functions: the shorter-reaching')
         write(s,'(a,es13.6,a)') '#   of silicate and graphite stops at ', lam_lo, &
              ' um, and a shorter floor is'
         call add_note(trim(s))
         call add_note('#   refused because (n, k) would freeze at the boundary value there.')
      else
         call header_wavelength_range('The grid is the T-matrix Q table''s own: ' // &
              'no EUV extension was requested,')
         call add_note('#   but that grid by itself already covers the ionizing band.')
      end if
      call add_note('#')
      call add_note('# How far to trust the EUV band of THIS model:')
      call add_note('#   the optics are Mie on the D03 dielectric functions at EVERY')
      call add_note('#   wavelength, so nothing changes character across the band, and the')
      call add_note('#   result can be checked directly against Draine''s own published')
      call add_note('#   table kext_albedo_WD_MW_3.1_60_D03.all (2003 vintage), which')
      call add_note('#   spans 1e-4 - 1e4 um.  Running this program prints that comparison')
      call add_note('#   decade by decade.')
      call header_common_format_and_scope('Densities: 3.5 g/cm^3 astrosilicate, ' // &
           '2.2 g/cm^3 graphite (Draine & Li 2007 sec. 2).')
      call add_note('#')
   end subroutine header_dl07


   subroutine header_mrn()
      character(len=200) :: s
      if (euv) then
         call add_note('# model = mrn_euv')
      else
         call add_note('# model = mrn')
      end if
      call add_note('# Size-integrated extinction, albedo and scattering asymmetry per H')
      call add_note('# for the Mathis, Rumpl & Nordsieck (1977) graphite + silicate model.')
      call add_note('#')
      call add_note('# Model definition: one power law per material,')
      call add_note('#   dn_i/da = A_i n_H a^-3.5 over 0.005 um <= a <= 0.25 um, cut sharply')
      call add_note('#   at both ends (Draine & Lee 1984 eq. 5.1; the cutoffs are MRN''s own')
      call add_note('#   estimates, which DL84 held fixed).  There are no PAHs in it.')
      call add_note('#   Normalization: log10 A = -25.16 (graphite), -25.11 (silicate)')
      call add_note('#   in cm^2.5/H -- the abundances Draine & Lee (1984) sec. Va adopted')
      call add_note('#   after fitting the Savage & Mathis average extinction curve, and')
      call add_note('#   the pair Draine''s own kext_albedo_MRN is the size integral of.')
      call add_note('# Size grid: 70-point log grid, the two cutoffs on its ends.')
      call add_note('# Optics:')
      call add_note('#   silicate : Mie on the D03 astrosilicate dielectric function')
      call add_note('#     ' // '../data/dielectric/index_silD03')
      call add_note('#   graphite : Mie on the D03 graphite dielectric functions,')
      call add_note('#     random orientation (1/3 parallel + 2/3 perpendicular)')
      call add_note('#     ' // '../data/dielectric/index_CpaD03, index_CpeD03')
      call add_note('#   MRN themselves used Wickramasinghe''s optical constants; DL84')
      call add_note('#   recomputed this size distribution on theirs, and D03 is the')
      call add_note('#   current revision of those.')
      if (euv) then
         call header_wavelength_range('The floor is set by the D03 dielectric ' // &
              'functions: the shorter-reaching')
         write(s,'(a,es13.6,a)') '#   of silicate and graphite stops at ', lam_lo, &
              ' um, and a shorter floor is'
         call add_note(trim(s))
         call add_note('#   refused because (n, k) would freeze at the boundary value there.')
      else
         call header_wavelength_range('The grid is the T-matrix Q table''s own: ' // &
              'no EUV extension was requested,')
         call add_note('#   but that grid by itself already covers the ionizing band.')
      end if
      call header_common_format_and_scope('Densities: 3.5 g/cm^3 astrosilicate, ' // &
           '2.2 g/cm^3 graphite, the values the D03 optics and the DL01 heat ' // &
           'capacities are defined with.  Draine''s own MRN table normalizes ' // &
           'K_abs by 3.3 and 2.24 instead, so his K_abs is 3.1% larger than ' // &
           'this one while C_ext/H, which the size distribution fixes, agrees.')
      call add_note('#')
   end subroutine header_mrn


   subroutine header_zubko()
      if (euv) then
         call add_note('# model = zubko_BARE_GR_S_euv')
      else
         call add_note('# model = zubko_BARE_GR_S')
      end if
      call add_note('# Size-integrated extinction, albedo and scattering asymmetry per H')
      call add_note('# for the Zubko, Dwek & Arendt (2004) BARE-GR-S dust model.')
      call add_note('#')
      call add_note('# Model definition: ZDA 2004, ApJS 152, 211, bare grains, PAH +')
      call add_note('#   graphite + silicate; the files are the copies distributed with the')
      call add_note('#   Camps et al. (2015) radiative-transfer benchmark.')
      call add_note('# Optics (all three components): the ZDA optics tables')
      call add_note('#   ' // D_ZUBKO // 'PAH_28_1201_neu.dat, Gra_121_1201.dat, suvSil_121_1201.dat')
      call add_note('#   A HOMOGENEOUS-SPHERE Mie calculation: ZDA 2004 sec. 3 says so for')
      call add_note('#   the bare grains (that paper''s effective-medium step belongs to its')
      call add_note('#   COMPOSITE models), and each file names the dielectric function it')
      call add_note('#   came from -- Draine''s eps_Sil and eps_Gra.  The PAH file is not Mie:')
      call add_note('#   it names the Li & Draine (2001) / Draine & Li (2007) cross sections.')
      call add_note('#   Checked against Bohren-Huffman Mie on the published eps_Sil: at')
      call add_note('#   fixed wavelength the Q ratio is radius-independent to 4 digits over')
      call add_note('#   3.55e-4 - 100 um and <cos> agrees to ~1e-4, leaving a 0.5-2% infrared')
      call add_note('#   residual that lies in the refractive index (the published file is on')
      call add_note('#   a different wavelength grid).  So Q_sca and <cos> are read from the')
      call add_note('#   same file as Q_abs rather than recomputed.')
      if (zubko_formula) then
         call add_note('# Size distribution: the ZDA log-polynomial FORMULA evaluated on the')
         call add_note('#   optics-table radii, with the coefficients read from')
         call add_note('#   ' // F_ZDA_CFG)
      else
         call add_note('# Size distribution: the tabulated SzDist files interpolated onto the')
         call add_note('#   optics-table radii, as listed by')
         call add_note('#   ' // F_ZDA_DESC)
      end if
      if (euv) then
         call header_wavelength_range('The range is the ZDA optics table itself, ' // &
              '1.24 keV down to 1.24e-4 eV.')
         call add_note('#   This model needs no EUV EXTENSION and could not have one: its')
         call add_note('#   tabulated optics already cover the whole ionizing band, and those')
         call add_note('#   tables ARE the model definition, so there is no dielectric')
         call add_note('#   function to extend them with.')
      else
         call header_wavelength_range('The ZDA optics table cut at the Lyman limit, ' // &
              '13.6 eV down to 1.24e-4 eV.')
         call add_note('#   The tables reach 1.24 keV; the wavelengths shortward of')
         call add_note('#   0.0912 um are dropped, table row by table row, because an')
         call add_note('#   interstellar radiation field carries no photon there.  Nothing is')
         call add_note('#   recomputed: every row here is the row kext_zubko_BARE_GR_S_euv.dat')
         call add_note('#   carries.  Use that file for a host that transports ionizing')
         call add_note('#   radiation.')
      end if
      call header_common_format_and_scope('Densities: the Density entry in the ' // &
           'header of each component''s own ZDA optics table.')
      call add_note('#')
   end subroutine header_zubko


   subroutine header_from_files()
      call add_note('# model = ' // trim(m%name))
      call add_note('# Size-integrated extinction, albedo and scattering asymmetry per H')
      call add_note('# for a file-defined dust model.')
      call add_note('#')
      call add_note('# Descriptor: ' // trim(desc))
      call add_note('# Data directory: ' // trim(ddir))
      call add_note('# The descriptor names, for each population, its ZDA optics table')
      call add_note('# (optics), its dn/da table (size distribution) and its calorimetry')
      call add_note('# table; it is reproduced verbatim below.')
      call add_note('#')
      call echo_descriptor()
      call header_wavelength_range('The range is the first population''s optics ' // &
           'table, which every other population must share.')
      call header_common_format_and_scope('Densities: the rho field of each ' // &
           'pop: line, or the optics table''s own Density where that field is 0.')
      call add_note('#')
   end subroutine header_from_files


   subroutine echo_descriptor()
      ! Copy the descriptor into the header so the product records exactly which
      ! optics / size / calorimetry files it was made from.
      integer :: u, ios2
      character(len=190) :: line
      open(newunit=u, file=trim(desc), status='old', action='read', iostat=ios2)
      if (ios2 /= 0) return
      do
         read(u,'(a)',iostat=ios2) line;  if (ios2 /= 0) exit
         if (len_trim(line) == 0) cycle
         call add_note('#   | ' // trim(line))
      end do
      close(u)
   end subroutine echo_descriptor


   ! ===================================================================
   subroutine check_internal_consistency()
      ! C_ext = C_abs + C_sca is an identity of the size integral, albedo is a
      ! fraction and <cos> a direction cosine; report the worst violation of
      ! each rather than assert, so a bad optics table shows up as a number.
      real(wp) :: dmax, r
      integer  :: jw, jbad
      dmax = 0.0_wp;  jbad = 1
      do jw = 1, nlam_out
         if (Cext(jw) > 0.0_wp) then
            r = abs(Cext(jw) - (Cabs(jw) + Csca(jw))) / Cext(jw)
            if (r > dmax) then;  dmax = r;  jbad = jw;  end if
         end if
      end do
      write(*,'(a)') ' internal consistency:'
      write(*,'(a,es10.3,a,es12.5,a)') '   max |Cext - (Cabs+Csca)|/Cext = ', dmax, &
           '  at ', m%lam(jbad), ' um'
      write(*,'(a,f12.8,a,f12.8)') '   albedo range  ', minval(alb),  ' .. ', maxval(alb)
      write(*,'(a,f12.8,a,f12.8)') '   <cos>  range  ', minval(gbar), ' .. ', maxval(gbar)
      write(*,'(a,es12.5,a,es12.5)') '   C_ext/H range ', minval(Cext), ' .. ', maxval(Cext)
   end subroutine check_internal_consistency


   ! ===================================================================
   subroutine compare_hd23_release()
      ! C_ext/H and C_sca/H against the HD23 release total columns
      ! (extinction.dat / scattering.dat: lambda, tau_Ad, tau_PAH, tau_tot).
      real(wp), allocatable :: wr(:), ext_tot(:), sca_tot(:)
      real(wp) :: a(4)
      integer  :: ur, ios2, n, i
      character(len=512) :: line
      real(wp) :: ce, cs, dext, dsca
      real(wp), parameter :: bands(5) = (/0.15_wp, 0.55_wp, 2.2_wp, 12.0_wp, 100.0_wp/)
      n = 0
      open(newunit=ur, file=F_EXT, status='old', action='read', iostat=ios2)
      if (ios2 /= 0) then; write(*,'(a)') ' (HD23 release not found; skipping)'; return; end if
      do
         read(ur,'(a)',iostat=ios2) line; if (ios2 /= 0) exit
         line = adjustl(line); if (len_trim(line)==0 .or. line(1:1)=='#') cycle
         n = n + 1
      end do
      rewind(ur)
      allocate(wr(n), ext_tot(n), sca_tot(n))
      i = 0
      do
         read(ur,'(a)',iostat=ios2) line; if (ios2 /= 0) exit
         line = adjustl(line); if (len_trim(line)==0 .or. line(1:1)=='#') cycle
         i = i + 1; read(line,*) a; wr(i) = a(1); ext_tot(i) = a(4)
      end do
      close(ur)
      open(newunit=ur, file=F_SCA, status='old', action='read', iostat=ios2)
      i = 0
      do
         read(ur,'(a)',iostat=ios2) line; if (ios2 /= 0) exit
         line = adjustl(line); if (len_trim(line)==0 .or. line(1:1)=='#') cycle
         i = i + 1; read(line,*) a; sca_tot(i) = a(4)
      end do
      close(ur)
      write(*,'(a)') ' validation vs HD23 release (total Ad+PAH):'
      write(*,'(a)') '   lam[um]   C_ext ours/HD23   C_sca ours/HD23'
      do i = 1, 5
         ce   = loginterp(m%lam, Cext, nlam_out, bands(i))
         cs   = loginterp(m%lam, Csca, nlam_out, bands(i))
         dext = loginterp(wr, ext_tot, n, bands(i))
         dsca = loginterp(wr, sca_tot, n, bands(i))
         write(*,'(f10.3,2(7x,f10.4))') bands(i), ce/dext, cs/dsca
      end do
   end subroutine compare_hd23_release


   ! ===================================================================
   subroutine compare_draine_kext(path, label)
      ! Point-by-point against Draine's published table for the SAME model,
      ! whichever model that is: kext_albedo_WD_MW_3.1_60_D03.all (Dec 2003
      ! vintage, the one his IDL dust_cross() reads) for DL07, kext_albedo_MRN
      ! for MRN.  For DL07 a later 2009 recomputation revised the FIR opacity
      ! up ~12% and the 2175 A bump down ~16%; our D03 optics match the 2003
      ! table.  Reported decade by decade in lambda, on the reference's own
      ! wavelengths (interpolating OUR curve to them).
      !
      ! In the X-ray the reference brackets every absorption edge (Fe K, Si K,
      ! Mg K, Fe L, O K, C K) with wavelength pairs a few 1e-3 apart, far
      ! finer than this model's log grid (dln lam = 0.0116).  At such a point
      ! the comparison measures our GRID, not our optics, so the table reports
      ! two counts: all overlapping reference points, and only those where the
      ! reference is locally no denser than our grid ("resolvable").
      character(len=*), intent(in) :: path, label
      integer,  parameter :: MAXREF = 4000
      real(wp) :: lr(MAXREF), ar(MAXREF), gr(MAXREF), cr(MAXREF), kr(MAXREF), dum
      real(wp) :: ce, alb_i, g_i, dc, da, dg, dref, dmod
      real(wp) :: sc, sa, sg, mc, ma, mg, lam1
      integer  :: u, ios2, nref, k, id, cnt, cnt_all, idec, j
      real(wp) :: sc_all, mc_all
      integer,  parameter :: NDEC = 8
      real(wp) :: dec_lo(NDEC)
      character(len=512) :: line
      write(*,'(a,a,a)') ' validation vs ', trim(label), ':'
      open(newunit=u, file=path, status='old', action='read', iostat=ios2)
      if (ios2 /= 0) then; write(*,'(a)') '   (reference not found; skipping)'; return; end if
      ! Rows are  lambda albedo <cos> C_ext/H K_abs [<cos^2>] [+ a trailing
      ! "NNNN eV" comment].  The sixth column is Draine's later vintage and is
      ! read when it is there and not required when it is not -- the MRN table
      ! carries five -- so ONE reader serves both.  The prose header fails the
      ! numeric read either way and is skipped by it.
      nref = 0
      do
         read(u,'(a)',iostat=ios2) line
         if (ios2 /= 0) exit
         if (len_trim(line) == 0) cycle
         read(line,*,iostat=ios2) lam1, ar(nref+1), gr(nref+1), cr(nref+1), kr(nref+1), dum
         if (ios2 /= 0) &
            read(line,*,iostat=ios2) lam1, ar(nref+1), gr(nref+1), cr(nref+1), kr(nref+1)
         if (ios2 /= 0) cycle
         nref = nref + 1;  lr(nref) = lam1
         if (nref >= MAXREF) exit
      end do
      close(u)
      write(*,'(a,i0,a)') '   read ', nref, ' reference rows'
      if (nref < 2) return

      do id = 1, NDEC
         dec_lo(id) = 10.0_wp**(id - 5)          ! 1e-4 .. 1e3, each decade wide
      end do
      write(*,'(a)') '   decade [um]     N_all  <|dCext|>  max|dCext|    N_res' // &
                     '  <|dCext|>  max|dCext|   <|d alb|>    <|d g|>'
      do id = 1, NDEC
         cnt = 0;  cnt_all = 0
         sc = 0.0_wp;  sa = 0.0_wp;  sg = 0.0_wp
         mc = 0.0_wp;  ma = 0.0_wp;  mg = 0.0_wp
         sc_all = 0.0_wp;  mc_all = 0.0_wp
         do k = 1, nref
            if (lr(k) < dec_lo(id) .or. lr(k) >= 10.0_wp*dec_lo(id)) cycle
            if (lr(k) < m%lam(1) .or. lr(k) > m%lam(nlam_out)) cycle
            if (cr(k) <= 0.0_wp) cycle
            ce    = loginterp(m%lam, Cext, nlam_out, lr(k))
            alb_i = lininterp(m%lam, alb,  nlam_out, lr(k))
            g_i   = lininterp(m%lam, gbar, nlam_out, lr(k))
            dc = abs(ce/cr(k) - 1.0_wp)
            da = abs(alb_i - ar(k))
            dg = abs(g_i   - gr(k))
            cnt_all = cnt_all + 1
            sc_all  = sc_all + dc;  mc_all = max(mc_all, dc)

            ! Is the reference locally coarser than our grid here?  Spacings are
            ! compared as |dln lambda| so that the test does not care which way
            ! the reference file runs (Draine's tables descend in lambda).
            dref = huge(1.0_wp)
            if (k > 1)    dref = min(dref, abs(log(lr(k)/lr(k-1))))
            if (k < nref) dref = min(dref, abs(log(lr(k+1)/lr(k))))
            j = bracket(m%lam, nlam_out, lr(k))
            dmod = log(m%lam(j)/m%lam(j-1))
            if (dref < dmod) cycle

            cnt = cnt + 1
            sc = sc + dc;  sa = sa + da;  sg = sg + dg
            mc = max(mc, dc);  ma = max(ma, da);  mg = max(mg, dg)
         end do
         if (cnt_all == 0) cycle
         idec = id - 5
         if (cnt == 0) cycle
         write(*,'(a,i3,a,i3,i8,2(1x,f10.6),i8,4(1x,f10.6))') '   1e', idec, ' - 1e', idec+1, &
            cnt_all, sc_all/real(cnt_all,wp), mc_all, &
            cnt, sc/real(cnt,wp), mc, sa/real(cnt,wp), sg/real(cnt,wp)
      end do
      write(*,'(a)') '   (N_res = reference points where the reference grid is no finer'
      write(*,'(a)') '    than ours; the rest straddle X-ray absorption edges our log grid'
      write(*,'(a)') '    cannot resolve, so there they measure the grid, not the optics.)'
   end subroutine compare_draine_kext


   ! ===================================================================
   function loginterp(x, y, n, x0) result(y0)
      ! Log-log interpolation: cross sections are positive and power-law-like.
      integer,  intent(in) :: n
      real(wp), intent(in) :: x(n), y(n), x0
      real(wp) :: y0
      integer  :: j
      j = bracket(x, n, x0)
      y0 = exp( log(max(y(j-1),tiny(0.0_wp))) + &
           (log(x0)-log(x(j-1)))/(log(x(j))-log(x(j-1))) * &
           (log(max(y(j),tiny(0.0_wp)))-log(max(y(j-1),tiny(0.0_wp)))) )
   end function loginterp

   function lininterp(x, y, n, x0) result(y0)
      ! Linear in log(lambda): albedo and <cos> are bounded, signed ratios,
      ! not power laws, so a log interpolant is neither defined nor right.
      integer,  intent(in) :: n
      real(wp), intent(in) :: x(n), y(n), x0
      real(wp) :: y0, t
      integer  :: j
      j = bracket(x, n, x0)
      t = (log(x0)-log(x(j-1))) / (log(x(j))-log(x(j-1)))
      y0 = (1.0_wp - t)*y(j-1) + t*y(j)
   end function lininterp

   function bracket(x, n, x0) result(j)
      ! Index j with x(j-1) <= x0 <= x(j) on an ascending grid.
      integer,  intent(in) :: n
      real(wp), intent(in) :: x(n), x0
      integer :: j, lo, hi, mid
      lo = 1;  hi = n
      do while (hi - lo > 1)
         mid = (lo + hi) / 2
         if (x(mid) <= x0) then;  lo = mid;  else;  hi = mid;  end if
      end do
      j = hi
   end function bracket

end program calc_kext
