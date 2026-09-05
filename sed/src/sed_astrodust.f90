module sed_astrodust_mod
   ! Astrodust SED solver. Two-stage API (init + solve) tailored for use
   ! inside a 3D radiative-transfer driver: dust optical and Planck-
   ! integral data are computed once at startup and reused per cell;
   ! only the local mean intensity J_lam varies per call.
   !
   ! API:
   !   call sed_init(qtable_path, NT_in, T_lo, T_hi)
   !   ...
   !   do icell = 1, ncells
   !      ... compute J_lam in this cell ...
   !      call sed_solve(J_lam, 'S1', lamI_lam)   ! or 'S2'
   !      ... use lamI_lam(:) ...
   !   end do
   !
   ! The driver `calc_sed.f90` calls this for the single
   ! Mathis-ISRF cell at U_mathis = 1.585 to compare with HD23.
   !
   ! Restructured into init/solve and
   ! limited to a single dust species (astrodust). The stochastic-vs-
   ! equilibrium decision and P(T) solver are unchanged in algorithm.

   use, intrinsic :: iso_fortran_env, only: real64, error_unit
   use constants,             only: wp
   use sed_mathlib,               only: interp, first_location, last_location
   use radfield,              only: bbody, calc_bbody, hardest_photon_energy, &
                                    cmb_temperature
   use p_sub,                 only: p_sub_setup, calc_Teq, calc_P
   use q_table_mod,           only: load_q_table, load_q_table_h5, &
                                    qt_n_lam=>n_lam, qt_n_aeff=>n_aeff, &
                                    qt_lam=>lam_t, qt_aeff=>aeff_t, &
                                    qt_qext=>qext, qt_qabs=>qabs, qt_qsca=>qsca, &
                                    qt_gpar=>gpar
   use q_table_jori_mod,      only: load_q_table_jori, falign_hd23, &
                                    qj_n_lam=>nj_lam, qj_lam=>lam_j, &
                                    qj_aeff=>aeff_j, qj_qpol_abs=>qpol_abs, &
                                    qj_qpol_ext=>qpol_ext, &
                                    qj_qbir_ext=>qbir_ext, qj_has_bir=>has_bir, &
                                    load_q_table_jori_euv, &
                                    qj_n_lam_euv=>nj_lam_euv, qj_lam_euv=>lam_j_euv, &
                                    qj_aeff_euv=>aeff_j_euv, &
                                    qj_qpol_abs_euv=>qpol_abs_euv, &
                                    qj_qpol_ext_euv=>qpol_ext_euv, &
                                    qj_qbir_ext_euv=>qbir_ext_euv, &
                                    qj_has_bir_euv=>has_bir_euv, &
                                    free_q_table_jori_euv
   use size_dist_mod,         only: hd23_size_distribution, sd_n=>n_size, &
                                    sd_aeff=>a_dist, sd_dn=>dn_ad, &
                                    sd_dn_pah=>dn_pah, sd_fion=>f_ion
   use enthalpy_astrodust_mod, only: enthalpy_S1, enthalpy_S2, RHO_AD, RHO_PAH
   use enthalpy,              only: enthalpy_DL01
   use qpah,                  only: qpah_dl07, qpah_ld01, qpah_abs, qpah_sca, &
                                    qpah_xsec_vintage, &
                                    nc_coeff, nc_integer, qpah_graphite_source
   use pah_ld01_mod,          only: q_pah_ld01, use_ld01_pah_xsec
   use stoch_qm_mod,          only: qm_solve_grain, qm_verbose
   ! DL07 (silicate + carbonaceous) model support
   use grain_dist_mod,        only: grain_dist_dl07, gd_apply_d03_reduction
   use q_silicate_mod,        only: q_silicate_abs, q_silicate_full, &
                                    silicate_index_lambda_range
   use q_graphite_mod,        only: q_graphite_full, &
                                    graphite_index_lambda_range
   ! Astrodust optics from the DH21 dielectric function, for the EUV band
   ! below the T-matrix Q table's 0.0912 um (13.6 eV) short-wavelength end.
   ! astrodust_index_at is the refractive index there, which the spheroid
   ! route (see euv_band_optics_i below) is handed; q_astrodust_full is the
   ! volume-equivalent-sphere Mie approximation, selected by
   ! euv_tmatrix = .false.
   use q_astrodust_mod,       only: astrodust_index_at, q_astrodust_full, &
                                    load_astrodust_index, &
                                    set_astrodust_index_path, &
                                    get_astrodust_index_path, &
                                    astrodust_index_lambda_range
   use pah_ioniz_mod,         only: pah_ionfrac
   use dust_model_mod,        only: dust_model_t, grain_pop_t, free_dust_model, &
                                    dust_set_alignment, dust_set_alignment_profile
   ! Size-integrated extinction tables (calc_kext.x products under data/),
   ! which the builders attach to the model and dust_extinction serves from.
   use kext_table_mod,        only: load_kext_table, tabulated_extinction_on_grid
   use scatmat_aligned_mod,   only: load_scatmat_aligned
   use q_component_mod, only: load_q_component
   use sed_paths,             only: sed_set_data_root, sed_get_data_root, sed_data_path, &
                                    SED_PATHLEN
   use sedust_product_mod,    only: sedust_dir, sedust_h5_file, read_sedust_grid, &
                                    lyman_index, &
                                    read_sedust_qtable, read_sedust_kext
   use zubko_io,              only: zda_comp_t, read_zda_config, zda_gofa, &
                                    read_zubko_optics, read_zubko_calor, &
                                    read_dnda_table, ZDA_MAXCOMP
   ! The DustEM file formats and the DustEM size distribution, for the two
   ! models this tree carries as DustEM definitions (THEMIS, G18 Model D).
   use dustem_io,             only: dustem_pop_t, DUSTEM_MAXPOP, &
                                    read_dustem_grain, dustem_size_distribution, &
                                    read_dustem_wavelengths, read_dustem_qtable, &
                                    read_dustem_gtable, read_dustem_heat_capacity, &
                                    optics_at_radii, grain_enthalpy_from_heat_capacity, &
                                    dustem_population_name
   implicit none
   private
   public :: sed_init, sed_solve, sed_solve_pah, sed_solve_qm_batch
   public :: sed_init_dl07, sed_solve_dl07
   ! Injection point for the spheroid (T-matrix) optics of the astrodust EUV
   ! band; see the abstract interface below.
   public :: euv_band_optics_i
   public :: sed_register_euv_band_optics, sed_forget_euv_band_optics
   ! Model-agnostic library API (path B: wraps the untouched solver core).
   ! One entry point for every coded model; the builders below stay as
   ! they are, so a caller that names its own files keeps working.
   public :: build_dust
   public :: dust_model_t, build_astrodust, build_dl07, build_zubko, dust_emission
   public :: build_mrn
   public :: build_from_files, dust_emission_single_teq, dust_extinction
   ! Models defined by DustEM input files: THEMIS and Guillet et al. (2018)
   ! Model D.  grain_pop_t comes with it so that a driver can ask a population
   ! whether it carries an asymmetry parameter.
   public :: build_dustem, grain_pop_t
   ! First-principles size integral over the model's own optics -- what the
   ! standalone calculators use, and what writes the tables dust_extinction
   ! then reads back.
   public :: size_integrated_extinction
   public :: dust_mass_per_H
   ! Solid densities of the two DL07 materials, exposed so that a program
   ! writing the model's optics products records the density the model itself
   ! integrates with rather than a second copy of the number.
   public :: RHO_ASTROSIL, RHO_GRAPHITE
   ! um -> cm, so that a program building a table on this module's grids uses
   ! the same conversion the module does rather than writing its own.
   public :: UM2CM
   ! Shortest wavelength each model's optics can be SOLVED at, and the stand-off
   ! factor behind them.  Every program that writes an EUV product of a model
   ! -- its Q tables, its extinction curve -- takes the floor from here, so
   ! that they land on ONE wavelength grid instead of on two that agree only to
   ! a few parts in 10^7 and cannot then share an axis on disk.
   public :: d03_euv_lambda_floor, astrodust_euv_lambda_floor, LAM_LO_MARGIN
   ! .true. iff the active model was built with polarized optics loaded; lets a
   ! host tell an intentionally scalar model from a polarized one.
   public :: dust_has_polarized_optics
   public :: dust_set_alignment, dust_set_alignment_profile
   public :: NLAM, NA, NT, lam, aeff, T_first, dn_ad, dn_pah, initialized
   ! Exposed so that external drivers can cross-check the optics:
   public :: Cabs, Csca, Cabs_pah, kappB_first, kappB_pah_first
   ! Polarized absorption cross section and alignment efficiency of the
   ! astrodust population (zero when the orientation-resolved table is absent).
   public :: Cpol, Cpol_ext, Cbir_ext, falign_ad, gsca_ad
   ! Charge-resolved PAH cross sections and number densities (neutral/cation),
   ! exposed so the MC SED builder can reproduce the same charge blend as
   ! the production sed_solve_pah (it loops both charge states).
   public :: Cabs_cneu, Cabs_cion, dn_cneu, dn_cion
   ! Exposed for mc/ single-grain P(T) cross-check (calc_P needs H, kappCMB):
   public :: H_first, H_pah_first, kappCMB, kappCMB_pah
   ! Toggle for the induced-emission (1 + J_lam/B_envelope) factor.
   ! Draine's emission kernel applies it, producing GROSS emission
   ! Cabs*B(T)*(1+J/B). But the
   ! published reference SEDs we compare to -- HD23 astrodust_irem.dat AND
   ! the DL07spec files -- are NET (gross minus the Cabs*J absorption),
   ! i.e. they do NOT carry this factor (verified 2026-05-17: enabling it
   ! pushed the >3000um band to +136% vs the HD23 release). The factor is
   ! (1 + n_gamma), n_gamma = J_lam*lam^5/(2 h c^2); negligible at FIR/
   ! sub-mm, ~0.7% at 1 mm, x2.4 at 1 cm (CMB occupation). Default .false.
   ! so our output is NET and matches the references. See [[induced-emission-factor]].
   public :: use_induced_emission
   ! Runtime stochastic-heating method selector. Values:
   !   'draine'    - grain-by-grain iterative T-window narrowing (Draine's method)
   !   'heuristic' - look-ahead narrowing (narrow_T_window)
   !   'qm'        - energy-space transition-matrix solver
   public :: stoch_method
   public :: gd_photon_cutoff
   ! Diagnostic-output toggle for the shared grain loop. Default .true. so the
   ! CLI drivers keep their solver diagnostics; dust_emission sets it from the
   ! model's `verbose` field so the library path stays silent by default.
   public :: sed_verbose

   real(wp), parameter :: PI    = 3.141592653589793238462643383279502884197d0
   real(wp), parameter :: UM2CM = 1.0e-4_wp

   ! Solid mass densities of the two DL07 / WD01 materials [g/cm^3], used to
   ! turn that model's size distribution into a dust mass per H.
   !
   ! Both are stated by the model's own paper, Draine & Li (2007) sec. 2, in the
   ! paragraph defining the effective radius a = (3M/4 pi rho)^(1/3):
   !
   !   astrosilicate  3.5   "amorphous silicate is assumed to have a mass
   !                        density rho = 3.5 g cm^-3".  The same value the DL01
   !                        astrosilicate enthalpy is built on (see
   !                        enthalpy_astrodust_mod, where the S1 stage rescales
   !                        N_atom by RHO_AD/3.5).
   !   graphite       2.2   "carbonaceous grains are assumed to have a mass
   !                        density due to graphitic carbon alone of
   !                        rho = 2.2 g cm^-3".  NOT the 2.24 of WD01 / LD01,
   !                        which is the density implied by the carbon-atom
   !                        count N_C = 470 (a/nm)^3 this builder selects.  DL07
   !                        kept that count and quoted the rounder density; the
   !                        model being built here is DL07's, so its number is
   !                        the one used.
   !
   ! The DL07 carbonaceous grains are ONE material sequence -- PAH-like at small
   ! a, graphite at large a, joined by the DL07 xi(a) weight -- with no size at
   ! which the solid changes, so the graphite density is used across the whole
   ! sequence rather than split at some radius.  (The astrodust model's PAHs take
   ! RHO_PAH = 2.0 instead.  That is HD23's own convention, tied to its
   ! N_C = 417 (a/nm)^3; the two numbers differ because the models differ, not
   ! because one corrects the other.)
   !
   ! CHECK against DL07 Table 3, model j_M = 7 (MW3.1_60, the size distribution
   ! this builder selects): the paper gives M_dust/M_H = 0.0104.  These densities
   ! reproduce it to about half a percent.  Draine's own kext_albedo_WD file
   ! header instead states M_dust/N_H = 1.870e-26 g/H, i.e. M_dust/M_H = 0.0112,
   ! which is 7% away from his Table 3; the header also normalizes the silicate
   ! distribution by 0.93 where the paper's Fig. 11 caption says 0.92.  The
   ! optics here follow the FILE (they agree with it to 0.1-0.3%), the mass
   ! follows the PAPER, and the two disagree by that ~1%.
   real(wp), parameter :: RHO_ASTROSIL = 3.5_wp
   real(wp), parameter :: RHO_GRAPHITE = 2.2_wp
   ! A model refuses a grid floor shorter than its own dielectric function's
   ! tabulation, because past it (n, k) would freeze at the boundary value.
   ! Asking for exactly that shortest wavelength puts the request on the
   ! rounding boundary of the refusal, so a floor is stood off it by this
   ! factor.
   real(wp), parameter :: LAM_LO_MARGIN = 1.001_wp
   ! The Lyman limit, 13.6 eV: where an interstellar radiation field stops, and
   ! the short-wavelength end of the non-ionizing products.
   real(wp), parameter :: LAM_LYMAN_UM = 0.0912_wp
   ! SI constants for the induced-emission factor (h*c^2 in J*m^2/s).
   real(wp), parameter :: H_SI    = 6.62606957e-34_wp
   real(wp), parameter :: C_SI    = 2.99792458e8_wp
   real(wp), parameter :: TWO_HCC = 2.0_wp * H_SI * C_SI**2

   ! ---- optics of the EUV band below the Q table ------------------------
   ! sed_init fills the wavelengths it prepends below the Q table's own
   ! short-wavelength end from one of two particles: the volume-equivalent
   ! SPHERE (Mie, q_astrodust_full, selected by euv_tmatrix = .false.), or the
   ! b/a = 1.400 oblate SPHEROID the table itself is made of. The spheroid is
   ! the continuation of the table -- same material, same shape, same
   ! random-orientation average -- and it is a T-matrix calculation.
   !
   ! That calculation is INJECTED rather than compiled in: this module holds a
   ! procedure pointer to it, so the SED library carries no reference to
   ! libtmatrix.a and links without it. euv_astrodust_tmatrix.f90 implements
   ! the interface below and registers it in one call; a build that leaves
   ! that file out has no spheroid route, and euv_tmatrix = .true. is then
   ! REFUSED (status 11) instead of being quietly answered with the sphere,
   ! which is a different particle.
   !
   ! The shipped Q table now reaches 1.0e-4 um (12398 eV), so a host whose
   ! lam_min lies on or above that end prepends no wavelengths at all and
   ! never enters this route.  It stays because lam_min is the host's to
   ! choose, and a band asked for below the table has to be answered with the
   ! table's own particle.
   abstract interface
      subroutine euv_band_optics_i(lam, a_um, nr, ki, qabs, qsca, gsca, status)
         !! Efficiencies of the astrodust grain at every (wavelength, radius)
         !! pair of the EUV band, normalized to pi a_eff^2 as the Q table is.
         !!   lam(n_euv)                   wavelengths [um], ascending, all
         !!                                below the Q table's own first point
         !!   a_um(NA)                     effective radii [um]
         !!   nr(n_euv), ki(n_euv)         refractive index m = nr + i*ki there
         !!   qabs/qsca/gsca(n_euv, NA)    absorption and scattering
         !!                                efficiencies and the scattering
         !!                                asymmetry <cos>
         !!   status   0  every point was produced and passed its physical
         !!               bounds
         !!         >  0  that many points failed a bound; the values are
         !!               still returned and the caller only warns
         !!         <  0  the arrays do not agree with the two grids, so
         !!               nothing was computed
         import :: wp
         real(wp), intent(in)  :: lam(:)
         real(wp), intent(in)  :: a_um(:)
         real(wp), intent(in)  :: nr(:), ki(:)
         real(wp), intent(out) :: qabs(:,:)
         real(wp), intent(out) :: qsca(:,:)
         real(wp), intent(out) :: gsca(:,:)
         integer,  intent(out) :: status
      end subroutine euv_band_optics_i
   end interface

   procedure(euv_band_optics_i), pointer :: euv_band_optics => null()

   ! Which build filled the module grids.  Bumped once per successful build and
   ! stamped into the model; dust_emission, which reads those grids through
   ! sed_grain_loop, compares against it.  See dust_model_t%build_id.
   integer, save :: active_build_id = 0
   logical :: initialized = .false.
   logical, save :: use_induced_emission = .false.

   ! True after the most recent astrodust build actually loaded the
   ! orientation-resolved polarized Q table (build_Cpol succeeded). False in an
   ! intentionally scalar build (load_polarized_optics = .false.), on a failed
   ! implicit-default load, and for models with no polarized optics (DL07,
   ! Zubko). build_astrodust reads it to decide whether the populations carry
   ! polarized optics; dust_has_polarized_optics then reports it off the model.
   logical, save :: polarized_optics_loaded = .false.
   ! dbdis-style photon-energy cutoff in the heuristic GD emission sum:
   ! a grain in a T bin cannot emit a photon more energetic than the bin
   ! enthalpy (hc/lambda <= H). Default OFF: on the coarse GD T-grid this
   ! over-suppresses the small-PAH FIR/submm by ~0.5-1% relative to dbdis
   ! (whose cutoff acts on its own fine adaptive energy bins); the window
   ! fixes alone reproduce dbdis to <1% per band. Kept as a toggle for
   ! experiments.
   logical, save :: gd_photon_cutoff = .false.
   real(wp), parameter :: HC_ERG_UM = 6.62606957e-27_wp * 2.99792458e10_wp * 1.0e4_wp
   ! Production default solver. 'heuristic' = the GD look-ahead-narrowing variant
   ! (faster); 'draine' = Draine's original GD (Guhathakurta & Draine 1989 +
   ! grain-by-grain narrowing); 'qm' = energy-space transition matrix. Selectable via the
   ! 'draine'/'qm' CLI toggles on the drivers.
   character(len=16), save :: stoch_method = 'heuristic'

   ! Guards the shared grain loop's solver diagnostics (see the public note).
   logical, save :: sed_verbose = .true.

   ! Module state set by sed_init
   integer  :: NLAM = 0, NA = 0, NT = 0
   ! Number of wavelength points the model grid carries into the extreme
   ! ultraviolet, below the T-matrix Q table's 0.0912 um (13.6 eV) end; 0 for
   ! the plain table grid. The table block is therefore lam(n_lam_euv+1:), and
   ! everything read off the Q tables -- the scalar Q's here, the polarized Q's
   ! in build_Cpol -- lands there. Every init path sets it, so it can never
   ! describe a grid other than the current one. The Planck integrals also read
   ! it, to anchor their internal ln(lambda) grid on the table range so that
   ! neither the step nor the sample points move as the grid widens -- see
   ! planck_integration_grid.
   integer, save :: n_lam_euv = 0
   real(wp), allocatable :: lam(:)              ! [um] (NLAM)
   real(wp), allocatable :: aeff(:)             ! [um] (NA), the size-dist grid
   real(wp), allocatable :: T_first(:)          ! [K]  (NT), log-spaced full range
   real(wp), allocatable :: dn_ad(:)            ! [1/H per bin] (NA)
   real(wp), allocatable :: Cabs(:,:)           ! [cm^2] (NLAM, NA)
   real(wp), allocatable :: Csca(:,:)           ! [cm^2] (NLAM, NA)
   ! Polarized absorption cross section of an aligned astrodust spheroid whose
   ! symmetry axis lies in the plane of the sky, and the HD23 alignment
   ! efficiency. Both are zero when the orientation-resolved table cannot be
   ! read, which makes the polarized emission vanish while everything else
   ! stays intact. See build_Cpol.
   ! NOTE ON NAMES: Cpol is the polarized ABSORPTION cross section, the one
   ! that drives the polarized emission. Cpol_ext is its EXTINCTION twin, the
   ! one that drives dichroic (polarized) extinction. They are built the same
   ! way from the same table and differ only in which Q the difference is
   ! taken from (Q_abs vs Q_ext). Cpol keeps its short name because the
   ! emission path (grain_pop_t%Cpol, set_pop, sed_grain_loop) is written
   ! around it.
   real(wp), allocatable :: Cpol(:,:)           ! [cm^2] (NLAM, NA)
   real(wp), allocatable :: Cpol_ext(:,:)       ! [cm^2] (NLAM, NA)
   ! Birefringent extinction cross section 0.5*(Cre_ext(3)-Cre_ext(2)) of the
   ! aligned astrodust spheroid, the UV-block optic of the extinction matrix.
   ! Built from the optional 4th (forward-amplitude real-part) block of the
   ! orientation-resolved table; stays zero when that block is absent. See
   ! build_Cpol.
   real(wp), allocatable :: Cbir_ext(:,:)       ! [cm^2] (NLAM, NA)
   real(wp), allocatable :: falign_ad(:)        ! (NA)
   ! Scattering asymmetry <cos> of the astrodust grains, taken from the same
   ! T-matrix Q table as Cabs/Csca and interpolated onto the size grid the same
   ! way. Only the extinction path (dust_extinction) uses it; the emission
   ! solver never needs it. Allocated by sed_init, dropped by sed_init_dl07.
   real(wp), allocatable :: gsca_ad(:,:)        ! (NLAM, NA)
   real(wp), allocatable :: kappB_first(:,:)    ! integral C_abs * B_lam dlam (NT, NA), wide grid
   real(wp), allocatable :: H_first(:,:,:)      ! enthalpy U(T, a, stage) (NT, NA, 2), wide grid
   ! CMB Planck integral (NA), at radfield's cmb_temperature() -- 2.725 K
   ! unless use_mathis_corrected is off, which restores Mathis (1983)'s 2.9 K.
   real(wp), allocatable :: kappCMB(:)

   ! Cached log copies for log-interpolation when narrowing T per grain
   real(wp), allocatable :: log_T_first(:)
   real(wp), allocatable :: log_H_first(:,:,:)
   real(wp), allocatable :: log_kappB_first(:,:)

   ! PAH-population data (set by sed_init alongside Astrodust). PAH cross
   ! sections from DL07 (qpah_dl07), mixed neutral + cation per the
   ! ionization fraction f_ion(a) of HD23 eq. (19) (size_dist_mod).
   ! Enthalpy uses DL01 'Car0' (carbonaceous) per HD23 §3.2.
   real(wp), allocatable :: dn_pah(:)             ! [1/H per bin] (NA)
   real(wp), allocatable :: Cabs_pah(:,:)         ! [cm^2] (NLAM, NA)
   real(wp), allocatable :: kappB_pah_first(:,:)  ! (NT, NA)
   real(wp), allocatable :: H_pah_first(:,:)      ! (NT, NA)  -- single 'Car0' enthalpy
   real(wp), allocatable :: kappCMB_pah(:)        ! (NA)
   real(wp), allocatable :: log_H_pah_first(:,:)
   real(wp), allocatable :: log_kappB_pah_first(:,:)

   ! DL07 carbonaceous charge states, solved as SEPARATE stochastically-heated
   ! populations (NOT pre-blended) because Teq/P(T)/emission are nonlinear in
   ! Cabs. Cabs_cneu/Cabs_cion are pure neutral/cation cross sections (cm^2);
   ! dn_cneu/dn_cion = full carbonaceous dn weighted by (1-fion)/fion. Both
   ! share the 'Car0' enthalpy (H_pah_first). The Cabs_pah/kappB_pah/kappCMB_pah
   ! arrays above are reused as charge-resolved scratch inside sed_solve_dl07.
   real(wp), allocatable :: Cabs_cneu(:,:), Cabs_cion(:,:)   ! [cm^2] (NLAM, NA)
   ! Scattering of the same carbonaceous grains.  ONE pair for both charge
   ! states: charge moves the PAH absorption features and nothing else, and
   ! what scatters is the graphite fraction xi_gra(a) of HD23 eq. 15, which
   ! does not know about charge either.
   real(wp), allocatable :: Csca_pah(:,:), gsca_pah(:,:)     ! [cm^2], - (NLAM, NA)
   ! Scattering of the DL07 carbonaceous grains.  Charge changes the PAH
   ! absorption features only, and a PAH is a molecule whose Rayleigh
   ! scattering is negligible, so the two charge states share one scattering
   ! description.  The xi_gra blend that splits the ABSORPTION between PAH and
   ! graphite therefore does not enter here: the graphite Q_sca carries the
   ! whole scattering -- the standard DL07 treatment.  Computed from
   ! the graphite dielectric function by Mie theory (q_graphite_full), on the
   ! same random-orientation average as the absorption, rather than
   ! interpolated off a precomputed Q table.
   real(wp), allocatable :: Csca_car(:,:), gsca_car(:,:)     ! [cm^2], - (NLAM, NA)
   real(wp), allocatable :: dn_cneu(:), dn_cion(:)           ! [1/H per bin] (NA)
   ! Charge-resolved kappB / kappCMB for the astrodust path (sed_init / sed_solve_pah /
   ! sed_solve_qm_batch), where both charge states must be available at once
   ! (the QM batch runs them in one parallel region, so they cannot share the
   ! Cabs_pah scratch the way sed_solve_dl07 does).
   real(wp), allocatable :: kappB_cneu(:,:), kappB_cion(:,:)       ! (NT, NA)
   real(wp), allocatable :: log_kappB_cneu(:,:), log_kappB_cion(:,:)
   real(wp), allocatable :: kappCMB_cneu(:), kappCMB_cion(:)       ! (NA)

   ! T-window narrowing: two algorithms available, switchable at runtime
   ! via the public `stoch_method` variable (see above). Both are kept
   ! in the source so they can be compared and so future debugging can
   ! fall back to either:
   !
   !   'draine'    -> grain-by-grain iterative refinement following Draine's
   !                  method
   !                  (initial guess from EEQ vs EEQSS, threshold PMIN=1e-13
   !                  tail trimming with 0.8/0.2 damping, 1.2x expansion).
   !   'heuristic' -> simple look-ahead heuristic: for each
   !                  grain set the next-grain window to
   !                  [Teq * exp(-NARROW_FAC*del), Teq * exp(+NARROW_FAC*del)]
   !                  with del = 2 * the previous grain's lnP > lnP_crit
   !                  half-width. Cheaper but does not threshold-target.
   !   'qm'        -> quantum-mechanical solver (not yet implemented).

   ! --- Draine iterative constants ---
   real(wp), parameter :: EV_TO_ERG      = 1.60218e-12_wp
   real(wp), parameter :: EEQSS_ERG      = 150.0_wp * EV_TO_ERG    ! steady-state threshold
   real(wp), parameter :: UMAXMIN_ERG    = 13.65_wp * EV_TO_ERG    ! single-UV-photon floor
   real(wp), parameter :: U_UV1_ERG      = 13.6_wp  * EV_TO_ERG    ! 1 hydrogen-ionizing photon
   real(wp), parameter :: HC_CGS_PER_CM  = 1.98645e-16_wp          ! erg per cm^-1 photon
   real(wp), parameter :: PMIN_LO        = 1.0e-13_wp              ! Draine v7 tail thresholds
   real(wp), parameter :: PMIN_UP        = 1.0e-13_wp
   integer,  parameter :: MAX_ITER_NARROW = 10

   ! --- Heuristic look-ahead constants ---
   real(wp), parameter :: P_crit         = 1.0e-15_wp
   real(wp), parameter :: dlnT_crit      = 0.5_wp
   real(wp), parameter :: lnP_crit       = log(P_crit)
   real(wp), parameter :: NARROW_FAC     = 0.60_wp

   integer,  parameter :: NSTAGE         = 2                       ! S1, S2

   ! Default location of the orientation-resolved DH21 spheroid table and its
   ! grid axes, resolved against the data root (sed_paths), so that they follow
   ! build_dust's data_dir like every other default the library opens.
   ! Overridable through sed_init / build_astrodust for a host naming its own.
   character(len=*), parameter :: QPOL_Q_DEF = &
        'astrodust/q_DH21Ad_P0.20_Fe0.00_1.400.dat.gz'
   character(len=*), parameter :: QPOL_W_DEF = 'astrodust/DH21_wave'
   character(len=*), parameter :: QPOL_A_DEF = 'astrodust/DH21_aeff'
   ! Companion table for the extreme ultraviolet, 0.0124 um (100 eV) up to the
   ! 0.0912 um (13.6 eV) node where the table above starts. Computed by
   ! tmatrix/driver/run_q_jori.f90 in its `euv` mode from the same DH21
   ! dielectric function and the same b/a = 1.4 spheroid, so it is the same
   ! material and shape as the table above, only at shorter wavelengths. It
   ! shares the DH21_aeff size axis (QPOL_A_DEF).
   character(len=*), parameter :: QPOL_EUV_Q_DEF = &
        'astrodust/q_astrodust_jori_euv_P0.20_Fe0.00_1.400.dat.gz'
   character(len=*), parameter :: QPOL_EUV_W_DEF = 'astrodust/DH21_wave_euv'

   ! Default size-integrated extinction table of each coded model, relative to
   ! the sed/ working directory -- what dust_extinction serves when the builder
   ! was given no kext_path.
   !
   ! The astrodust and DL07 defaults are the EUV products, the WIDEST grid each
   ! model has: they run from ~1e-4 um up to 39810 um, so one file covers a host
   ! that transports ionizing radiation and a host that does not, the latter's
   ! grid being a subset of the former's. Both were written on the same
   ! wavelength nodes as the T-matrix Q table above 0.0912 um, so an unextended
   ! model reads its own nodes straight out of them. The narrower non-EUV
   ! products (kext_astrodust_MW.dat, kext_dl07_MW.dat) remain available, but a
   ! host that wants one has to name it.
   !
   ! Zubko has the same choice for the opposite reason: its optics table IS the
   ! model definition and already reaches 0.001 um, so its _euv product is that
   ! table's own range and the narrower one is that range cut at the Lyman
   ! limit.  The default is the wider one here too.
   ! build_from_files has no default at all -- a file-defined model's product is
   ! named after the model, which is only known once the descriptor is read.
   character(len=*), parameter :: KEXT_ASTRODUST = 'astrodust/kext_astrodust_MW_euv.dat'
   character(len=*), parameter :: KEXT_DL07      = 'dl07/kext_dl07_MW_euv.dat'
   character(len=*), parameter :: KEXT_ZUBKO     = 'zubko/kext_zubko_BARE_GR_S_euv.dat'
   ! The same curve inside the model's own HDF5 product, which is tried FIRST.
   ! /kext is read on the whole wavelength axis, so it covers every grid the
   ! model can be built on -- which is what lets a host leave kext_path blank
   ! whatever grid it asked for, instead of having to pair a narrow model with
   ! a narrow table and a wide one with the _euv file.  A tree that has no
   ! product, or a build made without HDF5, falls back to the text defaults
   ! above and nothing changes for it.
   character(len=*), parameter :: KEXT_H5_ASTRODUST = 'astrodust/sedust_astrodust.h5'
   character(len=*), parameter :: KEXT_H5_DL07      = 'dl07/sedust_dl07.h5'
   character(len=*), parameter :: KEXT_H5_ZUBKO     = 'zubko/sedust_zubko.h5'
   character(len=*), parameter :: KEXT_H5_MRN       = 'mrn/sedust_mrn.h5'
   character(len=*), parameter :: KEXT_MRN          = 'mrn/kext_mrn_euv.dat'

   ! ---- the MRN (1977) power law -----------------------------------------
   ! dn_i/da = A_i n_H a^-3.5 over a sharp [a_min, a_max], one power law per
   ! material (Draine & Lee 1984, eq. 5.1).  Index 1 is graphite, 2 silicate,
   ! which is the order Draine's own parameter file for this model lists them
   ! in and the order the channels below carry.
   !
   ! The abundances are the ones DL84 sec. Va adopted after fitting the Savage
   ! & Mathis average extinction curve: "we adopted A_sil = 10^-25.11
   ! cm^2.5/H, and A_C = 10^-25.16 cm^2.5/H".  They are what Draine's own MRN
   ! model is built on -- his parameter file states the grain volumes per H,
   ! 2.49e-27 and 2.79e-27 cm^3/H, and those ARE these two A_i -- so the
   ! reference curve shipped beside this tree is the size integral of this
   ! model and nothing here has to be selected to match it.
   !
   ! log10 A_i, in cm^2.5 per H, as the paper prints it.  The A_i themselves
   ! are formed where they are used: a real exponent is not a constant
   ! expression, and rounding the powers by hand here would put a number in
   ! the code that the paper does not contain.
   real(wp), parameter :: MRN_ALPHA = -3.5_wp
   real(wp), parameter :: MRN_AMIN  = 5.0e-3_wp     ! [um]
   real(wp), parameter :: MRN_AMAX  = 0.25_wp       ! [um]
   real(wp), parameter :: LOG_A_DL84(2) = [-25.16_wp, -25.11_wp]
   ! Radii, log-spaced over that range: 70 points is 0.025 dex, twice the
   ! resolution of the DL07 size grid, which a power law cut sharply at both
   ! ends needs and the stochastic solve can afford over this size range.
   integer,  parameter :: MRN_NA    = 70

contains

   ! =====================================================================
   subroutine sed_register_euv_band_optics(proc)
      !! Name the calculation sed_init is to use for the astrodust EUV band
      !! when euv_tmatrix = .true.  One call before the model is built;
      !! euv_astrodust_tmatrix.f90 wraps it as use_tmatrix_euv_band_optics().
      procedure(euv_band_optics_i) :: proc
      euv_band_optics => proc
   end subroutine sed_register_euv_band_optics

   subroutine sed_forget_euv_band_optics()
      !! Undo the above: euv_tmatrix = .true. is then refused (status 11) and
      !! only the volume-equivalent sphere remains.
      euv_band_optics => null()
   end subroutine sed_forget_euv_band_optics

   ! =====================================================================
   subroutine sed_init(qtable_path, NT_in, T_lo, T_hi, status, &
                       qpol_path, qpol_wave_path, qpol_aeff_path, scatmat_path, &
                       load_polarized_optics, lam_min, astrodust_index_path, &
                       qpol_euv_path, qpol_euv_wave_path, euv_tmatrix, include_euv)
      character(len=*), intent(in) :: qtable_path
      integer,          intent(in) :: NT_in
      real(wp),         intent(in) :: T_lo, T_hi
      ! Optional status (0 = success). When present, a failed input read is
      ! reported through it instead of stopping the process; when absent the
      ! readers keep their message + stop behavior (as the CLI drivers expect).
      !   status = 1  Q-table load failed
      !   status = 3  aligned scattering table load failed (only when
      !               scatmat_path is supplied)
      !   status = 4  a polarized Q table explicitly requested via qpol_path
      !               could not be read (an implicit-default table degrades
      !               gracefully instead -- an explicit request must not vanish)
      !   status = 5  load_polarized_optics = .false. was combined with an
      !               explicit polarized-optics path (qpol_path / qpol_wave_path
      !               / qpol_aeff_path / scatmat_path) -- a contradictory request
      !   status = 6  astrodust dielectric function load failed (EUV band only)
      !   status = 7  lam_min below the astrodust dielectric function's own
      !               shortest wavelength (EUV band only)
      !   status = 8  the polarized table's wavelength grid does not match the
      !               model's grid (a missing EUV companion table is NOT an
      !               error -- that band degrades to a reported zero)
      !   status = 9  lam_min runs outside the EUV companion table's coverage
      !   status = 11 euv_tmatrix = .true. but the spheroid optics of the EUV
      !               band are not available: no implementation of
      !               euv_band_optics_i is registered, or the registered one
      !               reported that it could not compute the band
      integer, optional, intent(out) :: status
      ! Orientation-resolved DH21 table and its grid axes, supplying the
      ! polarized optics. Default to QPOL_*_DEF. An implicit-default table that
      ! cannot be read is NOT an error: Cpol and falign_ad stay zero and the run
      ! continues without polarized emission. A table supplied EXPLICITLY through
      ! qpol_path that cannot be read IS an error (status 4).
      character(len=*), optional, intent(in) :: qpol_path, qpol_wave_path, &
                                                qpol_aeff_path
      ! Aligned scattering table (run_scatmat_aligned.x product) for a polarized
      ! RT host. When present it is parsed into scatmat_aligned_mod here at init;
      ! its failure IS an error (status 3), unlike the implicit polarized-optics
      ! table, because a host that asked for it depends on it.
      character(len=*), optional, intent(in) :: scatmat_path
      ! Scalar-only switch. Absent or .true.: load the polarized optics as above.
      ! .false.: the polarized Q table is NEVER opened -- no default-path
      ! substitution, no decompression scratch file -- and Cpol / Cpol_ext /
      ! Cbir_ext / falign_ad stay allocated and zero, so the scalar Cext / Cabs /
      ! Csca / gbar and the SED are identical to a polarized build with zero
      ! alignment. Combining .false. with an explicit qpol_*/scatmat path is a
      ! caller contradiction (status 5).
      logical, optional, intent(in) :: load_polarized_optics
      ! Optional shortest wavelength [um] the model must cover. When it is
      ! shorter than the Q table's 0.0912 um the grid is carried down to it and
      ! the astrodust optics there come from the DH21 dielectric function
      ! (q_astrodust_mod) instead of the table. Absent = the table grid alone.
      ! The polarized optics of that EUV block come from the EUV companion
      ! table (QPOL_EUV_*_DEF); see build_Cpol.
      real(wp), optional, intent(in) :: lam_min
      ! Optional dielectric function for that EUV band. It must be the file the
      ! Q table was computed from -- same porosity, iron fraction and axial
      ! ratio -- or the model changes material at the seam. Omitted = the
      ! q_astrodust_mod default, which pairs with
      ! q_astrodust_P0.20_Fe0.00_1.400.dat.
      character(len=*), optional, intent(in) :: astrodust_index_path
      ! Optional EUV companion polarized table and its wavelength axis, used
      ! only when lam_min carries the grid below the main table. Omitted = the
      ! QPOL_EUV_*_DEF defaults, which pair with the same material. Unlike the
      ! main polarized table, a failure here is ALWAYS an error (status 8/9):
      ! the whole point of the companion table is that this band no longer
      ! returns an unannounced zero.
      character(len=*), optional, intent(in) :: qpol_euv_path, qpol_euv_wave_path
      ! How the EUV band's optics are computed. Default .true.: the T-matrix
      ! on the b/a = 1.400 oblate spheroid, i.e. the same particle and the same
      ! random-orientation average as the Q table, so the grain does not change
      ! shape at the seam. That route is the registered euv_band_optics_i, and
      ! asking for it without one registered is an error (status 11), not a
      ! silent substitution. .false. substitutes the volume-equivalent-sphere
      ! Mie approximation of q_astrodust_mod, which is ~2% low in the
      ! geometric-optics limit (that module measures it size by size) but costs
      ! milliseconds where the T-matrix costs minutes. Ignored when there is no
      ! EUV band.
      logical, optional, intent(in) :: euv_tmatrix
      ! Which wavelength axis to take out of an HDF5 product; ignored when
      ! qtable_path names a text table, which carries one axis and no choice.
      ! Default .false., the non-ionizing part -- the same grid the narrow text
      ! table gives, so naming the HDF5 file changes the source and not the
      ! model.
      logical, optional, intent(in) :: include_euv
      integer  :: i, ja, jw, jt, is, scstat, euv_stat
      real(wp) :: a_um, x, t, Q_neu, Q_ion, Q_sca_c, g_sca_c
      real(wp) :: qext1, qsca1, qabs1, gsca1
      real(wp) :: ad_lam_lo, ad_lam_hi
      real(wp), allocatable :: lam_grid(:)
      logical  :: rok, want_pol, pol_explicit, pol_loaded
      character(len=512) :: pol_q, pol_w, pol_a, pol_eq, pol_ew
      ! EUV band: route selector and the band's optics.
      logical  :: use_tm
      integer  :: euv_optics_stat, n_euv_bad
      real(wp), allocatable :: nr_euv(:), ki_euv(:)
      real(wp), allocatable :: q_euv_abs(:,:), q_euv_sca(:,:), q_euv_g(:,:)

      if (present(status)) status = 0
      use_tm = .true.
      if (present(euv_tmatrix)) use_tm = euv_tmatrix
      if (present(astrodust_index_path)) &
         call set_astrodust_index_path(astrodust_index_path)

      ! Scalar-only vs polarized: load the polarized optics unless told not to.
      want_pol = .true.
      if (present(load_polarized_optics)) want_pol = load_polarized_optics
      ! A scalar-only build must not also be handed polarized inputs: which does
      ! the caller mean? Report the contradiction rather than silently guessing.
      pol_explicit = present(qpol_path) .or. present(qpol_wave_path) .or. &
                     present(qpol_aeff_path) .or. present(scatmat_path) .or. &
                     present(qpol_euv_path) .or. present(qpol_euv_wave_path)
      if (.not. want_pol .and. pol_explicit) then
         if (present(status)) then
            status = 5;  return
         else
            write(error_unit,'(a)') ' sed_init: load_polarized_optics=.false. '// &
               'contradicts an explicit polarized-optics path (qpol_*/scatmat)'
            stop 1
         end if
      end if

      ! ---- Load Q table and size dist (modules cache their own state) ----
      ! An HDF5 path takes the table out of /qtable/astrodust; anything else is
      ! the text table run_tmatrix.x writes, read exactly as before.
      if (is_hdf5_path(qtable_path)) then
         call load_q_table_h5(qtable_path, euv_asked(include_euv), rok)
         if (.not. rok) then
            if (present(status)) then
               status = 1;  return
            else
               write(*,'(a,a)') ' sed_init: cannot read /qtable/astrodust from ', &
                                trim(qtable_path)
               stop 1
            end if
         end if
      else if (present(status)) then
         call load_q_table(qtable_path, ok=rok)
         if (.not. rok) then;  status = 1;  return;  end if
      else
         call load_q_table(qtable_path)
      end if
      ! HD23 size distributions, alignment and ionization fractions: analytic
      ! (size_dist_mod), on the 167-radius grid of the HD23 release table.
      call hd23_size_distribution()

      call euv_extended_lambda_grid(lam_grid, lam_min, n_extra=n_lam_euv)

      ! ---- EUV band below the Q table -----------------------------------
      ! Its optics come from the DH21 dielectric function, so load that file
      ! HERE, before anything is built, rather than letting the first optics
      ! call load it lazily inside the size loop: a missing file or a lam_min
      ! the file cannot cover has to reach the caller through `status`, which
      ! is what sed_init promises, and not stop the process out of an RT host.
      ! Also name the file, so its pairing with the Q table is visible.
      if (n_lam_euv > 0) then
         ! The spheroid route is an injected calculation, so a build that does
         ! not carry one cannot honor euv_tmatrix = .true.  Say so before
         ! anything is built rather than answering with the sphere, which is a
         ! different particle and would leave a step at the seam.
         if (use_tm .and. .not. associated(euv_band_optics)) then
            write(*,'(a)') ' sed_init: euv_tmatrix = .true., but no spheroid optics are'
            write(*,'(a)') '           registered for the EUV band.  Build with the'
            write(*,'(a)') '           T-matrix (euv_astrodust_tmatrix.f90) and call'
            write(*,'(a)') '           use_tmatrix_euv_band_optics() once beforehand,'
            write(*,'(a)') '           or pass euv_tmatrix = .false. for the'
            write(*,'(a)') '           volume-equivalent-sphere approximation.'
            if (present(status)) then
               status = 11;  return
            else
               stop 1
            end if
         end if
         if (present(status)) then
            call load_astrodust_index(ok=rok)
            if (.not. rok) then;  status = 6;  return;  end if
         else
            call load_astrodust_index()
         end if
         call astrodust_index_lambda_range(ad_lam_lo, ad_lam_hi)
         if (lam_grid(1) < ad_lam_lo) then
            write(*,'(a,es10.3,a)') ' sed_init: lam_min =', lam_grid(1), &
               ' um is shorter than the astrodust dielectric function,'
            write(*,'(a,es10.3,a)') '           which stops at', ad_lam_lo, &
               ' um; (n, k) below it would be frozen at the boundary value.'
            if (present(status)) then
               status = 7;  return
            else
               stop 1
            end if
         end if
         write(*,'(a,i0,a)') ' sed_init: EUV extension active, ', n_lam_euv, &
            ' wavelengths below the Q table.'
         write(*,'(a,a)')    '           astrodust dielectric function: ', &
            trim(get_astrodust_index_path())
         if (use_tm) then
            ! The axial ratio is stated by whoever implements the spheroid
            ! route, since that is where the particle is defined; naming a
            ! number here could only repeat it, or contradict it.
            write(*,'(a)')      '           EUV optics: T-matrix, oblate spheroid'// &
               ' (same particle as the table).'
         else
            write(*,'(a)')      '           EUV optics: Mie, volume-equivalent'// &
               ' sphere (shape approximation).'
         end if
      end if

      NLAM      = size(lam_grid)
      NA        = sd_n
      NT        = NT_in

      call free_shared_model_arrays()
      if (allocated(Cpol))     deallocate(Cpol, Cpol_ext, Cbir_ext, falign_ad)
      if (allocated(gsca_ad))  deallocate(gsca_ad)
      allocate(Cpol(NLAM, NA), Cpol_ext(NLAM, NA), Cbir_ext(NLAM, NA), falign_ad(NA))
      allocate(gsca_ad(NLAM, NA))
      allocate(lam(NLAM), aeff(NA), T_first(NT), dn_ad(NA))
      allocate(Cabs(NLAM, NA), Csca(NLAM, NA), kappB_first(NT, NA), &
               H_first(NT, NA, NSTAGE), kappCMB(NA))
      allocate(log_T_first(NT), log_H_first(NT, NA, NSTAGE), &
               log_kappB_first(NT, NA))
      allocate(dn_pah(NA), Cabs_pah(NLAM, NA), kappB_pah_first(NT, NA), &
               H_pah_first(NT, NA), kappCMB_pah(NA))
      allocate(log_H_pah_first(NT, NA), log_kappB_pah_first(NT, NA))
      ! Charge-resolved carbonaceous arrays (neutral + cation kept separate).
      if (allocated(Cabs_cneu)) deallocate(Cabs_cneu, Cabs_cion, dn_cneu, dn_cion, &
               kappB_cneu, kappB_cion, log_kappB_cneu, log_kappB_cion, &
               kappCMB_cneu, kappCMB_cion)
      if (allocated(Csca_pah)) deallocate(Csca_pah, gsca_pah)
      allocate(Cabs_cneu(NLAM, NA), Cabs_cion(NLAM, NA), dn_cneu(NA), dn_cion(NA), &
               kappB_cneu(NT, NA), kappB_cion(NT, NA), &
               log_kappB_cneu(NT, NA), log_kappB_cion(NT, NA), &
               kappCMB_cneu(NA), kappCMB_cion(NA))
      allocate(Csca_pah(NLAM, NA), gsca_pah(NLAM, NA))

      lam    = lam_grid
      aeff   = sd_aeff
      dn_ad  = sd_dn
      dn_pah = sd_dn_pah
      deallocate(lam_grid)

      ! ---- Build wide T grid (log-spaced) used for the smallest grain
      ! and as the source for narrowed T windows for subsequent grains ----
      do i = 1, NT
         t = log(T_lo) + (log(T_hi) - log(T_lo)) * real(i-1, wp) / real(NT-1, wp)
         T_first(i) = exp(t)
      end do
      log_T_first = log(T_first)

      ! ---- Setup p_sub's lambda-integration weights ----
      call p_sub_setup(lam)

      ! ---- EUV band, spheroid route -------------------------------------
      ! The registered euv_band_optics_i computes the whole band -- every
      ! wavelength on every radius -- in one call, on the DH21 astrodust index.
      ! Same material, same oblate spheroid and same random-orientation average
      ! as the Q table, so the extension continues the table's own calculation
      ! and the grain does not change shape at the seam.  The refractive index
      ! is a function of wavelength alone, so it is evaluated once here, before
      ! any threading inside that routine, rather than again for every radius.
      n_euv_bad = 0
      if (n_lam_euv > 0 .and. use_tm) then
         allocate(nr_euv(n_lam_euv), ki_euv(n_lam_euv))
         do jw = 1, n_lam_euv
            call astrodust_index_at(lam(jw), nr_euv(jw), ki_euv(jw))
         end do
         allocate(q_euv_abs(n_lam_euv, NA), q_euv_sca(n_lam_euv, NA), &
                  q_euv_g(n_lam_euv, NA))
         call euv_band_optics(lam(1:n_lam_euv), aeff, nr_euv, ki_euv, &
                              q_euv_abs, q_euv_sca, q_euv_g, euv_optics_stat)
         deallocate(nr_euv, ki_euv)
         if (euv_optics_stat < 0) then
            write(*,'(a)') ' sed_init: the registered EUV band optics could not'
            write(*,'(a)') '           compute the band.'
            if (present(status)) then
               status = 11;  return
            else
               stop 1
            end if
         end if
         ! Points that failed a physical bound: warned about by the routine
         ! itself, counted here for the summary below.
         n_euv_bad = euv_optics_stat
      end if

      ! ---- Cabs(NLAM, NA) and Csca(NLAM, NA) by interpolating Q in log(a) ----
      ! Q table grid (qt_aeff) is denser than dist grid (aeff); interpolate to
      ! the dist grid where the size sum lives.
      do ja = 1, NA
         a_um = aeff(ja)
         x = log(a_um)
         call interp_q_grid(x, qt_aeff, qt_qabs, Cabs(n_lam_euv+1:, ja))
         call interp_q_grid(x, qt_aeff, qt_qsca, Csca(n_lam_euv+1:, ja))
         ! Asymmetry <cos> comes from the same table on the same grid, but it
         ! is already dimensionless -- no pi a^2 conversion.
         call interp_q_grid(x, qt_aeff, qt_gpar, gsca_ad(n_lam_euv+1:, ja))
         ! EUV band below the Q table (n_lam_euv = 0 unless lam_min asked for
         ! it). The spheroid route was computed just above and is only copied
         ! in here; euv_tmatrix = .false. instead substitutes the
         ! volume-equivalent SPHERE (Mie) on the same DH21 dielectric function,
         ! and q_astrodust_mod carries the domain of validity of that shape
         ! approximation.
         do jw = 1, n_lam_euv
            if (use_tm) then
               Cabs(jw, ja)    = q_euv_abs(jw, ja)
               Csca(jw, ja)    = q_euv_sca(jw, ja)
               gsca_ad(jw, ja) = q_euv_g(jw, ja)
            else
               call q_astrodust_full(a_um, lam(jw), qext1, qsca1, qabs1, gsca1)
               Cabs(jw, ja)    = qabs1
               Csca(jw, ja)    = qsca1
               gsca_ad(jw, ja) = gsca1
            end if
         end do
         ! Convert Q -> C: C = pi * (a_cm)^2 * Q
         Cabs(:, ja) = Cabs(:, ja) * PI * (a_um * UM2CM)**2
         Csca(:, ja) = Csca(:, ja) * PI * (a_um * UM2CM)**2
      end do
      if (allocated(q_euv_abs)) deallocate(q_euv_abs, q_euv_sca, q_euv_g)
      if (n_euv_bad > 0) &
         write(*,'(a,i0,a)') ' sed_init: EUV optics: ', n_euv_bad, &
            ' point(s) failed a physical bound (values kept; see the warnings above).'

      ! ---- Cpol(NLAM, NA), falign_ad(NA) from the DH21 spheroid table ----
      ! Scalar-only build: skip the table entirely (no default-path substitution,
      ! no decompression scratch file) and leave the polarized arrays allocated
      ! and zero, exactly as a failed read would.
      polarized_optics_loaded = .false.
      if (want_pol) then
         pol_q = sed_data_path(QPOL_Q_DEF);  if (present(qpol_path))      pol_q = qpol_path
         pol_w = sed_data_path(QPOL_W_DEF);  if (present(qpol_wave_path)) pol_w = qpol_wave_path
         pol_a = sed_data_path(QPOL_A_DEF);  if (present(qpol_aeff_path)) pol_a = qpol_aeff_path
         pol_eq = sed_data_path(QPOL_EUV_Q_DEF)
         if (present(qpol_euv_path))      pol_eq = qpol_euv_path
         pol_ew = sed_data_path(QPOL_EUV_W_DEF)
         if (present(qpol_euv_wave_path)) pol_ew = qpol_euv_wave_path
         call build_Cpol(trim(pol_q), trim(pol_w), trim(pol_a), &
                         trim(pol_eq), trim(pol_ew), pol_loaded, euv_stat)
         polarized_optics_loaded = pol_loaded
         ! The EUV band never degrades silently: whichever way the table was
         ! named, a host that asked for those wavelengths gets physics or an
         ! error, so this is checked before the explicit/implicit rule below.
         if (euv_stat /= 0) then
            if (present(status)) then
               status = euv_stat;  return
            else
               write(error_unit,'(a,i0,a)') ' sed_init: EUV polarized optics '// &
                  'unavailable (code ', euv_stat, '); see the message above.'
               stop 1
            end if
         end if
         if (.not. pol_loaded .and. present(qpol_path)) then
            ! An explicitly requested table that cannot be read is an error: the
            ! capability the caller asked for must not disappear silently. (An
            ! implicit-default table degrades gracefully in build_Cpol above.)
            if (present(status)) then
               status = 4;  return
            else
               write(error_unit,'(a)') ' sed_init: explicitly requested '// &
                  'polarized table could not be read: '//trim(pol_q)
               stop 1
            end if
         end if
      else
         Cpol = 0.0_wp;  Cpol_ext = 0.0_wp;  Cbir_ext = 0.0_wp;  falign_ad = 0.0_wp
      end if

      ! ---- aligned scattering table for a polarized RT host (optional) ----
      if (present(scatmat_path)) then
         call load_scatmat_aligned(scatmat_path, scstat)
         if (scstat /= 0) then
            if (present(status)) then
               status = 3;  return
            else
               write(*,'(a,i0,a)') ' sed_init: cannot load aligned scattering table (code ', &
                    scstat, '): '//trim(scatmat_path)
               stop 1
            end if
         end if
      end if

      ! ---- kappB_first(NT, NA) = integral of Cabs * B_lambda over lambda ----
      call build_kappB()
      log_kappB_first = log(max(kappB_first, tiny(0.0_wp)))

      ! ---- kappCMB(NA) = CMB Planck integral (used in calc_P) ----
      call build_kappCMB()

      ! ---- H_first(NT, NA, 2) for the two astrodust enthalpy stages ----
      do is = 1, NSTAGE
         do ja = 1, NA
            do jt = 1, NT
               select case (is)
               case (1);  H_first(jt, ja, is) = enthalpy_S1(T_first(jt), aeff(ja))
               case (2);  H_first(jt, ja, is) = enthalpy_S2(T_first(jt), aeff(ja))
               end select
            end do
         end do
      end do
      log_H_first = log(max(H_first, tiny(0.0_wp)))

      ! ---- PAH population (HD23 §3.2 / Draine & Li 2007) ------------
      ! Neutral and cation are kept as SEPARATE stochastically-heated
      ! populations (NOT pre-blended): Teq, P(T) and the emission are all
      ! nonlinear in Cabs, so the correct mixture emission is
      !   (1-f_ion) E[C^neu] + f_ion E[C^ion],   each solved independently,
      ! not  E[(1-f_ion)C^neu + f_ion C^ion].  f_ion(a) from HD23
      ! eq. (19). C = pi*(a_cm)^2 * Q convention.
      do ja = 1, NA
         a_um = aeff(ja)
         do jw = 1, NLAM
            if (use_ld01_pah_xsec) then
               ! Use the Li & Draine (2001) PAH cross sections via the ported
               ! QPAH_LD01 routine (drude_ld01 + same xi_gra graphite blend),
               ! instead of our DL07 qpah_dl07.
               call qpah_ld01(0, a_um, lam(jw), Q_neu)
               call qpah_ld01(1, a_um, lam(jw), Q_ion)
            else
               call qpah_dl07(0, a_um, lam(jw), Q_neu)
               call qpah_dl07(1, a_um, lam(jw), Q_ion)
            end if
            Cabs_cneu(jw, ja) = Q_neu * PI * (a_um * UM2CM)**2
            Cabs_cion(jw, ja) = Q_ion * PI * (a_um * UM2CM)**2
            ! Scattering of the same grain: the graphite fraction of HD23
            ! eq. 15 scatters, the PAH part does not.  Its whole weight is in
            ! the far-UV -- 3.5% of this model's total tau_sca at 0.1 um in
            ! the HD23 release, and below 0.2% longward of 0.15 um.
            call qpah_sca(a_um, lam(jw), Q_sca_c, g_sca_c)
            Csca_pah(jw, ja) = Q_sca_c * PI * (a_um * UM2CM)**2
            gsca_pah(jw, ja) = g_sca_c
         end do
         dn_cneu(ja) = (1.0_wp - sd_fion(ja)) * dn_pah(ja)
         dn_cion(ja) =          sd_fion(ja)   * dn_pah(ja)
         ! Draine's charge cutoff:
         ! PAHs larger than ~100 A (a > 0.99999e-6 cm = 0.99999e-2 um) are
         ! treated as 100% ionized -- a single cation-only solve, no neutral
         ! channel. a_um is in microns.
         if (a_um > 0.99999e-2_wp) then
            dn_cneu(ja) = 0.0_wp
            dn_cion(ja) = dn_pah(ja)
         end if
      end do
      ! Charge-resolved kappB / kappCMB, built via the Cabs_pah scratch (build_*_pah
      ! read module Cabs_pah). Both stored so the QM batch can use both at once.
      Cabs_pah = Cabs_cneu
      call build_kappB_pah();   kappB_cneu   = kappB_pah_first
      call build_kappCMB_pah(); kappCMB_cneu = kappCMB_pah
      Cabs_pah = Cabs_cion
      call build_kappB_pah();   kappB_cion   = kappB_pah_first
      call build_kappCMB_pah(); kappCMB_cion = kappCMB_pah
      log_kappB_cneu = log(max(kappB_cneu, tiny(0.0_wp)))
      log_kappB_cion = log(max(kappB_cion, tiny(0.0_wp)))
      ! Restore the f_ion-blended Cabs_pah/kappB_pah/kappCMB_pah for any
      ! external driver that reads them.
      call build_Cabs_pah()
      call build_kappB_pah()
      call build_kappCMB_pah()
      log_kappB_pah_first = log(max(kappB_pah_first, tiny(0.0_wp)))
      ! Shared 'Car0' DL01 carbonaceous enthalpy for both charge states.
      do ja = 1, NA
         do jt = 1, NT
            H_pah_first(jt, ja) = enthalpy_DL01(T_first(jt), aeff(ja), 'Car0')
         end do
      end do
      log_H_pah_first = log(max(H_pah_first, tiny(0.0_wp)))

      initialized = .true.
   end subroutine sed_init

   ! =====================================================================
   subroutine sed_solve(J_lam, enthalpy_stage, lamI_lam_out)
      real(wp),         intent(in)  :: J_lam(:)        ! (NLAM)
      character(len=*), intent(in)  :: enthalpy_stage  ! 'S1' | 'S2'
      real(wp),         intent(out) :: lamI_lam_out(:) ! (NLAM)

      integer  :: is
      real(wp), allocatable :: Jout(:)

      if (.not. initialized) then
         write(*,'(a)') 'sed_solve: call sed_init first'
         stop 1
      end if

      select case (trim(enthalpy_stage))
      case ('S1');     is = 1
      case ('S2');     is = 2
      case default
         write(*,'(a,a)') 'sed_solve: unknown stage ', trim(enthalpy_stage)
         stop 1
      end select

      allocate(Jout(NLAM))

      call sed_grain_loop(NA, dn_ad, aeff, Cabs, kappB_first, H_first(:,:,is), &
                          log_H_first(:,:,is), log_kappB_first, kappCMB, &
                          J_lam, 'sil', Jout)

      if (use_induced_emission) call apply_induced_factor(J_lam, Jout)

      ! Unit conversion to HD23 convention (erg/s/sr/H):
      !   Jout = sum_a (dn_Ad/N_H per bin)(a) * Cabs(cm^2) * bbody(SI W/m^3/sr)
      ! Cabs * bbody:  cm^2 * W/m^3/sr
      !   = cm^2 * 10 erg/(s*cm^3*sr) per cm wavelength            [SI -> CGS]
      !   = 10 * erg/(s*cm*sr) per cm wavelength per grain
      ! After * (#grains/H per bin) and summing: 10 * erg/(s*cm*sr*H) per cm wave.
      ! lambda*I_lambda per H = lambda(cm) * Jout(CGS)
      !                       = (lambda_um * 1e-4) * (10 * Jout)
      !                       = lambda_um * Jout * 1e-3
      ! No 4*pi divisor (I_lam from a 1-H column in the optically-thin limit
      ! for an isotropic emitter is already integrated correctly without
      ! one -- the 4pi cancels between emission isotropy and the steradian
      ! denominator of B_lambda).
      lamI_lam_out = lam * Jout * 1.0e-3_wp

      deallocate(Jout)
   end subroutine sed_solve


   subroutine sed_solve_pah(J_lam, lamI_lam_out)
      ! PAH-population SED solve. Same dynamic-T algorithm as
      ! sed_solve(), but with PAH cross sections. Neutral and cation are
      ! solved as SEPARATE stochastically-heated populations (different
      ! absorption -> different T distribution) and summed, mirroring
      ! sed_solve_dl07 -- NOT pre-blended by f_ion. PAH size distribution
      ! of HD23 eq. (17), split into neutral/cation by
      ! f_ion(a); DL01 'Car0' carbonaceous enthalpy shared by both states.
      real(wp), intent(in)  :: J_lam(:)
      real(wp), intent(out) :: lamI_lam_out(:)

      real(wp), allocatable :: Jout(:), Jout_q(:)
      integer :: icharge

      if (.not. initialized) then
         write(*,'(a)') 'sed_solve_pah: call sed_init first'
         stop 1
      end if

      allocate(Jout(NLAM), Jout_q(NLAM))
      Jout = 0.0_wp
      do icharge = 0, 1
         if (icharge == 0) then
            call sed_grain_loop(NA, dn_cneu, aeff, Cabs_cneu, kappB_cneu, &
                                H_pah_first, log_H_pah_first, &
                                log_kappB_cneu, kappCMB_cneu, &
                                J_lam, 'pah', Jout_q)
         else
            call sed_grain_loop(NA, dn_cion, aeff, Cabs_cion, kappB_cion, &
                                H_pah_first, log_H_pah_first, &
                                log_kappB_cion, kappCMB_cion, &
                                J_lam, 'pah', Jout_q)
         end if
         Jout = Jout + Jout_q
      end do

      if (use_induced_emission) call apply_induced_factor(J_lam, Jout)

      ! Same unit conversion as sed_solve(): 1e-3 takes
      ! Cabs[cm^2] * bbody[SI W/m^3/sr] * dn[1/H] * lambda[um]
      ! to lambda*I_lambda in erg/s/sr/H (the HD23 convention).
      lamI_lam_out = lam * Jout * 1.0e-3_wp

      deallocate(Jout, Jout_q)
   end subroutine sed_solve_pah


   ! =====================================================================
   ! DL07 (Draine & Li 2007) model: amorphous silicate + carbonaceous
   ! (PAH + graphite) grains with WD01 size distributions (MW/LMC/SMC
   ! selected by sd_index). Reuses the shared SED-solver core; the
   ! astrodust path (sed_init / sed_solve) is untouched.
   !
   ! The two population slots are repurposed:
   !   dust slot  -> amorphous silicate  (q_silicate, enthalpy 'Sil')
   !   PAH slot   -> full carbonaceous   (qpah_dl07 blended optics,
   !                                      enthalpy 'Car0'),
   ! with the PAH ionization fraction computed directly from the WD01b
   ! grain-charging model (pah_ionfrac) at intensity u_isrf.
   ! =====================================================================
   subroutine sed_init_dl07(qtable_path, sd_index, u_isrf, &
                            NT_in, T_lo, T_hi, status, lam_min, lam_axis, include_euv, &
                            stored_q_dir)
      character(len=*), intent(in) :: qtable_path
      integer,          intent(in) :: sd_index, NT_in
      real(wp),         intent(in) :: u_isrf, T_lo, T_hi
      ! Optional status (0 = success). When present, a failed input read is
      ! reported through it instead of stopping the process; when absent the
      ! readers keep their message + stop behavior (as the CLI drivers expect).
      !   status = 1  Q-table load failed
      !   status = 7  lam_min below the D03 dielectric functions' own shortest
      !               wavelength (EUV band only)
      integer, optional, intent(out) :: status
      ! Optional shortest wavelength [um] the model must cover; see sed_init.
      ! The Q table supplies this model's wavelength grid only -- its silicate
      ! and carbonaceous optics are already Mie on the D03 dielectric
      ! functions, which run to 6.2e-5 um, so extending the grid is all the EUV
      ! needs here.
      real(wp), optional, intent(in) :: lam_min
      ! The model's wavelength axis, given outright.  This model takes only a
      ! grid from qtable_path -- its optics are Mie on the D03 dielectric
      ! functions at every wavelength -- so a caller that already has the axis
      ! passes it here and qtable_path is not read at all.  That is how the
      ! HDF5 product is used: data/dl07/sedust_dl07.h5 carries this model's own
      ! /grid/lambda, so a host running DL07 alone needs no astrodust file to
      ! borrow a grid from.  lam_min is then ignored: the axis is final.
      real(wp), optional, intent(in) :: lam_axis(:)
      ! Which wavelength axis to take when qtable_path names an HDF5 product.
      ! This model's own product carries its own /grid/lambda, so naming
      ! data/dl07/sedust_dl07.h5 here is the same as handing the axis through
      ! lam_axis.  Ignored for a text table, and by lam_axis when that is given.
      logical, optional, intent(in) :: include_euv
      ! Where this model's stored cross sections come from.  Omitted, the
      ! model's own directory under the data root, with the /qtable of an HDF5
      ! qtable_path tried ahead of it -- so one file supplies both the axis and
      ! the numbers, and no caller can pair the axis of one source with the
      ! optics of another.  Passed BLANK, NOTHING stored is read from anywhere
      ! and every optic is solved from the dielectric functions: that is what
      ! calc_qtable.x asks for, being the program that writes these tables, and
      ! what a test comparing two builds of one model asks for, so that both
      ! take the same route whatever grid they are on.
      character(len=*), optional, intent(in) :: stored_q_dir

      integer  :: i, ja, jw, jt
      real(wp) :: a_um, t, da, qabs1, Q_neu, Q_ion
      real(wp) :: geo, qext1, qsca1, gsca1
      real(wp) :: sil_lam_lo, sil_lam_hi, gra_lam_lo, gra_lam_hi, d03_lam_lo
      real(wp), allocatable :: fion(:), lna(:), lam_grid(:), lam_base(:)
      ! Stored cross sections of the four DL07 populations, when tables for
      ! this grid exist; got_* says whether each was found.
      real(wp), allocatable :: tQa(:,:), tQs(:,:), tGg(:,:)
      real(wp), allocatable :: uQa(:,:), uQs(:,:), uGg(:,:)
      real(wp), allocatable :: vQa(:,:), vQs(:,:), vGg(:,:)
      logical  :: got_sil, got_neu, got_ion, got_gra
      logical  :: rok
      character(len=512) :: q_h5, q_dir
      ! Draine's size grid: A(KA) = 1e-8*10^(0.55+(KA-1)*0.05) cm,
      ! NSIZE=84 (3.548 A .. 5.012 um, 0.05-dex log spacing). A(30)=100 A lands
      ! exactly on a node, so the 100 A charge cutoff sits on a grid point.
      integer,  parameter :: NSIZE_BD = 84
      real(wp), parameter :: A0_BD    = 0.55_wp   ! log10(a / 1e-8 cm) at KA=1
      real(wp), parameter :: DLGA_BD  = 0.05_wp   ! dex per step

      ! Lambda grid from the Q-table. The aeff grid is Draine's analytic 84-pt
      ! log grid built below -- NOT the size-dist file, whose dn columns are
      ! unused here (dn/da comes from grain_dist_dl07, the WD01 analytic model).
      if (present(status)) status = 0
      ! The two places stored optics can come from, resolved once.  When
      ! qtable_path names an HDF5 product, that product is also where the
      ! optics come from: one file, one set of numbers.
      q_h5 = ''
      if (is_hdf5_path(qtable_path)) q_h5 = qtable_path
      q_dir = sedust_dir(trim(sed_get_data_root()), 'dl07')
      if (present(stored_q_dir)) then
         q_dir = stored_q_dir
         ! Blank means no stored optics AT ALL, the product included.
         if (len_trim(stored_q_dir) == 0) q_h5 = ''
      end if
      if (present(lam_axis) .or. is_hdf5_path(qtable_path)) then
         ! The axis is given outright, or read from the product that carries
         ! it; either way no text Q table is opened.
         if (present(lam_axis)) then
            allocate(lam_grid(size(lam_axis)))
            lam_grid = lam_axis
         else
            call read_sedust_grid(qtable_path, euv_asked(include_euv), lam_grid, i, rok)
            if (.not. rok) then
               if (present(status)) then
                  status = 1;  return
               else
                  write(*,'(a,a)') ' sed_init_dl07: cannot read /grid/lambda from ', &
                                   trim(qtable_path)
                  stop 1
               end if
            end if
         end if
         ! lam_min still applies.  The product's axis is where this model's grid
         ! STARTS, not a ceiling on what the caller may ask for: a host that
         ! wants wavelengths shortward of the first node gets them prepended
         ! here, exactly as on the text route below.  Skipping this dropped
         ! lam_min on the floor whenever the axis came from HDF5.
         call move_alloc(lam_grid, lam_base)
         call euv_extended_lambda_grid(lam_grid, lam_min, base=lam_base)
         deallocate(lam_base)
         ! Everything shortward of the Q table's own 0.0912 um end is EUV here,
         ! and the coverage check below must see it as such.  Counted after the
         ! extension, since every point it prepends is shortward of that end.
         n_lam_euv = count(lam_grid < LAM_LYMAN_UM)
      else if (present(status)) then
         call load_q_table(qtable_path, ok=rok)
         if (.not. rok) then;  status = 1;  return;  end if
      else
         call load_q_table(qtable_path)
      end if
      ! Only the TEXT route needs the grid built here: the two branches above
      ! already have the axis, and neither loaded a Q table for this to read.
      if (.not. present(lam_axis) .and. .not. is_hdf5_path(qtable_path)) &
         call euv_extended_lambda_grid(lam_grid, lam_min, n_extra=n_lam_euv)

      ! This model's optics are Mie on the D03 dielectric functions throughout,
      ! so an EUV extension is a grid extension only -- as long as the grid
      ! stays inside those functions' own coverage. Past it `interp` freezes
      ! (n, k) at the boundary value, which would pass a constant index off as
      ! physics, so refuse instead. Silicate and graphite are both required.
      if (n_lam_euv > 0) then
         call silicate_index_lambda_range(sil_lam_lo, sil_lam_hi)
         call graphite_index_lambda_range(gra_lam_lo, gra_lam_hi)
         d03_lam_lo = max(sil_lam_lo, gra_lam_lo)
         if (lam_grid(1) < d03_lam_lo) then
            write(*,'(a,es10.3,a)') ' sed_init_dl07: lam_min =', lam_grid(1), &
               ' um is shorter than the D03 dielectric functions,'
            write(*,'(a,es10.3,a)') '                which stop at', d03_lam_lo, &
               ' um; (n, k) below it would be frozen at the boundary value.'
            if (present(status)) then
               status = 7;  deallocate(lam_grid);  return
            else
               stop 1
            end if
         end if
      end if

      NLAM = size(lam_grid)
      NA   = NSIZE_BD
      NT   = NT_in

      call free_shared_model_arrays()
      ! DL07 has no polarized optics; drop anything a previous astrodust
      ! init left behind rather than leave stale arrays on the wrong grid.
      polarized_optics_loaded = .false.
      if (allocated(Cpol)) deallocate(Cpol, Cpol_ext, Cbir_ext, falign_ad)
      ! Drop a stale astrodust asymmetry table; the DL07 scattering optics
      ! below are computed on this model's own size/wavelength grids.
      if (allocated(gsca_ad))  deallocate(gsca_ad)
      if (allocated(Csca_car)) deallocate(Csca_car, gsca_car)
      allocate(lam(NLAM), aeff(NA), T_first(NT), dn_ad(NA))
      allocate(Cabs(NLAM, NA), Csca(NLAM, NA), kappB_first(NT, NA), &
               H_first(NT, NA, NSTAGE), kappCMB(NA))
      allocate(gsca_ad(NLAM, NA), Csca_car(NLAM, NA), gsca_car(NLAM, NA))
      allocate(log_T_first(NT), log_H_first(NT, NA, NSTAGE), &
               log_kappB_first(NT, NA))
      allocate(dn_pah(NA), Cabs_pah(NLAM, NA), kappB_pah_first(NT, NA), &
               H_pah_first(NT, NA), kappCMB_pah(NA))
      allocate(log_H_pah_first(NT, NA), log_kappB_pah_first(NT, NA))
      if (allocated(Cabs_cneu)) deallocate(Cabs_cneu, Cabs_cion, dn_cneu, dn_cion)
      allocate(Cabs_cneu(NLAM, NA), Cabs_cion(NLAM, NA), dn_cneu(NA), dn_cion(NA))
      allocate(fion(NA), lna(NA))

      lam  = lam_grid
      deallocate(lam_grid)
      ! Draine's 84-pt log grid in microns (1e-8 cm = 1e-4 um).
      do ja = 1, NA
         aeff(ja) = 1.0e-4_wp * 10.0_wp**(A0_BD + real(ja-1,wp)*DLGA_BD)
      end do

      do i = 1, NT
         t = log(T_lo) + (log(T_hi)-log(T_lo))*real(i-1,wp)/real(NT-1,wp)
         T_first(i) = exp(t)
      end do
      log_T_first = log(T_first)
      call p_sub_setup(lam)

      do ja = 1, NA
         lna(ja) = log(aeff(ja))
      end do

      ! ---- Silicate population (dust slot) ----
      ! Mie on the D03 astrosilicate dielectric function keeps every return,
      ! so the scattering cross section and its asymmetry come from the same
      ! calculation as the absorption rather than from a tabulated albedo.
      call stored_q_on_model_grid(trim(q_h5), trim(q_dir), 'q_dl07_sil', 'sil', tQa, tQs, tGg, got_sil)
      do ja = 1, NA
         a_um = aeff(ja)
         geo  = PI * (a_um * UM2CM)**2
         do jw = 1, NLAM
            if (got_sil) then
               qabs1 = tQa(jw, ja);  qsca1 = tQs(jw, ja);  gsca1 = tGg(jw, ja)
            else
               call q_silicate_full(a_um, lam(jw), qext1, qsca1, qabs1, gsca1)
            end if
            Cabs(jw, ja)    = qabs1 * geo
            Csca(jw, ja)    = qsca1 * geo
            gsca_ad(jw, ja) = gsca1
         end do
         dn_ad(ja) = grain_dist_dl07(sd_index, 'sil', a_um) * bin_da(ja, lna)
      end do
      call build_kappB()
      log_kappB_first = log(max(kappB_first, tiny(0.0_wp)))
      call build_kappCMB()
      do ja = 1, NA
         do jt = 1, NT
            H_first(jt, ja, 1) = enthalpy_DL01(T_first(jt), aeff(ja), 'Sil ')
            H_first(jt, ja, 2) = H_first(jt, ja, 1)
         end do
      end do
      log_H_first = log(max(H_first, tiny(0.0_wp)))

      ! ---- Carbonaceous populations: neutral and cation kept SEPARATE ----
      ! Do NOT pre-blend the charge states: each is a distinct stochastically
      ! heated population (different absorption -> different T distribution),
      ! solved separately and summed in sed_solve_dl07 (matching DL07).
      do ja = 1, NA
         fion(ja) = pah_ionfrac(aeff(ja), u_isrf)   ! WD01b charging, direct
         ! Draine's charge cutoff: a > ~100 A treated as 100%
         ! ionized (single cation-only solve). aeff in microns.
         if (aeff(ja) > 0.99999e-2_wp) fion(ja) = 1.0_wp
      end do
      if (allocated(tQa)) deallocate(tQa, tQs, tGg)
      call stored_q_on_model_grid(trim(q_h5), trim(q_dir), 'q_dl07_pah_neu', 'pah_neu', tQa, tQs, tGg, got_neu)
      call stored_q_on_model_grid(trim(q_h5), trim(q_dir), 'q_dl07_pah_ion', 'pah_ion', uQa, uQs, uGg, got_ion)
      call stored_q_on_model_grid(trim(q_h5), trim(q_dir), 'q_dl07_gra', 'gra',     vQa, vQs, vGg, got_gra)
      if (sed_verbose) then
         if (got_sil .and. got_neu .and. got_ion .and. got_gra) then
            ! Name the source that was opened, not a directory written down
            ! here: this line was printed on the HDF5 route too, and named a
            ! directory that run had not touched.
            if (len_trim(q_h5) > 0) then
               write(*,'(a,a)') ' sed_init_dl07: optics read from ', trim(q_h5)
            else
               write(*,'(a,a)') ' sed_init_dl07: optics read from the stored' // &
                                ' tables under ', trim(q_dir)
            end if
         else
            write(*,'(a)') ' sed_init_dl07: optics solved from the dielectric' // &
                           ' functions (no stored table on this grid)'
         end if
      end if
      do ja = 1, NA
         a_um = aeff(ja)
         geo  = PI * (a_um * UM2CM)**2
         do jw = 1, NLAM
            if (got_neu) then
               Q_neu = tQa(jw, ja)
            else
               call qpah_abs(0, a_um, lam(jw), Q_neu)
            end if
            if (got_ion) then
               Q_ion = uQa(jw, ja)
            else
               call qpah_abs(1, a_um, lam(jw), Q_ion)
            end if
            Cabs_cneu(jw, ja) = Q_neu * geo
            Cabs_cion(jw, ja) = Q_ion * geo
            ! Scattering of the carbonaceous grain.  Charge shifts the PAH
            ! absorption features only, and a PAH is a molecule whose
            ! Rayleigh scattering is negligible, so both charge states
            ! scatter as the graphite sphere does -- random-orientation Mie
            ! (1/3 || + 2/3 perp) on the graphite dielectric function.
            if (got_gra) then
               qsca1 = vQs(jw, ja);  gsca1 = vGg(jw, ja)
            else
               call q_graphite_full(a_um, lam(jw), qext1, qsca1, qabs1, gsca1)
            end if
            Csca_car(jw, ja) = qsca1 * geo
            gsca_car(jw, ja) = gsca1
         end do
         ! full carbonaceous number per bin (graphite-split + PAH-split),
         ! partitioned into neutral / cation by the ionization fraction.
         da = ( grain_dist_dl07(sd_index, 'gra', a_um) &
              + grain_dist_dl07(sd_index, 'pah', a_um) ) * bin_da(ja, lna)
         dn_cneu(ja) = (1.0_wp - fion(ja)) * da
         dn_cion(ja) =          fion(ja)   * da
      end do
      ! Shared 'Car0' DL01 carbonaceous enthalpy for both charge states.
      do ja = 1, NA
         do jt = 1, NT
            H_pah_first(jt, ja) = enthalpy_DL01(T_first(jt), aeff(ja), 'Car0')
         end do
      end do
      log_H_pah_first = log(max(H_pah_first, tiny(0.0_wp)))
      ! kappB_pah_first / kappCMB_pah are built per charge state in
      ! sed_solve_dl07 (Cabs_pah is reused there as scratch).

      deallocate(fion, lna)
      initialized = .true.

   contains
      ! Log-spaced size-bin width da_i = a_i * dln(a), returned in CM
      ! (grain_dist_dl07 gives dn/da per cm of radius), central differences
      ! with one-sided ends. UM2CM converts the micron aeff axis to cm.
      pure function bin_da(j, lna_arr) result(da_out)
         integer,  intent(in) :: j
         real(wp), intent(in) :: lna_arr(:)
         real(wp) :: da_out
         ! Trapezoidal-in-log weights (following Draine's method):
         ! interior bins get the full central-difference dln(a); the two
         ! endpoints get half. dn/da is per cm of radius, so convert um->cm.
         if (j == 1) then
            da_out = aeff(j) * 0.5_wp * (lna_arr(2) - lna_arr(1))
         else if (j == NA) then
            da_out = aeff(j) * 0.5_wp * (lna_arr(NA) - lna_arr(NA-1))
         else
            da_out = aeff(j) * 0.5_wp * (lna_arr(j+1) - lna_arr(j-1))
         end if
         da_out = da_out * UM2CM       ! micron -> cm (dn/da is per cm)
      end function bin_da
   end subroutine sed_init_dl07


   subroutine sed_solve_dl07(J_lam, lamI_total, lamI_sil, lamI_carb)
      ! Total DL07 SED = silicate + carbonaceous, each via the shared
      ! grain loop (GD stochastic solver). Outputs lambda*I_lambda / N_H
      ! [erg s^-1 sr^-1 H^-1], same convention as sed_solve.
      real(wp), intent(in)  :: J_lam(:)
      real(wp), intent(out) :: lamI_total(:), lamI_sil(:), lamI_carb(:)
      real(wp), allocatable :: Jout_s(:), Jout_c(:), Jout_q(:)
      integer :: icharge

      if (.not. initialized) then
         write(*,'(a)') 'sed_solve_dl07: call sed_init_dl07 first'
         stop 1
      end if

      allocate(Jout_s(NLAM), Jout_c(NLAM), Jout_q(NLAM))

      ! Silicate population
      call sed_grain_loop(NA, dn_ad, aeff, Cabs, kappB_first, H_first(:,:,1), &
                          log_H_first(:,:,1), log_kappB_first, kappCMB, &
                          J_lam, 'sil', Jout_s)

      ! Carbonaceous: solve neutral and cation as SEPARATE stochastic
      ! populations (different absorption -> different T distribution) and
      ! sum. Cabs_pah/kappB_pah_first/kappCMB_pah are reused as scratch.
      Jout_c = 0.0_wp
      do icharge = 0, 1
         if (icharge == 0) then
            Cabs_pah = Cabs_cneu;  dn_pah = dn_cneu
         else
            Cabs_pah = Cabs_cion;  dn_pah = dn_cion
         end if
         call build_kappB_pah()                 ! reads Cabs_pah -> kappB_pah_first
         log_kappB_pah_first = log(max(kappB_pah_first, tiny(0.0_wp)))
         call build_kappCMB_pah()               ! reads Cabs_pah -> kappCMB_pah
         call sed_grain_loop(NA, dn_pah, aeff, Cabs_pah, kappB_pah_first, &
                             H_pah_first, log_H_pah_first, &
                             log_kappB_pah_first, kappCMB_pah, &
                             J_lam, 'pah', Jout_q)
         Jout_c = Jout_c + Jout_q
      end do

      if (use_induced_emission) then
         call apply_induced_factor(J_lam, Jout_s)
         call apply_induced_factor(J_lam, Jout_c)
      end if

      lamI_sil   = lam * Jout_s * 1.0e-3_wp
      lamI_carb  = lam * Jout_c * 1.0e-3_wp
      lamI_total = lamI_sil + lamI_carb

      deallocate(Jout_s, Jout_c, Jout_q)
   end subroutine sed_solve_dl07


   ! =====================================================================
   ! Shared grain-loop implementation for both astrodust and PAH
   ! populations. Dispatches to 'draine' iterative or 'heuristic'
   ! look-ahead T-window narrowing based on the module variable
   ! stoch_method.
   ! =====================================================================
   ! Polarized emission (optional): pass Cpol_pop, falign_pop and Jpol
   ! together to get the intrinsic polarized emissivity alongside Jout. The
   ! polarized accumulator is driven by exactly the same temperature weights
   ! as Jout -- every site that adds Cabs*B(T) adds Cpol*f_align*B(T) in the
   ! same breath -- so the P(T) solvers are untouched by this option.
   ! What comes back is the INTRINSIC rate: the geometric sin^2(gamma)
   ! projection and any turbulent depolarization belong to the radiative
   ! transfer, not here.
   subroutine sed_grain_loop(npop, dn_pop, aeff_pop, Cabs_pop, kappB_pop, H_pop, &
                              log_H_pop, log_kappB_pop, kappCMB_pop, &
                              J_lam, grain_type, Jout, Cpol_pop, falign_pop, Jpol)
      integer,          intent(in)  :: npop
      real(wp),         intent(in)  :: dn_pop(:)          ! (npop)
      real(wp),         intent(in)  :: aeff_pop(:)        ! (npop) [um] radii of this population
      real(wp),         intent(in)  :: Cabs_pop(:,:)      ! (NLAM, npop)
      real(wp),         intent(in)  :: kappB_pop(:,:)     ! (NT, npop)
      real(wp),         intent(in)  :: H_pop(:,:)         ! (NT, npop)
      real(wp),         intent(in)  :: log_H_pop(:,:)     ! (NT, npop)
      real(wp),         intent(in)  :: log_kappB_pop(:,:) ! (NT, npop)
      real(wp),         intent(in)  :: kappCMB_pop(:)     ! (npop)
      real(wp),         intent(in)  :: J_lam(:)           ! (NLAM)
      character(len=*), intent(in)  :: grain_type         ! 'sil' or 'pah'
      real(wp),         intent(out) :: Jout(:)            ! (NLAM)
      real(wp), optional, intent(in)  :: Cpol_pop(:,:)    ! (NLAM, npop) [cm^2]
      real(wp), optional, intent(in)  :: falign_pop(:)    ! (npop)
      real(wp), optional, intent(out) :: Jpol(:)          ! (NLAM)

      integer  :: ir, ii, loc1, loc2, iguard, n_guard_resolve
      integer  :: n_stoch, n_equil_eeq, n_equil_fail
      real(wp) :: Teq, EEQ, del, Tmin_n, Tmax_n, a_cm_qm, wpol
      real(wp), allocatable :: spec(:), P(:), lnP(:)
      real(wp), allocatable :: T(:), H(:), kappB(:)
      real(wp), allocatable :: Jout_local(:), emission_qm(:), Jpol_local(:)
      logical :: Equil, Equil_prev, converged, qm_ok, do_pol

      Jout = 0.0_wp
      ! Polarization is opt-in: all three arguments must be supplied together.
      do_pol = present(Cpol_pop) .and. present(falign_pop) .and. present(Jpol)
      if (present(Jpol)) Jpol = 0.0_wp

      select case (trim(stoch_method))

      case ('draine')
         ! ----- Grain-by-grain Draine iterative narrowing -----------------
         ! Each grain decides equilibrium independently (EEQ threshold),
         ! so the loop is OpenMP-parallelizable.
         n_stoch = 0; n_equil_eeq = 0; n_equil_fail = 0
         !$omp parallel default(none) &
         !$omp&   shared(npop, dn_pop, Cabs_pop, kappB_pop, H_pop, &
         !$omp&          log_H_pop, log_kappB_pop, kappCMB_pop, J_lam, &
         !$omp&          T_first, lam, Jout, NLAM, NT, &
         !$omp&          do_pol, Cpol_pop, falign_pop, Jpol) &
         !$omp&   private(ir, ii, Teq, EEQ, Equil, converged, wpol, &
         !$omp&           spec, P, lnP, T, H, kappB, Jout_local, Jpol_local) &
         !$omp&   reduction(+:n_stoch, n_equil_eeq, n_equil_fail)
         allocate(spec(NLAM), P(NT), lnP(NT), T(NT), H(NT), kappB(NT))
         allocate(Jout_local(NLAM), Jpol_local(NLAM))
         Jout_local = 0.0_wp
         Jpol_local = 0.0_wp
         wpol       = 0.0_wp
         !$omp do schedule(dynamic)
         do ir = 1, npop
            if (dn_pop(ir) <= 0.0_wp) cycle
            call calc_Teq(lam, Cabs_pop(:, ir), J_lam, T_first, &
                          kappB_pop(:, ir), Teq)
            call interp(T_first, H_pop(:, ir), Teq, EEQ)

            if (EEQ >= EEQSS_ERG) then
               Equil = .true.; n_equil_eeq = n_equil_eeq + 1
            else
               Equil = .false.
            end if

            if (.not. Equil) then
               call narrow_iterative(H_pop(:, ir), log_H_pop(:, ir), &
                                     kappB_pop(:, ir), log_kappB_pop(:, ir), &
                                     kappCMB_pop(ir), Cabs_pop(:, ir), J_lam, &
                                     Teq, EEQ, T, P, converged)
               if (.not. converged) then
                  Equil = .true.; n_equil_fail = n_equil_fail + 1
               else
                  n_stoch = n_stoch + 1
               end if
            end if

            if (do_pol) wpol = dn_pop(ir) * falign_pop(ir)

            if (Equil) then
               call calc_bbody(Teq, lam, spec)
               Jout_local = Jout_local + dn_pop(ir) * Cabs_pop(:, ir) * spec
               if (do_pol) Jpol_local = Jpol_local + wpol * Cpol_pop(:, ir) * spec
            else
               do ii = 1, NT
                  if (P(ii) > 0.0_wp) then
                     call calc_bbody(T(ii), lam, spec)
                     Jout_local = Jout_local + dn_pop(ir) * P(ii) * Cabs_pop(:, ir) * spec
                     if (do_pol) Jpol_local = Jpol_local + wpol * P(ii) * Cpol_pop(:, ir) * spec
                  end if
               end do
            end if
         end do
         !$omp end do
         !$omp critical
         Jout = Jout + Jout_local
         if (do_pol) Jpol = Jpol + Jpol_local
         !$omp end critical
         deallocate(spec, P, lnP, T, H, kappB, Jout_local, Jpol_local)
         !$omp end parallel
         if (sed_verbose) write(*,'(a,i4,a,i4,a,i4,a)') &
            '   [Draine narrowing: stoch=', n_stoch, ' eeq_gate=', n_equil_eeq, &
            ' fail_to_eq=', n_equil_fail, ']'

      case ('heuristic')
         ! ----- Heuristic look-ahead narrowing ---------------------------
         ! Serial: the window WIDTH (del) is inherited from the previous
         ! grain, but the window TOP is set per grain from the physical
         ! hot-tail bound (Draine's UMAX convention, vsg_td_emission_v7
         ! line 646):  T_max = H^{-1}(13.6 eV + 2 H(Teq)).
         ! No excursion can exceed "hardest single photon on top of twice
         ! the equilibrium enthalpy", so this places the top correctly
         ! without the draine-variant's iterative re-solves. The original
         ! look-ahead top (Teq*exp(+NARROW_FAC*del), inherited from the
         ! previous grain) could land inside the populated tail; calc_P's
         ! highest-bin correction then piles the clipped flux into the top
         ! bin, producing the spurious NIR hot-tail emission (+70% NIR for
         ! astrodust S1 vs the dbdis reference).
         ! A one-shot P(top) guard (re-solve with a raised top) backs up
         ! the analytic bound; it fires rarely.
         Equil_prev = .false.
         wpol       = 0.0_wp
         n_stoch = 0; n_guard_resolve = 0
         allocate(spec(NLAM), P(NT), lnP(NT), T(NT), H(NT), kappB(NT))
         del   = log(T_first(NT) / T_first(1))

         do ir = 1, npop
            if (dn_pop(ir) <= 0.0_wp) cycle
            call calc_Teq(lam, Cabs_pop(:, ir), J_lam, T_first, &
                          kappB_pop(:, ir), Teq)

            Equil = .false.
            if (Equil_prev) Equil = .true.

            if (.not. Equil) then
               ! --- window for THIS grain ---
               call interp(T_first, H_pop(:, ir), Teq, EEQ)
               ! The bound needs the hardest single photon the FIELD can
                ! The bound needs the hardest single photon the FIELD can
               ! deliver, hardest_photon_energy(lam, J_lam) -- NOT hc/lam(1),
               ! which is only the short end of the model's optics grid.  The
               ! two coincide for astrodust and DL07, whose Q tables stop at
               ! the Lyman limit, and there hc/0.0912 um = 13.595 eV < 13.6 eV
               ! so the max() returns U_UV1_ERG unchanged.  They diverge for
               ! Zubko/ZDA, whose ZDA optics tables start at 1.0e-3 um (1.24 keV)
               ! while the illuminating field still stops at the Lyman limit:
               ! reading the grid there would raise the top of the enthalpy bin
               ! set by a factor 91 with no photon behind it.  U_UV1_ERG
               ! (13.6 eV) is the correct bound only when the field itself
               ! stops at the Lyman limit; a field carried into the EUV raises
               ! it, and a top set from 13.6 eV would clip the excursions its
               ! hardest photons drive.
               call U_to_T(max(U_UV1_ERG, hardest_photon_energy(lam, J_lam)) &
                           + 2.0_wp*EEQ, &
                           H_pop(:, ir), log_H_pop(:, ir), Tmax_n)
               ! Pad the analytic top by one guard step (e^0.5 in T): the
               ! multi-photon tail at U ~ a few extends slightly past the
               ! single-photon bound at the 1e-12 level, which otherwise
               ! triggers the guard re-solve for ~40% of the grains
               ! (measured: 108 extra calc_P calls out of 271 -> +25%
               ! wall time). Padding costs only ~10% coarser bins over
               ! the window (irrelevant at NT=200; cf. the bin-count
               ! convergence test) and makes the guard fire rarely.
               Tmax_n = min(max(Tmax_n * 1.6487_wp, Teq * exp(0.05_wp)), &
                            T_first(NT))
               Tmin_n = max(Teq * exp(-NARROW_FAC * del), T_first(1))
               if (Tmin_n >= Tmax_n) Tmin_n = max(Tmax_n * exp(-0.5_wp), T_first(1))

               do iguard = 1, 3
                  call narrow_T_window(ir, log_H_pop, log_kappB_pop, &
                                       Tmin_n, Tmax_n, T, H, kappB)
                  call calc_P(lam, Cabs_pop(:, ir), J_lam, T, kappB, H, &
                              P, lnP, kappCMB_pop(ir))
                  ! P(top) guard: if the top bin is still populated above
                  ! the tail threshold, the window clipped the hot tail --
                  ! raise the top and re-solve (rare).
                  if (P(NT) <= 1.0e-12_wp * maxval(P) .or. &
                      Tmax_n >= T_first(NT)) exit
                  Tmax_n = min(Tmax_n * exp(0.5_wp), T_first(NT))
                  n_guard_resolve = n_guard_resolve + 1
               end do
               n_stoch = n_stoch + 1

               loc1 = first_location(lnP > lnP_crit)
               loc2 = last_location (lnP > lnP_crit)
               if (loc1 == 0 .or. loc2 == 0) then
                  Equil = .true.
               else
                  if (T(loc2) < Teq) Equil = .true.
                  if (ir < npop) then
                     del = max(log(T(loc2)/Teq), log(Teq/T(loc1))) * 2.0_wp
                     if (del < dlnT_crit) Equil = .true.
                  end if
               end if
            end if

            if (do_pol) wpol = dn_pop(ir) * falign_pop(ir)

            if (Equil) then
               call calc_bbody(Teq, lam, spec)
               Jout = Jout + dn_pop(ir) * Cabs_pop(:, ir) * spec
               if (do_pol) Jpol = Jpol + wpol * Cpol_pop(:, ir) * spec
            else
               do ii = 1, NT
                  if (P(ii) > 0.0_wp) then
                     call calc_bbody(T(ii), lam, spec)
                     if (gd_photon_cutoff) then
                        ! dbdis photon cutoff: a grain in bin ii cannot
                        ! emit a photon more energetic than its enthalpy.
                        where (lam < HC_ERG_UM / H(ii)) spec = 0.0_wp
                     end if
                     Jout = Jout + dn_pop(ir) * P(ii) * Cabs_pop(:, ir) * spec
                     if (do_pol) Jpol = Jpol + wpol * P(ii) * Cpol_pop(:, ir) * spec
                  end if
               end do
            end if
            Equil_prev = Equil
         end do
         deallocate(spec, P, lnP, T, H, kappB)
         if (sed_verbose) write(*,'(a,i4,a,i4,a)') '   [heuristic: stoch=', n_stoch, &
            ' guard_resolves=', n_guard_resolve, ']'

      case ('qm')
         ! ----- Quantum-mechanical (dbdis) solver -----------------------
         ! Each grain decides equilibrium independently (EEQ threshold),
         ! so the loop is OpenMP-parallelizable.
         n_stoch = 0; n_equil_eeq = 0; n_equil_fail = 0
         !$omp parallel default(none) &
         !$omp&   shared(npop, dn_pop, aeff_pop, Cabs_pop, kappB_pop, H_pop, &
         !$omp&          log_H_pop, log_kappB_pop, kappCMB_pop, J_lam, &
         !$omp&          T_first, lam, grain_type, Jout, NLAM, NT, &
         !$omp&          do_pol, Cpol_pop, falign_pop, Jpol) &
         !$omp&   private(ir, ii, Teq, EEQ, Equil, converged, a_cm_qm, qm_ok, wpol, &
         !$omp&           spec, P, lnP, T, H, kappB, Jout_local, emission_qm, Jpol_local) &
         !$omp&   reduction(+:n_stoch, n_equil_eeq, n_equil_fail)
         allocate(spec(NLAM), P(NT), lnP(NT), T(NT), H(NT), kappB(NT))
         allocate(Jout_local(NLAM), emission_qm(NLAM), Jpol_local(NLAM))
         Jout_local = 0.0_wp
         Jpol_local = 0.0_wp
         wpol       = 0.0_wp
         !$omp do schedule(dynamic)
         do ir = 1, npop
            block
               if (dn_pop(ir) <= 0.0_wp) cycle
               call calc_Teq(lam, Cabs_pop(:, ir), J_lam, T_first, &
                             kappB_pop(:, ir), Teq)
               call interp(T_first, H_pop(:, ir), Teq, EEQ)

               if (EEQ >= EEQSS_ERG) then
                  Equil = .true.; n_equil_eeq = n_equil_eeq + 1
               else
                  Equil = .false.
               end if

               if (do_pol) wpol = dn_pop(ir) * falign_pop(ir)

               if (.not. Equil) then
                  a_cm_qm = aeff_pop(ir) * UM2CM
                  call qm_solve_grain(NLAM, lam, Cabs_pop(:,ir), J_lam, &
                                      NT, T_first, H_pop(:,ir), &
                                      Teq, EEQ, EEQSS_ERG, &
                                      a_cm_qm, grain_type, &
                                      emission_qm, qm_ok)
                  if (qm_ok) then
                     do ii = 1, NLAM
                        Jout_local(ii) = Jout_local(ii) + dn_pop(ir) * emission_qm(ii) / &
                                         (4.0_wp * PI * lam(ii) * 1.0e-3_wp)
                     end do
                     ! qm_solve_grain returns the emitted spectrum with Cabs
                     ! already folded in, so the polarized counterpart is the
                     ! same spectrum rescaled by Cpol/Cabs. Wavelengths with
                     ! zero Cabs emit nothing and contribute nothing.
                     if (do_pol) then
                        do ii = 1, NLAM
                           if (Cabs_pop(ii, ir) > 0.0_wp) then
                              Jpol_local(ii) = Jpol_local(ii) + &
                                 wpol * emission_qm(ii) * &
                                 (Cpol_pop(ii, ir) / Cabs_pop(ii, ir)) / &
                                 (4.0_wp * PI * lam(ii) * 1.0e-3_wp)
                           end if
                        end do
                     end if
                     n_stoch = n_stoch + 1
                  else
                     ! QM failed: fall back to GD for this grain only.
                     call narrow_iterative(H_pop(:, ir), log_H_pop(:, ir), &
                                           kappB_pop(:, ir), log_kappB_pop(:, ir), &
                                           kappCMB_pop(ir), Cabs_pop(:, ir), J_lam, &
                                           Teq, EEQ, T, P, converged)
                     if (converged) then
                        do ii = 1, NT
                           if (P(ii) > 0.0_wp) then
                              call calc_bbody(T(ii), lam, spec)
                              Jout_local = Jout_local + dn_pop(ir) * P(ii) * Cabs_pop(:, ir) * spec
                              if (do_pol) Jpol_local = Jpol_local + &
                                 wpol * P(ii) * Cpol_pop(:, ir) * spec
                           end if
                        end do
                     else
                        Equil = .true.
                     end if
                     n_equil_fail = n_equil_fail + 1
                  end if
               end if

               if (Equil) then
                  call calc_bbody(Teq, lam, spec)
                  Jout_local = Jout_local + dn_pop(ir) * Cabs_pop(:, ir) * spec
                  if (do_pol) Jpol_local = Jpol_local + wpol * Cpol_pop(:, ir) * spec
               end if
            end block
         end do
         !$omp end do
         !$omp critical
         Jout = Jout + Jout_local
         if (do_pol) Jpol = Jpol + Jpol_local
         !$omp end critical
         deallocate(spec, P, lnP, T, H, kappB, Jout_local, Jpol_local)
         !$omp end parallel
         if (sed_verbose) write(*,'(a,i4,a,i4,a,i4,a)') &
            '   [QM solver: stoch=', n_stoch, ' eeq_gate=', n_equil_eeq, &
            ' fail_to_eq=', n_equil_fail, ']'

      case ('equil')
         ! Force EQUILIBRIUM for ALL grains (no stochastic heating): each
         ! grain emits B_lam(T_eq)*Cabs at its OWN size/composition-dependent
         ! equilibrium temperature, computed from that grain's absorption
         ! cross section. (Option 1: per-(type,size) equilibrium temperature.)
         allocate(spec(NLAM))
         n_equil_eeq = 0
         do ir = 1, npop
            if (dn_pop(ir) <= 0.0_wp) cycle
            call calc_Teq(lam, Cabs_pop(:, ir), J_lam, T_first, kappB_pop(:, ir), Teq)
            call calc_bbody(Teq, lam, spec)
            Jout = Jout + dn_pop(ir) * Cabs_pop(:, ir) * spec
            if (do_pol) Jpol = Jpol + dn_pop(ir) * falign_pop(ir) * Cpol_pop(:, ir) * spec
            n_equil_eeq = n_equil_eeq + 1
         end do
         deallocate(spec)
         if (sed_verbose) write(*,'(a,i4,a)') &
            '   [equil (single-grain Teq): ', n_equil_eeq, ' grains]'

      case default
         write(*,'(a,a)') 'sed_grain_loop: unknown stoch_method: ', trim(stoch_method)
         stop 1

      end select
   end subroutine sed_grain_loop


   ! =====================================================================
   ! Batch QM solver: process all 3 grain types (S1, S2, PAH)
   ! in ONE OpenMP parallel region.  Total work = 3 × NA grains ≈ 486,
   ! with ~210 stochastic grains, keeping all threads busy.
   ! =====================================================================
   subroutine sed_solve_qm_batch(J_lam, lamI_stages, lamI_pah)
      real(wp), intent(in)  :: J_lam(:)            ! (NLAM)
      real(wp), intent(out) :: lamI_stages(:,:)    ! (NLAM, 2)
      real(wp), intent(out) :: lamI_pah(:)         ! (NLAM)

      integer  :: total_grains, iw, itype, ir, ii, out_idx
      integer  :: n_stoch, n_equil
      real(wp) :: Teq, EEQ, a_cm_qm, dn_ir, kCMBg
      logical  :: Equil, qm_ok, converged
      character(len=3) :: gtype

      ! Thread-private arrays
      real(wp), allocatable :: spec(:), P(:), lnP(:), T(:), H_w(:), kappB_w(:)
      real(wp), allocatable :: Jout_local(:,:)    ! (NLAM, 3)
      real(wp), allocatable :: emission_qm(:), Cabs_g(:), kappB_g(:)
      real(wp), allocatable :: Hg(:), logHg(:), logkBg(:)
      real(wp) :: Jout_all(NLAM, 3)

      if (.not. initialized) then
         write(*,'(a)') 'sed_solve_qm_batch: call sed_init first'; stop 1
      end if

      ! 4 grain "types": 1=S1, 2=S2 (astrodust), 3=PAH-neutral, 4=PAH-cation.
      ! Neutral and cation are distinct stochastic populations (nonlinear in
      ! Cabs); both accumulate into the PAH output slot (out_idx=3).
      total_grains = 4 * NA
      Jout_all = 0.0_wp
      n_stoch = 0; n_equil = 0

      !$omp parallel default(none) &
      !$omp&   shared(total_grains, NA, NLAM, NT, J_lam, Jout_all, &
      !$omp&          lam, aeff, T_first, &
      !$omp&          dn_ad, Cabs, kappB_first, H_first, &
      !$omp&          log_H_first, log_kappB_first, kappCMB, &
      !$omp&          dn_cneu, dn_cion, Cabs_cneu, Cabs_cion, H_pah_first, &
      !$omp&          kappB_cneu, kappB_cion, log_H_pah_first, &
      !$omp&          log_kappB_cneu, log_kappB_cion, kappCMB_cneu, kappCMB_cion) &
      !$omp&   private(iw, itype, ir, ii, Teq, EEQ, Equil, converged, &
      !$omp&           a_cm_qm, dn_ir, gtype, qm_ok, out_idx, kCMBg, &
      !$omp&           spec, P, lnP, T, H_w, kappB_w, Jout_local, &
      !$omp&           emission_qm, Cabs_g, kappB_g, Hg, logHg, logkBg) &
      !$omp&   reduction(+:n_stoch, n_equil)
      allocate(spec(NLAM), P(NT), lnP(NT), T(NT), H_w(NT), kappB_w(NT))
      allocate(Jout_local(NLAM, 3))
      allocate(emission_qm(NLAM), Cabs_g(NLAM), kappB_g(NT))
      allocate(Hg(NT), logHg(NT), logkBg(NT))
      Jout_local = 0.0_wp

      !$omp do schedule(dynamic)
      do iw = 1, total_grains
         block
            itype = (iw - 1) / NA + 1    ! 1=S1, 2=S2, 3=PAH-neutral, 4=PAH-cation
            ir    = mod(iw - 1, NA) + 1   ! grain index within type

            ! Select the single-grain cross section, enthalpy, kappa and output slot.
            ! Neutral (itype=3) and cation (itype=4) are separate populations,
            ! both summed into the PAH output (out_idx=3).
            select case (itype)
            case (1, 2)
               dn_ir = dn_ad(ir);   gtype = 'sil';  out_idx = itype
               Cabs_g = Cabs(:,ir);          kappB_g = kappB_first(:,ir)
               Hg     = H_first(:,ir,itype); logHg   = log_H_first(:,ir,itype)
               logkBg = log_kappB_first(:,ir); kCMBg = kappCMB(ir)
            case (3)
               dn_ir = dn_cneu(ir);  gtype = 'pah';  out_idx = 3
               Cabs_g = Cabs_cneu(:,ir);  kappB_g = kappB_cneu(:,ir)
               Hg     = H_pah_first(:,ir);  logHg = log_H_pah_first(:,ir)
               logkBg = log_kappB_cneu(:,ir);  kCMBg = kappCMB_cneu(ir)
            case default   ! 4 = cation
               dn_ir = dn_cion(ir);  gtype = 'pah';  out_idx = 3
               Cabs_g = Cabs_cion(:,ir);  kappB_g = kappB_cion(:,ir)
               Hg     = H_pah_first(:,ir);  logHg = log_H_pah_first(:,ir)
               logkBg = log_kappB_cion(:,ir);  kCMBg = kappCMB_cion(ir)
            end select
            if (dn_ir <= 0.0_wp) cycle

            call calc_Teq(lam, Cabs_g, J_lam, T_first, kappB_g, Teq)
            call interp(T_first, Hg, Teq, EEQ)

            if (EEQ >= EEQSS_ERG) then
               Equil = .true.; n_equil = n_equil + 1
            else
               Equil = .false.
            end if

            if (.not. Equil) then
               a_cm_qm = aeff(ir) * UM2CM
               call qm_solve_grain(NLAM, lam, Cabs_g, J_lam, &
                                   NT, T_first, Hg, &
                                   Teq, EEQ, EEQSS_ERG, &
                                   a_cm_qm, gtype, emission_qm, qm_ok)

               if (qm_ok) then
                  do ii = 1, NLAM
                     Jout_local(ii, out_idx) = Jout_local(ii, out_idx) + &
                        dn_ir * emission_qm(ii) / (4.0_wp * PI * lam(ii) * 1.0e-3_wp)
                  end do
                  n_stoch = n_stoch + 1
               else
                  ! QM failed: fall back to GD
                  call narrow_iterative(Hg, logHg, kappB_g, logkBg, &
                                        kCMBg, Cabs_g, J_lam, &
                                        Teq, EEQ, T, P, converged)
                  if (converged) then
                     do ii = 1, NT
                        if (P(ii) > 0.0_wp) then
                           call calc_bbody(T(ii), lam, spec)
                           Jout_local(:,out_idx) = Jout_local(:,out_idx) + &
                              dn_ir * P(ii) * Cabs_g * spec
                        end if
                     end do
                  else
                     Equil = .true.
                  end if
               end if
            end if

            if (Equil) then
               call calc_bbody(Teq, lam, spec)
               Jout_local(:,out_idx) = Jout_local(:,out_idx) + dn_ir * Cabs_g * spec
            end if
         end block
      end do
      !$omp end do

      !$omp critical
      Jout_all = Jout_all + Jout_local
      !$omp end critical
      deallocate(spec, P, lnP, T, H_w, kappB_w, Jout_local)
      !$omp end parallel

      if (sed_verbose) write(*,'(a,i4,a,i4,a)') &
         '   [QM batch: stoch=', n_stoch, ' equil=', n_equil, ']'

      ! Unit conversion: Jout → lamI_lam
      do ii = 1, NSTAGE
         lamI_stages(:, ii) = lam * Jout_all(:, ii) * 1.0e-3_wp
      end do
      lamI_pah = lam * Jout_all(:, 3) * 1.0e-3_wp

      if (use_induced_emission) then
         do ii = 1, NSTAGE
            call apply_induced_factor(J_lam, lamI_stages(:, ii))
         end do
         call apply_induced_factor(J_lam, lamI_pah)
      end if
   end subroutine sed_solve_qm_batch


   subroutine narrow_T_window(ir, log_H_p, log_kappB_p, Tmin_n, Tmax_n, &
                              T_out, H_out, kappB_out)
      ! Heuristic narrowing helper. Build a log-spaced T grid of
      ! NT points in [Tmin_n, Tmax_n] and interpolate H_p(:, ir) and
      ! kappB_p(:, ir) (passed in log) onto it via log-log interpolation
      ! against log_T_first. Used by the stoch_method='heuristic'
      ! branch in sed_grain_loop.
      integer,  intent(in)  :: ir
      real(wp), intent(in)  :: log_H_p(:,:), log_kappB_p(:,:)
      real(wp), intent(in)  :: Tmin_n, Tmax_n
      real(wp), intent(out) :: T_out(:), H_out(:), kappB_out(:)
      integer  :: i
      real(wp) :: lT, lT_lo, lT_hi, lH, lk
      lT_lo = log(Tmin_n)
      lT_hi = log(Tmax_n)
      do i = 1, NT
         lT = lT_lo + (lT_hi - lT_lo) * real(i-1, wp) / real(NT-1, wp)
         T_out(i) = exp(lT)
         call interp(log_T_first, log_H_p    (:, ir), lT, lH)
         call interp(log_T_first, log_kappB_p(:, ir), lT, lk)
         H_out(i)     = exp(lH)
         kappB_out(i) = exp(lk)
      end do
   end subroutine narrow_T_window


   subroutine build_TgrigGHk(TMIN, TMAX, log_H_wide, log_kappB_wide, &
                             T_out, H_out, kappB_out)
      ! Helper: log-spaced T grid of NT points in [TMIN, TMAX], with
      ! H and kappB interpolated (log-log) from the wide tables.
      real(wp), intent(in)  :: TMIN, TMAX
      real(wp), intent(in)  :: log_H_wide(:), log_kappB_wide(:)
      real(wp), intent(out) :: T_out(:), H_out(:), kappB_out(:)
      integer  :: i
      real(wp) :: lT, lT_lo, lT_hi, lH, lk
      lT_lo = log(TMIN)
      lT_hi = log(TMAX)
      do i = 1, NT
         lT = lT_lo + (lT_hi - lT_lo) * real(i-1, wp) / real(NT-1, wp)
         T_out(i) = exp(lT)
         call interp(log_T_first, log_H_wide,     lT, lH)
         call interp(log_T_first, log_kappB_wide, lT, lk)
         H_out(i)     = exp(lH)
         kappB_out(i) = exp(lk)
      end do
   end subroutine build_TgrigGHk


   subroutine U_to_T(U, H_wide, log_H_wide, T_out)
      ! Invert enthalpy: T such that H_wide(T) = U.
      ! Uses log-log interp on (log_H_wide, log_T_first) against log(U).
      real(wp), intent(in)  :: U
      real(wp), intent(in)  :: H_wide(:), log_H_wide(:)
      real(wp), intent(out) :: T_out
      real(wp) :: lT

      if (U <= H_wide(1)) then
         T_out = T_first(1)
         return
      end if
      if (U >= H_wide(NT)) then
         T_out = T_first(NT)
         return
      end if
      call interp(log_H_wide, log_T_first, log(U), lT)
      T_out = exp(lT)
   end subroutine U_to_T


   subroutine narrow_iterative(H_wide, log_H_wide, kappB_wide, log_kappB_wide, &
                               kappCMB_r, Cabs_r, J_lam, &
                               Teq, EEQ, T_out, P_out, converged)
      ! Iterative T-window selection following Draine's method:
      !   - Initial guess from EEQ (UMAX = max(13.6eV + 2·EEQ, 13.65eV);
      !     UMIN = 0 or EEQ/5).
      !   - Build log-T grid in [TMIN, TMAX] = [H^{-1}(UMIN), H^{-1}(UMAX)],
      !     interpolate H, kappB, solve P (calc_P).
      !   - Adjust UMAX (UMIN) based on P(NT)/Pmax vs PMIN_UP (P(2)/Pmax
      !     vs PMIN_LO): shrink toward last bin where P > threshold (with
      !     0.8/0.2 damping), or expand ×1.2 / bisect to UMAXHI / UMINLO.
      !   - Iterate until no UMIN/UMAX change (converged) or MAX_ITER.
      ! Returns converged P(NT), T(NT). If degenerate (Pmax=0 or window
      ! collapses) converged = .false. -> caller uses equilibrium.
      real(wp), intent(in)  :: H_wide(:), log_H_wide(:)
      real(wp), intent(in)  :: kappB_wide(:), log_kappB_wide(:)
      real(wp), intent(in)  :: kappCMB_r, Cabs_r(:), J_lam(:)
      real(wp), intent(in)  :: Teq, EEQ
      real(wp), intent(out) :: T_out(:), P_out(:)
      logical,  intent(out) :: converged

      real(wp) :: UMIN, UMAX, UMINHI, UMINLO, UMAXHI, UMAXLO
      real(wp) :: TMIN, TMAX, Pmax
      real(wp), allocatable :: H(:), kappB(:), lnP(:), U(:)
      integer  :: i, iter, JCUT
      logical  :: refine
      real(wp), parameter :: BIG = 1.0e70_wp

      allocate(H(NT), kappB(NT), lnP(NT), U(NT))

      ! Initial guesses (Draine v7 lines 637-651).  The single-photon term is
      ! the hardest photon the FIELD carries, hardest_photon_energy(lam,J_lam),
      ! not hc/lam(1): the latter is the short end of the model's optics grid,
      ! which coincides with the field only when the grid stops where the
      ! illumination does.  It does for astrodust and DL07 (Lyman limit,
      ! 13.595 eV, so U_UV1_ERG stays selected); it does not for Zubko/ZDA,
      ! whose ZDA optics tables reach 1.0e-3 um and would hand back 1.24 keV for a
      ! field that carries nothing below 0.0912 um.  A field genuinely carried
      ! into the EUV does raise the bound, which is the point.  This is only
      ! the starting window -- the loop below expands UMAX until the tail is
      ! resolved -- but starting it away from the true single-photon bound
      ! wastes iterations, and starting it too high leaves the bins coarse
      ! because the contraction is guarded and iteration-capped.
      UMAX = max(max(U_UV1_ERG, hardest_photon_energy(lam, J_lam)) &
                 + 2.0_wp*EEQ, UMAXMIN_ERG)
      if (EEQ < 0.1_wp * EEQSS_ERG) then
         UMIN = 0.0_wp
      else
         UMIN = EEQ / 5.0_wp
      end if
      UMAXHI = BIG
      UMAXLO = 0.0_wp
      UMINHI = BIG
      UMINLO = 0.0_wp

      ! defensive; the loop assigns both before any reachable read
      Pmax  = 0.0_wp
      P_out = 0.0_wp
      converged = .false.
      do iter = 1, MAX_ITER_NARROW
         call U_to_T(UMIN, H_wide, log_H_wide, TMIN)
         call U_to_T(UMAX, H_wide, log_H_wide, TMAX)
         if (TMAX <= TMIN) exit                              ! degenerate

         call build_TgrigGHk(TMIN, TMAX, log_H_wide, log_kappB_wide, &
                             T_out, H, kappB)
         U = H

         call calc_P(lam, Cabs_r, J_lam, T_out, kappB, H, P_out, lnP, kappCMB_r)
         Pmax = maxval(P_out)
         if (Pmax <= 0.0_wp) exit                            ! degenerate

         refine = .false.

         ! ---- Adjust UMAX (Draine v7 lines 1276-1341) ----
         if (P_out(NT)/Pmax <= PMIN_UP .and. UMAX > UMAXMIN_ERG) then
            ! tail at top falls below threshold -> shrink UMAX
            UMAXHI = UMAX
            JCUT = NT
            do i = NT, 1, -1
               JCUT = i
               if (U(i) < UMAXMIN_ERG) exit
               if (P_out(i)/Pmax > PMIN_UP) exit
            end do
            if (UMAX > 1.02_wp*UMAXLO .and. UMAX > 1.01_wp*U(JCUT)) then
               UMAX = 0.8_wp * U(JCUT) + 0.2_wp * UMAX
               if (UMAX < UMAXLO)     UMAX = 1.01_wp * UMAXLO
               if (UMAX < UMAXMIN_ERG) UMAX = UMAXMIN_ERG
               refine = .true.
            end if
         elseif (P_out(NT)/Pmax > PMIN_UP) then
            ! tail at top still high -> expand UMAX
            UMAXLO = UMAX
            if (1.2_wp*UMAX < UMAXHI) then
               UMAX = 1.2_wp * UMAX
               refine = .true.
            elseif (UMAX/UMAXHI - 1.0_wp > 0.01_wp) then
               UMAX = 0.5_wp * (UMAXHI + UMAX)
               refine = .true.
            end if
         end if

         ! ---- Adjust UMIN (Draine v7 lines 1345-1381) ----
         if (P_out(2)/Pmax < PMIN_LO) then
            ! tail at bottom too low -> increase UMIN
            UMINLO = UMIN
            JCUT = 1
            do i = 1, NT
               JCUT = i
               if (P_out(i)/Pmax > PMIN_LO) exit
            end do
            if (UMIN < 0.95_wp * U(JCUT)) then
               UMIN = 0.2_wp * UMIN + 0.8_wp * U(JCUT)
               refine = .true.
            end if
         elseif (UMIN > HC_CGS_PER_CM .and. P_out(1)/Pmax > PMIN_LO) then
            ! tail at bottom still elevated and UMIN > 1 cm^-1 -> reduce
            UMINHI = UMIN
            if (0.8_wp * UMIN > UMINLO) then
               UMIN = max(HC_CGS_PER_CM, 0.8_wp * UMIN)
               refine = .true.
            elseif ((UMIN - UMINLO) > 0.01_wp*UMIN .and. UMIN/HC_CGS_PER_CM > 20.0_wp) then
               UMIN = max(0.5_wp*(UMINLO + UMIN), HC_CGS_PER_CM)
               refine = .true.
            end if
         end if

         if (.not. refine) then
            converged = .true.
            exit
         end if
      end do

      ! If iteration exited without ever computing a usable P, fall back.
      if (.not. converged) then
         if (Pmax > 0.0_wp) then
            ! best-effort: keep the last P even if formally unconverged
            converged = .true.
         end if
      end if

      deallocate(H, kappB, lnP, U)
   end subroutine narrow_iterative


   subroutine build_Cabs_pah()
      ! Cabs_pah(NLAM, NA) by mixing neutral and cation per f_ion.
      integer  :: ja, iw
      real(wp) :: Q_neu, Q_ion, ksi
      do ja = 1, NA
         ksi = sd_fion(ja)
         do iw = 1, NLAM
            call qpah_dl07(0, aeff(ja), lam(iw), Q_neu)
            call qpah_dl07(1, aeff(ja), lam(iw), Q_ion)
            Cabs_pah(iw, ja) = (1.0_wp - ksi)*Q_neu + ksi*Q_ion
         end do
         Cabs_pah(:, ja) = Cabs_pah(:, ja) * PI * (aeff(ja) * UM2CM)**2
      end do
   end subroutine build_Cabs_pah


   subroutine build_kappB_pah()
      ! Same algorithm as build_kappB() but using Cabs_pah → kappB_pah_first.
      integer  :: NW_INT, n_below
      real(wp), allocatable :: w(:), lnw(:), Cross(:)
      real(wp), allocatable :: Bt(:,:)
      real(wp) :: lnlam(NLAM), w1, dlnw
      integer  :: jt, ja, iw

      call planck_integration_grid(lam, NW_INT, n_below, w1, dlnw)
      allocate(w(NW_INT), lnw(NW_INT), Cross(NW_INT))

      do iw = 1, NLAM
         lnlam(iw) = log(lam(iw))
      end do
      do iw = 1, NW_INT
         w(iw)   = w1 * exp(real(iw-1-n_below, wp) * dlnw)
         lnw(iw) = log(w(iw))
      end do

      ! The Planck factor depends only on (T, w): evaluate it once instead of
      ! once for every size.
      allocate(Bt(NW_INT, NT))
      do jt = 1, NT
         do iw = 1, NW_INT
            Bt(iw, jt) = bbody(T_first(jt), w(iw))
         end do
      end do

      kappB_pah_first = 0.0_wp
      do ja = 1, NA
         do iw = 1, NW_INT
            call interp(lnlam, Cabs_pah(:, ja), lnw(iw), Cross(iw))
         end do
         do jt = 1, NT
            kappB_pah_first(jt, ja) = sum(Cross * Bt(:, jt) * w) * dlnw
         end do
      end do
      deallocate(Bt, w, lnw, Cross)
   end subroutine build_kappB_pah


   subroutine build_kappCMB_pah()
      ! Same as build_kappCMB() but using Cabs_pah.
      ! calc_P subtracts this from the grain's own emission, so kappB - kappCMB
      ! is the NET cooling rate and a grain cannot cool below its surroundings.
      ! T_CMB must therefore be the temperature of the CMB the FIELD carries;
      ! radfield's cmb_temperature() is the single place that value is written
      ! down.  It was hard-coded 2.9 K here while J_Mathis had moved to
      ! 2.725 K, which held grains up against photons the field never supplied.
      real(wp)            :: T_CMB
      real(wp), parameter :: lam_min  = 1000.0_wp
      integer,  parameter :: NW_INT   = 101
      real(wp) :: w(NW_INT), spec(NW_INT), Cabs_w(NW_INT), lam_max, dlnw
      integer  :: ja, iw

      T_CMB   = cmb_temperature()
      lam_max = maxval(lam)
      kappCMB_pah = 0.0_wp
      if (lam_max <= lam_min) return

      dlnw = log(lam_max/lam_min) / real(NW_INT-1, wp)
      do iw = 1, NW_INT
         w(iw)    = lam_min * exp(real(iw-1, wp) * dlnw)
         spec(iw) = bbody(T_CMB, w(iw))
      end do
      do ja = 1, NA
         do iw = 1, NW_INT
            call interp(lam, Cabs_pah(:, ja), w(iw), Cabs_w(iw))
         end do
         kappCMB_pah(ja) = sum(Cabs_w * spec * w) * dlnw
      end do
   end subroutine build_kappCMB_pah


   ! =====================================================================
   ! Internal helpers
   ! =====================================================================

   subroutine stored_q_on_model_grid(h5, dir, base, comp, Qa, Qs, Gg, found, rho)
      ! Load the stored cross-section table of one population whose axes ARE
      ! the model grid this build just fixed.  Two wavelength sets ship for each
      ! population, so the one whose lambda and a_eff match is the one this
      ! model was built on.  found = .false. when neither does -- a caller on a
      ! grid we stored no table for -- and the caller then solves the optics as
      ! it always did.
      !
      ! WHERE THE NUMBERS COME FROM IS AN ARGUMENT, not module state.  It used
      ! to be a module variable that only build_dust ever assigned, so every
      ! other caller silently took the text route whatever product it had named
      ! -- which is how calc_kext.x came to take a model's wavelength axis from
      ! its HDF5 product and that model's optics from the seven-digit text
      ! tables, and write the two into one file as /kext and /qtable.
      !
      !   h5   the model's HDF5 product, or blank.  Tried first; its two
      !        wavelength sets are the two halves of one axis, not two files.
      !   dir  a directory of text tables ('.../dl07/'), or blank.  The same
      !        numbers to their seven written digits.
      !
      ! Both blank means "no stored optics": the caller solves from the
      ! dielectric functions.  That is what calc_qtable.x asks for, being the
      ! program that writes these tables.
      !
      ! Every node must agree to a relative 1e-6.  That is set by the TEXT
      ! tables' own written precision, not by any physics -- and it is still
      ! four orders tighter than it needs to be to tell one node from the next,
      ! since these grids step by about 1% in wavelength.
      character(len=*),      intent(in)  :: h5, dir
      character(len=*),      intent(in)  :: base
      ! Group name inside the HDF5 product ('sil', 'gra', 'pah_neu', ...).
      character(len=*),      intent(in)  :: comp
      real(wp), allocatable, intent(out) :: Qa(:,:), Qs(:,:), Gg(:,:)
      logical,               intent(out) :: found
      real(wp), optional,    intent(inout) :: rho
      character(len=*), parameter :: SUF(2) = [character(len=4) :: '_euv', '    ']
      logical,          parameter :: WIDE(2) = [.true., .false.]
      real(wp), allocatable :: tl(:), ta(:), Qe(:,:)
      real(wp) :: rho_h5
      integer :: k, nw, na_t, j, i_lyman
      logical :: ok

      found = .false.

      if (len_trim(h5) > 0) then
         do k = 1, 2
            call read_sedust_grid(trim(h5), WIDE(k), tl, i_lyman, ok)
            if (.not. ok) exit          ! no such file: fall through to the text
            found = grid_matches(tl)
            deallocate(tl)
            if (.not. found) cycle
            call read_sedust_qtable(trim(h5), comp, WIDE(k), ta, Qe, Qa, Qs, &
                                    Gg, rho_h5, ok)
            if (.not. ok) then;  found = .false.;  cycle;  end if
            found = size(ta) == NA
            if (found) found = axis_matches(ta, aeff)
            if (found .and. present(rho) .and. rho_h5 > 0.0_wp) rho = rho_h5
            deallocate(ta, Qe)
            if (found) return
            deallocate(Qa, Qs, Gg)
         end do
         found = .false.
      end if

      if (len_trim(dir) == 0) return

      do k = 1, 2
         call load_q_component(trim(dir)//base//trim(SUF(k))//'.dat', nw, na_t, &
                               tl, ta, Qa, Qs, Gg, ok, rho=rho)
         if (.not. ok) cycle
         if (nw == NLAM .and. na_t == NA) then
            found = .true.
            do j = 1, NLAM
               if (abs(tl(j)/lam(j) - 1.0_wp) > 1.0e-6_wp) then;  found = .false.;  exit;  end if
            end do
            if (found) then
               do j = 1, NA
                  if (abs(ta(j)/aeff(j) - 1.0_wp) > 1.0e-6_wp) then;  found = .false.;  exit;  end if
               end do
            end if
         end if
         deallocate(tl, ta)
         if (found) return
         deallocate(Qa, Qs, Gg)
      end do

   contains

      logical function grid_matches(t)
         real(wp), intent(in) :: t(:)
         grid_matches = size(t) == NLAM
         if (grid_matches) grid_matches = axis_matches(t, lam)
      end function grid_matches

      logical function axis_matches(t, x)
         real(wp), intent(in) :: t(:), x(:)
         integer :: i
         axis_matches = .true.
         do i = 1, size(t)
            if (abs(t(i)/x(i) - 1.0_wp) > 1.0e-6_wp) then
               axis_matches = .false.;  return
            end if
         end do
      end function axis_matches

   end subroutine stored_q_on_model_grid

   subroutine free_shared_model_arrays()
      ! Drop the module arrays a model build fills, each under its own guard.
      !
      ! They used to be dropped as one list guarded on `lam` alone, which holds
      ! only while every builder allocates the same set.  build_zubko does not:
      ! it takes its grid from its own optics tables and allocates `lam` without
      ! `aeff`, so a process that built the Zubko model and then a DL07 or
      ! astrodust one hit "Attempt to DEALLOCATE unallocated 'aeff'".  One
      ! builder per process never saw it; anything that builds two did.
      if (allocated(lam))                 deallocate(lam)
      if (allocated(aeff))                deallocate(aeff)
      if (allocated(T_first))             deallocate(T_first)
      if (allocated(dn_ad))               deallocate(dn_ad)
      if (allocated(Cabs))                deallocate(Cabs)
      if (allocated(Csca))                deallocate(Csca)
      if (allocated(kappB_first))         deallocate(kappB_first)
      if (allocated(H_first))             deallocate(H_first)
      if (allocated(kappCMB))             deallocate(kappCMB)
      if (allocated(log_T_first))         deallocate(log_T_first)
      if (allocated(log_H_first))         deallocate(log_H_first)
      if (allocated(log_kappB_first))     deallocate(log_kappB_first)
      if (allocated(dn_pah))              deallocate(dn_pah)
      if (allocated(Cabs_pah))            deallocate(Cabs_pah)
      if (allocated(kappB_pah_first))     deallocate(kappB_pah_first)
      if (allocated(H_pah_first))         deallocate(H_pah_first)
      if (allocated(kappCMB_pah))         deallocate(kappCMB_pah)
      if (allocated(log_H_pah_first))     deallocate(log_H_pah_first)
      if (allocated(log_kappB_pah_first)) deallocate(log_kappB_pah_first)
   end subroutine free_shared_model_arrays


   pure logical function is_hdf5_path(path)
      ! An HDF5 product is named by its suffix, so one path argument serves
      ! both sources and no caller needs a second flag to say which it meant.
      character(len=*), intent(in) :: path
      integer :: k
      k = len_trim(path)
      is_hdf5_path = .false.
      if (k > 3) is_hdf5_path = (path(k-2:k) == '.h5')
   end function is_hdf5_path


   pure logical function euv_asked(include_euv)
      ! Default .false.: the non-ionizing part of a product's axis, which is
      ! what an interstellar radiation field illuminates and what the narrow
      ! text products carry.
      logical, optional, intent(in) :: include_euv
      euv_asked = .false.
      if (present(include_euv)) euv_asked = include_euv
   end function euv_asked


   real(wp) function d03_euv_lambda_floor() result(lam_min)
      ! Shortest wavelength Mie on the D03 dielectric functions can be solved
      ! at, which is the floor of every model whose optics come from them --
      ! DL07 and MRN both.  Both materials are required at every wavelength, so
      ! it is the LONGER of the two functions' short-wavelength ends -- the same
      ! max(silicate, graphite) those builders refuse below -- stood off it so
      ! that asking for this value is not on the rounding boundary of that
      ! refusal.
      real(wp) :: sil_lo, sil_hi, gra_lo, gra_hi
      call silicate_index_lambda_range(sil_lo, sil_hi)
      call graphite_index_lambda_range(gra_lo, gra_hi)
      lam_min = LAM_LO_MARGIN * max(sil_lo, gra_lo)
   end function d03_euv_lambda_floor


   real(wp) function astrodust_euv_lambda_floor() result(lam_min)
      ! The same for the astrodust model, whose EUV band is solved on the DH21
      ! dielectric function of the composite grain.
      real(wp) :: ad_lo, ad_hi
      call astrodust_index_lambda_range(ad_lo, ad_hi)
      lam_min = LAM_LO_MARGIN * ad_lo
   end function astrodust_euv_lambda_floor


   subroutine euv_extended_lambda_grid(lam_out, lam_min, n_extra, base)
      ! Model wavelength grid = the T-matrix Q table's grid, optionally carried
      ! below its short-wavelength end.  That end is now 1.0e-4 um (12398 eV):
      ! the ionizing band a photoionization RT host transports is INSIDE the
      ! astrodust table, so nothing is prepended for it any more.  For that
      ! model nothing can be: the DH21 dielectric function stops at
      ! 1.000032e-4 um, longward of the table's own first wavelength, so every
      ! lam_min a caller may legally ask for lands inside the table.  DL07 is
      ! where this still does something -- its D03 optical constants reach
      ! 6.205e-5 um, past the table.  (The orientation-resolved polarized table
      ! is a separate product and still ends at 0.0912 um; build_Cpol aligns the
      ! two on their long-wavelength end and reports the block below it.)
      !
      ! When lam_min IS shorter than the table's first wavelength, log-spaced
      ! points are prepended from lam_min up to just below it.  Their spacing
      ! is at most the table's own spacing at its short-wavelength end, taken
      ! from the table itself (dln lam = 0.00794 on the current axis, where the
      ! nodes below 0.0912 um are the dielectric function's own), so the
      ! extension is never coarser than the grid it joins.  lam_out(1) is set
      ! to lam_min exactly, so the caller's requested floor is covered rather
      ! than approached.  n_extra = 0 (and the plain table grid) when lam_min is
      ! absent, non-positive, or not shorter than the table -- which is what
      ! keeps the unextended model bit-identical.
      real(wp), allocatable, intent(out) :: lam_out(:)
      real(wp), optional,    intent(in)  :: lam_min
      ! Number of points prepended; 0 when the grid is the plain table grid.
      integer,  optional,    intent(out) :: n_extra
      ! Axis to prepend onto.  Absent means the loaded text Q table's own, the
      ! only axis there is on that route; a model whose axis came from the HDF5
      ! product passes it here, no text table having been opened for qt_lam to
      ! hold.
      real(wp), optional,    intent(in)  :: base(:)
      real(wp), allocatable :: b(:)
      real(wp) :: dln_qt, dln_ext, span
      integer  :: j, nx

      if (present(base)) then
         allocate(b(size(base)));  b = base
      else
         allocate(b(qt_n_lam));    b = qt_lam
      end if

      nx = 0
      if (present(lam_min)) then
         if (lam_min > 0.0_wp .and. lam_min < b(1)) then
            span   = log(b(1) / lam_min)
            dln_qt = log(b(2) / b(1))
            ! span > 0 inside this branch, so the ceiling is already >= 1.
            nx     = ceiling(span / dln_qt)
         end if
      end if
      if (present(n_extra)) n_extra = nx

      allocate(lam_out(nx + size(b)))
      if (nx > 0) then
         dln_ext = log(b(1) / lam_min) / real(nx, wp)
         do j = 1, nx
            lam_out(j) = b(1) * exp(-real(nx - j + 1, wp) * dln_ext)
         end do
         lam_out(1) = lam_min
      end if
      lam_out(nx+1:) = b
   end subroutine euv_extended_lambda_grid


   subroutine interp_q_grid(loga_target, aeff_in, q_in, q_out)
      ! Interpolate q_in(NLAM, NA_in) at log(a_target) -> q_out(NLAM)
      ! using log-linear interpolation in a. Clamps to grid edges.
      real(wp), intent(in)  :: loga_target
      real(wp), intent(in)  :: aeff_in(:)               ! (NA_in)
      real(wp), intent(in)  :: q_in(:,:)                ! (NLAM, NA_in)
      real(wp), intent(out) :: q_out(:)                 ! (NLAM)
      integer  :: NA_in, lo, hi, mid
      real(wp) :: x_lo, x_hi, t

      NA_in = size(aeff_in)
      if (loga_target <= log(aeff_in(1))) then
         q_out = q_in(:, 1); return
      end if
      if (loga_target >= log(aeff_in(NA_in))) then
         q_out = q_in(:, NA_in); return
      end if
      lo = 1; hi = NA_in
      do while (hi - lo > 1)
         mid = (lo + hi) / 2
         if (log(aeff_in(mid)) <= loga_target) then
            lo = mid
         else
            hi = mid
         end if
      end do
      x_lo = log(aeff_in(lo))
      x_hi = log(aeff_in(hi))
      t = (loga_target - x_lo) / (x_hi - x_lo)
      q_out = (1.0_wp - t) * q_in(:, lo) + t * q_in(:, hi)
   end subroutine interp_q_grid


   subroutine build_Cpol(q_file, wave_file, aeff_file, euv_q_file, euv_wave_file, &
                         loaded, euv_status)
      ! Fill Cpol(NLAM, NA), Cpol_ext(NLAM, NA), Cbir_ext(NLAM, NA) and
      ! falign_ad(NA) from the orientation-resolved DH21 spheroid table:
      !
      !   Q_pol,abs = 0.5 * (Q_abs(k perp a, E perp a) - Q_abs(k perp a, E || a))
      !   C_pol     = Q_pol,abs * pi * a_eff^2
      !   Q_pol,ext = 0.5 * (Q_ext(k perp a, E perp a) - Q_ext(k perp a, E || a))
      !   C_pol_ext = Q_pol,ext * pi * a_eff^2
      !   Q_bir,ext = 0.5 * (Q_re(k perp a, E perp a) - Q_re(k perp a, E || a))
      !   C_bir_ext = Q_bir,ext * pi * a_eff^2
      !
      ! C_pol drives the polarized emission; C_pol_ext drives the dichroic
      ! (polarized) extinction and C_bir_ext the birefringence a radiative
      ! transfer host needs for the extinction matrix. C_bir_ext is built only
      ! from a 4-block table (has_bir); it stays zero for an older 3-block table.
      !
      ! i.e. the absorption difference a perfectly aligned grain with its
      ! symmetry axis in the plane of the sky presents to the two linear
      ! polarizations. The grain loop weights this by dn * f_align(a), so the
      ! result is the INTRINSIC polarized emission; the geometric projection
      ! onto the sky is the radiative transfer's job.
      !
      ! MIXED RANDOM-ORIENTATION AVERAGES. Cabs above comes from our own
      ! T-matrix run, whose random-orientation average is exact, whereas the
      ! release table's average is the 1/3 trace (Q1+Q2+Q3)/3. The two are
      ! different approximations, so this routine deliberately does NOT touch
      ! Cabs -- it only adds the polarized channel. Measured on the shared
      ! (lambda, a_eff) grid, |trace-average / exact - 1| for Q_abs has a
      ! median of 0.022% for lambda > 30 um (0.023% over 30-100 um, 0.025%
      ! over 100-1000 um, 0.020% beyond), with a worst case of 2.0% among the
      ! grains that carry 99.999% of the geometric cross section (a <= 0.5 um).
      ! Polarized emission is a far-infrared phenomenon and, at a << lambda,
      ! both averages approach the same Rayleigh limit -- so the mixture is
      ! harmless here. It would NOT be in the ultraviolet, where the median
      ! rises to 0.21% and the worst case to 9.5%.
      character(len=*), intent(in)  :: q_file, wave_file, aeff_file
      ! The extreme-ultraviolet companion table and its wavelength axis, read
      ! only when the model grid actually reaches below the main table
      ! (n_pol_euv > 0). It shares aeff_file as its size axis.
      character(len=*), intent(in)  :: euv_q_file, euv_wave_file
      ! .true. iff the orientation-resolved table was read successfully. .false.
      ! leaves the polarized arrays allocated and zero (graceful degradation);
      ! the caller decides whether that is an error (see sed_init status 4).
      logical,          intent(out) :: loaded
      ! Extreme-ultraviolet report, 0 when the band was covered (or not asked
      ! for). Unlike the main table, this one never degrades silently: a host
      ! that asked for EUV wavelengths AND polarized optics is never left with
      ! an unannounced zero: a missing companion table is reported on stderr and
      ! the band stays zero (a documented deficit), while a table that does not
      ! reach the requested wavelengths is an error.
      !   euv_status = 8  the polarized table's wavelength grid does not match
      !                   the model's table block (count or node positions)
      !   euv_status = 9  the requested grid runs outside the wavelengths the
      !                   EUV companion table covers
      integer,          intent(out) :: euv_status
      integer  :: ja, jw
      logical  :: rok
      logical  :: euv_pol      ! an EUV companion table was read for this grid
      ! Model wavelengths lying SHORTWARD of the polarized table's own first
      ! node -- the block this routine has to fill from somewhere other than
      ! that table.  It is derived from the table's coverage, not from
      ! n_lam_euv: the scalar Q table and the polarized one are separate
      ! products and no longer start at the same wavelength (the scalar table
      ! reaches 1.0e-4 um, the polarized ones 0.0912 um), so the count of
      ! wavelengths the scalar table does not cover says nothing about the
      ! block the polarized table leaves open.  The two agree, and this reduces
      ! to n_lam_euv, whenever the two tables do share a short-wavelength end.
      integer  :: n_pol_euv
      real(wp), allocatable :: qeuv_a(:)      ! one size, all EUV table lambdas

      Cpol       = 0.0_wp
      Cpol_ext   = 0.0_wp
      Cbir_ext   = 0.0_wp
      falign_ad  = 0.0_wp
      loaded     = .false.
      euv_status = 0

      call load_q_table_jori(q_file, wave_file, aeff_file, ok=rok)
      if (.not. rok) then
         ! Reported unconditionally (not verbose-gated): a missing polarized
         ! table silently changes the model's capability, so always say so.
         write(error_unit,'(a,a)') &
            ' sed_init: no polarized optics (cannot read ', trim(q_file)//')'
         return
      end if

      ! The polarized table has to sit on the LONG-wavelength end of the model
      ! grid, node for node: both are the DH21 axis, and the model grid is that
      ! axis possibly extended shortward. Anything else is inconsistent input,
      ! not a recoverable condition. The n_pol_euv points below the table's
      ! first node are the block it does not reach; they are treated separately.
      n_pol_euv = NLAM - qj_n_lam
      if (n_pol_euv < 0) then
         write(error_unit,'(a,i0,a,i0)') ' build_Cpol: polarized table has ', &
            qj_n_lam, ' wavelengths but the model grid only has ', NLAM
         euv_status = 8;  return
      end if
      do jw = 1, qj_n_lam
         if (abs(qj_lam(jw) - lam(n_pol_euv+jw)) > 1.0e-10_wp * abs(lam(n_pol_euv+jw))) then
            write(error_unit,'(a,i0)') &
               ' build_Cpol: polarized and Q wavelength grids differ at jw=', jw
            euv_status = 8;  return
         end if
      end do

      ! EUV EXTENSION (n_pol_euv > 0). The points shortward of the polarized
      ! table's first node -- which is where the model grid runs below it,
      ! whether because lam_min extended the grid or because the scalar Q table
      ! itself reaches further than the polarized one -- have no entry in that
      ! table, and their polarized optics come from the companion
      ! table computed on the DH21_wave_euv axis by run_q_jori.f90's `euv`
      ! mode: the SAME first-principles core (Rayleigh dipole / Mishchenko
      ! T-matrix / geometric optics), the same DH21 dielectric function and the
      ! same b/a = 1.4 spheroid, only at shorter wavelengths.
      !
      ! This band is NOT a neighborhood of small error, which is why it is
      ! computed rather than set to zero: the alignment-weighted,
      ! size-integrated |C_pol,ext| / C_ext RISES from 1.3e-3 at the join to
      ! 3.6e-3 near 20.6 eV -- a factor 2.9 -- before decaying to 1.6e-4 at
      ! 100 eV, and it changes sign relative to the optical band between 0.1072
      ! and 0.1059 um. |C_bir,ext| / C_ext falls monotonically from 7.2e-3 at
      ! the join and changes sign near 51 eV. Zero would be the right asymptote
      ! only for x >> 1 AND |m - 1| << 1, and |m - 1| is 0.765 at the join,
      ! 0.564 at 20 eV, and below 0.1 only above ~50 eV, while the grains that
      ! carry the signal sit at x = 5-50. A sphere has exactly zero dichroic
      ! extinction and exactly zero birefringence, so the scalar EUV optics
      ! (q_astrodust_mod, Mie on the volume-equivalent sphere) could not have
      ! supplied these entries either.
      !
      ! WHAT REMAINS UNRESOLVED. The companion table keeps the T-matrix only
      ! where two independent convergence settings agree (run_q_jori.f90,
      ! cross_sections_large_x), which in practice certifies it to x ~ 56;
      ! above that the geometric-optics limit is used, and it carries
      ! Q_pol,ext = Q_bir,ext = 0 exactly. No extrapolation is put in their
      ! place. The dichroic extinction lost that way is nothing at the join,
      ! ~8% of the band total at 0.031 um, ~20% at 0.022 um and ~46% at
      ! 0.0124 um -- but by then |C_pol,ext| / C_ext is itself below 6e-4,
      ! about 1% of the 5.5e-2 the V band reaches, so the residual is a small
      ! fraction of an already small number. C_pol (the ABSORPTION dichroism,
      ! which drives polarized emission) has no such gap: the geometric-optics
      ! limit gets it from the opaque-grain Fresnel surface integral.
      !
      ! The table stops at 0.0124 um (100 eV) because that opaque limit needs
      ! the chord optical depth 4 Im(m) x to be large, and for astrodust it is
      ! 3.7 at x = 50 there and falls below 1 by 200 eV. A grid reaching
      ! shortward of the table is rejected (euv_status = 9) rather than filled
      ! with a limit known to be invalid.
      !
      ! GRIDS. The size axes DO differ (169 table nodes vs the
      ! size-distribution grid), so both blocks interpolate in log(a) exactly
      ! as Cabs/Csca do. The EUV wavelengths differ too -- the model's EUV
      ! spacing is set at run time by lam_min (euv_extended_lambda_grid), so it
      ! cannot coincide with the precomputed axis -- and are therefore
      ! interpolated in log(lambda) as well.
      euv_pol = .false.
      if (n_pol_euv > 0) then
         call load_q_table_jori_euv(euv_q_file, euv_wave_file, aeff_file, ok=euv_pol)
         if (.not. euv_pol) then
            ! No EUV companion table. This is the DEFAULT state: the polarized
            ! optics are generated band by band, on demand, rather than swept
            ! over a full EUV axis, so most models simply have none here. Say so
            ! and leave the EUV block at the zero set on entry -- a KNOWN
            ! DEFICIT, not the correct value. Measured first-principles size
            ! integrals (tmatrix/driver/euv_polarized_optics.f90) put the true
            ! |C_pol,ext| / C_ext at 1.26e-3 at the 0.0912 um join, rising to
            ! 3.64e-3 near 20.6 eV before decaying to 1.56e-4 at 100 eV, with
            ! the sign opposite to the optical band. Generate the band with
            ! tmatrix/driver/run_q_jori.f90 (euv mode) if the host needs it.
            write(error_unit,'(a,a)') &
               ' build_Cpol: no EUV polarized table (', trim(euv_q_file)//')'
            write(error_unit,'(a,es10.3,a)') &
               '             dichroic extinction and birefringence are zero below', &
               lam(n_pol_euv+1), ' um (zero by omission, not by physics).'
         end if
      end if
      if (euv_pol) then
         if (lam(1) < qj_lam_euv(1) * (1.0_wp - 1.0e-12_wp)) then
            write(error_unit,'(a,es10.3,a,es10.3,a)') &
               ' sed_init: lam_min =', lam(1), &
               ' um is shorter than the EUV polarized table, which starts at', &
               qj_lam_euv(1), ' um.'
            write(error_unit,'(a)') &
               '           Below it the opaque geometric-optics limit that '// &
               'table relies on is invalid.'
            euv_status = 9
            return
         end if
         if (lam(n_pol_euv) > qj_lam_euv(qj_n_lam_euv) * (1.0_wp + 1.0e-12_wp)) then
            write(error_unit,'(a,es10.3,a,es10.3,a)') &
               ' sed_init: the EUV block reaches', lam(n_pol_euv), &
               ' um but the EUV polarized table stops at', &
               qj_lam_euv(qj_n_lam_euv), ' um.'
            euv_status = 9
            return
         end if
      else
         ! No EUV block on this grid: drop any companion table a previous build
         ! left loaded, so the module state always describes the current model.
         call free_q_table_jori_euv()
      end if

      allocate(qeuv_a(max(qj_n_lam_euv, 1)))
      do ja = 1, NA
         call interp_q_grid(log(aeff(ja)), qj_aeff, qj_qpol_abs, Cpol(n_pol_euv+1:, ja))
         call interp_q_grid(log(aeff(ja)), qj_aeff, qj_qpol_ext, Cpol_ext(n_pol_euv+1:, ja))
         if (euv_pol) then
            call interp_q_grid(log(aeff(ja)), qj_aeff_euv, qj_qpol_abs_euv, qeuv_a)
            call interp_loglam_grid(qj_lam_euv, qeuv_a, lam(1:n_pol_euv), Cpol(1:n_pol_euv, ja))
            call interp_q_grid(log(aeff(ja)), qj_aeff_euv, qj_qpol_ext_euv, qeuv_a)
            call interp_loglam_grid(qj_lam_euv, qeuv_a, lam(1:n_pol_euv), Cpol_ext(1:n_pol_euv, ja))
         end if
         Cpol(:, ja)     = Cpol(:, ja)     * PI * (aeff(ja) * UM2CM)**2
         Cpol_ext(:, ja) = Cpol_ext(:, ja) * PI * (aeff(ja) * UM2CM)**2
         falign_ad(ja)   = falign_hd23(aeff(ja))
      end do

      ! Birefringent extinction. Each block is filled from whichever table
      ! carries the 4th (forward-amplitude real-part) block; the block of a
      ! 3-block table keeps the zero set at entry. The two are independent --
      ! withholding a birefringence one table does carry, because the other
      ! does not, would discard a measured optic.
      if (qj_has_bir .or. (euv_pol .and. qj_has_bir_euv)) then
         do ja = 1, NA
            if (qj_has_bir) &
               call interp_q_grid(log(aeff(ja)), qj_aeff, qj_qbir_ext, &
                                  Cbir_ext(n_pol_euv+1:, ja))
            if (euv_pol .and. qj_has_bir_euv) then
               call interp_q_grid(log(aeff(ja)), qj_aeff_euv, qj_qbir_ext_euv, qeuv_a)
               call interp_loglam_grid(qj_lam_euv, qeuv_a, lam(1:n_pol_euv), &
                                       Cbir_ext(1:n_pol_euv, ja))
            end if
            Cbir_ext(:, ja) = Cbir_ext(:, ja) * PI * (aeff(ja) * UM2CM)**2
         end do
      end if
      deallocate(qeuv_a)

      ! One table with a Q_re block and one without leaves the birefringence
      ! zero on one side of the join and not the other -- a step the host
      ! should know about, so say it rather than let it pass as physics.
      if (euv_pol .and. (qj_has_bir .neqv. qj_has_bir_euv)) then
         if (qj_has_bir) then
            write(error_unit,'(a,es10.3,a)') &
               ' build_Cpol: the EUV polarized table has no Q_re block, so the '// &
               'birefringence is zero below', lam(n_pol_euv+1), ' um.'
         else
            write(error_unit,'(a,es10.3,a)') &
               ' build_Cpol: the main polarized table has no Q_re block, so the '// &
               'birefringence is zero above', lam(n_pol_euv+1), ' um.'
         end if
      end if

      loaded = .true.
   end subroutine build_Cpol


   subroutine interp_loglam_grid(lam_in, q_in, lam_out, q_out)
      ! Linear interpolation of q_in(lam_in) onto lam_out, in log(lambda).
      ! Both axes ascending; values outside lam_in are clamped to its ends,
      ! which the caller has already excluded (build_Cpol rejects a grid that
      ! runs outside the table). Used for the EUV block, whose model
      ! wavelengths are set at run time by lam_min and so never coincide with
      ! the precomputed axis.
      real(wp), intent(in)  :: lam_in(:), q_in(:)     ! (N_in)
      real(wp), intent(in)  :: lam_out(:)             ! (N_out)
      real(wp), intent(out) :: q_out(:)               ! (N_out)
      integer  :: n_in, k, lo, hi, mid
      real(wp) :: x, t

      n_in = size(lam_in)
      do k = 1, size(lam_out)
         x = log(lam_out(k))
         if (x <= log(lam_in(1))) then
            q_out(k) = q_in(1);  cycle
         end if
         if (x >= log(lam_in(n_in))) then
            q_out(k) = q_in(n_in);  cycle
         end if
         lo = 1;  hi = n_in
         do while (hi - lo > 1)
            mid = (lo + hi) / 2
            if (log(lam_in(mid)) <= x) then
               lo = mid
            else
               hi = mid
            end if
         end do
         t = (x - log(lam_in(lo))) / (log(lam_in(hi)) - log(lam_in(lo)))
         q_out(k) = (1.0_wp - t) * q_in(lo) + t * q_in(hi)
      end do
   end subroutine interp_loglam_grid


   subroutine build_kappB()
      ! kappB_first(jt, ja) = integral_lam Cabs(lam, ja) * B_lam(T_first(jt), lam) dlam
      ! Uses a denser internal log-lam grid over [min(lam), max(lam)] and
      ! trapezoidal log-integration, matching setup_kappB1's algorithm. The
      ! grid comes from planck_integration_grid, so neither the step nor the
      ! sample points depend on how far the model grid reaches into the EUV.
      integer  :: NW_INT, n_below
      real(wp), allocatable :: w(:), lnw(:), Cross(:)
      real(wp), allocatable :: Bt(:,:)
      real(wp) :: lnlam(NLAM), w1, dlnw
      integer  :: jt, ja, iw

      call planck_integration_grid(lam, NW_INT, n_below, w1, dlnw)
      allocate(w(NW_INT), lnw(NW_INT), Cross(NW_INT))

      do iw = 1, NLAM
         lnlam(iw) = log(lam(iw))
      end do
      do iw = 1, NW_INT
         w(iw)   = w1 * exp(real(iw-1-n_below, wp) * dlnw)
         lnw(iw) = log(w(iw))
      end do

      ! The Planck factor depends only on (T, w): evaluate it once instead of
      ! once for every size.
      allocate(Bt(NW_INT, NT))
      do jt = 1, NT
         do iw = 1, NW_INT
            Bt(iw, jt) = bbody(T_first(jt), w(iw))
         end do
      end do

      kappB_first = 0.0_wp
      do ja = 1, NA
         do iw = 1, NW_INT
            call interp(lnlam, Cabs(:, ja), lnw(iw), Cross(iw))
         end do
         do jt = 1, NT
            kappB_first(jt, ja) = sum(Cross * Bt(:, jt) * w) * dlnw
         end do
      end do
      deallocate(Bt, w, lnw, Cross)
   end subroutine build_kappB


   subroutine planck_integration_grid(lam_in, nw, n_below, w1, dlnw)
      ! Internal log-lambda grid for the Planck integrals over lam_in.
      !
      ! NW_TABLE points span the model's own optics-table wavelength range,
      ! [lam_in(n_below+1), lam_in(size(lam_in))], at the step that range
      ! implies.  An EUV extension only widens the interval downward, and it
      ! carries no Planck signal: below 0.0912 um B_lam is ~1e-231 of the peak
      ! at 288 K and still only ~1e-10 at 5000 K.  Spending a fixed budget of
      ! points on the widened interval would therefore coarsen the sampling of
      ! the range that does carry the integrand, moving kappB by ~1e-4 relative
      ! for a purely numerical reason.  Instead the extra points are PREPENDED
      ! at the same step, so both the step and the sample points over the table
      ! range are independent of how far the model reaches into the EUV -- which
      ! is what makes kappB, an infrared quantity, invariant under the
      ! extension.
      !
      ! With no extension n_below is 0 and w1 is lam_in(1), so the grid is
      ! exactly the historical NW_TABLE points and an unextended model
      ! integrates bit for bit as before.
      real(wp), intent(in)  :: lam_in(:)
      integer,  intent(out) :: nw        ! total number of points
      integer,  intent(out) :: n_below   ! points prepended below the table range
      real(wp), intent(out) :: w1        ! anchor: shortest table wavelength [um]
      real(wp), intent(out) :: dlnw      ! step in ln(lambda)
      integer, parameter :: NW_TABLE = 1001
      integer :: nl, nx

      nl = size(lam_in)
      nx = 0
      ! n_lam_euv counts the prepended points of the ACTIVE model's grid, so it
      ! applies only to a grid of that model's length.
      if (n_lam_euv > 0 .and. nl == NLAM) nx = n_lam_euv

      w1      = lam_in(nx+1)
      dlnw    = log(lam_in(nl) / w1) / real(NW_TABLE-1, wp)
      n_below = 0
      if (nx > 0) n_below = ceiling(log(w1 / lam_in(1)) / dlnw)
      nw      = NW_TABLE + n_below
   end subroutine planck_integration_grid


   subroutine build_kappCMB()
      ! kappCMB(ja) = integral_(lam>1mm) Cabs(lam, ja) * B_lam(T_CMB, lam) dlam
      ! See setup_kappCMB.
      ! calc_P subtracts this from the grain's own emission, so kappB - kappCMB
      ! is the NET cooling rate and a grain cannot cool below its surroundings.
      ! T_CMB must therefore be the temperature of the CMB the FIELD carries;
      ! radfield's cmb_temperature() is the single place that value is written
      ! down.  It was hard-coded 2.9 K here while J_Mathis had moved to
      ! 2.725 K, which held grains up against photons the field never supplied.
      real(wp)            :: T_CMB
      real(wp), parameter :: lam_min  = 1000.0_wp     ! [um]
      integer,  parameter :: NW_INT   = 101
      real(wp) :: w(NW_INT), spec(NW_INT), Cabs_w(NW_INT), lam_max, dlnw
      integer  :: ja, iw

      T_CMB   = cmb_temperature()
      lam_max = maxval(lam)
      kappCMB = 0.0_wp
      if (lam_max <= lam_min) return

      dlnw = log(lam_max/lam_min) / real(NW_INT-1, wp)
      do iw = 1, NW_INT
         w(iw)    = lam_min * exp(real(iw-1, wp) * dlnw)
         spec(iw) = bbody(T_CMB, w(iw))
      end do
      do ja = 1, NA
         do iw = 1, NW_INT
            call interp(lam, Cabs(:, ja), w(iw), Cabs_w(iw))
         end do
         kappCMB(ja) = sum(Cabs_w * spec * w) * dlnw
      end do
   end subroutine build_kappCMB


   subroutine apply_induced_factor(J_lam, Jout)
      ! Multiply Jout by the induced-emission factor (1 + J_lam/B_env)
      ! where B_env(lambda) = 2*h*c^2 / lambda^5 is the Planck envelope.
      ! Mirrors Draine's method.
      ! Pulled out of the grain loop because it varies only with
      ! wavelength (no T or grain dependence).
      real(wp), intent(in)    :: J_lam(:)
      real(wp), intent(inout) :: Jout(:)
      real(wp) :: lam_m
      integer :: k
      do k = 1, NLAM
         lam_m = lam(k) * 1.0e-6_wp
         Jout(k) = Jout(k) * (1.0_wp + J_lam(k) * lam_m**5 / TWO_HCC)
      end do
   end subroutine apply_induced_factor


   ! =====================================================================
   ! Model-agnostic library layer (path B). The validated solver core
   ! (sed_grain_loop & helpers) is UNTOUCHED; these routines package a
   ! model's populations into a dust_model_t and run the core per
   ! population. The module-global grids (lam, aeff, T_first, NLAM, NT)
   ! are the *active model's* working set -- build_<model> sets them via
   ! sed_init*, so dust_emission(m,...) is correct as long as m is the
   ! model most recently built (one active model at a time).
   ! =====================================================================

   ! Copy one population's arrays (from the module globals) into a grain_pop_t.
   subroutine set_pop(p, gtype, chan, dn_in, Cabs_in, kappB_in, H_in, &
                      log_H_in, log_kappB_in, kappCMB_in, Cpol_in, falign_in, &
                      Csca_in, Cpol_ext_in, gsca_in, Cbir_ext_in, rho_bulk)
      type(grain_pop_t), intent(inout) :: p
      character(len=*),  intent(in)    :: gtype
      integer,           intent(in)    :: chan
      real(wp),          intent(in)    :: dn_in(:), kappCMB_in(:)
      real(wp),          intent(in)    :: Cabs_in(:,:), kappB_in(:,:), H_in(:,:)
      real(wp),          intent(in)    :: log_H_in(:,:), log_kappB_in(:,:)
      ! Polarized optics. Supply BOTH to make this population contribute to
      ! the polarized emission; leave them out and p%Cpol / p%falign stay
      ! unallocated, which is how dust_emission recognizes an unpolarized
      ! population.
      real(wp), optional, intent(in)   :: Cpol_in(:,:), falign_in(:)
      ! Extinction-side optics, read only by dust_extinction. Each is optional
      ! and independent: a population that does not scatter (the PAHs) simply
      ! leaves them out and contributes zero to those terms of the size
      ! integral. The emission path never touches them.
      real(wp), optional, intent(in)   :: Csca_in(:,:), Cpol_ext_in(:,:), gsca_in(:,:)
      ! Birefringent extinction, read only by dust_extinction. Optional and
      ! independent: a population without it (the PAHs, or an astrodust model
      ! built from a 3-block table) leaves it out and contributes zero.
      real(wp), optional, intent(in)   :: Cbir_ext_in(:,:)
      ! Solid mass density of the material [g/cm^3], read only by
      ! dust_mass_per_H. Omitted, the population states none and contributes no
      ! mass; see grain_pop_t.
      real(wp), optional, intent(in)   :: rho_bulk
      p%grain_type = gtype
      p%out_channel = chan
      p%aeff      = aeff          ! [um] module-global size grid (set by sed_init)
      p%dn        = dn_in
      p%Cabs      = Cabs_in
      p%kappB     = kappB_in
      p%H         = H_in
      p%log_H     = log_H_in
      p%log_kappB = log_kappB_in
      p%kappCMB   = kappCMB_in
      if (present(rho_bulk)) p%rho_bulk = rho_bulk
      if (present(Cpol_in) .and. present(falign_in)) then
         p%Cpol   = Cpol_in
         p%falign = falign_in
      end if
      if (present(Csca_in))     p%Csca     = Csca_in
      if (present(Cpol_ext_in)) p%Cpol_ext = Cpol_ext_in
      if (present(gsca_in))     p%gsca     = gsca_in
      ! Store C_bir_ext only when the table supplied it (all-zero for a 3-block
      ! table); an all-zero array would contribute nothing anyway, but leaving
      ! it unallocated keeps the has-birefringence test honest.
      if (present(Cbir_ext_in)) then
         if (any(Cbir_ext_in /= 0.0_wp)) p%Cbir_ext = Cbir_ext_in
      end if
   end subroutine set_pop


   ! Build the HD23 astrodust model into m. Channels: AD_S1, AD_S2, PAH
   ! (PAH = neutral + cation populations summed into one channel).
   subroutine build_astrodust(m, qtable_path, NT_in, T_lo, T_hi, status, &
                              qpol_path, qpol_wave_path, qpol_aeff_path, scatmat_path, &
                              load_polarized_optics, lam_min, astrodust_index_path, &
                              qpol_euv_path, qpol_euv_wave_path, kext_path, &
                              euv_tmatrix, include_euv)
      type(dust_model_t), intent(out) :: m
      character(len=*),   intent(in)  :: qtable_path
      integer,            intent(in)  :: NT_in
      real(wp),           intent(in)  :: T_lo, T_hi
      ! Optional status (0 = success, non-zero = model build failed). When
      ! present, a failed input read is reported through it instead of stopping
      ! the process; when absent the build stops on error (CLI behavior).
      !   status = 1  Q-table load failed
      !   status = 3  aligned scattering table load failed (only when
      !               scatmat_path is supplied)
      !   status = 4  a polarized Q table explicitly requested via qpol_path
      !               could not be read (an implicit-default table degrades
      !               gracefully instead)
      !   status = 5  load_polarized_optics = .false. combined with an explicit
      !               polarized-optics path (qpol_*/scatmat) -- a contradiction
      !   status = 6  astrodust dielectric function load failed (EUV band only)
      !   status = 7  lam_min below the astrodust dielectric function's own
      !               shortest wavelength (EUV band only)
      !   status = 8  the polarized table's wavelength grid does not match the
      !               model's grid (a missing EUV companion table is NOT an
      !               error -- that band degrades to a reported zero)
      !   status = 9  the EUV band of the grid runs outside the wavelengths the
      !               EUV companion table covers
      !   status = 10 an explicitly named extinction table (kext_path) could
      !               not be read
      !   status = 11 euv_tmatrix = .true. but the spheroid optics of the EUV
      !               band are not available; see sed_init
      integer, optional,  intent(out) :: status
      ! Orientation-resolved DH21 table + grid axes for the polarized optics,
      ! forwarded to sed_init. Omit to use the defaults; an implicit-default
      ! table that cannot be read leaves the model unpolarized without failing
      ! the build, but an explicit qpol_path that cannot be read fails it
      ! (status 4).
      character(len=*), optional, intent(in) :: qpol_path, qpol_wave_path, &
                                                qpol_aeff_path
      ! Aligned scattering table (run_scatmat_aligned.x product) for a polarized
      ! RT host, forwarded to sed_init. Omit to skip it; when supplied its
      ! failure fails the build (status 3), unlike the implicit polarized table.
      character(len=*), optional, intent(in) :: scatmat_path
      ! Scalar-only switch, forwarded to sed_init. .false. builds a model with
      ! no polarized optics (the polarized Q table is never opened); combining
      ! it with an explicit qpol_*/scatmat path is a contradiction (status 5).
      logical, optional, intent(in) :: load_polarized_optics
      ! Optional shortest wavelength [um] the model must cover, for a host that
      ! transports shortward of the T-matrix Q table's 0.0912 um (13.6 eV) end
      ! -- a photoionization RT spanning 6-100 eV, say. Omitting it gives the
      ! table grid, unchanged. In the extended band Cpol_ext and Cbir_ext come
      ! from the EUV companion table, computed from the same dielectric
      ! function and the same spheroid; build_Cpol states what that table
      ! resolves and what it leaves at the geometric-optics zero. lam_min below
      ! 0.0124 um (100 eV) is refused (status 9) rather than answered with an
      ! invalid limit.
      real(wp), optional, intent(in)  :: lam_min
      ! Optional dielectric function for that EUV band; must be the file
      ! qtable_path was computed from. See sed_init.
      character(len=*), optional, intent(in) :: astrodust_index_path
      ! Optional EUV companion polarized table and its wavelength axis,
      ! forwarded to sed_init. Omit for the shipped defaults.
      character(len=*), optional, intent(in) :: qpol_euv_path, qpol_euv_wave_path
      ! Size-integrated extinction table dust_extinction serves this model's
      ! scalar optics from. Omitting it takes KEXT_ASTRODUST, and a default
      ! that cannot be read is not an error -- the model is still built, and
      ! only dust_extinction is left with nothing to return (status 2 there).
      ! A kext_path that cannot be read FAILS the build (status 10): a host
      ! naming a file that is not there is a configuration error.
      character(len=*), optional, intent(in) :: kext_path
      ! How the EUV band's optics are computed; see sed_init. Default .true.
      ! = the T-matrix on the oblate spheroid of the Q table, which requires a
      ! registered euv_band_optics_i and fails with status 11 without one.
      ! .false. = the volume-equivalent-sphere Mie approximation, far cheaper
      ! and ~2% low in the geometric-optics limit.
      logical, optional, intent(in) :: euv_tmatrix
      ! Which wavelength axis to take when qtable_path names an HDF5 product;
      ! see sed_init.  Ignored for a text table.
      logical, optional, intent(in) :: include_euv
      logical :: kext_ok

      if (present(status)) status = 0

      ! Astrodust/HD23 optics: Nc=417 (rho=2.0), D16 turbostratic graphite.
      nc_coeff = 417.0d0;  nc_integer = .false.;  qpah_graphite_source = 'd16_sphere'
      ! Forward the optional paths straight through -- an absent optional stays
      ! absent in sed_init -- so sed_init substitutes the defaults and decides
      ! the explicit-vs-implicit and scalar-vs-polarized behavior from presence.
      call sed_init(qtable_path, NT_in, T_lo, T_hi, status=status, &
                    qpol_path=qpol_path, qpol_wave_path=qpol_wave_path, &
                    qpol_aeff_path=qpol_aeff_path, scatmat_path=scatmat_path, &
                    load_polarized_optics=load_polarized_optics, &
                    lam_min=lam_min, &
                    astrodust_index_path=astrodust_index_path, &
                    qpol_euv_path=qpol_euv_path, &
                    qpol_euv_wave_path=qpol_euv_wave_path, &
                    euv_tmatrix=euv_tmatrix, include_euv=include_euv)  ! sets globals
      if (present(status)) then
         if (status /= 0) return
      end if

      m%name = 'astrodust'
      active_build_id = active_build_id + 1;  m%build_id = active_build_id
      m%NA = NA;  m%NLAM = NLAM;  m%NT = NT
      m%lam = lam;  m%aeff = aeff;  m%T_first = T_first;  m%log_T_first = log_T_first
      m%use_induced_emission = use_induced_emission
      m%stoch_method = stoch_method
      ! Channels: AD (astrodust grains, production S2 enthalpy) + PAH
      ! (neutral+cation summed). S1 is an alternative diagnostic enthalpy stage
      ! -- it is NOT a separate population, so it is excluded from the model to
      ! avoid double-counting the astrodust silicate in the total SED.
      m%n_channel = 2
      allocate(m%channel_name(2))
      m%channel_name = [character(len=16):: 'AD', 'PAH']

      allocate(m%pops(3))
      ! Only the astrodust grains are aligned. HD23 take the PAHs to be
      ! unaligned, so the two PAH populations get no polarized optics and
      ! contribute nothing to the polarized emission. They do scatter, though:
      ! what scatters in a PAH population is the graphite fraction xi_gra(a)
      ! of HD23 eq. 15 -- there is no PAH scattering cross section, only a PAH
      ! absorption one -- so Csca_pah carries xi_gra(a) times the graphite
      ! sphere and both charge states share it. Leaving it out costs 3.2% of
      ! tau_sca at 0.1 um against the HD23 release and nothing longward of
      ! 0.3 um.
      !
      ! The scattering optics (Csca / gsca) are attached in both cases, so
      ! Cext / Csca / gbar and the total SED are identical between a polarized
      ! and a scalar build. The polarized optics (Cpol / Cpol_ext / Cbir_ext /
      ! falign) are attached only when the orientation-resolved table was
      ! actually loaded; a scalar-only build or a failed implicit-default load
      ! leaves the population unpolarized, which dust_has_polarized_optics reports.
      !
      ! The astrodust grains carry the porosity-corrected astrodust density and
      ! the two PAH charge states the HD23 PAH density, so dust_mass_per_H sums
      ! this model's mass per H from the same size distribution the optics use.
      if (polarized_optics_loaded) then
         call set_pop(m%pops(1), 'sil', 1, dn_ad, Cabs, kappB_first, H_first(:,:,2), &
                      log_H_first(:,:,2), log_kappB_first, kappCMB, &
                      Cpol_in=Cpol, falign_in=falign_ad, &
                      Csca_in=Csca, Cpol_ext_in=Cpol_ext, gsca_in=gsca_ad, &
                      Cbir_ext_in=Cbir_ext, rho_bulk=RHO_AD)
      else
         call set_pop(m%pops(1), 'sil', 1, dn_ad, Cabs, kappB_first, H_first(:,:,2), &
                      log_H_first(:,:,2), log_kappB_first, kappCMB, &
                      Csca_in=Csca, gsca_in=gsca_ad, rho_bulk=RHO_AD)
      end if
      call set_pop(m%pops(2), 'pah', 2, dn_cneu, Cabs_cneu, kappB_cneu, H_pah_first, &
                   log_H_pah_first, log_kappB_cneu, kappCMB_cneu, &
                   Csca_in=Csca_pah, gsca_in=gsca_pah, rho_bulk=RHO_PAH)
      call set_pop(m%pops(3), 'pah', 2, dn_cion, Cabs_cion, kappB_cion, H_pah_first, &
                   log_H_pah_first, log_kappB_cion, kappCMB_cion, &
                   Csca_in=Csca_pah, gsca_in=gsca_pah, rho_bulk=RHO_PAH)

      call load_model_extinction_table(m, sed_data_path(KEXT_ASTRODUST), kext_path, kext_ok, &
                                       default_h5=sed_data_path(KEXT_H5_ASTRODUST))
      if (.not. kext_ok) then
         if (present(status)) then
            status = 10;  return
         else if (present(kext_path)) then
            write(*,'(a)') ' build_astrodust: cannot read the extinction table '// &
                 trim(kext_path)
            stop 1
         end if
      end if
   end subroutine build_astrodust


   ! Build the DL07 model into m. Channels: SIL, CARB (carbonaceous =
   ! neutral + cation summed). Reuses sed_init_dl07 to set the globals.
   subroutine build_dl07(m, qtable_path, sd_index, u_isrf, &
                         NT_in, T_lo, T_hi, status, lam_min, kext_path, lam_axis, &
                         include_euv, stored_q_dir, pah_xsec)
      type(dust_model_t), intent(out) :: m
      character(len=*),   intent(in)  :: qtable_path
      integer,            intent(in)  :: sd_index, NT_in
      real(wp),           intent(in)  :: u_isrf, T_lo, T_hi
      ! Optional status (0 = success, non-zero = model build failed). When
      ! present, a failed input read is reported through it instead of stopping
      ! the process; when absent the build stops on error (CLI behavior).
      !   status = 1  Q-table load failed
      !   status = 7  lam_min below the D03 dielectric functions' own shortest
      !               wavelength (forwarded from sed_init_dl07; EUV band only)
      !   status = 5  an explicitly named extinction table (kext_path) could
      !               not be read
      !   status = 8  pah_xsec is not one of 'dl07' | 'ld01'
      integer, optional,  intent(out) :: status
      ! Optional shortest wavelength [um] the model must cover; see
      ! build_astrodust. This model's optics are dielectric-function Mie
      ! throughout, so the extension is a grid extension only.
      real(wp), optional, intent(in)  :: lam_min
      ! Size-integrated extinction table dust_extinction serves this model's
      ! scalar optics from; see build_astrodust. Omitting it takes KEXT_DL07.
      character(len=*), optional, intent(in) :: kext_path
      ! The model's wavelength axis, given outright instead of taken from
      ! qtable_path; see sed_init_dl07.  This is how the HDF5 product supplies
      ! this model's own grid.
      real(wp), optional, intent(in) :: lam_axis(:)
      ! Which axis to take when qtable_path names an HDF5 product; see
      ! sed_init_dl07.
      logical, optional, intent(in) :: include_euv
      ! Where this model's stored cross sections come from, blank for none at
      ! all; see sed_init_dl07.
      character(len=*), optional, intent(in) :: stored_q_dir
      ! Which published PAH absorption cross section the carbonaceous blend
      ! takes: 'dl07' (default, the model's own -- Draine & Li 2007) or 'ld01'
      ! (Li & Draine 2001, the earlier vintage Draine's 2003 kext_albedo table
      ! was computed with).  The two differ only in the carbonaceous ABSORPTION;
      ! the silicate optics, the graphite scattering and the charge mixing are
      ! the same, so a pair of builds isolates that one cross section.
      character(len=*), optional, intent(in) :: pah_xsec
      logical :: kext_ok
      character(len=8) :: vintage

      if (present(status)) status = 0

      vintage = 'dl07';  if (present(pah_xsec)) vintage = pah_xsec
      if (trim(vintage) /= 'dl07' .and. trim(vintage) /= 'ld01') then
         if (present(status)) then
            status = 8;  return
         else
            write(*,'(a,a,a)') ' build_dl07: pah_xsec = ''', trim(vintage), &
               ''' is not one of dl07 | ld01'
            stop 1
         end if
      end if

      ! DL07 carbonaceous optics (matching Draine): Nc=470 (rho~2.2, NINT),
      ! D03 graphite, and the Draine-2003a 0.93 abundance reduction.
      nc_coeff = 470.0d0;  nc_integer = .true.;  qpah_graphite_source = 'd03_sphere'
      gd_apply_d03_reduction = .true.
      qpah_xsec_vintage = trim(vintage)
      ! The stored cross sections are the DL07 vintage -- that is what
      ! calc_qtable.x computed -- so reading them back would leave an 'ld01'
      ! request with no effect at all.  That vintage is solved from the
      ! dielectric functions instead.
      if (trim(vintage) == 'dl07') then
         call sed_init_dl07(qtable_path, sd_index, u_isrf, NT_in, T_lo, T_hi, &
                            status=status, lam_min=lam_min, lam_axis=lam_axis, &
                            include_euv=include_euv, stored_q_dir=stored_q_dir)
      else
         call sed_init_dl07(qtable_path, sd_index, u_isrf, NT_in, T_lo, T_hi, &
                            status=status, lam_min=lam_min, lam_axis=lam_axis, &
                            include_euv=include_euv, stored_q_dir='')
      end if
      if (present(status)) then
         if (status /= 0) return
      end if

      m%name = 'dl07'
      active_build_id = active_build_id + 1;  m%build_id = active_build_id
      m%NA = NA;  m%NLAM = NLAM;  m%NT = NT
      m%lam = lam;  m%aeff = aeff;  m%T_first = T_first;  m%log_T_first = log_T_first
      m%use_induced_emission = use_induced_emission
      m%stoch_method = stoch_method
      m%n_channel = 2
      allocate(m%channel_name(2))
      m%channel_name = [character(len=16):: 'SIL', 'CARB']

      ! sed_init_dl07 stores silicate in dn_ad/Cabs/H_first(:,:,1) and the
      ! carbonaceous charge states in dn_cneu/cion / Cabs_cneu/cion, but it
      ! does NOT build the charge-resolved kappB/kappCMB -- the production solver
      ! sed_solve_dl07 builds those on the fly.  The library copies them into
      ! the populations, so build them here from the DL07 Cabs, mirroring
      ! sed_solve_dl07 (Cabs_pah / kappB_pah_first are reused as scratch).
      ! (Without this, build_dl07 only works after a prior build_astrodust
      !  has allocated these arrays, and would then reuse stale astrodust
      !  kappB for the DL07 carbonaceous grains.)
      if (allocated(kappB_cneu)) deallocate(kappB_cneu, kappB_cion, &
            log_kappB_cneu, log_kappB_cion, kappCMB_cneu, kappCMB_cion)
      allocate(kappB_cneu(NT, NA), kappB_cion(NT, NA), &
               log_kappB_cneu(NT, NA), log_kappB_cion(NT, NA), &
               kappCMB_cneu(NA), kappCMB_cion(NA))
      Cabs_pah = Cabs_cneu
      call build_kappB_pah();    kappB_cneu   = kappB_pah_first
      call build_kappCMB_pah();  kappCMB_cneu = kappCMB_pah
      Cabs_pah = Cabs_cion
      call build_kappB_pah();    kappB_cion   = kappB_pah_first
      call build_kappCMB_pah();  kappCMB_cion = kappCMB_pah
      log_kappB_cneu = log(max(kappB_cneu, tiny(0.0_wp)))
      log_kappB_cion = log(max(kappB_cion, tiny(0.0_wp)))

      allocate(m%pops(3))
      ! All three populations scatter, so all three carry their scattering
      ! optics into dust_extinction.  The two charge states share one
      ! scattering description and differ in dn and absorption only.
      ! Densities: the DL07 astrosilicate and graphite bulk values.  Both
      ! carbonaceous charge states are the same continuous material sequence, so
      ! both take the graphite density (see RHO_GRAPHITE).
      call set_pop(m%pops(1), 'sil', 1, dn_ad, Cabs, kappB_first, H_first(:,:,1), &
                   log_H_first(:,:,1), log_kappB_first, kappCMB, &
                   Csca_in=Csca, gsca_in=gsca_ad, rho_bulk=RHO_ASTROSIL)
      call set_pop(m%pops(2), 'pah', 2, dn_cneu, Cabs_cneu, kappB_cneu, H_pah_first, &
                   log_H_pah_first, log_kappB_cneu, kappCMB_cneu, &
                   Csca_in=Csca_car, gsca_in=gsca_car, rho_bulk=RHO_GRAPHITE)
      call set_pop(m%pops(3), 'pah', 2, dn_cion, Cabs_cion, kappB_cion, H_pah_first, &
                   log_H_pah_first, log_kappB_cion, kappCMB_cion, &
                   Csca_in=Csca_car, gsca_in=gsca_car, rho_bulk=RHO_GRAPHITE)

      ! Each vintage's own curve: what dust_extinction serves is the size
      ! integral of the very cross sections the model was built on, not of the
      ! other vintage's.  The DL07 one is unmarked, so the shipped names and
      ! the /kext group do not move.
      call load_model_extinction_table(m, sed_data_path(dl07_kext_default(vintage)), &
                                       kext_path, kext_ok, &
                                       default_h5=sed_data_path(KEXT_H5_DL07), &
                                       h5_group='kext'//trim(dl07_kext_tag(vintage)))
      if (.not. kext_ok) then
         if (present(status)) then
            status = 5;  return
         else if (present(kext_path)) then
            write(*,'(a)') ' build_dl07: cannot read the extinction table '//trim(kext_path)
            stop 1
         end if
      end if
   end subroutine build_dl07

   ! Build the Mathis, Rumpl & Nordsieck (1977) graphite + silicate model into
   ! m.  Channels: GRA, SIL.
   !
   ! WHAT THE MODEL IS.  Two materials, each with the power law
   !
   !     dn_i/da = A_i n_H a^-3.5 ,     0.005 um <= a <= 0.25 um
   !
   ! cut sharply at both ends (Draine & Lee 1984, eq. 5.1; the cutoffs are
   ! MRN's own estimates, which DL84 held fixed).  There are NO PAHs in it:
   ! the smallest grain is a 50 A graphite sphere and the 2175 A feature is
   ! graphite's, so the emergent SED carries no aromatic features at all.
   ! That is the model, not a gap in the solve.
   !
   ! OPTICS.  Mie on the Draine (2003) dielectric functions -- the very
   ! q_silicate_full and q_graphite_full the DL07 model's silicate and
   ! graphite come from, graphite as 1/3 E||c + 2/3 E-perp-c.  MRN worked with
   ! Wickramasinghe's optical constants; DL84 recomputed the same size
   ! distribution on theirs, and D03 is the current revision of those.  So the
   ! grid is free: like DL07, this model takes only a wavelength AXIS from
   ! qtable_path and solves every optic on it.
   !
   ! ENTHALPY AND DENSITY.  The Draine & Li (2001) heat capacities, 'Sil' and
   ! 'Car0', and RHO_ASTROSIL / RHO_GRAPHITE -- the densities those capacities
   ! and the D03 optics are both defined with.  They are not the 1984-vintage
   ! 3.3 and 2.24 behind Draine's own MRN table, so the dust mass per H here
   ! is larger than his by 3.1% and K_abs = C_abs/M_dust with it.  C_ext/H,
   ! which is what the size distribution fixes, is unaffected and is what the
   ! reference comparison tests.
   subroutine build_mrn(m, qtable_path, NT_in, T_lo, T_hi, status, lam_min, &
                        kext_path, lam_axis, include_euv, stored_q_dir)
      type(dust_model_t), intent(out) :: m
      character(len=*),   intent(in)  :: qtable_path
      integer,            intent(in)  :: NT_in
      real(wp),           intent(in)  :: T_lo, T_hi
      ! Optional status (0 = success, non-zero = model build failed). When
      ! present, a failed input read is reported through it instead of stopping
      ! the process; when absent the build stops on error (CLI behavior).
      !   status = 1  wavelength axis could not be read
      !   status = 2  lam_min below the D03 dielectric functions' own shortest
      !               wavelength (EUV band only)
      !   status = 3  an explicitly named extinction table (kext_path) could
      !               not be read
      integer, optional,  intent(out) :: status
      ! Optional shortest wavelength [um] the model must cover; see
      ! build_astrodust.  This model's optics are dielectric-function Mie
      ! throughout, so the extension is a grid extension only.
      real(wp), optional, intent(in)  :: lam_min
      ! Size-integrated extinction table dust_extinction serves this model's
      ! scalar optics from; see build_astrodust.  Omitting it takes this
      ! model's own product, then KEXT_MRN.
      character(len=*), optional, intent(in) :: kext_path
      ! The model's wavelength axis, given outright instead of taken from
      ! qtable_path.  This is how the HDF5 product supplies this model's own
      ! grid; qtable_path is then not read at all.
      real(wp), optional, intent(in) :: lam_axis(:)
      ! Which axis to take when qtable_path names an HDF5 product.
      logical, optional, intent(in) :: include_euv
      ! Where this model's stored cross sections come from, blank for none at
      ! all -- which is what calc_qtable.x asks for, being the program that
      ! writes them.  Omitted, the model's own directory under the data root,
      ! with the /qtable of an HDF5 qtable_path tried ahead of it.
      character(len=*), optional, intent(in) :: stored_q_dir

      integer  :: ja, jw, jt, k
      real(wp) :: a_um, geo, t, dlga, qext1, qsca1, qabs1, gsca1
      real(wp) :: sil_lam_lo, sil_lam_hi, gra_lam_lo, gra_lam_hi, d03_lam_lo
      real(wp) :: A_norm(2), acm
      real(wp), allocatable :: lam_grid(:), lam_base(:), lna(:)
      real(wp), allocatable :: Cabs_gra(:,:), Csca_gra(:,:), gsca_gra(:,:)
      real(wp), allocatable :: Cabs_sil(:,:), Csca_sil(:,:), gsca_sil(:,:)
      real(wp), allocatable :: kappB_gra(:,:), kappB_sil(:,:)
      real(wp), allocatable :: kappCMB_gra(:), kappCMB_sil(:)
      real(wp), allocatable :: H_gra(:,:), H_sil(:,:), dn_gra(:), dn_sil(:)
      real(wp), allocatable :: tQa(:,:), tQs(:,:), tGg(:,:)
      real(wp), allocatable :: uQa(:,:), uQs(:,:), uGg(:,:)
      logical  :: got_gra, got_sil, rok, kext_ok
      character(len=512) :: q_h5, q_dir

      if (present(status)) status = 0
      A_norm = 10.0_wp ** LOG_A_DL84

      ! ---- wavelength axis -------------------------------------------------
      ! The two places stored optics can come from, resolved once.  When
      ! qtable_path names an HDF5 product, that product is also where the
      ! optics come from: one file, one set of numbers.
      q_h5 = ''
      if (is_hdf5_path(qtable_path)) q_h5 = qtable_path
      q_dir = sedust_dir(trim(sed_get_data_root()), 'mrn')
      if (present(stored_q_dir)) then
         q_dir = stored_q_dir
         if (len_trim(stored_q_dir) == 0) q_h5 = ''
      end if

      if (present(lam_axis) .or. is_hdf5_path(qtable_path)) then
         if (present(lam_axis)) then
            allocate(lam_grid(size(lam_axis)))
            lam_grid = lam_axis
         else
            call read_sedust_grid(qtable_path, euv_asked(include_euv), lam_grid, k, rok)
            if (.not. rok) then
               if (present(status)) then
                  status = 1;  return
               else
                  write(*,'(a,a)') ' build_mrn: cannot read /grid/lambda from ', &
                                   trim(qtable_path)
                  stop 1
               end if
            end if
         end if
         ! lam_min still applies: the product's axis is where this model's grid
         ! STARTS, not a ceiling on what the caller may ask for.
         call move_alloc(lam_grid, lam_base)
         call euv_extended_lambda_grid(lam_grid, lam_min, base=lam_base)
         deallocate(lam_base)
         n_lam_euv = count(lam_grid < LAM_LYMAN_UM)
      else
         if (present(status)) then
            call load_q_table(qtable_path, ok=rok)
            if (.not. rok) then;  status = 1;  return;  end if
         else
            call load_q_table(qtable_path)
         end if
         call euv_extended_lambda_grid(lam_grid, lam_min, n_extra=n_lam_euv)
      end if

      ! Both materials are Mie on the D03 dielectric functions, so an EUV
      ! extension is a grid extension only -- as long as the grid stays inside
      ! those functions' own coverage.  Past it `interp` freezes (n, k) at the
      ! boundary value, which would pass a constant index off as physics.
      if (n_lam_euv > 0) then
         call silicate_index_lambda_range(sil_lam_lo, sil_lam_hi)
         call graphite_index_lambda_range(gra_lam_lo, gra_lam_hi)
         d03_lam_lo = max(sil_lam_lo, gra_lam_lo)
         if (lam_grid(1) < d03_lam_lo) then
            write(*,'(a,es10.3,a)') ' build_mrn: lam_min =', lam_grid(1), &
               ' um is shorter than the D03 dielectric functions,'
            write(*,'(a,es10.3,a)') '            which stop at', d03_lam_lo, &
               ' um; (n, k) below it would be frozen at the boundary value.'
            if (present(status)) then
               status = 2;  deallocate(lam_grid);  return
            else
               stop 1
            end if
         end if
      end if

      ! ---- shared grids ----------------------------------------------------
      NLAM = size(lam_grid)
      NA   = MRN_NA
      NT   = NT_in

      call free_shared_model_arrays()
      ! This model has no polarized optics; drop anything a previous astrodust
      ! build left behind rather than leave stale arrays on the wrong grid.
      polarized_optics_loaded = .false.
      if (allocated(Cpol)) deallocate(Cpol, Cpol_ext, Cbir_ext, falign_ad)
      if (allocated(gsca_ad))  deallocate(gsca_ad)
      if (allocated(Csca_car)) deallocate(Csca_car, gsca_car)
      allocate(lam(NLAM), aeff(NA), T_first(NT), log_T_first(NT))
      allocate(Cabs(NLAM, NA), Csca(NLAM, NA), gsca_ad(NLAM, NA))
      allocate(kappB_first(NT, NA), kappCMB(NA))
      lam = lam_grid
      deallocate(lam_grid)

      ! Log-spaced radii with the two cutoffs ON the grid, so the sharp ends of
      ! the power law are grid points and not interpolated across.
      dlga = log(MRN_AMAX / MRN_AMIN) / real(NA - 1, wp)
      do ja = 1, NA
         aeff(ja) = MRN_AMIN * exp(dlga * real(ja - 1, wp))
      end do
      allocate(lna(NA))
      lna = log(aeff)

      do jt = 1, NT
         t = log(T_lo) + (log(T_hi) - log(T_lo)) * real(jt-1, wp) / real(NT-1, wp)
         T_first(jt) = exp(t)
      end do
      log_T_first = log(T_first)
      call p_sub_setup(lam)

      ! ---- optics ----------------------------------------------------------
      call stored_q_on_model_grid(trim(q_h5), trim(q_dir), 'q_mrn_gra', 'gra', &
                                  tQa, tQs, tGg, got_gra)
      call stored_q_on_model_grid(trim(q_h5), trim(q_dir), 'q_mrn_sil', 'sil', &
                                  uQa, uQs, uGg, got_sil)
      if (sed_verbose) then
         if (got_gra .and. got_sil) then
            if (len_trim(q_h5) > 0) then
               write(*,'(a,a)') ' build_mrn: optics read from ', trim(q_h5)
            else
               write(*,'(a,a)') ' build_mrn: optics read from the stored'// &
                                ' tables under ', trim(q_dir)
            end if
         else
            write(*,'(a)') ' build_mrn: optics solved from the dielectric'// &
                           ' functions (no stored table on this grid)'
         end if
      end if

      allocate(Cabs_gra(NLAM, NA), Csca_gra(NLAM, NA), gsca_gra(NLAM, NA))
      allocate(Cabs_sil(NLAM, NA), Csca_sil(NLAM, NA), gsca_sil(NLAM, NA))
      allocate(dn_gra(NA), dn_sil(NA))
      do ja = 1, NA
         a_um = aeff(ja)
         geo  = PI * (a_um * UM2CM)**2
         do jw = 1, NLAM
            if (got_gra) then
               qabs1 = tQa(jw, ja);  qsca1 = tQs(jw, ja);  gsca1 = tGg(jw, ja)
            else
               call q_graphite_full(a_um, lam(jw), qext1, qsca1, qabs1, gsca1)
            end if
            Cabs_gra(jw, ja) = qabs1 * geo
            Csca_gra(jw, ja) = qsca1 * geo
            gsca_gra(jw, ja) = gsca1
            if (got_sil) then
               qabs1 = uQa(jw, ja);  qsca1 = uQs(jw, ja);  gsca1 = uGg(jw, ja)
            else
               call q_silicate_full(a_um, lam(jw), qext1, qsca1, qabs1, gsca1)
            end if
            Cabs_sil(jw, ja) = qabs1 * geo
            Csca_sil(jw, ja) = qsca1 * geo
            gsca_sil(jw, ja) = gsca1
         end do
         ! dn per bin [1/H] = A_i a^-3.5 da, with da the trapezoidal-in-log
         ! width in CM (the power law is per cm of radius); the two endpoints
         ! carry half a bin, which is what makes the sum the integral over the
         ! closed interval [a_min, a_max].
         acm = a_um * UM2CM
         dn_gra(ja) = A_norm(1) * acm**MRN_ALPHA * bin_da(ja)
         dn_sil(ja) = A_norm(2) * acm**MRN_ALPHA * bin_da(ja)
      end do
      if (allocated(tQa)) deallocate(tQa, tQs, tGg)
      if (allocated(uQa)) deallocate(uQa, uQs, uGg)

      ! ---- Planck-averaged opacities and enthalpies, material by material ---
      ! build_kappB / build_kappCMB read the module Cabs on the module grids,
      ! so each material passes through them in turn.
      allocate(kappB_gra(NT, NA), kappB_sil(NT, NA))
      allocate(kappCMB_gra(NA), kappCMB_sil(NA))
      Cabs = Cabs_gra
      call build_kappB();    kappB_gra   = kappB_first
      call build_kappCMB();  kappCMB_gra = kappCMB
      Cabs = Cabs_sil
      call build_kappB();    kappB_sil   = kappB_first
      call build_kappCMB();  kappCMB_sil = kappCMB

      allocate(H_gra(NT, NA), H_sil(NT, NA))
      do ja = 1, NA
         do jt = 1, NT
            H_gra(jt, ja) = enthalpy_DL01(T_first(jt), aeff(ja), 'Car0')
            H_sil(jt, ja) = enthalpy_DL01(T_first(jt), aeff(ja), 'Sil ')
         end do
      end do

      ! ---- assemble --------------------------------------------------------
      m%name = 'mrn'
      active_build_id = active_build_id + 1;  m%build_id = active_build_id
      m%NA = NA;  m%NLAM = NLAM;  m%NT = NT
      m%lam = lam;  m%aeff = aeff;  m%T_first = T_first;  m%log_T_first = log_T_first
      m%use_induced_emission = use_induced_emission
      m%stoch_method = stoch_method
      m%n_channel = 2
      allocate(m%channel_name(2))
      m%channel_name = [character(len=16):: 'GRA', 'SIL']

      allocate(m%pops(2))
      ! Both materials scatter, so both carry their scattering optics into
      ! dust_extinction.
      call set_pop(m%pops(1), 'gra', 1, dn_gra, Cabs_gra, kappB_gra, H_gra, &
                   log(max(H_gra, tiny(0.0_wp))), &
                   log(max(kappB_gra, tiny(0.0_wp))), kappCMB_gra, &
                   Csca_in=Csca_gra, gsca_in=gsca_gra, rho_bulk=RHO_GRAPHITE)
      call set_pop(m%pops(2), 'sil', 2, dn_sil, Cabs_sil, kappB_sil, H_sil, &
                   log(max(H_sil, tiny(0.0_wp))), &
                   log(max(kappB_sil, tiny(0.0_wp))), kappCMB_sil, &
                   Csca_in=Csca_sil, gsca_in=gsca_sil, rho_bulk=RHO_ASTROSIL)

      call load_model_extinction_table(m, sed_data_path(KEXT_MRN), &
                                       kext_path, kext_ok, &
                                       default_h5=sed_data_path(KEXT_H5_MRN))
      if (.not. kext_ok) then
         if (present(status)) then
            status = 3;  return
         else if (present(kext_path)) then
            write(*,'(a)') ' build_mrn: cannot read the extinction table '//trim(kext_path)
            stop 1
         end if
      end if

      deallocate(lna, Cabs_gra, Csca_gra, gsca_gra, Cabs_sil, Csca_sil, gsca_sil)
      deallocate(kappB_gra, kappB_sil, kappCMB_gra, kappCMB_sil)
      deallocate(H_gra, H_sil, dn_gra, dn_sil)

   contains
      ! Size-bin width da_i in CM on the log radius grid: the full central
      ! difference inside, half of it at each end.
      pure function bin_da(j) result(da_out)
         integer, intent(in) :: j
         real(wp) :: da_out
         if (j == 1) then
            da_out = aeff(j) * 0.5_wp * (lna(2) - lna(1))
         else if (j == NA) then
            da_out = aeff(j) * 0.5_wp * (lna(NA) - lna(NA-1))
         else
            da_out = aeff(j) * 0.5_wp * (lna(j+1) - lna(j-1))
         end if
         da_out = da_out * UM2CM
      end function bin_da
   end subroutine build_mrn



   ! Build the Zubko (ZDA 2004) BARE-GR-S model into m. Three components
   ! (PAH, Graphite, Silicate), each with its OWN size grid (the component's
   ! ZDA optics-table radii), so the populations carry component-by-component grids.
   ! Size distribution from the ZDA formula; optics from the ZDA optics tables
   ! (Cabs = Qabs*pi*a^2); enthalpy from the specific-heat calorimetry tables
   ! (H = u_spec(T)*rho*(4pi/3)a^3). The shared lambda grid is the optics
   ! grid (all 3 components share 1201 wavelengths). Channels: PAH, GRA, SIL.
   subroutine build_zubko(m, config_path, data_dir, NT_in, T_lo, T_hi, status, kext_path, &
                          lam_min, include_euv, qtable_path, optics)
      type(dust_model_t), intent(out) :: m
      character(len=*),   intent(in)  :: config_path, data_dir
      integer,            intent(in)  :: NT_in
      real(wp),           intent(in)  :: T_lo, T_hi
      ! Optional status (0 = success, non-zero = model build failed). When
      ! present, a bad input is reported through it instead of stopping; when
      ! absent the build stops on error (CLI behavior).
      !   status = 1  config read failed
      !   status = 2  fewer than 3 components in the config
      !   status = 3  a component's optics read failed
      !   status = 4  a component's size/wavelength grid is inconsistent
      !   status = 5  a component's calorimetry read failed
      !   status = 6  an explicitly named extinction table (kext_path) could
      !               not be read
      !   status = 7  lam_min shorter than this model's own optics table, which
      !               it cannot extend
      !   status = 8  optics is not one of 'zda' | 'mie_d03'
      integer, optional,  intent(out) :: status
      ! Size-integrated extinction table dust_extinction serves this model's
      ! scalar optics from; see build_astrodust. Omitting it takes KEXT_ZUBKO.
      character(len=*), optional, intent(in) :: kext_path
      ! Shortest wavelength [um] the model must COVER.  One meaning for every
      ! model: astrodust and DL07 meet it by extending the grid on the
      ! dielectric function their optics came from, and this model meets it
      ! when its own table already reaches -- the ZDA tables ARE the model and
      ! nothing here can solve a new wavelength for them -- or refuses it
      ! (status 7) when they do not.  It never NARROWS the grid.  It used to,
      ! which is how a host with one physically motivated floor and one
      ! build_dust call silently truncated this model alone while extending the
      ! other two.  Narrowing is what include_euv is for.
      real(wp), optional, intent(in) :: lam_min
      ! Whether to keep the ionizing part of the model's own optics grid.
      ! Default .true., the whole table.  .false. cuts at lyman_index, the same
      ! index cut the HDF5 products carry for astrodust and DL07, so one
      ! argument means one thing across the models.  Pure row selection
      ! on the tables: no value is recomputed.  The ZDA tables start at
      ! 1.0e-3 um (1.24 keV), 91 times harder than a field illuminated to the
      ! Lyman limit, and every wavelength kept costs solver time in every cell.
      logical, optional, intent(in) :: include_euv
      ! The model's HDF5 product, when the caller has one.  This is where the
      ! stored cross sections come from; blank or omitted leaves the text
      ! tables under data_dir and then the ZDA optics tables themselves.  It is
      ! an argument rather than module state so that a program writing /kext
      ! into a product takes its optics from the /qtable of that same product.
      character(len=*), optional, intent(in) :: qtable_path
      ! WHICH optics inside that product.  Two sets are stored:
      !   'zda'      (default) the tables the Camps et al. (2015) benchmark
      !              distributes, which is what the seven codes it compares
      !              against read, so a run against that benchmark measures the
      !              stochastic-heating solver and not an optics difference.
      !              Their own headers: Zubko's multilayer-sphere code
      !              (1997-2002) for the silicate and graphite, Misselt's 2009
      !              implementation of Li & Draine (2001) / Draine & Li (2007)
      !              for the PAHs.
      !   'mie_d03'  this tree's own recomputation, Mie on the Draine (2003)
      !              optical constants.  It reproduces the distributed silicate
      !              to a mean relative 2.4e-6 and the graphite to 0.45% in
      !              Q_sca and 8% in Q_abs; the PAH component is a different
      !              implementation of the same LD01/DL07 prescription.
      ! Both live in /qtable of one file, the second under names suffixed
      ! _mie_d03, so the choice is a group name and not a second product to
      ! keep in step.  A tree without the HDF5 product falls back to the text
      ! q_zubko_*.dat, which are the recomputation, and then to the
      ! distributed tables themselves.
      character(len=*), optional, intent(in) :: optics
      logical               :: kext_ok

      type(zda_comp_t)      :: comps(ZDA_MAXCOMP)
      integer               :: ncomp, ic, jt, ja, jw, nsize, nwave, ntc, k_lo
      real(wp)              :: rho, vol_fac, mass, dlna, uspec, t, wdev
      real(wp)              :: rho_h5
      real(wp), allocatable :: qe_h5(:,:)
      real(wp), allocatable :: a_opt(:), lam_opt(:), qa(:,:), qs(:,:), gg(:,:)
      real(wp), allocatable :: Tcal(:), Ucal(:), Ccal(:), Hcol(:)
      logical               :: rok, wide
      character(len=512)    :: q_h5
      character(len=16)     :: qset
      character(len=32)     :: h5comp
      character(len=16)     :: cn(3)
      character(len=8)      :: gt(3)
      character(len=64)     :: optf

      if (present(status)) status = 0
      wide = .true.;  if (present(include_euv)) wide = include_euv
      q_h5 = '';      if (present(qtable_path)) q_h5 = qtable_path
      qset = 'zda';   if (present(optics)) qset = optics
      if (trim(qset) /= 'zda' .and. trim(qset) /= 'mie_d03') then
         if (present(status)) then
            status = 8;  return
         else
            write(*,'(a,a,a)') ' build_zubko: optics = ''', trim(qset), &
               ''' is not one of zda | mie_d03'
            stop 1
         end if
      end if

      cn = [character(len=16):: 'PAH', 'GRA', 'SIL']
      gt = [character(len=8) :: 'pah', 'gra', 'sil']

      if (present(status)) then
         call read_zda_config(config_path, ncomp, comps, ok=rok)
         if (.not. rok) then;  status = 1;  return;  end if
      else
         call read_zda_config(config_path, ncomp, comps)
      end if
      if (ncomp < 3) then
         if (present(status)) then
            status = 2;  return
         else
            write(*,'(a)') ' build_zubko: expected 3 components'; stop 1
         end if
      end if

      m%name = 'zubko'
      active_build_id = active_build_id + 1;  m%build_id = active_build_id
      m%use_induced_emission = use_induced_emission
      m%stoch_method = stoch_method
      m%n_channel = 3
      allocate(m%channel_name(3));  m%channel_name = cn
      allocate(m%pops(3))

      do ic = 1, 3
         ! Cross-section file name from the config ('Cross Sections=...').
         optf = trim(comps(ic)%xsec)//'.dat'
         ! The ZDA model is defined by these tables, so the scattering side is
         ! read from the same file as the absorption (Q_sca and <cos> columns)
         ! rather than recomputed here.  They are a HOMOGENEOUS-SPHERE Mie
         ! calculation -- ZDA 2004 sec. 3 says so for the bare grains, the
         ! effective-medium step of that paper belonging to its COMPOSITE models --
         ! and this tree reproduces the silicate one to a mean relative 2.4e-6 over
         ! all 121 x 1201 cells, by Mie on the Draine (2003) optical constants under
         ! data/dielectric/.  That fixes what they are made of: their headers name
         ! Draine's older eps_Sil / eps_Gra, but the D03 revision is what reproduces
         ! them, and only D03 reaches the 1e4 um end of their grid at all.  The
         ! graphite agrees in Q_sca to 0.45% and in Q_abs to 8% out to 1e3 um.
         !
         ! The PAH file is a COMPOSITE, and measurably so -- its header names the
         ! Li & Draine (2001) / Draine & Li (2007) cross sections and a 2009
         ! implementation by Misselt, not Zubko's code, and comparing it cell by
         ! cell against Gra_121_1201.dat shows what it is made of:
         !   Q_sca and <cos>  IDENTICAL to graphite, every size, all 1201
         !                    wavelengths.  The PAH population scatters as the
         !                    graphite sphere of the same radius; there is no
         !                    PAH scattering calculation in it.
         !   Q_abs            identical to graphite shortward of 0.0585 um
         !                    (21.2 eV, the DL07 PAH-to-graphite transition) and
         !                    the PAH prescription longward of it.
         ! Its 28 radii are graphite's first 28 except the smallest, which the
         ! file states as 3.50e-4 um while its own x = 2 pi a / lambda column was
         ! computed with graphite's 3.55e-4.
         !
         ! This tree recomputes the PAH absorption with qpah (a different
         ! implementation of the same LD01/DL07 prescription) and stores NO
         ! scattering for it, so the model built from the stored tables has a
         ! non-scattering PAH population.  Measured on the distributed tables,
         ! the PAH share of the model's C_sca is at most 0.17% (at 0.036 um) and
         ! below 0.06% everywhere else, so this is a small approximation -- but
         ! it is one, and it is not what the benchmark's tables do.
         !
         ! This tree's OWN table first: calc_qtable.x recomputed each component
         ! on the ZDA grids by the routines above.  The wide table is the one to
         ! load -- include_euv below cuts it, exactly as it cuts the shipped
         ! grid.  Falling back to the distributed file keeps a tree without
         ! the stored tables working, and changes which optics the model is.
         rok = .false.
         if (len_trim(q_h5) > 0) then
            h5comp = trim(gt(ic))
            if (trim(qset) == 'mie_d03') h5comp = trim(gt(ic))//'_mie_d03'
            call read_sedust_grid(trim(q_h5), .true., lam_opt, k_lo, rok)
            if (rok) then
               call read_sedust_qtable(trim(q_h5), trim(h5comp), .true., a_opt, &
                                       qe_h5, qa, qs, gg, rho_h5, rok)
               if (rok) then
                  nwave = size(lam_opt);  nsize = size(a_opt)
                  if (rho_h5 > 0.0_wp) rho = rho_h5
                  deallocate(qe_h5)
               end if
            end if
            if (.not. rok) then
               if (allocated(lam_opt)) deallocate(lam_opt)
               if (allocated(a_opt))   deallocate(a_opt)
            end if
            if (rok .and. ic == 1 .and. sed_verbose) write(*,'(a,a,a,a)') &
               ' build_zubko: ', trim(qset), ' optics read from ', trim(q_h5)
         end if
         if (.not. rok) then
            call load_q_component(trim(data_dir)//'q_zubko_'//trim(gt(ic))//'_euv.dat', &
                                  nwave, nsize, lam_opt, a_opt, qa, qs, gg, rok, rho=rho)
            ! Only on the route that was actually taken: this line was printed
            ! after a successful HDF5 read as well, and named a directory that
            ! run had not opened.
            if (rok .and. sed_verbose) write(*,'(a,a,a,a)') &
               ' build_zubko: ', trim(gt(ic)), &
               ' optics read from the stored tables under ', trim(data_dir)
         end if
         if (.not. rok) then
            ! Per COMPONENT, not once for the first: a table that is missing,
            ! stale or truncated takes only its own component down the fallback,
            ! and announcing that for ic = 1 alone hid exactly such a case --
            ! all three q_zubko_*_euv.dat in this tree were truncated at ~8% of
            ! their rows and the run said nothing.
            if (sed_verbose) write(*,'(a,a,a)') &
               ' build_zubko: ', trim(gt(ic)), &
               ' optics read from the distributed ZDA optics tables'
            if (present(status)) then
               call read_zubko_optics(trim(data_dir)//trim(optf), nsize, nwave, &
                                      a_opt, lam_opt, qa, qs, rho, ok=rok, gpar=gg)
               if (.not. rok) then;  status = 3;  return;  end if
            else
               call read_zubko_optics(trim(data_dir)//trim(optf), nsize, nwave, &
                                      a_opt, lam_opt, qa, qs, rho, gpar=gg)
            end if
         end if

         ! lam_min is a COVERAGE requirement, not a cut: this model cannot
         ! solve a wavelength its tables do not carry, so a floor below them is
         ! refused rather than answered with a frozen or extrapolated optic.
         if (present(lam_min)) then
            if (lam_min > 0.0_wp .and. lam_min < lam_opt(1)) then
               if (present(status)) then
                  status = 7;  return
               else
                  write(*,'(a,es12.5,a,es12.5)') ' build_zubko: lam_min ', lam_min, &
                     ' um is shorter than the optics table, which starts at ', lam_opt(1)
                  stop 1
               end if
            end if
         end if

         ! Narrow the grid to the requested coverage before anything downstream
         ! sees it, so the shared-grid checks, the cross sections and the
         ! Planck-averaged opacities all follow from the same axis.  Pure row
         ! selection on the tables: no value is recomputed.
         if (.not. wide) then
            k_lo = lyman_index(lam_opt)
            if (k_lo > 1) then
               nwave   = nwave - k_lo + 1
               lam_opt = lam_opt(k_lo:)
               qa      = qa(k_lo:, :)
               qs      = qs(k_lo:, :)
               gg      = gg(k_lo:, :)
            end if
         end if

         ! The endpoint dln(a) below reads a_opt(2), so demand at least 2 radii.
         if (nsize < 2) then
            if (present(status)) then
               status = 4;  return
            else
               write(*,'(a,i0,a,i0)') ' build_zubko: component ', ic, &
                  ' needs >= 2 radii, got ', nsize
               stop 1
            end if
         end if

         ! On the first component, fix the shared lambda + T grids (globals)
         ! and the calc_P setup. All three components share the lambda grid.
         if (ic == 1) then
            NLAM = nwave;  NT = NT_in
            ! This model's grid is the optics file's own, narrowed by lam_min if
            ! one was given; nothing is ever prepended below it.
            n_lam_euv = 0
            if (allocated(lam)) deallocate(lam, T_first, log_T_first)
            allocate(lam(NLAM), T_first(NT), log_T_first(NT))
            lam = lam_opt
            do jt = 1, NT
               t = log(T_lo) + (log(T_hi)-log(T_lo))*real(jt-1,wp)/real(NT-1,wp)
               T_first(jt) = exp(t)
            end do
            log_T_first = log(T_first)
            call p_sub_setup(lam)
            m%NLAM = NLAM;  m%NT = NT;  m%NA = nsize
            m%lam = lam;  m%T_first = T_first;  m%log_T_first = log_T_first
            allocate(m%aeff(0))           ! grids held per population; model aeff unused
         else
            ! All components must share the lambda grid fixed on component 1.
            if (nwave /= NLAM) then
               if (present(status)) then
                  status = 4;  return
               else
                  write(*,'(a,i0,a,i0,a,i0)') ' build_zubko: component ', ic, &
                     ' wavelength count ', nwave, ' /= ', NLAM
                  stop 1
               end if
            end if
            wdev = maxval(abs(lam_opt - lam) / lam)
            if (wdev > 1.0e-6_wp) then
               if (present(status)) then
                  status = 4;  return
               else
                  write(*,'(a,i0,a,es12.4)') ' build_zubko: component ', ic, &
                     ' wavelength grid mismatch, max rel dev = ', wdev
                  stop 1
               end if
            end if
         end if

         ! --- component-by-component working set in the module globals (scratch) ---
         NA = nsize
         if (allocated(Cabs)) deallocate(Cabs, kappB_first, kappCMB)
         if (allocated(Csca)) deallocate(Csca, gsca_ad)
         allocate(Cabs(NLAM, nsize), kappB_first(NT, nsize), kappCMB(nsize))
         allocate(Csca(NLAM, nsize), gsca_ad(NLAM, nsize))
         do ja = 1, nsize
            do jw = 1, NLAM
               Cabs(jw, ja)    = qa(jw, ja) * PI * (a_opt(ja)*UM2CM)**2   ! cm^2
               Csca(jw, ja)    = qs(jw, ja) * PI * (a_opt(ja)*UM2CM)**2   ! cm^2
               gsca_ad(jw, ja) = gg(jw, ja)
            end do
         end do
         call build_kappB()         ! Cabs, lam, T_first, NA -> kappB_first
         call build_kappCMB()       ! -> kappCMB

         ! --- enthalpy H(T,a) = u_spec(T) * rho * (4pi/3) a_cm^3 ---
         if (present(status)) then
            call read_zubko_calor(trim(data_dir)//trim(comps(ic)%calor), ntc, Tcal, Ucal, Ccal, ok=rok)
            if (.not. rok) then;  status = 5;  return;  end if
         else
            call read_zubko_calor(trim(data_dir)//trim(comps(ic)%calor), ntc, Tcal, Ucal, Ccal)
         end if
         vol_fac = (4.0_wp/3.0_wp) * PI
         allocate(Hcol(NT))
         block
            real(wp), allocatable :: Hmat(:,:)
            allocate(Hmat(NT, nsize))
            do ja = 1, nsize
               mass = rho * vol_fac * (a_opt(ja)*UM2CM)**3       ! g
               do jt = 1, NT
                  uspec = uspec_interp(Tcal, Ucal, ntc, T_first(jt))
                  Hmat(jt, ja) = uspec * mass                    ! erg
               end do
            end do
            ! --- assemble the population ---
            m%pops(ic)%grain_type = gt(ic)
            m%pops(ic)%out_channel = ic
            ! Bulk density of this component, as its own ZDA optics file
            ! declares it ("Density in gr/cm^3") -- the same number the
            ! enthalpy above turns into a grain mass, so the model's dust mass
            ! per H and its heat capacity refer to one and the same material.
            m%pops(ic)%rho_bulk = rho
            m%pops(ic)%Cabs    = Cabs
            m%pops(ic)%Csca    = Csca
            m%pops(ic)%gsca    = gsca_ad
            m%pops(ic)%kappB   = kappB_first
            m%pops(ic)%log_kappB = log(max(kappB_first, tiny(0.0_wp)))
            m%pops(ic)%H       = Hmat
            m%pops(ic)%log_H   = log(max(Hmat, tiny(0.0_wp)))
            m%pops(ic)%kappCMB = kappCMB
            deallocate(Hmat)
         end block
         deallocate(Hcol)

         ! --- size distribution: dn per bin from the ZDA formula ---
         ! dn_bin[1/H] = (dn/da) * a * dln(a) = f_formula(a) * a_um * dln(a)
         ! (the cm<->um unit factors cancel in (dn/da)*da on a log grid).
         block
            real(wp), allocatable :: dn(:)
            allocate(dn(nsize))
            do ja = 1, nsize
               ! Trapezoid in ln a: the interior bins get the full central
               ! difference and the two ENDPOINTS GET HALF, the same weights
               ! sed_init_dl07's bin_da uses.  A full step at an endpoint
               ! would place half a bin BEYOND the end of the grid, and where
               ! the grid begins at the size distribution's own lower cutoff
               ! -- which is where the ZDA tables begin -- that half bin lies
               ! outside the distribution's support and the steeply rising
               ! integrand there makes it count.
               if (ja == 1) then
                  dlna = 0.5_wp * log(a_opt(2)/a_opt(1))
               else if (ja == nsize) then
                  dlna = 0.5_wp * log(a_opt(nsize)/a_opt(nsize-1))
               else
                  dlna = 0.5_wp * log(a_opt(ja+1)/a_opt(ja-1))
               end if
               dn(ja) = zda_gofa(comps(ic), a_opt(ja)) * a_opt(ja) * dlna
            end do
            m%pops(ic)%dn = dn
            deallocate(dn)
         end block

         m%pops(ic)%aeff = a_opt        ! [um] radii of this component (needed by 'qm')

         deallocate(a_opt, lam_opt, qa, qs, gg, Tcal, Ucal, Ccal)
      end do

      call load_model_extinction_table(m, sed_data_path(zubko_kext_default(qset)), kext_path, kext_ok, &
                                       default_h5=sed_data_path(KEXT_H5_ZUBKO), &
                                       h5_group='kext'//trim(zubko_kext_tag(qset)))
      if (.not. kext_ok) then
         if (present(status)) then
            status = 6;  return
         else if (present(kext_path)) then
            write(*,'(a)') ' build_zubko: cannot read the extinction table '//trim(kext_path)
            stop 1
         end if
      end if
   end subroutine build_zubko


   ! Linear interpolation of specific enthalpy u_spec(T) [erg/gm], with linear
   ! extrapolation above the table's Tmax (the ZDA convention "extrapolate high
   ! with T") and clamping below Tmin.
   pure function uspec_interp(Tt, Ut, n, T) result(u)
      real(wp), intent(in) :: Tt(:), Ut(:)
      integer,  intent(in) :: n
      real(wp), intent(in) :: T
      real(wp) :: u, f
      integer  :: lo, hi, mid
      if (T <= Tt(1)) then
         u = Ut(1) * (T / Tt(1))            ! ~ low-T clamp (avoids <0)
         return
      end if
      if (T >= Tt(n)) then
         ! extrapolate high with u ~ T (slope from the last interval)
         f = (Ut(n) - Ut(n-1)) / (Tt(n) - Tt(n-1))
         u = Ut(n) + f * (T - Tt(n))
         return
      end if
      lo = 1;  hi = n
      do while (hi - lo > 1)
         mid = (lo + hi) / 2
         if (Tt(mid) <= T) then;  lo = mid;  else;  hi = mid;  end if
      end do
      f = (T - Tt(lo)) / (Tt(hi) - Tt(lo))
      u = (1.0_wp - f) * Ut(lo) + f * Ut(hi)
   end function uspec_interp


   ! Generic file-defined model loader. Reads a small descriptor:
   !   name = <model name>          (optional)
   !   pop: <grain_type> <channel> <optics_file> <dnda_file> <calor_file> <rho>
   !   pop: ...                     (one line per population)
   ! Each population's optics is a ZDA optics table, the size distribution is a
   ! 2-column a[um] dn/da[cm^-1 H^-1] table, and the enthalpy a specific-heat
   ! calorimetry table; all files are sought under data_dir. This is the
   ! data-driven path (build_astrodust/dl07/zubko are the coded builders).
   subroutine build_from_files(m, descriptor_path, data_dir, NT_in, T_lo, T_hi, status, &
                               kext_path, lam_min, include_euv)
      type(dust_model_t), intent(out) :: m
      character(len=*),   intent(in)  :: descriptor_path, data_dir
      integer,            intent(in)  :: NT_in
      real(wp),           intent(in)  :: T_lo, T_hi
      ! Optional status (0 = success, non-zero = model build failed). When
      ! present, a bad input is reported through it instead of stopping; when
      ! absent the build stops on error (CLI behavior).
      !   status = 1  descriptor open failed
      !   status = 2  too many pop: lines (MAXP exceeded)
      !   status = 3  a population has an invalid channel
      !   status = 4  no pop: lines found
      !   status = 5  a population's optics read failed
      !   status = 6  a population's size/wavelength grid is inconsistent
      !   status = 7  a population's size-distribution read failed
      !   status = 8  a population's calorimetry read failed
      !   status = 9  the extinction table named by kext_path could not be read
      !   status = 10 lam_min shorter than the model's own optics tables
      integer, optional,  intent(out) :: status
      ! Size-integrated extinction table dust_extinction serves this model's
      ! scalar optics from. A file-defined model has NO default: its product is
      ! named after the model, which only the descriptor knows, so a host that
      ! wants extinction out of such a model must name the file. Omitting it
      ! leaves m%kext_n = 0 and dust_extinction with nothing to return.
      character(len=*), optional, intent(in) :: kext_path
      ! Shortest wavelength [um] the model must COVER.  A file-defined model is
      ! its tables, with no dielectric function behind them, so it meets this
      ! when they already reach and refuses it (status 10) when they do not.
      ! It never narrows the grid; see build_zubko, which takes it the same way.
      real(wp), optional, intent(in) :: lam_min
      ! Whether to keep the ionizing part of the tables' own grid; default
      ! .true.  .false. cuts at lyman_index, as it does for every other model.
      logical, optional, intent(in) :: include_euv
      logical               :: kext_ok, wide

      integer, parameter :: MAXP = 16
      character(len=8)   :: p_gt(MAXP)
      integer            :: p_ch(MAXP)
      character(len=64)  :: p_opt(MAXP), p_dn(MAXP), p_cal(MAXP)
      real(wp)           :: p_rho(MAXP)
      integer            :: npop, u, ios, ip, jt, ja, jw, nsize, nwave, ntc, ndn, nchan, ic, nline
      integer            :: k_lo
      real(wp)           :: t, rho, mass, vf, dlna, uspec, fa, loga, wdev
      logical            :: rok
      character(len=256) :: line
      real(wp), allocatable :: a_opt(:), lam_opt(:), qa(:,:), qs(:,:), gg(:,:)
      real(wp), allocatable :: a_dn(:), f_dn(:), la_dn(:), lf_dn(:), Tc(:), Uc(:), Cc(:)

      if (present(status)) status = 0
      wide = .true.;  if (present(include_euv)) wide = include_euv

      npop = 0;  nline = 0;  m%name = 'file_model'
      active_build_id = active_build_id + 1;  m%build_id = active_build_id
      open(newunit=u, file=trim(descriptor_path), status='old', action='read', iostat=ios)
      if (ios /= 0) then
         if (present(status)) then
            status = 1;  return
         else
            write(*,'(a,a)') ' build_from_files: cannot open ', trim(descriptor_path); stop 1
         end if
      end if
      do
         read(u,'(a)', iostat=ios) line;  if (ios /= 0) exit
         nline = nline + 1
         line = adjustl(line)
         if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
         if (line(1:4) == 'pop:') then
            if (npop >= MAXP) then
               if (present(status)) then
                  close(u);  status = 2;  return
               else
                  write(*,'(a,i0,a,i0)') ' build_from_files: too many pop: lines (max ', &
                     MAXP, ') at input line ', nline
                  stop 1
               end if
            end if
            npop = npop + 1
            read(line(5:), *) p_gt(npop), p_ch(npop), p_opt(npop), p_dn(npop), &
                              p_cal(npop), p_rho(npop)
            if (p_ch(npop) < 1) then
               if (present(status)) then
                  close(u);  status = 3;  return
               else
                  write(*,'(a,i0,a,i0)') ' build_from_files: population ', npop, &
                     ' has invalid channel ', p_ch(npop)
                  stop 1
               end if
            end if
         else if (index(line,'name') > 0 .and. index(line,'=') > 0) then
            m%name = trim(adjustl(line(index(line,'=')+1:)))
         end if
      end do
      close(u)
      if (npop == 0) then
         if (present(status)) then
            status = 4;  return
         else
            write(*,'(a)') ' build_from_files: no pop: lines found'; stop 1
         end if
      end if
      nchan = maxval(p_ch(1:npop))

      m%use_induced_emission = use_induced_emission
      m%stoch_method = stoch_method
      m%n_channel = nchan
      allocate(m%channel_name(nchan))
      do ic = 1, nchan
         write(m%channel_name(ic), '(a,i0)') 'CH', ic
      end do
      allocate(m%pops(npop))

      do ip = 1, npop
         if (present(status)) then
            call read_zubko_optics(trim(data_dir)//trim(p_opt(ip)), nsize, nwave, &
                                   a_opt, lam_opt, qa, qs, rho, ok=rok, gpar=gg)
            if (.not. rok) then;  status = 5;  return;  end if
         else
            call read_zubko_optics(trim(data_dir)//trim(p_opt(ip)), nsize, nwave, &
                                   a_opt, lam_opt, qa, qs, rho, gpar=gg)
         end if
         if (p_rho(ip) > 0.0_wp) rho = p_rho(ip)         ! descriptor rho overrides file

         ! Narrow the grid to the requested coverage before anything downstream
         ! sees it.  Pure row selection on the tables: no value is recomputed.
         ! lam_min is a COVERAGE requirement, not a cut; see build_zubko.
         if (present(lam_min)) then
            if (lam_min > 0.0_wp .and. lam_min < lam_opt(1)) then
               if (present(status)) then
                  status = 10;  return
               else
                  write(*,'(a,es12.5,a,es12.5)') ' build_from_files: lam_min ', lam_min, &
                     ' um is shorter than the optics table, which starts at ', lam_opt(1)
                  stop 1
               end if
            end if
         end if

         if (.not. wide) then
            k_lo = lyman_index(lam_opt)
            if (k_lo > 1) then
               nwave   = nwave - k_lo + 1
               lam_opt = lam_opt(k_lo:)
               qa      = qa(k_lo:, :)
               qs      = qs(k_lo:, :)
               gg      = gg(k_lo:, :)
            end if
         end if

         ! The endpoint dln(a) below reads a_opt(2), so demand at least 2 radii.
         if (nsize < 2) then
            if (present(status)) then
               status = 6;  return
            else
               write(*,'(a,i0,a,i0)') ' build_from_files: population ', ip, &
                  ' needs >= 2 radii, got ', nsize
               stop 1
            end if
         end if

         if (ip == 1) then
            NLAM = nwave;  NT = NT_in
            ! This model's grid is the optics files' own, narrowed by include_euv if
            ! one was given; nothing is ever prepended below them.
            n_lam_euv = 0
            if (allocated(lam)) deallocate(lam, T_first, log_T_first)
            allocate(lam(NLAM), T_first(NT), log_T_first(NT))
            lam = lam_opt
            do jt = 1, NT
               t = log(T_lo) + (log(T_hi)-log(T_lo))*real(jt-1,wp)/real(NT-1,wp)
               T_first(jt) = exp(t)
            end do
            log_T_first = log(T_first)
            call p_sub_setup(lam)
            m%NLAM = NLAM;  m%NT = NT;  m%NA = nsize
            m%lam = lam;  m%T_first = T_first;  m%log_T_first = log_T_first
            allocate(m%aeff(0))
         else
            ! Every population must share the lambda grid fixed on population 1.
            if (nwave /= NLAM) then
               if (present(status)) then
                  status = 6;  return
               else
                  write(*,'(a,i0,a,i0,a,i0)') ' build_from_files: population ', ip, &
                     ' wavelength count ', nwave, ' /= ', NLAM
                  stop 1
               end if
            end if
            wdev = maxval(abs(lam_opt - lam) / lam)
            if (wdev > 1.0e-6_wp) then
               if (present(status)) then
                  status = 6;  return
               else
                  write(*,'(a,i0,a,es12.4)') ' build_from_files: population ', ip, &
                     ' wavelength grid mismatch, max rel dev = ', wdev
                  stop 1
               end if
            end if
         end if

         NA = nsize
         if (allocated(Cabs)) deallocate(Cabs, kappB_first, kappCMB)
         if (allocated(Csca)) deallocate(Csca, gsca_ad)
         allocate(Cabs(NLAM, nsize), kappB_first(NT, nsize), kappCMB(nsize))
         allocate(Csca(NLAM, nsize), gsca_ad(NLAM, nsize))
         do ja = 1, nsize
            do jw = 1, NLAM
               Cabs(jw, ja)    = qa(jw, ja) * PI * (a_opt(ja)*UM2CM)**2
               Csca(jw, ja)    = qs(jw, ja) * PI * (a_opt(ja)*UM2CM)**2
               gsca_ad(jw, ja) = gg(jw, ja)
            end do
         end do
         call build_kappB();  call build_kappCMB()

         if (present(status)) then
            call read_dnda_table(trim(data_dir)//trim(p_dn(ip)), ndn, a_dn, f_dn, ok=rok)
            if (.not. rok) then;  status = 7;  return;  end if
         else
            call read_dnda_table(trim(data_dir)//trim(p_dn(ip)), ndn, a_dn, f_dn)
         end if
         allocate(la_dn(ndn), lf_dn(ndn))
         la_dn = log(a_dn);  lf_dn = log(max(f_dn, tiny(0.0_wp)))
         if (present(status)) then
            call read_zubko_calor(trim(data_dir)//trim(p_cal(ip)), ntc, Tc, Uc, Cc, ok=rok)
            if (.not. rok) then;  deallocate(la_dn, lf_dn);  status = 8;  return;  end if
         else
            call read_zubko_calor(trim(data_dir)//trim(p_cal(ip)), ntc, Tc, Uc, Cc)
         end if
         vf = (4.0_wp/3.0_wp) * PI

         block
            real(wp), allocatable :: dn(:), Hmat(:,:)
            allocate(dn(nsize), Hmat(NT, nsize))
            do ja = 1, nsize
               loga = log(a_opt(ja))
               if (a_opt(ja) < a_dn(1) .or. a_opt(ja) > a_dn(ndn)) then
                  fa = 0.0_wp                              ! outside size-dist range
               else
                  call interp(la_dn, lf_dn, loga, fa);  fa = exp(fa)
               end if
               ! Trapezoid in ln a: the interior bins get the full central
               ! difference and the two ENDPOINTS GET HALF, the same weights
               ! sed_init_dl07's bin_da uses.  A full step at an endpoint
               ! would place half a bin BEYOND the end of the grid, and where
               ! the grid begins at the size distribution's own lower cutoff
               ! -- which is where the ZDA tables begin -- that half bin lies
               ! outside the distribution's support and the steeply rising
               ! integrand there makes it count.
               if (ja == 1) then
                  dlna = 0.5_wp * log(a_opt(2)/a_opt(1))
               else if (ja == nsize) then
                  dlna = 0.5_wp * log(a_opt(nsize)/a_opt(nsize-1))
               else
                  dlna = 0.5_wp * log(a_opt(ja+1)/a_opt(ja-1))
               end if
               ! dn[1/H] = (dn/da)[cm^-1] * da[cm] = f * (a_cm) * dln(a)
               dn(ja) = fa * (a_opt(ja)*UM2CM) * dlna
               mass = rho * vf * (a_opt(ja)*UM2CM)**3
               do jt = 1, NT
                  uspec = uspec_interp(Tc, Uc, ntc, T_first(jt))
                  Hmat(jt, ja) = uspec * mass
               end do
            end do
            m%pops(ip)%grain_type = p_gt(ip)
            m%pops(ip)%out_channel = p_ch(ip)
            ! Bulk density: the descriptor's rho if it gave one, otherwise the
            ! optics file's own declared density -- whichever the enthalpy just
            ! above used, so mass per H and heat capacity stay consistent.
            m%pops(ip)%rho_bulk = rho
            m%pops(ip)%dn = dn
            m%pops(ip)%Cabs = Cabs
            m%pops(ip)%Csca = Csca
            m%pops(ip)%gsca = gsca_ad
            m%pops(ip)%kappB = kappB_first
            m%pops(ip)%log_kappB = log(max(kappB_first, tiny(0.0_wp)))
            m%pops(ip)%H = Hmat
            m%pops(ip)%log_H = log(max(Hmat, tiny(0.0_wp)))
            m%pops(ip)%kappCMB = kappCMB
            deallocate(dn, Hmat)
         end block

         m%pops(ip)%aeff = a_opt        ! [um] radii of this population (needed by 'qm')

         deallocate(a_opt, lam_opt, qa, qs, gg, a_dn, f_dn, la_dn, lf_dn, Tc, Uc, Cc)
      end do

      call load_model_extinction_table(m, '', kext_path, kext_ok)
      if (.not. kext_ok) then
         if (present(status)) then
            status = 9;  return
         else if (present(kext_path)) then
            write(*,'(a)') ' build_from_files: cannot read the extinction table '// &
                 trim(kext_path)
            stop 1
         end if
      end if
   end subroutine build_from_files

   ! =====================================================================
   ! A model defined by DustEM input files -- THEMIS (Jones et al. 2013 for
   ! the carbon grains, Koehler et al. 2014 for the silicates) and Guillet
   ! et al. (2018) Model D, the two models Hensley & Draine (2023), ApJ 948,
   ! 55, Sec. 6.2.2 and Fig. 16 compare their astrodust model against.
   !
   ! grain_path names the GRAIN_*.DAT that IS the model definition; data_dir
   ! is the directory the files it implies hang under,
   !
   !   <data_dir>oprop/LAMBDA.DAT           the common wavelength grid
   !   <data_dir>oprop/Q_<gtype>.DAT        Qabs, Qsca
   !   <data_dir>oprop/G_<gtype>.DAT        <cos theta>, where one is shipped
   !   <data_dir>hcap/C_<gtype>.DAT         heat capacity per unit volume
   !
   ! keeping the DustEM subdirectory names so that the paths a GRAIN_*.DAT
   ! implies resolve the way they do inside DustEM itself.
   !
   ! Two things differ from a table-defined model such as build_from_files.
   ! First, the size grid is the MODEL's own -- nsize points spaced evenly in
   ! ln a between the amin and amax of the GRAIN line -- and the optics tables,
   ! which have a radius grid of their own, are interpolated onto it.  Second,
   ! the optics tables are shared by every model that names the same gtype, so
   ! GRAIN_*.DAT is the only file here that belongs to one model alone.
   !
   ! Each population becomes one output channel, named by its gtype, which is
   ! the layout of DustEM's own EXT_/SED_ products.  Channel names are 16
   ! characters, so a longer gtype is carried truncated there; the population's
   ! own optics and size distribution are unaffected.
   subroutine build_dustem(m, grain_path, data_dir, NT_in, T_lo, T_hi, status, &
                           kext_path, lam_min, include_euv, gsca_missing, qtable_path)
      type(dust_model_t), intent(out) :: m
      character(len=*),   intent(in)  :: grain_path, data_dir
      integer,            intent(in)  :: NT_in
      real(wp),           intent(in)  :: T_lo, T_hi
      ! Optional status (0 = success, non-zero = model build failed). When
      ! present, a bad input is reported through it instead of stopping; when
      ! absent the build stops on error (CLI behavior).
      !   status = 1  GRAIN_*.DAT could not be read, or names a size-distribution
      !               keyword this code does not implement
      !   status = 2  oprop/LAMBDA.DAT could not be read
      !   status = 3  a population's Q_<gtype>.DAT could not be read
      !   status = 4  a population's size distribution could not be formed
      !   status = 5  a population's model radii fall outside its optics tables
      !   status = 6  a population's hcap/C_<gtype>.DAT could not be read, or
      !               does not cover that population's model radii
      !   status = 7  a population's G_<gtype>.DAT exists but could not be read
      !   status = 8  the requested extinction table failed to load
      !   status = 9  lam_min is shorter than the DustEM wavelength grid
      integer, optional,  intent(out) :: status
      ! Optional precomputed size-integrated extinction curve for
      ! dust_extinction to serve; see build_astrodust.  Omitted, the model's
      ! own default: /kext of the product below, then <data_dir>kext_<model>.dat.
      ! Neither being there leaves the model with no extinction to serve and
      ! does NOT fail the build, which is what every other model does too.
      character(len=*), optional, intent(in) :: kext_path
      ! Shortest wavelength [um] the model must COVER.  This model IS its
      ! tables, with no dielectric function behind them, so it meets the
      ! requirement when they already reach and refuses it (status 9) when they
      ! do not.  It never narrows the grid; build_zubko takes it the same way.
      real(wp), optional, intent(in) :: lam_min
      ! Whether to keep the ionizing part of the tables' own grid; default
      ! .true.  .false. cuts at lyman_index, as it does for every other model.
      logical, optional, intent(in) :: include_euv
      ! The gtypes, comma separated, for which the distribution ships no
      ! G_<gtype>.DAT, so that a caller can say which populations are the
      ! reason m%gsca_complete came back .false.  Blank when every scattering
      ! population carries an asymmetry parameter.
      character(len=*), optional, intent(out) :: gsca_missing
      ! The model's HDF5 product, when the caller has one.  It holds this
      ! model's own wavelength axis and, for each population, the cross
      ! sections ALREADY interpolated onto that population's model radii --
      ! the very numbers the text route below computes, so the two routes must
      ! agree, and test_dustem_product.x measures that they do.  Blank or
      ! omitted reads the DustEM text tables under data_dir, and so does a
      ! product whose grids do not match this GRAIN_*.DAT.  It is an argument
      ! rather than module state so that a program writing /kext into a product
      ! takes its optics from the /qtable of that same product.
      character(len=*), optional, intent(in) :: qtable_path

      type(dustem_pop_t) :: pops(DUSTEM_MAXPOP)
      real(wp) :: G0
      integer  :: npop, ip, ja, jt, jw, nwave, nwave_full, nsize, ntab, k_lo, nchan
      integer  :: nsize_c, ntemp_c
      real(wp) :: t, dlna
      logical  :: rok, kok, wide, has_g, use_h5, qfrom_h5, hg
      integer  :: k0, k1, k_h5
      real(wp) :: rho_h5
      character(len=512) :: qpath, gpath, cpath, q_h5, mname
      character(len=140) :: pname
      real(wp), allocatable :: lam_full(:), a_tab(:), qa_tab(:,:), qs_tab(:,:), gg_tab(:,:)
      real(wp), allocatable :: a_cm(:), lna(:), dnda(:), a_um(:)
      real(wp), allocatable :: ac_tab(:), logT_c(:), logC_c(:,:)
      real(wp), allocatable :: qa(:,:), qs(:,:), gfac(:,:)
      real(wp), allocatable :: lam_h5(:), a_h5(:), qe_h5(:,:), qa_h5(:,:), qs_h5(:,:)
      real(wp), allocatable :: gg_h5(:,:)

      if (present(status)) status = 0
      if (present(gsca_missing)) gsca_missing = ''
      wide = .true.;  if (present(include_euv)) wide = include_euv
      q_h5 = '';      if (present(qtable_path)) q_h5 = qtable_path

      ! The model's own directory names it.  Every model in this tree lives in
      ! data/<model>/, so the last component of data_dir IS <model>, and the
      ! default extinction curve and the product are named from it the same way
      ! every other model's are.  A blank data_dir leaves the generic name.
      mname = 'dustem'
      k1 = len_trim(data_dir)
      if (k1 > 0) then
         if (data_dir(k1:k1) == '/') k1 = k1 - 1
         k0 = index(data_dir(1:max(k1,1)), '/', back=.true.) + 1
         if (k1 >= k0) mname = data_dir(k0:k1)
      end if
      m%name = trim(mname)
      active_build_id = active_build_id + 1;  m%build_id = active_build_id

      ! ---- the model definition -----------------------------------------
      call read_dustem_grain(trim(grain_path), G0, npop, pops, ok=rok)
      if (.not. rok) then
         if (present(status)) then
            status = 1;  return
         else
            write(*,'(a,a)') ' build_dustem: cannot use ', trim(grain_path);  stop 1
         end if
      end if

      ! ---- the common wavelength grid ------------------------------------
      call read_dustem_wavelengths(trim(data_dir)//'oprop/LAMBDA.DAT', nwave_full, &
                                   lam_full, ok=rok)
      if (.not. rok) then
         if (present(status)) then
            status = 2;  return
         else
            write(*,'(a,a)') ' build_dustem: cannot read ', &
               trim(data_dir)//'oprop/LAMBDA.DAT'
            stop 1
         end if
      end if

      ! ---- the stored optics, when the caller named the product -----------
      ! Usable only if the product is THIS model's: it carries one wavelength
      ! axis, and a set of cross sections filed against a different one is not
      ! this model's optics whatever the group is called.  A mismatch here
      ! reads the DustEM text tables instead, which is what a tree without the
      ! product does anyway.
      use_h5 = .false.
      if (len_trim(q_h5) > 0) then
         call read_sedust_grid(trim(q_h5), .true., lam_h5, k_h5, rok)
         if (rok) then
            if (size(lam_h5) == nwave_full) &
               use_h5 = maxval(abs(lam_h5/lam_full - 1.0_wp)) <= 1.0e-12_wp
            deallocate(lam_h5)
         end if
         if (use_h5 .and. sed_verbose) write(*,'(a,a)') &
            ' build_dustem: optics read from ', trim(q_h5)
      end if

      ! lam_min is a COVERAGE requirement, not a cut; see build_zubko.
      if (present(lam_min)) then
         if (lam_min > 0.0_wp .and. lam_min < lam_full(1)) then
            if (present(status)) then
               status = 9;  return
            else
               write(*,'(a,es12.5,a,es12.5)') ' build_dustem: lam_min ', lam_min, &
                  ' um is shorter than the DustEM grid, which starts at ', lam_full(1)
               stop 1
            end if
         end if
      end if

      ! Narrow the grid to the requested coverage before anything downstream
      ! sees it.  Pure row selection on the tables: no value is recomputed.
      k_lo = 1
      if (.not. wide) k_lo = lyman_index(lam_full)

      NLAM = nwave_full - k_lo + 1
      nwave = NLAM
      NT    = NT_in
      ! This model's grid is its own wavelength table, narrowed by include_euv
      ! if one was asked for.  Nothing is ever prepended below it.
      n_lam_euv = 0
      if (allocated(lam)) deallocate(lam, T_first, log_T_first)
      allocate(lam(NLAM), T_first(NT), log_T_first(NT))
      lam = lam_full(k_lo:nwave_full)
      do jt = 1, NT
         t = log(T_lo) + (log(T_hi)-log(T_lo))*real(jt-1,wp)/real(NT-1,wp)
         T_first(jt) = exp(t)
      end do
      log_T_first = log(T_first)
      call p_sub_setup(lam)

      m%use_induced_emission = use_induced_emission
      m%stoch_method = stoch_method
      m%NLAM = NLAM;  m%NT = NT;  m%NA = 0
      m%lam = lam;  m%T_first = T_first;  m%log_T_first = log_T_first
      allocate(m%aeff(0))

      ! One channel per population, as DustEM's own products are laid out.
      nchan = npop
      m%n_channel = nchan
      allocate(m%channel_name(nchan))
      do ip = 1, npop
         m%channel_name(ip) = pops(ip)%gtype
      end do
      allocate(m%pops(npop))
      m%gsca_complete = .true.

      ! ---- population by population --------------------------------------
      do ip = 1, npop
         nsize = pops(ip)%nsize

         ! --- size distribution, the DustEM formula ---
         call dustem_size_distribution(pops(ip), a_cm, lna, dnda, ok=rok)
         if (.not. rok) then
            if (present(status)) then
               status = 4;  return
            else
               write(*,'(a,a)') ' build_dustem: bad size distribution for ', &
                  trim(pops(ip)%gtype)
               stop 1
            end if
         end if
         allocate(a_um(nsize))
         a_um = a_cm / UM2CM

         ! --- optics, from the product when it carries this population -------
         ! The stored table is already on THIS population's model radii, so
         ! there is nothing to interpolate: the interpolation happened once,
         ! in calc_qtable.x, by the same optics_at_radii the text route below
         ! calls.  That interpolation is in radius alone, so cutting the
         ! wavelength axis and interpolating commute and the two routes give
         ! the same numbers bit for bit.  A group whose radii are not this
         ! model's is not this population's optics, and takes only itself down
         ! to the text tables.
         qfrom_h5 = .false.;  has_g = .false.
         if (use_h5) then
            pname = dustem_population_name(pops, npop, ip)
            call read_sedust_qtable(trim(q_h5), trim(pname), wide, a_h5, qe_h5, &
                                    qa_h5, qs_h5, gg_h5, rho_h5, rok, has_g=hg)
            if (rok) then
               if (size(a_h5) == nsize .and. size(qa_h5,1) == nwave) &
                  qfrom_h5 = maxval(abs(a_h5/a_um - 1.0_wp)) <= 1.0e-12_wp
            end if
            if (.not. qfrom_h5) then
               if (allocated(a_h5))  deallocate(a_h5)
               if (allocated(qe_h5)) deallocate(qe_h5)
               if (allocated(qa_h5)) deallocate(qa_h5)
               if (allocated(qs_h5)) deallocate(qs_h5)
               if (allocated(gg_h5)) deallocate(gg_h5)
            end if
         end if

         if (qfrom_h5) then
            allocate(qa(nwave, nsize), qs(nwave, nsize))
            qa = qa_h5;  qs = qs_h5
            has_g = hg
            if (has_g) then
               allocate(gfac(nwave, nsize));  gfac = gg_h5
            end if
            deallocate(a_h5, qe_h5, qa_h5, qs_h5, gg_h5)
         else
            ! --- optics on the table's own radius grid ---
            qpath = trim(data_dir)//'oprop/Q_'//trim(pops(ip)%gtype)//'.DAT'
            call read_dustem_qtable(trim(qpath), nwave_full, ntab, a_tab, qa_tab, qs_tab, ok=rok)
            if (.not. rok) then
               if (present(status)) then
                  status = 3;  return
               else
                  write(*,'(a,a)') ' build_dustem: cannot read ', trim(qpath);  stop 1
               end if
            end if

            ! --- onto the model's radius grid, linearly in radius as DustEM does ---
            allocate(qa(nwave, nsize), qs(nwave, nsize))
            call optics_at_radii(a_tab, qa_tab(k_lo:nwave_full,:), a_um, qa, ok=rok, &
                                 what='Q_abs')
            if (rok) call optics_at_radii(a_tab, qs_tab(k_lo:nwave_full,:), a_um, qs, ok=rok, &
                                          what='Q_sca')
            if (.not. rok) then
               if (present(status)) then
                  status = 5;  return
               else
                  write(*,'(a,a)') ' build_dustem: optics do not cover the sizes of ', &
                     trim(pops(ip)%gtype)
                  stop 1
               end if
            end if

            ! --- the asymmetry parameter, where the distribution ships one ---
            gpath = trim(data_dir)//'oprop/G_'//trim(pops(ip)%gtype)//'.DAT'
            inquire(file=trim(gpath), exist=has_g)
            if (has_g) then
               call read_dustem_gtable(trim(gpath), nwave_full, ntab, ac_tab, gg_tab, ok=rok)
               if (.not. rok) then
                  if (present(status)) then
                     status = 7;  return
                  else
                     write(*,'(a,a)') ' build_dustem: cannot read ', trim(gpath);  stop 1
                  end if
               end if
               allocate(gfac(nwave, nsize))
               call optics_at_radii(ac_tab, gg_tab(k_lo:nwave_full,:), a_um, gfac, ok=rok, &
                                    what='<cos theta>')
               if (.not. rok) then
                  if (present(status)) then
                     status = 5;  return
                  else
                     write(*,'(a,a)') ' build_dustem: the g table does not cover the sizes of ', &
                        trim(pops(ip)%gtype)
                     stop 1
                  end if
               end if
               deallocate(ac_tab, gg_tab)
            end if
         end if

         if (.not. has_g) then
            ! No G_ file in the distribution, and so no g dataset in the
            ! product either.  The population's gsca is left UNALLOCATED,
            ! which is the honest state: this model has no asymmetry parameter
            ! for it.  size_integrated_extinction reads m%gsca_complete and
            ! returns <cos theta> = 0 for the whole model rather than averaging
            ! over the populations that do carry g.
            m%gsca_complete = .false.
            if (present(gsca_missing)) then
               if (len_trim(gsca_missing) == 0) then
                  gsca_missing = trim(pops(ip)%gtype)
               else
                  gsca_missing = trim(gsca_missing)//', '//trim(pops(ip)%gtype)
               end if
            end if
         end if

         ! --- cross sections and the Planck integrals -----------------------
         NA = nsize
         if (allocated(Cabs)) deallocate(Cabs, kappB_first, kappCMB)
         if (allocated(Csca)) deallocate(Csca, gsca_ad)
         allocate(Cabs(NLAM, nsize), kappB_first(NT, nsize), kappCMB(nsize))
         allocate(Csca(NLAM, nsize), gsca_ad(NLAM, nsize))
         do ja = 1, nsize
            do jw = 1, NLAM
               Cabs(jw, ja) = qa(jw, ja) * PI * a_cm(ja)**2
               Csca(jw, ja) = qs(jw, ja) * PI * a_cm(ja)**2
            end do
         end do
         if (has_g) gsca_ad = gfac
         call build_kappB();  call build_kappCMB()

         ! --- enthalpy from the tabulated heat capacity --------------------
         cpath = trim(data_dir)//'hcap/C_'//trim(pops(ip)%gtype)//'.DAT'
         call read_dustem_heat_capacity(trim(cpath), nsize_c, ntemp_c, ac_tab, logT_c, logC_c, &
                                        ok=rok)
         if (.not. rok) then
            if (present(status)) then
               status = 6;  return
            else
               write(*,'(a,a)') ' build_dustem: cannot read ', trim(cpath);  stop 1
            end if
         end if

         block
            real(wp), allocatable :: dn(:), Hmat(:,:)
            allocate(dn(nsize), Hmat(NT, nsize))
            do ja = 1, nsize
               ! Bin width in ln a: the interior bins get the full central
               ! difference and the two ENDPOINTS GET HALF, which is the same
               ! convention sed_init_dl07's bin_da uses.  That half weight is
               ! what makes sum_a Cabs(a) dn(a) the very trapezoid in ln a that
               ! DustEM's own XINTEG2 performs over this same grid, so a
               ! comparison against a DustEM product is a comparison of the two
               ! models and not of two quadrature rules.  It also makes
               ! sum_a (4/3) pi a^3 rho dn(a) come out at m_p * (Mdust/MH)
               ! exactly, the mass per H the GRAIN line states.
               if (ja == 1) then
                  dlna = 0.5_wp * (lna(2) - lna(1))
               else if (ja == nsize) then
                  dlna = 0.5_wp * (lna(nsize) - lna(nsize-1))
               else
                  dlna = 0.5_wp * (lna(ja+1) - lna(ja-1))
               end if
               ! dn[1/H] = (dn/da)[cm^-1 H^-1] * da[cm] = (dn/da) * a_cm * dln(a)
               dn(ja) = dnda(ja) * a_cm(ja) * dlna
               call grain_enthalpy_from_heat_capacity(a_um(ja), ac_tab, logT_c, logC_c, &
                                                      T_first, Hmat(:, ja), ok=rok)
               if (.not. rok) then
                  deallocate(dn, Hmat)
                  if (present(status)) then
                     ! The CALORIMETRY does not cover this radius, which is a
                     ! different failure from the optics not covering it, and
                     ! is reported as one so that build_dust can name the
                     ! stage.  It used to share the optics code.
                     status = 6;  return
                  else
                     write(*,'(a,a)') ' build_dustem: the C table does not cover the ' // &
                        'sizes of ', trim(pops(ip)%gtype)
                     stop 1
                  end if
               end if
            end do
            ! The grain type the stochastic-heating solver's quantum branch
            ! keys on, from the material each DustEM gtype names.  It selects
            ! an atom count and a mode spectrum there and nothing else; the
            ! model's own optics and enthalpy, which the default solver uses,
            ! come from the tables above and do not consult it.
            m%pops(ip)%grain_type  = solver_grain_type(pops(ip)%gtype)
            m%pops(ip)%out_channel = ip
            m%pops(ip)%rho_bulk    = pops(ip)%rho
            m%pops(ip)%dn          = dn
            m%pops(ip)%Cabs        = Cabs
            m%pops(ip)%Csca        = Csca
            if (has_g) m%pops(ip)%gsca = gsca_ad
            m%pops(ip)%kappB     = kappB_first
            m%pops(ip)%log_kappB = log(max(kappB_first, tiny(0.0_wp)))
            m%pops(ip)%H         = Hmat
            m%pops(ip)%log_H     = log(max(Hmat, tiny(0.0_wp)))
            m%pops(ip)%kappCMB   = kappCMB
            deallocate(dn, Hmat)
         end block

         m%pops(ip)%aeff = a_um          ! [um] radii of this population
         ! Descriptive only: the populations of a DustEM model need not share a
         ! size grid (G18D's PAHs have 10 nodes and its large grains 25), so
         ! the model-level count is the largest of them.
         m%NA = max(m%NA, nsize)

         deallocate(a_cm, lna, dnda, a_um, qa, qs)
         if (allocated(a_tab))  deallocate(a_tab)
         if (allocated(qa_tab)) deallocate(qa_tab)
         if (allocated(qs_tab)) deallocate(qs_tab)
         if (allocated(gfac))   deallocate(gfac)
         deallocate(ac_tab, logT_c, logC_c)
      end do
      deallocate(lam_full)

      ! The model's own curve, in the model's own directory, with /kext of its
      ! product ahead of the text file -- the same two defaults, in the same
      ! order, that every other model has.  Neither being there is not a build
      ! failure: the model is complete for emission, and dust_extinction
      ! reports that it has nothing to serve.
      call load_model_extinction_table(m, trim(data_dir)//'kext_'//trim(mname)//'.dat', &
                                       kext_path, kok, &
                                       default_h5=trim(data_dir)//'sedust_'// &
                                                  trim(mname)//'.h5')
      if (.not. kok) then
         if (present(status)) then
            status = 8;  return
         else if (present(kext_path)) then
            write(*,'(a,a)') ' build_dustem: cannot read the extinction table ', &
                 trim(kext_path)
            stop 1
         end if
      end if
   end subroutine build_dustem


   ! Which of the solver's three grain materials a DustEM gtype is.  It is read
   ! only by the quantum stochastic-heating branch, which needs an atom count
   ! and a vibrational mode spectrum for the grain; the tabulated enthalpy that
   ! the default solver integrates comes from the model's own C_ file.
   !   CM20, amC*   amorphous carbon        -> 'gra'
   !   PAH*         polycyclic aromatics    -> 'pah'
   !   aPyM5, aOlM5, aSil*  amorphous silicates -> 'sil'
   ! A gtype outside this list gets 'sil', the solver's silicate mode spectrum.
   pure function solver_grain_type(gtype) result(gt)
      character(len=*), intent(in) :: gtype
      character(len=8) :: gt
      character(len=len(gtype)) :: g
      integer :: i, k
      g = gtype
      do i = 1, len(g)
         k = iachar(g(i:i))
         if (k >= iachar('A') .and. k <= iachar('Z')) g(i:i) = achar(k + 32)
      end do
      if (index(g, 'pah') == 1) then
         gt = 'pah'
      else if (index(g, 'cm20') == 1 .or. index(g, 'amc') == 1 .or. index(g, 'gra') == 1) then
         gt = 'gra'
      else
         gt = 'sil'
      end if
   end function solver_grain_type


   ! =====================================================================
   subroutine build_dust(m, model, data_dir, NT_in, T_lo, T_hi, include_euv, status, &
                         lam_min, kext_path, sd_index, u_isrf, &
                         config_path, astrodust_index_path, euv_tmatrix, &
                         load_polarized_optics, scatmat_path, message, zubko_optics)
      ! One entry point for every model this library codes, built from ONE
      ! directory and ONE flag:
      !
      !   call build_dust(m, 'astrodust', '../data', 100, 1.0d0, 3000.0d0, &
      !                   include_euv = .false., status = st)
      !
      ! WHAT IT RESOLVES.  data_dir is the SEDust data directory.  The optics
      ! come from data_dir/sedust_<model>.h5 -- the wavelength axis, the
      ! cross-section tables and the extinction curve alike -- so a host names
      ! one directory instead of a Q table and an extinction table
      ! separately, and ships only the models it runs.  What is NOT in that
      ! file, because it is not optics, still comes from beside it: the ZDA
      ! config and the calorimetry.  The HD23 and WD01 size distributions are
      ! analytic and need no file.
      !
      ! WHAT include_euv DECIDES.  The grid, and it decides it the same way for
      ! every model: .false. (the default) builds on the non-ionizing part of
      ! the model's own axis, .true. on the whole of it.  The cut is an INDEX
      ! cut at lyman_index -- the last node at or below 0.0912 um -- so the
      ! grid a host gets COVERS the Lyman limit whatever model it named.  It
      ! used to be an index cut for two models and a wavelength cut for the
      ! third, and the third then began 1.2e-5 of itself INSIDE the limit,
      ! which a host on its own Lyman-limit grid was refused for.
      !
      ! Because the extinction curve is read from the same file and covers the
      ! wider axis, it covers either grid -- which is what retires the old
      ! pairing of a narrow model with an _euv-only kext table, a mismatch that
      ! used to fail at dust_extinction rather than at the build.
      !
      ! WHAT lam_min DECIDES.  The shortest wavelength the model must COVER,
      ! and nothing else.  astrodust and DL07 meet it by extending their grid
      ! on the dielectric function their optics come from; zubko, the two
      ! DustEM-defined models and a file-defined one, which have no dielectric
      ! function behind their tables, meet it when those tables already reach
      ! and REFUSE it when they do not.  No model truncates for it.  It used to
      ! truncate for two of them, so one host with one physically motivated
      ! floor got a longer grid from two models and a shorter one from the
      ! others.
      !
      ! WHERE THE DATA IS.  For the length of the build, data_dir is the data
      ! root: the dielectric functions, the default extinction curves and the
      ! stored cross-section tables all resolve inside it, not against the
      ! working directory.  The previous root is restored on the way out.  A
      ! host can therefore put the data anywhere and stop changing directory
      ! around the call.
      !
      ! FALLBACK.  A tree with no HDF5 product, or a build made with HDF5=0,
      ! gets the text products through the same builders, and each says on
      ! stdout which source it used.  Nothing here requires the library to have
      ! been compiled against HDF5.
      type(dust_model_t), intent(out) :: m
      ! 'astrodust' | 'dl07' | 'mrn' | 'zubko' | 'themis' | 'g18d' | 'from_files'
      character(len=*),   intent(in)  :: model
      character(len=*),   intent(in)  :: data_dir
      ! The internal temperature grid: NT points, log-spaced over [T_lo, T_hi].
      ! It is what the EMISSION side is solved on -- the enthalpy H(T, a) and
      ! the Planck-averaged opacity kappB(T, a) are tabulated on it, calc_Teq
      ! interpolates T_eq in that table, and the stochastic solver's window is
      ! clamped to [T_lo, T_hi] rather than extrapolated past it.  So the range
      ! must BRACKET the grain temperatures the field produces, and NT sets how
      ! finely T_eq is resolved.  The extinction never reads any of it: a host
      ! that only wants dust_extinction can leave all three out.
      integer,  optional, intent(in)  :: NT_in
      real(wp), optional, intent(in)  :: T_lo, T_hi
      ! Carry the ionizing band.  Default .false.
      logical,  optional, intent(in)  :: include_euv
      ! 0 = success.  ONE vocabulary, whatever the model: the builders
      ! number their own stages differently (an unreadable extinction table is
      ! 10 for astrodust, 5 for DL07, 3 for MRN, 6 for zubko, 8 for a
      ! DustEM-defined model and 9 for from_files), and a
      ! host using the single entry point should not have to branch on the
      ! model name to read a code.  Each builder's code is mapped onto this
      ! list:
      !   status = 0   built
      !   status = 1   optics table (Q table, or a component's optics)
      !   status = 2   size distribution (retired: none is read from a file)
      !   status = 3   dielectric function
      !   status = 4   lam_min not coverable by this model
      !   status = 5   extinction table
      !   status = 6   calorimetry
      !   status = 7   grid inconsistent between components
      !   status = 8   EUV spheroid optics unavailable (no T-matrix registered)
      !   status = 9   model definition (the ZDA config, or the descriptor)
      !   status = 10  polarized optics (the aligned or oriented tables)
      !   status = 90  model name not one of the seven
      !   status = 92  zubko_optics not one of 'zda' | 'mie_d03'
      !   status = 91  'from_files' without config_path (the descriptor)
      integer,  optional, intent(out) :: status
      ! What went wrong, in words, for a host that has to print one line before
      ! a collective abort.  Blank on success.
      character(len=*), optional, intent(out) :: message
      real(wp), optional, intent(in)  :: lam_min
      ! Extinction curve to serve.  Omitted, /kext of the model's own product.
      character(len=*), optional, intent(in) :: kext_path
      ! DL07 only: which of the WD01 size distributions, and the field strength
      ! its PAH ionization balance is computed at.  Defaults 7 and 1.0, the
      ! reference MW R_V = 3.1, b_C = 6e-5 model at U = 1.
      integer,  optional, intent(in)  :: sd_index
      real(wp), optional, intent(in)  :: u_isrf
      ! The file that DEFINES the model, where the model is a set of files:
      ! the ZDA config for zubko, the GRAIN_*.DAT for themis and g18d, the
      ! descriptor for from_files -- and there it is required.  Omitted, each
      ! of those models takes the file its own directory ships.
      character(len=*), optional, intent(in) :: config_path
      character(len=*), optional, intent(in) :: astrodust_index_path
      logical,  optional, intent(in)  :: euv_tmatrix
      ! astrodust only.  The orientation-resolved table and its axes are
      ! resolved from data_dir like everything else, so a host names neither;
      ! these two are what it may still want to say.  load_polarized_optics =
      ! .false. builds the model scalar-only.  scatmat_path names the aligned
      ! scattering matrix and is what asks for it to be loaded at all -- it is
      ! a large table and nothing loads it unasked.
      logical,  optional, intent(in)  :: load_polarized_optics
      character(len=*), optional, intent(in) :: scatmat_path
      ! zubko only: which of the two stored optics sets to build on, 'zda'
      ! (default, the benchmark's own tables) or 'mie_d03' (this tree's
      ! recomputation).  See build_zubko.
      character(len=*), optional, intent(in) :: zubko_optics

      character(len=512) :: h5, cfg, ddir, adir
      character(len=SED_PATHLEN) :: saved_root
      integer  :: nt, st
      real(wp) :: tlo, thi
      real(wp), allocatable :: lam_h5(:)
      integer  :: sdi, i_lyman
      real(wp) :: uisrf
      logical  :: wide, got

      if (present(status)) status = 0
      if (present(message)) message = ''
      st = 0
      wide = euv_asked(include_euv)
      ! Defaults: the grid the manual quotes as typical, wide enough for the
      ! CMB floor at one end and a stochastically heated small grain at the
      ! other.
      nt   = 200;        if (present(NT_in)) nt   = NT_in
      tlo  = 2.7_wp;     if (present(T_lo))  tlo  = T_lo
      thi  = 5.0e3_wp;   if (present(T_hi))  thi  = T_hi
      ddir = data_dir
      if (len_trim(ddir) == 0) ddir = '../data'
      ! Point the whole library at this directory for the length of the build,
      ! so the dielectric functions and the default extinction curves resolve
      ! inside it and not against the working directory.  Restored on the way
      ! out, at every exit.
      saved_root = sed_get_data_root()
      call sed_set_data_root(trim(ddir))

      h5  = sedust_h5_file(trim(ddir), model)
      ! kext_path is NOT given a value here.  Naming a file is a contract --
      ! the builder fails if it cannot be read -- and passing the product path
      ! unconditionally turned the model's own SOFT default into a hard one,
      ! so a tree without the curve could no longer build a model for emission
      ! alone.  Absent, it propagates as absent, and each builder falls to its
      ! own default: /kext of the product, then the model's text curve, then no
      ! extinction and a status from dust_extinction rather than a failed
      ! build.  Both defaults resolve against the data root set above.

      sdi   = 7;        if (present(sd_index)) sdi   = sd_index
      uisrf = 1.0_wp;   if (present(u_isrf))   uisrf = u_isrf

      select case (trim(model))

      case ('astrodust')
         ! The Q table IS the optics here, so the HDF5 path goes straight in;
         ! sed_init reads /qtable/astrodust from it and cuts at i_lyman.  The
         ! orientation-resolved table and its two axes come from the same
         ! directory, so a host that moved data_dir moves them with it.
         adir = sedust_dir(trim(ddir), 'astrodust')
         call build_astrodust(m, trim(h5), nt, tlo, thi, status=st, &
                              lam_min=lam_min, astrodust_index_path=astrodust_index_path, &
                              kext_path=kext_path, euv_tmatrix=euv_tmatrix, &
                              include_euv=wide, &
                              load_polarized_optics=load_polarized_optics, &
                              qpol_path=trim(adir)//'q_DH21Ad_P0.20_Fe0.00_1.400.dat.gz', &
                              qpol_wave_path=trim(adir)//'DH21_wave', &
                              qpol_aeff_path=trim(adir)//'DH21_aeff', &
                              scatmat_path=scatmat_path)
         call report(st, [1, 2, 10, 10, 10, 3, 4, 10, 10, 5, 8])

      case ('dl07')
         ! This model takes only a grid from a Q table -- its optics are Mie on
         ! the D03 dielectric functions -- so hand it its own axis and it needs
         ! no astrodust product to borrow one from.
         call read_sedust_grid(trim(h5), wide, lam_h5, i_lyman, got)
         if (got) then
            call build_dl07(m, trim(h5), sdi, uisrf, nt, tlo, thi, &
                            status=st, lam_min=lam_min, kext_path=kext_path, &
                            lam_axis=lam_h5)
            deallocate(lam_h5)
         else
            ! No product to read: the text route, on the astrodust Q table the
            ! model's grid has always come from.
            call build_dl07(m, trim(ddir)//dl07_text_qtable(wide), sdi, uisrf, &
                            nt, tlo, thi, status=st, lam_min=lam_min, &
                            kext_path=kext_path)
         end if
         call report(st, [1, 2, 0, 0, 5, 0, 4])

      case ('mrn')
         ! Like DL07, this model takes only a wavelength axis from a product --
         ! its optics are Mie on the D03 dielectric functions -- so hand it its
         ! own axis and nothing else has to be read for a grid.
         call read_sedust_grid(trim(h5), wide, lam_h5, i_lyman, got)
         if (got) then
            call build_mrn(m, trim(h5), nt, tlo, thi, status=st, lam_min=lam_min, &
                           kext_path=kext_path, lam_axis=lam_h5)
            deallocate(lam_h5)
         else
            ! No product to read: the text route, on the astrodust Q table the
            ! DL07 grid has always come from, so the two models that share a
            ! calculation share an axis.
            call build_mrn(m, trim(ddir)//dl07_text_qtable(wide), nt, tlo, thi, &
                           status=st, lam_min=lam_min, kext_path=kext_path)
         end if
         call report(st, [1, 4, 5])

      case ('zubko')
         cfg = trim(sedust_dir(trim(ddir), 'zubko'))//'ZDA_BARE_GR_S_Config.dat'
         if (present(config_path)) cfg = config_path
         ! nt/tlo/thi, not the optional dummies: passing NT_in/T_lo/T_hi on
         ! through meant that a caller who omitted them handed absent
         ! optionals to non-optional arguments.
         call build_zubko(m, trim(cfg), trim(sedust_dir(trim(ddir), 'zubko')), &
                          nt, tlo, thi, &
                          status=st, kext_path=kext_path, lam_min=lam_min, &
                          include_euv=wide, qtable_path=trim(h5), &
                          optics=zubko_optics)
         call report(st, [9, 9, 1, 7, 6, 5, 4, 92])

      case ('themis', 'g18d')
         ! Two published models carried as DustEM input files -- THEMIS (Jones
         ! et al. 2013 for the carbon grains, Koehler et al. 2014 for the
         ! silicates) and Guillet et al. (2018) Model D.  Both are SCALAR: the
         ! distribution ships orientation-averaged Q tables, so these models
         ! carry no polarized optics and none can be formed here; a host that
         ! asks dust_has_polarized_optics about one gets .false.
         !
         ! The model definition is one file inside the model's own directory,
         ! and everything it implies -- oprop/, hcap/, the extinction curve,
         ! the product -- hangs under that same directory, so naming the model
         ! names all of it.  G18D ships no G_ file for its two large
         ! populations, so that model carries no <cos theta>; the model object
         ! says so through m%gsca_complete and the size integral returns zero
         ! rather than the average of the populations that do carry one.
         if (trim(model) == 'themis') then
            cfg = trim(sedust_dir(trim(ddir), 'themis'))//'GRAIN_J13.DAT'
         else
            cfg = trim(sedust_dir(trim(ddir), 'g18d'))//'GRAIN_G17_ModelD.DAT'
         end if
         if (present(config_path)) cfg = config_path
         call build_dustem(m, trim(cfg), trim(sedust_dir(trim(ddir), model)), &
                           nt, tlo, thi, status=st, kext_path=kext_path, &
                           lam_min=lam_min, include_euv=wide, qtable_path=trim(h5))
         call report(st, [9, 1, 1, 2, 1, 6, 1, 5, 4])

      case ('from_files')
         if (.not. present(config_path)) then
            call finish(91, 'build_dust: from_files needs config_path (the descriptor)')
            return
         end if
         call build_from_files(m, config_path, trim(ddir), nt, tlo, thi, &
                               status=st, kext_path=kext_path, lam_min=lam_min, &
                               include_euv=wide)
         call report(st, [9, 9, 9, 9, 1, 7, 2, 6, 5, 4])

      case default
         call finish(90, 'build_dust: unknown model '''//trim(model)// &
                         ''' (astrodust | dl07 | mrn | zubko | themis | g18d | from_files)')
         return
      end select

      call sed_set_data_root(trim(saved_root))

   contains

      subroutine report(code, map)
         ! Map a builder's own stage number onto build_dust's single
         ! vocabulary.  map(k) is what that builder's code k means here; 0
         ! marks a code the builder does not use.
         integer, intent(in) :: code, map(:)
         if (code <= 0) return
         if (code <= size(map)) then
            if (map(code) > 0) then
               call finish(map(code), '');  return
            end if
         end if
         call finish(code, '')       ! 90/91 and anything unmapped, as given
      end subroutine report

      subroutine finish(code, msg)
         ! One exit for every failure: restore the data root, then report
         ! through status, or stop if the caller asked for no status.
         integer,          intent(in) :: code
         character(len=*), intent(in) :: msg
         call sed_set_data_root(trim(saved_root))
         if (present(message)) then
            if (len_trim(msg) > 0) then
               message = msg
            else
               message = 'build_dust: '//trim(model)//' failed at '//trim(stage_name(code))
            end if
         end if
         if (present(status)) then
            status = code
         else
            if (len_trim(msg) > 0) then
               write(*,'(a)') trim(msg)
            else
               write(*,'(a,a,a,a)') ' build_dust: ', trim(model), ' failed at ', &
                                    trim(stage_name(code))
            end if
            stop 1
         end if
      end subroutine finish

      function stage_name(code) result(nm)
         integer, intent(in) :: code
         character(len=48)   :: nm
         select case (code)
         case (1);   nm = 'the optics table'
         case (2);   nm = 'the size distribution'
         case (3);   nm = 'the dielectric function'
         case (4);   nm = 'lam_min, which this model cannot cover'
         case (5);   nm = 'the extinction table'
         case (6);   nm = 'the calorimetry'
         case (7);   nm = 'the grid, which its components disagree on'
         case (8);   nm = 'the EUV spheroid optics (no T-matrix registered)'
         case (9);   nm = 'the model definition'
         case (10);  nm = 'the polarized optics'
         case (90);  nm = 'the model name'
         case (91);  nm = 'the missing descriptor'
         case (92);  nm = 'the zubko_optics name'
         case default;  write(nm,'(a,i0)') 'stage ', code
         end select
      end function stage_name

      function dl07_text_qtable(w) result(p)
         ! The astrodust T-matrix table the DL07 grid comes from in the text
         ! route.  It is astrodust's, so it lives in that model's directory.
         logical, intent(in) :: w
         character(len=128)  :: p
         if (w) then
            p = '/astrodust/q_astrodust_P0.20_Fe0.00_1.400_euv.dat'
         else
            p = '/astrodust/q_astrodust_P0.20_Fe0.00_1.400.dat'
         end if
      end function dl07_text_qtable


   end subroutine build_dust


   ! Generic single-cell dust emission for the (active) model m. Loops the
   ! populations through the untouched sed_grain_loop, sums per output
   ! channel, applies the induced factor (if enabled) and the HD23 unit
   ! convention. lamI_total(NLAM) is the summed SED; optional
   ! lamI_chan(NLAM, n_channel) returns the SED of each channel.
   ! REQUIRES: m is the most recently built model.  This routine solves through
   ! sed_grain_loop, which reads the module grids (lam, NLAM, NT, T_first) that
   ! the last build filled, so it can answer for that model and no other; the
   ! guard below turns a stale model into a status instead of another model's
   ! numbers.
   subroutine dust_emission(m, J_lam, lamI_total, lamI_chan, status, lamI_pol)
      type(dust_model_t), intent(in)  :: m
      real(wp),           intent(in)  :: J_lam(:)              ! (NLAM)
      real(wp),           intent(out) :: lamI_total(:)         ! (NLAM)
      real(wp), optional, intent(out) :: lamI_chan(:,:)        ! (NLAM, n_channel)
      ! Optional error report (0 = success). When present, a bad model is
      ! reported through it instead of stopping the process; when absent the
      ! original stop-on-error behavior is kept (as the CLI drivers expect).
      !   status = 1  unknown stoch_method
      !   status = 2  'qm' selected but a population is missing its radii
      !   status = 3  m is not the most recently built model
      integer,  optional, intent(out) :: status
      ! Optional INTRINSIC polarized emission, same shape and units as
      ! lamI_total. Only populations carrying both Cpol and falign contribute;
      ! anything else (PAHs, any model built without polarized optics) adds
      ! zero. The geometric sin^2(gamma) projection onto the plane of the sky
      ! and any turbulent depolarization are the radiative transfer's job and
      ! are deliberately NOT applied here.
      real(wp), optional, intent(out) :: lamI_pol(:)           ! (NLAM)
      real(wp), allocatable :: Jout_pop(:), Jchan(:,:), Jpol_pop(:), Jpol_sum(:)
      integer :: ip, ic
      logical :: want_pol, pop_pol

      if (present(status)) status = 0
      want_pol = present(lamI_pol)

      ! This routine solves on the MODULE grids, which the last build filled,
      ! so it can only answer for the model that build produced.  Refuse any
      ! other rather than return another model's numbers under this one's name.
      if (m%build_id /= active_build_id .or. active_build_id == 0) then
         if (present(status)) then
            status = 3;  return
         else
            write(*,'(a)') 'dust_emission: m is not the most recently built model'
            write(*,'(a)') '   rebuild it, or query models one at a time'
            stop 1
         end if
      end if

      ! Validate the model's chosen solver before doing any work.
      select case (trim(m%stoch_method))
      case ('heuristic', 'draine', 'qm', 'equil')
         ! supported
      case default
         if (present(status)) then
            status = 1;  return
         else
            write(*,'(a,a)') 'dust_emission: unknown stoch_method: ', trim(m%stoch_method)
            stop 1
         end if
      end select

      ! The 'qm' solver reads each population's radii; refuse rather than read
      ! an unallocated array (all builders now fill them, so this is a guard).
      if (trim(m%stoch_method) == 'qm') then
         do ip = 1, size(m%pops)
            if (.not. allocated(m%pops(ip)%aeff)) then
               if (present(status)) then
                  status = 2;  return
               else
                  write(*,'(a,i0)') 'dust_emission: qm needs radii but pop is unset, ip=', ip
                  stop 1
               end if
            end if
         end do
      end if

      ! Honor the model's chosen solver ('heuristic'/'draine'/'qm'/'equil')
      ! and its diagnostic verbosity (library path stays silent by default).
      stoch_method = m%stoch_method
      sed_verbose  = m%verbose
      if (trim(m%stoch_method) == 'qm') qm_verbose = m%verbose

      allocate(Jout_pop(m%NLAM), Jchan(m%NLAM, m%n_channel))
      allocate(Jpol_pop(m%NLAM), Jpol_sum(m%NLAM))
      Jchan    = 0.0_wp
      Jpol_sum = 0.0_wp
      do ip = 1, size(m%pops)
         ! Ask for polarized emission only from a population that actually
         ! carries polarized optics; the rest go down the untouched path.
         pop_pol = want_pol .and. allocated(m%pops(ip)%Cpol) &
                             .and. allocated(m%pops(ip)%falign)
         ! size count for each population (Zubko-like models have
         ! component-by-component grids)
         if (pop_pol) then
            call sed_grain_loop(size(m%pops(ip)%dn), m%pops(ip)%dn, m%pops(ip)%aeff, &
                                m%pops(ip)%Cabs, &
                                m%pops(ip)%kappB, m%pops(ip)%H, m%pops(ip)%log_H, &
                                m%pops(ip)%log_kappB, m%pops(ip)%kappCMB, &
                                J_lam, trim(m%pops(ip)%grain_type), Jout_pop, &
                                m%pops(ip)%Cpol, m%pops(ip)%falign, Jpol_pop)
            Jpol_sum = Jpol_sum + Jpol_pop
         else
            call sed_grain_loop(size(m%pops(ip)%dn), m%pops(ip)%dn, m%pops(ip)%aeff, &
                                m%pops(ip)%Cabs, &
                                m%pops(ip)%kappB, m%pops(ip)%H, m%pops(ip)%log_H, &
                                m%pops(ip)%log_kappB, m%pops(ip)%kappCMB, &
                                J_lam, trim(m%pops(ip)%grain_type), Jout_pop)
         end if
         ic = m%pops(ip)%out_channel
         Jchan(:, ic) = Jchan(:, ic) + Jout_pop
      end do

      do ic = 1, m%n_channel
         if (m%use_induced_emission) call apply_induced_factor(J_lam, Jchan(:, ic))
         Jchan(:, ic) = m%lam * Jchan(:, ic) * 1.0e-3_wp
      end do

      lamI_total = sum(Jchan, dim=2)
      if (present(lamI_chan)) lamI_chan = Jchan
      if (want_pol) then
         ! Same unit conversion the total emission gets.
         if (m%use_induced_emission) call apply_induced_factor(J_lam, Jpol_sum)
         lamI_pol = m%lam * Jpol_sum * 1.0e-3_wp
      end if
      deallocate(Jout_pop, Jchan, Jpol_pop, Jpol_sum)
   end subroutine dust_emission


   ! Size-distribution-integrated extinction of the (active) model m, per H
   ! atom, computed from first principles out of the model's own optics.
   !
   ! This is what the standalone calculators use, and what calc_kext.x writes
   ! into the data/kext_*.dat tables; dust_extinction reads those tables back
   ! and does NOT come through here for its scalar optics.  The two therefore
   ! agree by construction as long as the table was made from the same model.
   !
   ! The size integral is the plain binned sum over each population, because
   ! dn(a) already carries the bin width:
   !   C_abs/H  = sum_pop sum_a dn(a) * Cabs(lambda, a)
   !   C_sca/H  = sum_pop sum_a dn(a) * Csca(lambda, a)
   !   C_ext/H  = C_abs/H + C_sca/H
   !   <cos>    = sum dn * Csca * g  /  sum dn * Csca      (scattering-weighted)
   !   C_polext = sum_pop sum_a dn(a) * Cpol_ext(lambda, a) * f_align(a)
   !   C_birext = sum_pop sum_a dn(a) * Cbir_ext(lambda, a) * f_align(a)
   ! A population whose Csca / gsca / Cpol_ext / Cbir_ext are unallocated
   ! contributes zero to those terms; the PAHs are that case for the polarized
   ! pair, being unaligned, but not for the scalar one -- what scatters in a
   ! PAH population is the graphite fraction xi_gra(a) of HD23 eq. 15, since a
   ! PAH has an absorption cross section and no scattering one. C_birext is
   ! also zero when the model was built from a 3-block table with no
   ! birefringence.
   !
   ! Units: all cross sections [cm^2/H]; gbar dimensionless. C_polext and
   ! C_birext are the MAXIMUM dichroic and birefringent extinction -- the size
   ! integral and the f_align weight are done here, but the sin^2(gamma) geometry
   ! factor and any turbulent depolarization are the radiative transfer's job and
   ! are NOT applied. C_polext is the IQ-block optic and C_birext the UV-block
   ! optic of the extinction matrix (see extinction_matrix_aligned).
   !
   ! This routine reads NOTHING but m: every array it sums is a component of
   ! the model argument, so unlike dust_emission it does NOT require m to be
   ! the most recently built model, and two models may be integrated one after
   ! the other.  It carried the opposite claim in this comment for a while, and
   ! at least one host wrote that restriction down as though it were real.
   subroutine size_integrated_extinction(m, Cext, Cabs, Csca, gbar, Cpol_ext, Cbir_ext, &
                                         albedo, status)
      type(dust_model_t), intent(in)  :: m
      real(wp),           intent(out) :: Cext(:), Cabs(:), Csca(:)   ! (NLAM) [cm^2/H]
      ! Scattering-weighted asymmetry; 0 where nothing scatters.
      real(wp), optional, intent(out) :: gbar(:)                     ! (NLAM)
      real(wp), optional, intent(out) :: Cpol_ext(:)                 ! (NLAM) [cm^2/H]
      ! Birefringent extinction; 0 for a 3-block model or where nothing aligns.
      real(wp), optional, intent(out) :: Cbir_ext(:)                 ! (NLAM) [cm^2/H]
      ! Scattering albedo C_sca/C_ext; 0 where the medium is transparent.
      ! Derived here so that every caller gets the same convention at the
      ! wavelengths where C_ext underflows to zero.
      real(wp), optional, intent(out) :: albedo(:)                   ! (NLAM)
      ! Optional error report (0 = success). When present, a size mismatch is
      ! reported through it instead of stopping the process; when absent such a
      ! call stops the run, matching dust_emission.
      !   status = 1  an output array is not of size m%NLAM
      integer,  optional, intent(out) :: status
      real(wp), allocatable :: gnum(:)
      integer :: ip, ja, jw, na_p
      logical :: bad, g_complete

      if (present(status)) status = 0

      bad = size(Cext) /= m%NLAM .or. size(Cabs) /= m%NLAM .or. size(Csca) /= m%NLAM
      if (present(gbar))     bad = bad .or. size(gbar)     /= m%NLAM
      if (present(Cpol_ext)) bad = bad .or. size(Cpol_ext) /= m%NLAM
      if (present(Cbir_ext)) bad = bad .or. size(Cbir_ext) /= m%NLAM
      if (present(albedo))   bad = bad .or. size(albedo)   /= m%NLAM
      if (bad) then
         if (present(status)) then
            status = 1;  return
         else
            write(*,'(a,i0)') 'size_integrated_extinction: output arrays must be of '// &
                 'size m%NLAM=', m%NLAM
            stop 1
         end if
      end if

      allocate(gnum(m%NLAM))
      Cabs = 0.0_wp;  Csca = 0.0_wp;  gnum = 0.0_wp

      ! Is <cos theta> defined for this model at all?  gnum sums over the
      ! populations that carry an asymmetry parameter while Csca sums over
      ! every population that scatters, so if one population scatters without
      ! a g the ratio below is a partial numerator over a complete
      ! denominator -- a number that is the asymmetry of nothing and that
      ! reads as a small g rather than as a missing one.  The model then gets
      ! <cos theta> = 0 everywhere, and the model object says why through
      ! gsca_complete.  This is the state of Guillet et al. (2018) Model D,
      ! whose DustEM distribution ships no G_ file for its two large
      ! populations; every model built from a Mie or T-matrix calculation
      ! carries g for each scattering population and is unaffected.
      g_complete = .true.
      do ip = 1, size(m%pops)
         if (allocated(m%pops(ip)%Csca) .and. .not. allocated(m%pops(ip)%gsca)) &
            g_complete = .false.
      end do

      do ip = 1, size(m%pops)
         na_p = size(m%pops(ip)%dn)
         do ja = 1, na_p
            do jw = 1, m%NLAM
               Cabs(jw) = Cabs(jw) + m%pops(ip)%dn(ja) * m%pops(ip)%Cabs(jw, ja)
            end do
         end do
         if (allocated(m%pops(ip)%Csca)) then
            do ja = 1, na_p
               do jw = 1, m%NLAM
                  Csca(jw) = Csca(jw) + m%pops(ip)%dn(ja) * m%pops(ip)%Csca(jw, ja)
               end do
            end do
            if (allocated(m%pops(ip)%gsca)) then
               do ja = 1, na_p
                  do jw = 1, m%NLAM
                     gnum(jw) = gnum(jw) + m%pops(ip)%dn(ja) &
                                * m%pops(ip)%Csca(jw, ja) * m%pops(ip)%gsca(jw, ja)
                  end do
               end do
            end if
         end if
      end do

      call size_integrated_polarized_extinction(m, Cpol_ext, Cbir_ext)

      Cext = Cabs + Csca
      if (present(gbar)) then
         gbar = 0.0_wp
         if (g_complete) then
            do jw = 1, m%NLAM
               if (Csca(jw) > 0.0_wp) gbar(jw) = gnum(jw) / Csca(jw)
            end do
         end if
      end if
      if (present(albedo)) then
         albedo = 0.0_wp
         do jw = 1, m%NLAM
            if (Cext(jw) > 0.0_wp) albedo(jw) = Csca(jw) / Cext(jw)
         end do
      end if
      deallocate(gnum)
   end subroutine size_integrated_extinction


   ! Dust mass per H nucleon [g/H] of the model m -- the model's own size
   ! distribution weighted by the solid density of each population's material:
   !
   !   M_dust/N_H = sum_pop rho_bulk(pop) * sum_a (4/3) pi a_cm^3 dn_pop(a)
   !
   ! dn(a) already carries the width of its size bin, so no da enters the sum,
   ! and a_eff is the population's own effective radius in um.  A population
   ! whose model states no density (rho_bulk = 0) contributes nothing.
   !
   ! This is the constant that converts between the two ways dust opacity is
   ! quoted: the cross section per H this library serves and the mass opacity of
   ! Draine's tables,
   !
   !   K_abs [cm^2/g] = (C_abs/H) / (M_dust/N_H) ,
   !
   ! which is the trailing column calc_kext.x writes into every data/kext_*.dat.
   ! An RT host that carries a dust mass density rather than an H column density
   ! needs the same number to turn one into the other.
   !
   ! The value is a property of the MODEL alone -- size distribution and grain
   ! densities -- so it is wavelength- and radiation-field-independent.
   ! REQUIRES: nothing beyond a built model; it reads m only.
   real(wp) function dust_mass_per_H(m) result(Mdust_H)
      type(dust_model_t), intent(in) :: m
      real(wp) :: acm3
      integer  :: ip, ja, na_p

      Mdust_H = 0.0_wp
      if (.not. allocated(m%pops)) return

      do ip = 1, size(m%pops)
         if (m%pops(ip)%rho_bulk <= 0.0_wp) cycle
         if (.not. allocated(m%pops(ip)%dn) .or. .not. allocated(m%pops(ip)%aeff)) cycle
         na_p = min(size(m%pops(ip)%dn), size(m%pops(ip)%aeff))
         do ja = 1, na_p
            acm3 = (m%pops(ip)%aeff(ja) * UM2CM)**3
            Mdust_H = Mdust_H + (4.0_wp/3.0_wp)*PI*acm3 * &
                      m%pops(ip)%dn(ja) * m%pops(ip)%rho_bulk
         end do
      end do
   end function dust_mass_per_H


   ! Dichroic and birefringent extinction per H of the (active) model m, the
   ! f_align-weighted size integrals
   !   C_polext = sum_pop sum_a dn(a) * Cpol_ext(lambda, a) * f_align(a)
   !   C_birext = sum_pop sum_a dn(a) * Cbir_ext(lambda, a) * f_align(a)
   ! computed from the model's own orientation-resolved optics.  Both outputs
   ! are optional; an absent one is not computed.  A population without
   ! polarized optics or without an alignment efficiency -- the PAHs, which
   ! HD23 take to be unaligned -- contributes zero, as does a model built from
   ! a 3-block table with no birefringence.
   !
   ! Shared by size_integrated_extinction and dust_extinction, which is why
   ! the two cannot disagree on the polarized optics even though they get
   ! their scalar optics from different places.
   subroutine size_integrated_polarized_extinction(m, Cpol_ext, Cbir_ext)
      type(dust_model_t), intent(in)  :: m
      real(wp), optional, intent(out) :: Cpol_ext(:)   ! (NLAM) [cm^2/H]
      real(wp), optional, intent(out) :: Cbir_ext(:)   ! (NLAM) [cm^2/H]
      integer :: ip, ja, jw, na_p

      if (.not. present(Cpol_ext) .and. .not. present(Cbir_ext)) return
      if (present(Cpol_ext)) Cpol_ext = 0.0_wp
      if (present(Cbir_ext)) Cbir_ext = 0.0_wp
      if (.not. allocated(m%pops)) return

      do ip = 1, size(m%pops)
         na_p = size(m%pops(ip)%dn)
         if (present(Cpol_ext)) then
            if (allocated(m%pops(ip)%Cpol_ext) .and. allocated(m%pops(ip)%falign)) then
               do ja = 1, na_p
                  do jw = 1, m%NLAM
                     Cpol_ext(jw) = Cpol_ext(jw) + m%pops(ip)%dn(ja) &
                                    * m%pops(ip)%Cpol_ext(jw, ja) * m%pops(ip)%falign(ja)
                  end do
               end do
            end if
         end if
         if (present(Cbir_ext)) then
            if (allocated(m%pops(ip)%Cbir_ext) .and. allocated(m%pops(ip)%falign)) then
               do ja = 1, na_p
                  do jw = 1, m%NLAM
                     Cbir_ext(jw) = Cbir_ext(jw) + m%pops(ip)%dn(ja) &
                                    * m%pops(ip)%Cbir_ext(jw, ja) * m%pops(ip)%falign(ja)
                  end do
               end do
            end if
         end if
      end do
   end subroutine size_integrated_polarized_extinction


   ! Transport optics of the model m per H atom, on the model's own wavelength
   ! grid -- the extinction counterpart of dust_emission, and the routine an RT
   ! host calls.
   !
   ! WHERE THE NUMBERS COME FROM.  The four SCALAR outputs -- Cext, Cabs, Csca
   ! and gbar -- are read from the size-integrated extinction table the builder
   ! attached to the model (m%kext_*, one of the data/kext_*.dat products of
   ! calc_kext.x), interpolated onto m%lam.  The transport optics of a dust
   ! model do not depend on the transport, so there is nothing for a host to
   ! gain by re-running the size integral in the middle of a run; the table is
   ! the same integral, done once, by calc_kext.x, and recorded.  The
   ! first-principles route is still available under its own name,
   ! size_integrated_extinction, and is what the standalone calculators use.
   !
   ! The two POLARIZED outputs -- Cpol_ext and Cbir_ext -- are NOT tabulated.
   ! They are computed here, exactly as before, from the model's own
   ! orientation-resolved optics weighted by the alignment efficiency
   ! f_align(a), because f_align is a runtime state an RT host resets cell by
   ! cell through dust_set_alignment: a table fixed at build time could not
   ! follow it.  (The astrodust product does carry a dichroic column, and the
   ! EUV products deliberately do not; neither is read here.)  This asymmetry
   ! -- scalar optics from a table, polarized optics from the model -- is
   ! deliberate.
   !
   ! INTERPOLATION.  Cext, Cabs and Csca are interpolated linearly in
   ! log(lambda)-log(C); gbar linearly in log(lambda), because it is signed.
   ! A model wavelength that coincides with a table node takes that node's
   ! value unchanged, so a model whose grid IS the grid the table was written
   ! on -- the usual case, since both come from the same Q table -- gets the
   ! tabulated numbers back exactly.  Nothing is extrapolated: a model
   ! wavelength outside the table is an error (status 3), not a frozen
   ! boundary value.
   !
   ! Units: all cross sections [cm^2/H]; gbar and albedo dimensionless.
   ! C_polext and C_birext are the MAXIMUM dichroic and birefringent
   ! extinction -- the size integral and the f_align weight are done here, but
   ! the sin^2(gamma) geometry factor and any turbulent depolarization are the
   ! radiative transfer's job and are NOT applied.  C_polext is the IQ-block
   ! optic and C_birext the UV-block optic of the extinction matrix (see
   ! extinction_matrix_aligned).
   !
   ! Reads only m, like size_integrated_extinction: it serves m%kext_* and no
   ! module grid, so it is not restricted to the most recently built model.
   subroutine dust_extinction(m, Cext, Cabs, Csca, gbar, Cpol_ext, Cbir_ext, albedo, status)
      type(dust_model_t), intent(in)  :: m
      real(wp),           intent(out) :: Cext(:), Cabs(:), Csca(:)   ! (NLAM) [cm^2/H]
      ! Scattering asymmetry <cos>, as tabulated.
      real(wp), optional, intent(out) :: gbar(:)                     ! (NLAM)
      real(wp), optional, intent(out) :: Cpol_ext(:)                 ! (NLAM) [cm^2/H]
      ! Birefringent extinction; 0 for a 3-block model or where nothing aligns.
      real(wp), optional, intent(out) :: Cbir_ext(:)                 ! (NLAM) [cm^2/H]
      ! Scattering albedo C_sca/C_ext; 0 where the medium is transparent.
      ! Derived here so that every caller gets the same convention at the
      ! wavelengths where C_ext underflows to zero.
      real(wp), optional, intent(out) :: albedo(:)                   ! (NLAM)
      ! Optional error report (0 = success). When present, a bad call is
      ! reported through it instead of stopping the process; when absent such a
      ! call stops the run, matching dust_emission.
      !   status = 1  an output array is not of size m%NLAM
      !   status = 2  no extinction table was loaded for this model, so there
      !               are no scalar optics to return (see the builders'
      !               kext_path argument)
      !   status = 3  the model's wavelength grid runs outside the table
      integer,  optional, intent(out) :: status
      integer :: jw
      logical :: bad, ok

      if (present(status)) status = 0

      bad = size(Cext) /= m%NLAM .or. size(Cabs) /= m%NLAM .or. size(Csca) /= m%NLAM
      if (present(gbar))     bad = bad .or. size(gbar)     /= m%NLAM
      if (present(Cpol_ext)) bad = bad .or. size(Cpol_ext) /= m%NLAM
      if (present(Cbir_ext)) bad = bad .or. size(Cbir_ext) /= m%NLAM
      if (present(albedo))   bad = bad .or. size(albedo)   /= m%NLAM
      if (bad) then
         if (present(status)) then
            status = 1;  return
         else
            write(*,'(a,i0)') 'dust_extinction: output arrays must be of size m%NLAM=', m%NLAM
            stop 1
         end if
      end if

      if (m%kext_n <= 0) then
         if (present(status)) then
            status = 2;  return
         else
            write(*,'(a)') 'dust_extinction: no extinction table loaded for model '// &
                 trim(m%name)//'; name one with the builder''s kext_path argument'
            stop 1
         end if
      end if

      ! gbar is forwarded as an optional, so this routine allocates nothing and
      ! stays callable from an RT host's photon threads.
      call tabulated_extinction_on_grid(m%kext_n, m%kext_lam, m%kext_Cext, m%kext_Cabs, &
                                        m%kext_Csca, m%kext_gbar, m%lam, &
                                        Cext, Cabs, Csca, gbar, ok)
      if (.not. ok) then
         if (present(status)) then
            status = 3;  return
         else
            write(*,'(a)') 'dust_extinction: the model grid runs outside the extinction '// &
                 'table '//trim(m%kext_path)
            stop 1
         end if
      end if

      if (present(albedo)) then
         albedo = 0.0_wp
         do jw = 1, m%NLAM
            if (Cext(jw) > 0.0_wp) albedo(jw) = Csca(jw) / Cext(jw)
         end do
      end if

      ! Polarized optics: still computed from the model, never tabulated.
      call size_integrated_polarized_extinction(m, Cpol_ext, Cbir_ext)
   end subroutine dust_extinction


   ! Attach the size-integrated extinction table this model will serve to an RT
   ! host through dust_extinction.  Called at the END of every builder, once the
   ! model stands, because the table paths are relative to sed/ and a host that
   ! changes directory around the build call would find them broken later.
   !
   ! kext_path names the file explicitly; omitting it falls back on
   ! default_path, and an empty default_path means this model has none to try.
   !
   ! ok = .false. ONLY when an EXPLICITLY named table could not be read -- a
   ! host naming a file that is not there is a configuration error the builder
   ! must report.  A missing or unreadable DEFAULT leaves m%kext_n = 0 and
   ! ok = .true.: the standalone drivers (calc_sed.x astrodust and the rest) build
   ! models to compute emission and never call dust_extinction, and calc_kext.x
   ! builds a model precisely in order to WRITE the table that is not there yet.
   ! Neither may be stopped by its absence.
   pure function dl07_kext_tag(vintage) result(t)
      ! The suffix the DL07-model products carry when the carbonaceous
      ! absorption is NOT this model's own DL07 vintage.  The model's own is
      ! unmarked, so the shipped file names and the /kext group stay put.
      character(len=*), intent(in) :: vintage
      character(len=16) :: t
      t = ''
      if (trim(vintage) /= 'dl07') t = '_'//trim(vintage)
   end function dl07_kext_tag


   pure function dl07_kext_default(vintage) result(p)
      ! The size-integrated extinction curve of that vintage, as calc_kext.x
      ! writes it.  The widest grid, for the reason load_model_extinction_table
      ! gives: the curve is interpolated onto whatever grid the model is built
      ! on, so the widest one covers every such grid.
      character(len=*), intent(in) :: vintage
      character(len=96) :: p
      if (trim(vintage) == 'dl07') then
         p = KEXT_DL07
      else
         p = 'dl07/kext_'//trim(vintage)//'_MW_euv.dat'
      end if
   end function dl07_kext_default


   pure function zubko_kext_tag(qset) result(t)
      ! The suffix a zubko extinction curve carries when it is the size
      ! integral of the non-default optics.  The default set is unmarked, so
      ! the shipped file and group names do not move.
      character(len=*), intent(in) :: qset
      character(len=16) :: t
      t = ''
      if (trim(qset) /= 'zda') t = '_'//trim(qset)
   end function zubko_kext_tag


   pure function zubko_kext_default(qset) result(p)
      ! The text curve behind the product, named the same way.
      character(len=*), intent(in) :: qset
      character(len=96) :: p
      p = 'zubko/kext_zubko_BARE_GR_S'//trim(zubko_kext_tag(qset))//'_euv.dat'
   end function zubko_kext_default


   subroutine load_model_extinction_table(m, default_path, kext_path, ok, default_h5, &
                                          h5_group)
      type(dust_model_t),         intent(inout) :: m
      character(len=*),           intent(in)    :: default_path
      character(len=*), optional, intent(in)    :: kext_path
      ! The model's HDF5 product.  When the caller names no table of its own,
      ! /kext of this file is tried before default_path, and default_path is
      ! what a tree without the product falls back to.
      character(len=*), optional, intent(in)    :: default_h5
      ! Which curve inside that product, '/kext' by default.  A model storing
      ! more than one set of optics stores the matching curve beside each, so
      ! that what dust_extinction serves is the size integral of the very
      ! optics the model was built on.
      character(len=*), optional, intent(in)    :: h5_group
      logical,                    intent(out)   :: ok

      character(len=512)    :: path
      real(wp), allocatable :: alb(:)
      integer               :: n
      logical               :: read_ok

      ok = .true.
      m%kext_n = 0;  m%kext_path = ''
      if (allocated(m%kext_lam))  deallocate(m%kext_lam)
      if (allocated(m%kext_Cext)) deallocate(m%kext_Cext)
      if (allocated(m%kext_Cabs)) deallocate(m%kext_Cabs)
      if (allocated(m%kext_Csca)) deallocate(m%kext_Csca)
      if (allocated(m%kext_gbar)) deallocate(m%kext_gbar)

      if (present(kext_path)) then
         path = kext_path
      else
         path = default_path
         if (present(default_h5)) then
            if (h5_kext_readable(default_h5)) path = default_h5
         end if
      end if
      if (len_trim(path) == 0) return

      ! The albedo column is read and dropped: dust_extinction derives the
      ! albedo from the interpolated Csca/Cext, so that it stays consistent with
      ! the cross sections returned alongside it and every caller gets the same
      ! convention where Cext underflows.
      if (is_hdf5_path(path)) then
         ! /kext of the model's own product.  The WHOLE axis is taken, never the
         ! non-ionizing slice: this curve is interpolated onto the model grid, so
         ! the widest one covers every grid the model can be built on and there
         ! is nothing left for a caller to select wrongly.  That is what retires
         ! the pairing of a narrow model with an _euv-only curve.
         if (present(h5_group)) then
            call read_sedust_kext(trim(path), .true., n, m%kext_lam, alb, m%kext_gbar, &
                                  m%kext_Cext, m%kext_Cabs, m%kext_Csca, read_ok, &
                                  group = h5_group)
         else
            call read_sedust_kext(trim(path), .true., n, m%kext_lam, alb, m%kext_gbar, &
                                  m%kext_Cext, m%kext_Cabs, m%kext_Csca, read_ok)
         end if
      else
         call load_kext_table(trim(path), n, m%kext_lam, alb, m%kext_gbar, &
                              m%kext_Cext, m%kext_Cabs, m%kext_Csca, read_ok)
      end if
      if (allocated(alb)) deallocate(alb)
      if (.not. read_ok) then
         ! A table named by the caller is required, and its failure reaches that
         ! caller as a build status. The default one is optional -- the emission
         ! model is complete without it -- but dropping it in silence would take
         ! dust_extinction away from a host with nothing to show for it.
         ! Announce the loss. The usual causes are a table that has not been
         ! generated yet and one left behind by an earlier wavelength grid,
         ! neither of which the reader can distinguish from a corrupt file.
         if (present(kext_path)) then
            ok = .false.
         else
            write(*,'(a,a)') ' WARNING: could not read the default extinction table ', trim(path)
            write(*,'(a)')   '          The model is built, but dust_extinction has nothing to serve.'
         end if
         return
      end if
      m%kext_n    = n
      m%kext_path = path
   end subroutine load_model_extinction_table


   logical function h5_kext_readable(path) result(yes)
      ! Does this product carry a /kext this build can read?  Asked before the
      ! default is chosen, so that a tree without the file, or a build made
      ! without HDF5, quietly takes the text table instead of announcing a loss
      ! it has not suffered.
      character(len=*), intent(in) :: path
      real(wp), allocatable :: l(:), a(:), g(:), ce(:), ca(:), cs(:)
      integer :: n
      call read_sedust_kext(trim(path), .true., n, l, a, g, ce, ca, cs, yes)
      if (yes) deallocate(l, a, g, ce, ca, cs)
   end function h5_kext_readable


   ! .true. iff the model m carries polarized optics, i.e. at least one
   ! population has an allocated polarized cross section (Cpol on the emission
   ! side or Cpol_ext on the extinction side). This distinguishes an
   ! intentionally scalar-only astrodust model -- one built with
   ! load_polarized_optics = .false., or an implicit-default build whose table
   ! was absent -- from a polarized one, and is .false. for models that never
   ! carry polarized optics (DL07, Zubko). It is independent of the
   ! aligned-scattering table, whose separate load state is visible as scm_loaded.
   pure logical function dust_has_polarized_optics(m) result(has)
      type(dust_model_t), intent(in) :: m
      integer :: ip
      has = .false.
      if (.not. allocated(m%pops)) return
      do ip = 1, size(m%pops)
         if (allocated(m%pops(ip)%Cpol) .or. allocated(m%pops(ip)%Cpol_ext)) then
            has = .true.
            return
         end if
      end do
   end function dust_has_polarized_optics


   ! Option 2: a SINGLE equilibrium temperature for the WHOLE model,
   ! regardless of grain type or size. The total (type- and size-integrated,
   ! dn-weighted) absorption cross section per H,
   !     Cabs_tot(lam) = sum_pop sum_a dn(a) * Cabs(lam, a),
   ! is heated by J_lam to one T_eq (energy balance
   ! int Cabs_tot*J dlam = int Cabs_tot*B(T_eq) dlam), and ALL dust emits
   ! lamI = lam * Cabs_tot * B_lam(T_eq) * 1e-3.  Optional Teq_out returns
   ! that single temperature.
   subroutine dust_emission_single_teq(m, J_lam, lamI_total, Teq_out)
      type(dust_model_t), intent(in)  :: m
      real(wp),           intent(in)  :: J_lam(:)          ! (NLAM)
      real(wp),           intent(out) :: lamI_total(:)     ! (NLAM)
      real(wp), optional, intent(out) :: Teq_out
      real(wp), allocatable :: Cabs_tot(:), kappB_tot(:), spec(:)
      real(wp) :: Teq
      integer  :: ip, ja

      allocate(Cabs_tot(m%NLAM), kappB_tot(m%NT), spec(m%NLAM))
      Cabs_tot = 0.0_wp
      do ip = 1, size(m%pops)
         do ja = 1, size(m%pops(ip)%dn)
            if (m%pops(ip)%dn(ja) <= 0.0_wp) cycle
            Cabs_tot = Cabs_tot + m%pops(ip)%dn(ja) * m%pops(ip)%Cabs(:, ja)
         end do
      end do

      call planck_integral_one(m%lam, Cabs_tot, m%T_first, m%NT, kappB_tot)
      call calc_Teq(m%lam, Cabs_tot, J_lam, m%T_first, kappB_tot, Teq)
      call calc_bbody(Teq, m%lam, spec)
      lamI_total = m%lam * (Cabs_tot * spec) * 1.0e-3_wp
      if (present(Teq_out)) Teq_out = Teq
      deallocate(Cabs_tot, kappB_tot, spec)
   end subroutine dust_emission_single_teq


   ! Planck integral kappB(jt) = int_lam Cabs1(lam) * B_lam(T_first(jt)) dlam
   ! for a SINGLE Cabs1(lam) array, on a denser log-lam grid (same algorithm
   ! as build_kappB, which does it for all sizes). Used by the single-Teq mode.
   subroutine planck_integral_one(lam_in, Cabs1, T_in, ntemp, kappB1)
      real(wp), intent(in)  :: lam_in(:), Cabs1(:), T_in(:)
      integer,  intent(in)  :: ntemp
      real(wp), intent(out) :: kappB1(:)
      integer  :: NW_INT, n_below
      real(wp), allocatable :: w(:), lnw(:), Cross(:), B(:)
      real(wp), allocatable :: lnlam(:)
      real(wp) :: w1, dlnw
      integer  :: nl, jt, iw

      call planck_integration_grid(lam_in, NW_INT, n_below, w1, dlnw)
      allocate(w(NW_INT), lnw(NW_INT), Cross(NW_INT), B(NW_INT))
      nl = size(lam_in)
      allocate(lnlam(nl));  lnlam = log(lam_in)
      do iw = 1, NW_INT
         w(iw)   = w1 * exp(real(iw-1-n_below, wp) * dlnw)
         lnw(iw) = log(w(iw))
      end do
      do iw = 1, NW_INT
         call interp(lnlam, Cabs1, lnw(iw), Cross(iw))
      end do
      do jt = 1, ntemp
         do iw = 1, NW_INT
            B(iw) = bbody(T_in(jt), w(iw))
         end do
         kappB1(jt) = sum(Cross * B * w) * dlnw
      end do
      deallocate(lnlam, w, lnw, Cross, B)
   end subroutine planck_integral_one

end module sed_astrodust_mod
