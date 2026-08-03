module sed_astrodust_mod
   ! Astrodust SED solver. Two-stage API (init + solve) tailored for use
   ! inside a 3D radiative-transfer driver: dust optical and Planck-
   ! integral data are computed once at startup and reused per cell;
   ! only the local mean intensity J_lam varies per call.
   !
   ! API:
   !   call sed_init(qtable_path, sizedist_path, NT_in, T_lo, T_hi)
   !   ...
   !   do icell = 1, ncells
   !      ... compute J_lam in this cell ...
   !      call sed_solve(J_lam, 'S1', lamI_lam)   ! or 'S2'
   !      ... use lamI_lam(:) ...
   !   end do
   !
   ! The driver `main_astrodust.f90` calls this for the single
   ! Mathis-ISRF cell at U_mathis = 1.585 to compare with HD23.
   !
   ! Restructured into init/solve and
   ! limited to a single dust species (astrodust). The stochastic-vs-
   ! equilibrium decision and P(T) solver are unchanged in algorithm.

   use, intrinsic :: iso_fortran_env, only: real64, int64
   use constants,             only: wp
   use sed_mathlib,               only: interp, first_location, last_location, locate
   use radfield,              only: bbody, calc_bbody, hardest_photon_energy, &
                                    cmb_temperature
   use p_sub,                 only: p_sub_setup, calc_Teq, calc_P
   use q_table_mod,           only: load_q_table, &
                                    qt_n_lam=>n_lam, qt_n_aeff=>n_aeff, &
                                    qt_lam=>lam_t, qt_aeff=>aeff_t, &
                                    qt_qext=>qext, qt_qabs=>qabs, qt_qsca=>qsca, &
                                    qt_gpar=>gpar
   use size_dist_mod,         only: load_size_dist, sd_n=>n_size, &
                                    sd_aeff=>a_dist, sd_dn=>dn_ad, &
                                    sd_dn_pah=>dn_pah, sd_fion=>f_ion
   use enthalpy_astrodust_mod, only: enthalpy_S1, enthalpy_S2, RHO_AD, RHO_PAH
   use enthalpy,              only: enthalpy_DL01
   use qpah,                  only: qpah_dl07, qpah_ld01, &
                                    nc_coeff, nc_integer, qpah_use_d03_graphite
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
   use dust_model_mod,        only: dust_model_t, grain_pop_t, free_dust_model
   use kext_table_mod,        only: load_kext_table
   use zubko_io,              only: zda_comp_t, read_zda_config, zda_gofa, &
                                    read_zubko_optics, read_zubko_calor, &
                                    read_dnda_table, ZDA_MAXCOMP
   !$ use omp_lib, only: omp_get_max_threads
   implicit none
   private
   public :: sed_init, sed_solve, sed_solve_pah, sed_solve_qm_batch
   public :: sed_init_dl07, sed_solve_dl07
   ! Injection point for the spheroid (T-matrix) optics of the astrodust EUV
   ! band; see the abstract interface below.
   public :: euv_band_optics_i
   public :: sed_register_euv_band_optics, sed_forget_euv_band_optics
   ! Model-agnostic library API (path B: wraps the untouched solver core).
   public :: dust_model_t, build_astrodust, build_dl07, build_zubko, dust_emission
   public :: build_from_files, dust_emission_single_teq, dust_extinction
   ! The size integral itself, which is what writes the data/kext_*.dat tables
   ! that dust_extinction then serves; see its own header.
   public :: size_integrated_extinction
   ! Dust mass per H of a built model, the constant that turns a cross section
   ! per H into a mass opacity; see its own header.
   public :: dust_mass_per_H
   public :: NLAM, NA, NT, lam, aeff, T_first, dn_ad, dn_pah, initialized
   ! Exposed so that external drivers can cross-check the optics:
   public :: Cabs, Csca, gsca_ad, Cabs_pah, kappB_first, kappB_pah_first
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

   ! Precomputed size-integrated extinction curve each coded model serves to an
   ! RT host through dust_extinction, when the caller names no other file.  The
   ! paths are relative to the sed/ directory the builders are run from.
   !
   ! The default is the EUV product wherever the model has one, because its
   ! wavelength range CONTAINS the plain grid's: both are the same size integral
   ! on the same optics, the extended one merely carried below the T-matrix Q
   ! table's 0.0912 um end, so its nodes coincide with the plain grid's over the
   ! overlap.  One default therefore covers a host that transports ionizing
   ! radiation and one that does not.  The narrower non-EUV product is the file
   ! to name explicitly, not the other way round.  The Zubko model has no EUV
   ! extension and needs none: its own optics table already reaches 1.24 keV.
   character(len=*), parameter :: KEXT_ASTRODUST = '../data/kext_astrodust_MW_euv.dat'
   character(len=*), parameter :: KEXT_DL07      = '../data/kext_dl07_MW_euv.dat'
   character(len=*), parameter :: KEXT_ZUBKO     = '../data/kext_zubko_BARE_GR_S.dat'
   ! SI constants for the induced-emission factor (h*c^2 in J*m^2/s).
   real(wp), parameter :: H_SI    = 6.62606957e-34_wp
   real(wp), parameter :: C_SI    = 2.99792458e8_wp
   real(wp), parameter :: TWO_HCC = 2.0_wp * H_SI * C_SI**2

   ! ---- optics of the EUV band below the Q table ------------------------
   ! sed_init fills the wavelengths it prepends below the Q table's 0.0912 um
   ! end from one of two particles: the volume-equivalent SPHERE (Mie,
   ! q_astrodust_full, selected by euv_tmatrix = .false.), or the b/a = 1.400
   ! oblate SPHEROID the table itself is made of. The spheroid is the
   ! continuation of the table -- same material, same shape, same
   ! random-orientation average -- and it is a T-matrix calculation.
   !
   ! That calculation is INJECTED rather than compiled in: this module holds a
   ! procedure pointer to it, so the SED library carries no reference to
   ! libtmatrix.a and links without it. euv_astrodust_tmatrix.f90 implements
   ! the interface below and registers it in one call; a build that leaves
   ! that file out has no spheroid route, and euv_tmatrix = .true. is then
   ! REFUSED (status 6) instead of being quietly answered with the sphere,
   ! which is a different particle.
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

   logical :: initialized = .false.
   logical, save :: use_induced_emission = .false.

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
   ! Number of wavelength points that were prepended BELOW the optics table's
   ! own short-wavelength end (the EUV extension, see euv_extended_lambda_grid);
   ! 0 for every model built on a plain table grid. The Planck integrals read it
   ! to anchor their internal ln(lambda) grid on the table range, so that neither
   ! the step nor the sample points move as the grid widens -- see
   ! planck_integration_grid.
   integer  :: n_lam_euv = 0
   real(wp), allocatable :: lam(:)              ! [um] (NLAM)
   real(wp), allocatable :: aeff(:)             ! [um] (NA), the size-dist grid
   real(wp), allocatable :: T_first(:)          ! [K]  (NT), log-spaced full range
   real(wp), allocatable :: dn_ad(:)            ! [1/H per bin] (NA)
   real(wp), allocatable :: Cabs(:,:)           ! [cm^2] (NLAM, NA)
   real(wp), allocatable :: Csca(:,:)           ! [cm^2] (NLAM, NA)
   ! Scattering asymmetry <cos> of the astrodust grains, taken from the same
   ! T-matrix Q table as Cabs/Csca and interpolated onto the size grid the same
   ! way. Only the size integral (size_integrated_extinction) uses it; the
   ! emission solver never needs it. Allocated by sed_init, dropped by
   ! sed_init_dl07.
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
   ! ionization fraction f_ion(a) read from size_distribution.dat
   ! column 4. Enthalpy uses DL01 'Car0' (carbonaceous) per HD23 §3.2.
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
   real(wp), allocatable :: dn_cneu(:), dn_cion(:)           ! [1/H per bin] (NA)
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
      !! Undo the above: euv_tmatrix = .true. is then refused (status 6) and
      !! only the volume-equivalent sphere remains.
      euv_band_optics => null()
   end subroutine sed_forget_euv_band_optics

   ! =====================================================================
   subroutine sed_init(qtable_path, sizedist_path, NT_in, T_lo, T_hi, status, lam_min, &
                       astrodust_index_path, euv_tmatrix)
      character(len=*), intent(in) :: qtable_path, sizedist_path
      integer,          intent(in) :: NT_in
      real(wp),         intent(in) :: T_lo, T_hi
      ! Optional status (0 = success). When present, a failed input read is
      ! reported through it instead of stopping the process; when absent the
      ! readers keep their message + stop behavior (as the CLI drivers expect).
      !   status = 1  Q-table load failed
      !   status = 2  size-distribution load failed
      !   status = 3  astrodust dielectric function load failed (EUV band only)
      !   status = 4  lam_min below the astrodust dielectric function's own
      !               shortest wavelength (EUV band only)
      !   status = 6  euv_tmatrix = .true. but the spheroid optics of the EUV
      !               band are not available: no implementation of
      !               euv_band_optics_i is registered, or the registered one
      !               reported that it could not compute the band
      integer, optional, intent(out) :: status
      ! Optional shortest wavelength [um] the model must cover. When it is
      ! shorter than the Q table's 0.0912 um the grid is carried down to it and
      ! the astrodust optics there are computed from the DH21 dielectric
      ! function instead of read off the table. Absent = the table grid alone.
      real(wp), optional, intent(in) :: lam_min
      ! Optional dielectric function for that EUV band. It must be the file the
      ! Q table was computed from -- same porosity, iron fraction and axial
      ! ratio -- or the model changes material at the seam. Omitted = the
      ! q_astrodust_mod default, which pairs with
      ! q_astrodust_P0.20_Fe0.00_1.400.dat.
      character(len=*), optional, intent(in) :: astrodust_index_path
      ! How the EUV band's optics are computed. Default .true.: the T-matrix
      ! on the b/a = 1.400 oblate spheroid, i.e. the same particle and the same
      ! random-orientation average as the Q table, so the grain does not change
      ! shape at the seam. That route is the registered euv_band_optics_i, and
      ! asking for it without one registered is an error (status 6), not a
      ! silent substitution. .false. substitutes the volume-equivalent-sphere
      ! Mie approximation of q_astrodust_mod, which is ~2% low in the
      ! geometric-optics limit (that module measures it size by size) but costs
      ! milliseconds where the T-matrix costs minutes. Ignored when there is no
      ! EUV band.
      logical, optional, intent(in) :: euv_tmatrix
      integer  :: i, ja, jw, jt, is, n_euv
      real(wp) :: a_um, x, t, Q_neu, Q_ion
      real(wp) :: qext1, qsca1, qabs1, gsca1
      real(wp) :: ad_lam_lo, ad_lam_hi
      real(wp), allocatable :: lam_grid(:)
      logical  :: rok
      ! EUV band: route selector, the band's optics, and the cost report.
      logical  :: use_tm
      integer  :: euv_stat, n_euv_bad, n_thread
      integer(int64) :: tick0, tick1, tick_rate
      real(wp), allocatable :: nr_euv(:), ki_euv(:)
      real(wp), allocatable :: q_euv_abs(:,:), q_euv_sca(:,:), q_euv_g(:,:)

      if (present(status)) status = 0
      if (present(astrodust_index_path)) &
         call set_astrodust_index_path(astrodust_index_path)
      use_tm = .true.
      if (present(euv_tmatrix)) use_tm = euv_tmatrix

      ! ---- Load Q table and size dist (modules cache their own state) ----
      if (present(status)) then
         call load_q_table(qtable_path, ok=rok)
         if (.not. rok) then;  status = 1;  return;  end if
         call load_size_dist(sizedist_path, ok=rok)
         if (.not. rok) then;  status = 2;  return;  end if
      else
         call load_q_table(qtable_path)
         call load_size_dist(sizedist_path)
      end if

      call euv_extended_lambda_grid(lam_grid, lam_min, n_extra=n_euv)

      ! ---- EUV band below the Q table -----------------------------------
      ! Its optics come from the DH21 dielectric function, so load that file
      ! HERE, before anything is built, rather than letting the first optics
      ! call load it lazily inside the size loop.  Two reasons: a missing file
      ! or a lam_min the file cannot cover has to reach the caller through
      ! `status`, which is what sed_init promises, and not stop the process out
      ! of an RT host; and the size loop below is threaded, where a lazy load
      ! would be several threads writing one cache at once.
      ! Also name the file, so its pairing with the Q table is visible.
      if (n_euv > 0) then
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
               status = 6;  return
            else
               stop 1
            end if
         end if
         if (present(status)) then
            call load_astrodust_index(ok=rok)
            if (.not. rok) then;  status = 3;  return;  end if
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
               status = 4;  return
            else
               stop 1
            end if
         end if
         write(*,'(a,i0,a)') ' sed_init: EUV extension active, ', n_euv, &
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
      n_lam_euv = n_euv
      NA        = sd_n
      NT        = NT_in

      if (allocated(lam))      deallocate(lam, aeff, T_first, dn_ad, &
                                          Cabs, Csca, kappB_first, H_first, kappCMB, &
                                          log_T_first, log_H_first, log_kappB_first, &
                                          dn_pah, Cabs_pah, kappB_pah_first, &
                                          H_pah_first, kappCMB_pah, &
                                          log_H_pah_first, log_kappB_pah_first)
      allocate(lam(NLAM), aeff(NA), T_first(NT), dn_ad(NA))
      allocate(Cabs(NLAM, NA), Csca(NLAM, NA), kappB_first(NT, NA), &
               H_first(NT, NA, NSTAGE), kappCMB(NA))
      ! Astrodust asymmetry table, allocated under its own guard because
      ! sed_init_dl07 drops it independently of the lam-guarded arrays above.
      if (allocated(gsca_ad)) deallocate(gsca_ad)
      allocate(gsca_ad(NLAM, NA))
      allocate(log_T_first(NT), log_H_first(NT, NA, NSTAGE), &
               log_kappB_first(NT, NA))
      allocate(dn_pah(NA), Cabs_pah(NLAM, NA), kappB_pah_first(NT, NA), &
               H_pah_first(NT, NA), kappCMB_pah(NA))
      allocate(log_H_pah_first(NT, NA), log_kappB_pah_first(NT, NA))
      ! Charge-resolved carbonaceous arrays (neutral + cation kept separate).
      if (allocated(Cabs_cneu)) deallocate(Cabs_cneu, Cabs_cion, dn_cneu, dn_cion, &
               kappB_cneu, kappB_cion, log_kappB_cneu, log_kappB_cion, &
               kappCMB_cneu, kappCMB_cion)
      allocate(Cabs_cneu(NLAM, NA), Cabs_cion(NLAM, NA), dn_cneu(NA), dn_cion(NA), &
               kappB_cneu(NT, NA), kappB_cion(NT, NA), &
               log_kappB_cneu(NT, NA), log_kappB_cion(NT, NA), &
               kappCMB_cneu(NA), kappCMB_cion(NA))

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

      ! ---- Cabs(NLAM, NA) and Csca(NLAM, NA) by interpolating Q in log(a) ----
      ! Q table grid (qt_aeff) is denser than dist grid (aeff); interpolate to
      ! the dist grid where the size sum lives.
      !
      ! Threading.  Each ja writes its own column of Cabs / Csca / gsca_ad and
      ! reads nothing another ja writes, so the loop carries no dependence and
      ! the result cannot depend on the schedule or the thread count.  The
      ! region is entered only when there IS an EUV band (if clause): without
      ! one the loop is three interpolations per size, and the old serial path
      ! is what runs.  With one, the Mie route computes the band here; the
      ! spheroid route was computed just above, threaded on its own, and the
      ! loop only copies it in.
      n_euv_bad = 0
      call system_clock(tick0, tick_rate)

      ! ---- EUV band, spheroid route -------------------------------------
      ! The registered euv_band_optics_i computes the whole band -- every
      ! wavelength on every radius -- in one call, on the DH21 astrodust index.
      ! Same material, same oblate spheroid and same random-orientation average
      ! as the Q table, so the extension continues the table's own calculation
      ! and the grain does not change shape at the seam.  The refractive index
      ! is a function of wavelength alone, so it is evaluated once here, before
      ! any threading, rather than again for every grain radius.
      if (n_euv > 0 .and. use_tm) then
         allocate(nr_euv(n_euv), ki_euv(n_euv))
         do jw = 1, n_euv
            call astrodust_index_at(lam(jw), nr_euv(jw), ki_euv(jw))
         end do
         allocate(q_euv_abs(n_euv, NA), q_euv_sca(n_euv, NA), q_euv_g(n_euv, NA))
         call euv_band_optics(lam(1:n_euv), aeff, nr_euv, ki_euv, &
                              q_euv_abs, q_euv_sca, q_euv_g, euv_stat)
         deallocate(nr_euv, ki_euv)
         if (euv_stat < 0) then
            write(*,'(a)') ' sed_init: the registered EUV band optics could not'
            write(*,'(a)') '           compute the band.'
            if (present(status)) then
               status = 6;  return
            else
               stop 1
            end if
         end if
         ! Points that failed a physical bound: warned about by the routine
         ! itself, counted here for the summary below.
         n_euv_bad = euv_stat
      end if

      !$omp parallel do if(n_euv > 0) default(none) &
      !$omp&   shared(NA, n_euv, use_tm, lam, aeff, Cabs, Csca, gsca_ad, &
      !$omp&          q_euv_abs, q_euv_sca, q_euv_g, &
      !$omp&          qt_aeff, qt_qabs, qt_qsca, qt_gpar) &
      !$omp&   private(ja, jw, a_um, x, qext1, qsca1, qabs1, gsca1) &
      !$omp&   schedule(dynamic)
      do ja = 1, NA
         a_um = aeff(ja)
         x = log(a_um)
         call interp_q_grid(x, qt_aeff, qt_qabs, Cabs(n_euv+1:, ja))
         call interp_q_grid(x, qt_aeff, qt_qsca, Csca(n_euv+1:, ja))
         ! Asymmetry <cos> comes from the same table on the same grid, but it
         ! is already dimensionless -- no pi a^2 conversion.
         call interp_q_grid(x, qt_aeff, qt_gpar, gsca_ad(n_euv+1:, ja))
         ! EUV band below the Q table (n_euv = 0 unless lam_min asked for it).
         ! euv_tmatrix = .false. substitutes the volume-equivalent SPHERE
         ! (Mie) for the spheroid; q_astrodust_mod measures what that costs in
         ! accuracy.
         do jw = 1, n_euv
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
      !$omp end parallel do
      if (allocated(q_euv_abs)) deallocate(q_euv_abs, q_euv_sca, q_euv_g)
      call system_clock(tick1)
      if (n_euv > 0) then
         n_thread = 1
         !$ n_thread = omp_get_max_threads()
         write(*,'(a,i0,a,f9.2,a,i0,a)') '           EUV band: ', n_euv * NA, &
            ' points in ', real(tick1 - tick0, wp) / real(tick_rate, wp), &
            ' s on ', n_thread, ' thread(s).'
         if (n_euv_bad > 0) write(*,'(a,i0,a)') '           ', n_euv_bad, &
            ' of them failed a physical bound (warnings above).'
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
      ! not  E[(1-f_ion)C^neu + f_ion C^ion].  f_ion(a) from
      ! size_distribution.dat column 4. C = pi*(a_cm)^2 * Q convention.
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
      ! from size_distribution.dat col 3, split into neutral/cation by
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
   subroutine sed_init_dl07(qtable_path, sizedist_path, sd_index, u_isrf, &
                            NT_in, T_lo, T_hi, status, lam_min)
      character(len=*), intent(in) :: qtable_path, sizedist_path
      integer,          intent(in) :: sd_index, NT_in
      real(wp),         intent(in) :: u_isrf, T_lo, T_hi
      ! Optional status (0 = success). When present, a failed input read is
      ! reported through it instead of stopping the process; when absent the
      ! readers keep their message + stop behavior (as the CLI drivers expect).
      !   status = 1  Q-table load failed
      !   status = 2  size-distribution load failed
      !   status = 4  lam_min below the D03 dielectric functions' own shortest
      !               wavelength (EUV band only).  There is no code 3 here:
      !               this model needs no separate EUV dielectric function, so
      !               nothing corresponds to sed_init's "index load failed".
      integer, optional, intent(out) :: status
      ! Optional shortest wavelength [um] the model must cover; see sed_init.
      ! The Q table supplies this model's wavelength grid only -- its silicate
      ! and carbonaceous optics are already Mie on the D03 dielectric
      ! functions, which run to 6.2e-5 um, so extending the grid is all the EUV
      ! needs here.
      real(wp), optional, intent(in) :: lam_min

      integer  :: i, ja, jw, jt, n_euv
      real(wp) :: a_um, t, da, qabs1, Q_neu, Q_ion
      real(wp) :: geo, qext1, qsca1, gsca1
      real(wp) :: sil_lam_lo, sil_lam_hi, gra_lam_lo, gra_lam_hi, d03_lam_lo
      real(wp), allocatable :: fion(:), lna(:), lam_grid(:)
      logical  :: rok
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
      if (present(status)) then
         call load_q_table(qtable_path, ok=rok)
         if (.not. rok) then;  status = 1;  return;  end if
         call load_size_dist(sizedist_path, ok=rok)
         if (.not. rok) then;  status = 2;  return;  end if
      else
         call load_q_table(qtable_path)
         call load_size_dist(sizedist_path)
      end if

      call euv_extended_lambda_grid(lam_grid, lam_min, n_extra=n_euv)

      ! This model's optics are Mie on the D03 dielectric functions throughout,
      ! so an EUV extension is a grid extension only -- as long as the grid
      ! stays inside those functions' own coverage. Past it `interp` freezes
      ! (n, k) at the boundary value, which would pass a constant index off as
      ! physics, so refuse instead. Silicate and graphite are both required.
      if (n_euv > 0) then
         call silicate_index_lambda_range(sil_lam_lo, sil_lam_hi)
         call graphite_index_lambda_range(gra_lam_lo, gra_lam_hi)
         d03_lam_lo = max(sil_lam_lo, gra_lam_lo)
         if (lam_grid(1) < d03_lam_lo) then
            write(*,'(a,es10.3,a)') ' sed_init_dl07: lam_min =', lam_grid(1), &
               ' um is shorter than the D03 dielectric functions,'
            write(*,'(a,es10.3,a)') '                which stop at', d03_lam_lo, &
               ' um; (n, k) below it would be frozen at the boundary value.'
            if (present(status)) then
               status = 4;  deallocate(lam_grid);  return
            else
               stop 1
            end if
         end if
      end if

      NLAM      = size(lam_grid)
      n_lam_euv = n_euv
      NA        = NSIZE_BD
      NT        = NT_in

      if (allocated(lam)) deallocate(lam, aeff, T_first, dn_ad, &
            Cabs, Csca, kappB_first, H_first, kappCMB, &
            log_T_first, log_H_first, log_kappB_first, &
            dn_pah, Cabs_pah, kappB_pah_first, H_pah_first, kappCMB_pah, &
            log_H_pah_first, log_kappB_pah_first)
      ! Drop anything a previous astrodust init left behind rather than leave a
      ! stale array on the wrong grid; the DL07 scattering optics below are
      ! computed on this model's own size/wavelength grids.
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
      ! Mie on the D03 astrosilicate dielectric function keeps every return, so
      ! the scattering cross section and its asymmetry come from the same
      ! calculation as the absorption rather than from a tabulated albedo.
      do ja = 1, NA
         a_um = aeff(ja)
         geo  = PI * (a_um * UM2CM)**2
         do jw = 1, NLAM
            call q_silicate_full(a_um, lam(jw), qext1, qsca1, qabs1, gsca1)
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
      do ja = 1, NA
         a_um = aeff(ja)
         geo  = PI * (a_um * UM2CM)**2
         do jw = 1, NLAM
            call qpah_dl07(0, a_um, lam(jw), Q_neu)
            call qpah_dl07(1, a_um, lam(jw), Q_ion)
            Cabs_cneu(jw, ja) = Q_neu * geo
            Cabs_cion(jw, ja) = Q_ion * geo
            ! Scattering of the carbonaceous grain.  Charge shifts the PAH
            ! absorption features only, and a PAH is a molecule whose Rayleigh
            ! scattering is negligible, so both charge states scatter as the
            ! graphite sphere does -- random-orientation Mie (1/3 || + 2/3 perp)
            ! on the graphite dielectric function -- the standard DL07
            ! treatment.
            call q_graphite_full(a_um, lam(jw), qext1, qsca1, qabs1, gsca1)
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
   subroutine sed_grain_loop(npop, dn_pop, aeff_pop, Cabs_pop, kappB_pop, H_pop, &
                              log_H_pop, log_kappB_pop, kappCMB_pop, &
                              J_lam, grain_type, Jout)
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

      integer  :: ir, ii, loc1, loc2, iguard, n_guard_resolve
      integer  :: n_stoch, n_equil_eeq, n_equil_fail
      real(wp) :: Teq, EEQ, del, Tmin_n, Tmax_n, a_cm_qm
      real(wp), allocatable :: spec(:), P(:), lnP(:)
      real(wp), allocatable :: T(:), H(:), kappB(:)
      real(wp), allocatable :: Jout_local(:), emission_qm(:)
      logical :: Equil, Equil_prev, converged, qm_ok

      Jout = 0.0_wp

      select case (trim(stoch_method))

      case ('draine')
         ! ----- Grain-by-grain Draine iterative narrowing -----------------
         ! Each grain decides equilibrium independently (EEQ threshold),
         ! so the loop is OpenMP-parallelizable.
         n_stoch = 0; n_equil_eeq = 0; n_equil_fail = 0
         !$omp parallel default(none) &
         !$omp&   shared(npop, dn_pop, Cabs_pop, kappB_pop, H_pop, &
         !$omp&          log_H_pop, log_kappB_pop, kappCMB_pop, J_lam, &
         !$omp&          T_first, lam, Jout, NLAM, NT) &
         !$omp&   private(ir, ii, Teq, EEQ, Equil, converged, &
         !$omp&           spec, P, lnP, T, H, kappB, Jout_local) &
         !$omp&   reduction(+:n_stoch, n_equil_eeq, n_equil_fail)
         allocate(spec(NLAM), P(NT), lnP(NT), T(NT), H(NT), kappB(NT))
         allocate(Jout_local(NLAM))
         Jout_local = 0.0_wp
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

            if (Equil) then
               call calc_bbody(Teq, lam, spec)
               Jout_local = Jout_local + dn_pop(ir) * Cabs_pop(:, ir) * spec
            else
               do ii = 1, NT
                  if (P(ii) > 0.0_wp) then
                     call calc_bbody(T(ii), lam, spec)
                     Jout_local = Jout_local + dn_pop(ir) * P(ii) * Cabs_pop(:, ir) * spec
                  end if
               end do
            end if
         end do
         !$omp end do
         !$omp critical
         Jout = Jout + Jout_local
         !$omp end critical
         deallocate(spec, P, lnP, T, H, kappB, Jout_local)
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
               ! deliver, hardest_photon_energy(lam, J_lam) -- NOT hc/lam(1),
               ! which is only the short end of the model's optics grid.  The
               ! two coincide for astrodust and DL07, whose Q tables stop at
               ! the Lyman limit, and there hc/0.0912 um = 13.595 eV < 13.6 eV
               ! so the max() returns U_UV1_ERG unchanged.  They diverge for
               ! Zubko/ZDA, whose DustEM tables start at 1.0e-3 um (1.24 keV)
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

            if (Equil) then
               call calc_bbody(Teq, lam, spec)
               Jout = Jout + dn_pop(ir) * Cabs_pop(:, ir) * spec
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
         !$omp&          T_first, lam, grain_type, Jout, NLAM, NT) &
         !$omp&   private(ir, ii, Teq, EEQ, Equil, converged, a_cm_qm, qm_ok, &
         !$omp&           spec, P, lnP, T, H, kappB, Jout_local, emission_qm) &
         !$omp&   reduction(+:n_stoch, n_equil_eeq, n_equil_fail)
         allocate(spec(NLAM), P(NT), lnP(NT), T(NT), H(NT), kappB(NT))
         allocate(Jout_local(NLAM), emission_qm(NLAM))
         Jout_local = 0.0_wp
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
               end if
            end block
         end do
         !$omp end do
         !$omp critical
         Jout = Jout + Jout_local
         !$omp end critical
         deallocate(spec, P, lnP, T, H, kappB, Jout_local)
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
      ! whose DustEM tables reach 1.0e-3 um and would hand back 1.24 keV for a
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

   subroutine euv_extended_lambda_grid(lam_out, lam_min, n_extra)
      ! Model wavelength grid = the T-matrix Q table's grid, optionally carried
      ! into the extreme ultraviolet below its short-wavelength end (0.0912 um
      ! = 13.6 eV, the Lyman limit).  A photoionization RT host transports
      ! 6-100 eV, so half of that band lies off the table.
      !
      ! When lam_min is shorter than the table's first wavelength, log-spaced
      ! points are prepended from lam_min up to just below it.  Their spacing
      ! is at most the table's own spacing at its short-wavelength end
      ! (dln lam = 0.01156), so the extension is never coarser than the grid it
      ! joins.  lam_out(1) is set to lam_min exactly, so the caller's requested
      ! floor is covered rather than approached.  n_extra = 0 (and the plain
      ! table grid) when lam_min is absent, non-positive, or not shorter than
      ! the table -- which is what keeps the unextended model bit-identical.
      real(wp), allocatable, intent(out) :: lam_out(:)
      real(wp), optional,    intent(in)  :: lam_min
      ! Number of points prepended; 0 when the grid is the plain table grid.
      integer,  optional,    intent(out) :: n_extra
      real(wp) :: dln_qt, dln_ext, span
      integer  :: j, nx

      nx = 0
      if (present(lam_min)) then
         if (lam_min > 0.0_wp .and. lam_min < qt_lam(1)) then
            span   = log(qt_lam(1) / lam_min)
            dln_qt = log(qt_lam(2) / qt_lam(1))
            ! span > 0 inside this branch, so the ceiling is already >= 1.
            nx     = ceiling(span / dln_qt)
         end if
      end if
      if (present(n_extra)) n_extra = nx

      allocate(lam_out(nx + qt_n_lam))
      if (nx > 0) then
         dln_ext = log(qt_lam(1) / lam_min) / real(nx, wp)
         do j = 1, nx
            lam_out(j) = qt_lam(1) * exp(-real(nx - j + 1, wp) * dln_ext)
         end do
         lam_out(1) = lam_min
      end if
      lam_out(nx+1:) = qt_lam
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
                      log_H_in, log_kappB_in, kappCMB_in, Csca_in, gsca_in, &
                      rho_bulk)
      type(grain_pop_t), intent(inout) :: p
      character(len=*),  intent(in)    :: gtype
      integer,           intent(in)    :: chan
      real(wp),          intent(in)    :: dn_in(:), kappCMB_in(:)
      real(wp),          intent(in)    :: Cabs_in(:,:), kappB_in(:,:), H_in(:,:)
      real(wp),          intent(in)    :: log_H_in(:,:), log_kappB_in(:,:)
      ! Extinction-side optics, read only by size_integrated_extinction. Each is
      ! optional and independent: a population that does not scatter (the PAHs)
      ! simply leaves them out and contributes zero to those terms of the size
      ! integral. The emission path never touches them.
      real(wp), optional, intent(in)   :: Csca_in(:,:), gsca_in(:,:)
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
      if (present(Csca_in))  p%Csca     = Csca_in
      if (present(gsca_in))  p%gsca     = gsca_in
      if (present(rho_bulk)) p%rho_bulk = rho_bulk
   end subroutine set_pop


   ! Build the HD23 astrodust model into m. Channels: AD_S1, AD_S2, PAH
   ! (PAH = neutral + cation populations summed into one channel).
   subroutine build_astrodust(m, qtable_path, sizedist_path, NT_in, T_lo, T_hi, status, &
                              lam_min, astrodust_index_path, kext_path, euv_tmatrix)
      type(dust_model_t), intent(out) :: m
      character(len=*),   intent(in)  :: qtable_path, sizedist_path
      integer,            intent(in)  :: NT_in
      real(wp),           intent(in)  :: T_lo, T_hi
      ! Optional status (0 = success, non-zero = model build failed). When
      ! present, a failed input read is reported through it instead of stopping
      ! the process; when absent the build stops on error (CLI behavior).
      !   status = 1  Q-table load failed
      !   status = 2  size-distribution load failed
      !   status = 3  astrodust dielectric function load failed (EUV band only)
      !   status = 4  lam_min below that dielectric function's shortest
      !               wavelength (EUV band only)
      !   status = 5  the requested extinction table failed to load
      !   status = 6  euv_tmatrix = .true. but the spheroid optics of the EUV
      !               band are not available; see sed_init
      integer, optional,  intent(out) :: status
      ! Optional shortest wavelength [um] the model must cover, for a host that
      ! transports shortward of the T-matrix Q table's 0.0912 um (13.6 eV) end
      ! -- a photoionization RT spanning 6-100 eV, say. Omitting it gives the
      ! table grid, unchanged.
      real(wp), optional, intent(in)  :: lam_min
      ! Optional dielectric function for that EUV band; must be the file
      ! qtable_path was computed from. See sed_init.
      character(len=*), optional, intent(in) :: astrodust_index_path
      ! Optional precomputed size-integrated extinction curve for
      ! dust_extinction to serve. Omitted, KEXT_ASTRODUST is tried and a model
      ! without it simply has no extinction to serve; named, a file that cannot
      ! be read fails the build.
      character(len=*), optional, intent(in) :: kext_path
      ! How the EUV band's optics are computed; see sed_init. Default .true.
      ! = the T-matrix on the oblate spheroid of the Q table, which requires a
      ! registered euv_band_optics_i and fails with status 6 without one.
      ! .false. = the volume-equivalent-sphere Mie approximation, far cheaper
      ! and ~2% low in the geometric-optics limit.
      logical, optional, intent(in) :: euv_tmatrix
      logical :: kok
      character(len=512) :: kpath

      if (present(status)) status = 0

      ! Astrodust/HD23 optics: Nc=417 (rho=2.0), D16 turbostratic graphite.
      nc_coeff = 417.0d0;  nc_integer = .false.;  qpah_use_d03_graphite = .false.
      call sed_init(qtable_path, sizedist_path, NT_in, T_lo, T_hi, status=status, &
                    lam_min=lam_min, astrodust_index_path=astrodust_index_path, &
                    euv_tmatrix=euv_tmatrix)  ! sets globals
      if (present(status)) then
         if (status /= 0) return
      end if

      m%name = 'astrodust'
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
      ! The astrodust grains carry all the scattering, so Csca / gsca go to this
      ! population alone; the PAHs enter the size integral through absorption only.
      call set_pop(m%pops(1), 'sil', 1, dn_ad, Cabs, kappB_first, H_first(:,:,2), &
                   log_H_first(:,:,2), log_kappB_first, kappCMB, &
                   Csca_in=Csca, gsca_in=gsca_ad, rho_bulk=RHO_AD)
      call set_pop(m%pops(2), 'pah', 2, dn_cneu, Cabs_cneu, kappB_cneu, H_pah_first, &
                   log_H_pah_first, log_kappB_cneu, kappCMB_cneu, rho_bulk=RHO_PAH)
      call set_pop(m%pops(3), 'pah', 2, dn_cion, Cabs_cion, kappB_cion, H_pah_first, &
                   log_H_pah_first, log_kappB_cion, kappCMB_cion, rho_bulk=RHO_PAH)

      call load_model_extinction_table(m, KEXT_ASTRODUST, kext_path, kok, kpath)
      if (.not. kok) then
         if (present(status)) then
            status = 5;  return
         else
            write(*,'(a,a)') ' build_astrodust: cannot read the extinction table ', trim(kpath)
            stop 1
         end if
      end if
   end subroutine build_astrodust


   ! Build the DL07 model into m. Channels: SIL, CARB (carbonaceous =
   ! neutral + cation summed). Reuses sed_init_dl07 to set the globals.
   subroutine build_dl07(m, qtable_path, sizedist_path, sd_index, u_isrf, &
                         NT_in, T_lo, T_hi, status, lam_min, kext_path)
      type(dust_model_t), intent(out) :: m
      character(len=*),   intent(in)  :: qtable_path, sizedist_path
      integer,            intent(in)  :: sd_index, NT_in
      real(wp),           intent(in)  :: u_isrf, T_lo, T_hi
      ! Optional status (0 = success, non-zero = model build failed). When
      ! present, a failed input read is reported through it instead of stopping
      ! the process; when absent the build stops on error (CLI behavior).
      !   status = 1  Q-table load failed
      !   status = 2  size-distribution load failed
      !   status = 4  lam_min below the D03 dielectric functions' shortest
      !               wavelength (EUV band only); see sed_init_dl07
      !   status = 5  the requested extinction table failed to load
      integer, optional,  intent(out) :: status
      ! Optional shortest wavelength [um] the model must cover; see
      ! build_astrodust. This model's optics are dielectric-function Mie
      ! throughout, so the extension is a grid extension only.
      real(wp), optional, intent(in)  :: lam_min
      ! Optional precomputed size-integrated extinction curve for
      ! dust_extinction to serve; see build_astrodust. Omitted, KEXT_DL07.
      character(len=*), optional, intent(in) :: kext_path
      logical :: kok
      character(len=512) :: kpath

      if (present(status)) status = 0

      ! DL07 carbonaceous optics (matching Draine): Nc=470 (rho~2.2, NINT),
      ! D03 graphite, and the Draine-2003a 0.93 abundance reduction.
      nc_coeff = 470.0d0;  nc_integer = .true.;  qpah_use_d03_graphite = .true.
      gd_apply_d03_reduction = .true.
      call sed_init_dl07(qtable_path, sizedist_path, sd_index, u_isrf, NT_in, T_lo, T_hi, &
                         status=status, lam_min=lam_min)
      if (present(status)) then
         if (status /= 0) return
      end if

      m%name = 'dl07'
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
      ! Both the silicate and the two carbonaceous charge states scatter, so all
      ! three populations carry their scattering optics into the size integral.
      ! The charge states share one scattering description (the graphite sphere);
      ! they differ in dn and in absorption only.
      call set_pop(m%pops(1), 'sil', 1, dn_ad, Cabs, kappB_first, H_first(:,:,1), &
                   log_H_first(:,:,1), log_kappB_first, kappCMB, &
                   Csca_in=Csca, gsca_in=gsca_ad, rho_bulk=RHO_ASTROSIL)
      call set_pop(m%pops(2), 'pah', 2, dn_cneu, Cabs_cneu, kappB_cneu, H_pah_first, &
                   log_H_pah_first, log_kappB_cneu, kappCMB_cneu, &
                   Csca_in=Csca_car, gsca_in=gsca_car, rho_bulk=RHO_GRAPHITE)
      call set_pop(m%pops(3), 'pah', 2, dn_cion, Cabs_cion, kappB_cion, H_pah_first, &
                   log_H_pah_first, log_kappB_cion, kappCMB_cion, &
                   Csca_in=Csca_car, gsca_in=gsca_car, rho_bulk=RHO_GRAPHITE)

      call load_model_extinction_table(m, KEXT_DL07, kext_path, kok, kpath)
      if (.not. kok) then
         if (present(status)) then
            status = 5;  return
         else
            write(*,'(a,a)') ' build_dl07: cannot read the extinction table ', trim(kpath)
            stop 1
         end if
      end if
   end subroutine build_dl07


   ! Build the Zubko (ZDA 2004) BARE-GR-S model into m. Three components
   ! (PAH, Graphite, Silicate), each with its OWN size grid (the component's
   ! DustEM Q-table radii), so the populations carry component-by-component grids.
   ! Size distribution from the ZDA formula; optics from the DustEM Q-tables
   ! (Cabs = Qabs*pi*a^2); enthalpy from the specific-heat calorimetry tables
   ! (H = u_spec(T)*rho*(4pi/3)a^3). The shared lambda grid is the optics
   ! grid (all 3 components share 1201 wavelengths). Channels: PAH, GRA, SIL.
   subroutine build_zubko(m, config_path, data_dir, NT_in, T_lo, T_hi, status, kext_path)
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
      !   status = 6  the requested extinction table failed to load
      integer, optional,  intent(out) :: status
      ! Optional precomputed size-integrated extinction curve for
      ! dust_extinction to serve; see build_astrodust. Omitted, KEXT_ZUBKO.
      character(len=*), optional, intent(in) :: kext_path

      type(zda_comp_t)      :: comps(ZDA_MAXCOMP)
      integer               :: ncomp, ic, jt, ja, jw, nsize, nwave, ntc
      real(wp)              :: rho, vol_fac, mass, dlna, uspec, t, wdev
      real(wp), allocatable :: a_opt(:), lam_opt(:), qa(:,:), qs(:,:), gg(:,:)
      real(wp), allocatable :: Tcal(:), Ucal(:), Ccal(:), Hcol(:)
      logical               :: rok, kok
      character(len=512)    :: kpath
      character(len=16)     :: cn(3)
      character(len=8)      :: gt(3)
      character(len=64)     :: optf

      if (present(status)) status = 0

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
      m%use_induced_emission = use_induced_emission
      m%stoch_method = stoch_method
      m%n_channel = 3
      allocate(m%channel_name(3));  m%channel_name = cn
      allocate(m%pops(3))

      do ic = 1, 3
         ! Cross-section file name from the config ('Cross Sections=...').
         optf = trim(comps(ic)%xsec)//'.dat'
         ! The ZDA model is defined by these tables -- they are Zubko's own
         ! multilayer-sphere solution, so the scattering side is read from the
         ! same file as the absorption (Q_sca and <cos> columns) rather than
         ! recomputed from a dielectric function this model does not supply.
         if (present(status)) then
            call read_zubko_optics(trim(data_dir)//trim(optf), nsize, nwave, &
                                   a_opt, lam_opt, qa, qs, rho, ok=rok, gpar=gg)
            if (.not. rok) then;  status = 3;  return;  end if
         else
            call read_zubko_optics(trim(data_dir)//trim(optf), nsize, nwave, &
                                   a_opt, lam_opt, qa, qs, rho, gpar=gg)
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
            ! This model's grid is its own optics table, with no EUV extension.
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
            ! rho is the component's own density from its DustEM optics file,
            ! and the next component overwrites the scalar, so record it here.
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
               if (ja == 1) then
                  dlna = log(a_opt(2)/a_opt(1))
               else if (ja == nsize) then
                  dlna = log(a_opt(nsize)/a_opt(nsize-1))
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

      call load_model_extinction_table(m, KEXT_ZUBKO, kext_path, kok, kpath)
      if (.not. kok) then
         if (present(status)) then
            status = 6;  return
         else
            write(*,'(a,a)') ' build_zubko: cannot read the extinction table ', trim(kpath)
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
   ! Each population's optics is a DustEM Q-table, the size distribution is a
   ! 2-column a[um] dn/da[cm^-1 H^-1] table, and the enthalpy a specific-heat
   ! calorimetry table; all files are sought under data_dir. This is the
   ! data-driven path (build_astrodust/dl07/zubko are the coded builders).
   subroutine build_from_files(m, descriptor_path, data_dir, NT_in, T_lo, T_hi, status, &
                               kext_path)
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
      !   status = 9  the requested extinction table failed to load
      integer, optional,  intent(out) :: status
      ! Optional precomputed size-integrated extinction curve for
      ! dust_extinction to serve. A file-defined model has no default: the
      ! descriptor names optics, sizes and calorimetry but not a precomputed
      ! extinction curve, so unless this is given the model has none to serve.
      character(len=*), optional, intent(in) :: kext_path

      integer, parameter :: MAXP = 16
      character(len=8)   :: p_gt(MAXP)
      integer            :: p_ch(MAXP)
      character(len=64)  :: p_opt(MAXP), p_dn(MAXP), p_cal(MAXP)
      real(wp)           :: p_rho(MAXP)
      integer            :: npop, u, ios, ip, jt, ja, jw, nsize, nwave, ntc, ndn, nchan, ic, nline
      real(wp)           :: t, rho, mass, vf, dlna, uspec, fa, loga, wdev
      logical            :: rok, kok
      character(len=256) :: line
      character(len=512) :: kpath
      real(wp), allocatable :: a_opt(:), lam_opt(:), qa(:,:), qs(:,:), gg(:,:)
      real(wp), allocatable :: a_dn(:), f_dn(:), la_dn(:), lf_dn(:), Tc(:), Uc(:), Cc(:)

      if (present(status)) status = 0

      npop = 0;  nline = 0;  m%name = 'file_model'
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
            ! This model's grid is its own optics table, with no EUV extension.
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
               if (ja == 1) then
                  dlna = log(a_opt(2)/a_opt(1))
               else if (ja == nsize) then
                  dlna = log(a_opt(nsize)/a_opt(nsize-1))
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
            ! Density of this population's material -- the descriptor's value
            ! where it gave one, otherwise the optics file's own header.
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

      call load_model_extinction_table(m, '', kext_path, kok, kpath)
      if (.not. kok) then
         if (present(status)) then
            status = 9;  return
         else
            write(*,'(a,a)') ' build_from_files: cannot read the extinction table ', trim(kpath)
            stop 1
         end if
      end if
   end subroutine build_from_files


   ! Generic single-cell dust emission for the (active) model m. Loops the
   ! populations through the untouched sed_grain_loop, sums per output
   ! channel, applies the induced factor (if enabled) and the HD23 unit
   ! convention. lamI_total(NLAM) is the summed SED; optional
   ! lamI_chan(NLAM, n_channel) returns the SED of each channel.
   ! REQUIRES: m is the most recently built model (its grids == the globals).
   subroutine dust_emission(m, J_lam, lamI_total, lamI_chan, status)
      type(dust_model_t), intent(in)  :: m
      real(wp),           intent(in)  :: J_lam(:)              ! (NLAM)
      real(wp),           intent(out) :: lamI_total(:)         ! (NLAM)
      real(wp), optional, intent(out) :: lamI_chan(:,:)        ! (NLAM, n_channel)
      ! Optional error report (0 = success). When present, a bad model is
      ! reported through it instead of stopping the process; when absent the
      ! original stop-on-error behavior is kept (as the CLI drivers expect).
      !   status = 1  unknown stoch_method
      !   status = 2  'qm' selected but a population is missing its radii
      integer,  optional, intent(out) :: status
      real(wp), allocatable :: Jout_pop(:), Jchan(:,:)
      integer :: ip, ic

      if (present(status)) status = 0

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
      Jchan = 0.0_wp
      do ip = 1, size(m%pops)
         ! size count for each population (Zubko-like models have
         ! component-by-component grids)
         call sed_grain_loop(size(m%pops(ip)%dn), m%pops(ip)%dn, m%pops(ip)%aeff, &
                             m%pops(ip)%Cabs, &
                             m%pops(ip)%kappB, m%pops(ip)%H, m%pops(ip)%log_H, &
                             m%pops(ip)%log_kappB, m%pops(ip)%kappCMB, &
                             J_lam, trim(m%pops(ip)%grain_type), Jout_pop)
         ic = m%pops(ip)%out_channel
         Jchan(:, ic) = Jchan(:, ic) + Jout_pop
      end do

      do ic = 1, m%n_channel
         if (m%use_induced_emission) call apply_induced_factor(J_lam, Jchan(:, ic))
         Jchan(:, ic) = m%lam * Jchan(:, ic) * 1.0e-3_wp
      end do

      lamI_total = sum(Jchan, dim=2)
      if (present(lamI_chan)) lamI_chan = Jchan
      deallocate(Jout_pop, Jchan)
   end subroutine dust_emission


   ! Transport optics an RT host takes from the model: the size-integrated
   ! extinction curve of data/kext_*.dat, returned on the model's own
   ! wavelength grid m%lam and per H atom.  This is the extinction twin of
   ! dust_emission -- one model object supplies both, so the light a ray loses
   ! and the light the cell puts back refer to the same grains.
   !
   ! The curve is READ FROM THE PRECOMPUTED TABLE the builder loaded
   ! (m%kext_path); it is not recomputed here.  The size integral itself lives
   ! in size_integrated_extinction, which is what wrote that table, so a host
   ! transports exactly the numbers that were checked against the HD23 release
   ! and Draine's own kext_albedo tables.  The table has to be read at BUILD
   ! time because its path is relative to the sed/ directory, and a host is free
   ! to change directory before it ever calls this routine.
   !
   ! Interpolation onto m%lam, needed only when the table grid and the model
   ! grid differ: log(lambda)-log(C) linear for Cext, Cabs and Csca, which are
   ! positive and locally power-law in lambda; linear in log(lambda) for <cos>,
   ! which changes sign and so has no logarithm.  Where a model wavelength
   ! coincides with a table wavelength the table value is returned unchanged.
   ! Outside the tabulated range there is no optics to serve, so the call is
   ! refused (status 3) rather than extrapolated.
   !
   ! albedo is re-derived as C_sca/C_ext instead of interpolating the table's
   ! own albedo column, because the cross-section columns are written with more
   ! digits; 0 where C_ext underflows.  <cos> is taken from the table as it
   ! stands -- it is already the scattering-weighted average, and re-weighting
   ! it would need the size-resolved optics the table does not carry.
   !
   ! Units: all cross sections [cm^2/H]; gbar and albedo dimensionless.
   subroutine dust_extinction(m, Cext, Cabs, Csca, gbar, albedo, status)
      type(dust_model_t), intent(in)  :: m
      real(wp),           intent(out) :: Cext(:), Cabs(:), Csca(:)   ! (NLAM) [cm^2/H]
      ! Scattering-weighted asymmetry <cos>.
      real(wp), optional, intent(out) :: gbar(:)                     ! (NLAM)
      ! Scattering albedo C_sca/C_ext; 0 where the medium is transparent.
      ! Derived here so that every caller gets the same convention at the
      ! wavelengths where C_ext underflows to zero.
      real(wp), optional, intent(out) :: albedo(:)                   ! (NLAM)
      ! Optional error report (0 = success). When present, a bad request is
      ! reported through it instead of stopping the process; when absent such a
      ! call stops the run, matching dust_emission.
      !   status = 1  an output array is not of size m%NLAM
      !   status = 2  no extinction table was loaded for this model
      !   status = 3  m%lam reaches outside the table's wavelength range
      integer,  optional, intent(out) :: status
      ! A model wavelength can miss a table wavelength by the precision the
      ! table's lambda column is written with -- six significant digits in the
      ! oldest product, so up to 5e-6 in relative terms -- and that must not
      ! read as "outside the table".  This allowance covers such rounding and is
      ! still three orders of magnitude below the grid spacing
      ! dln(lambda) = 0.0116, so it cannot mask a real gap.
      real(wp), parameter :: EDGE_TOL = 1.0e-5_wp
      real(wp) :: lam0, t, g_here
      integer  :: jw, j, n
      logical  :: bad

      if (present(status)) status = 0

      bad = size(Cext) /= m%NLAM .or. size(Cabs) /= m%NLAM .or. size(Csca) /= m%NLAM
      if (present(gbar))   bad = bad .or. size(gbar)   /= m%NLAM
      if (present(albedo)) bad = bad .or. size(albedo) /= m%NLAM
      if (bad) then
         if (present(status)) then
            status = 1;  return
         else
            write(*,'(a,i0)') 'dust_extinction: output arrays must be of size m%NLAM=', m%NLAM
            stop 1
         end if
      end if

      n = m%kext_n
      if (n <= 0) then
         if (present(status)) then
            status = 2;  return
         else
            write(*,'(a)') 'dust_extinction: no extinction table was loaded for this model'
            write(*,'(a)') '   name the file through the builder''s kext_path argument'
            stop 1
         end if
      end if

      do jw = 1, m%NLAM
         if (log(m%lam(jw)/m%kext_lam(1)) < -EDGE_TOL .or. &
             log(m%lam(jw)/m%kext_lam(n)) >  EDGE_TOL) then
            if (present(status)) then
               status = 3;  return
            else
               write(*,'(a,es12.5,a)') 'dust_extinction: lambda = ', m%lam(jw), &
                    ' um is outside the extinction table'
               write(*,'(a,es12.5,a,es12.5,a,a)') '   table covers ', m%kext_lam(1), ' - ', &
                    m%kext_lam(n), ' um: ', trim(m%kext_path)
               stop 1
            end if
         end if
      end do

      do jw = 1, m%NLAM
         lam0 = m%lam(jw)
         call locate(m%kext_lam, lam0, j)
         j = min(max(j, 1), n-1)               ! bracket [j, j+1]
         if (lam0 == m%kext_lam(j)) then
            Cext(jw) = m%kext_Cext(j);  Cabs(jw) = m%kext_Cabs(j)
            Csca(jw) = m%kext_Csca(j);  g_here   = m%kext_gbar(j)
         else if (lam0 == m%kext_lam(j+1)) then
            Cext(jw) = m%kext_Cext(j+1);  Cabs(jw) = m%kext_Cabs(j+1)
            Csca(jw) = m%kext_Csca(j+1);  g_here   = m%kext_gbar(j+1)
         else
            ! Fraction of the ln(lambda) interval, clamped so that a wavelength
            ! sitting within EDGE_TOL of an end returns that end's value rather
            ! than a hair of extrapolation.
            t = log(lam0/m%kext_lam(j)) / log(m%kext_lam(j+1)/m%kext_lam(j))
            t = min(max(t, 0.0_wp), 1.0_wp)
            Cext(jw) = power_law_in_lambda(m%kext_Cext(j), m%kext_Cext(j+1), t)
            Cabs(jw) = power_law_in_lambda(m%kext_Cabs(j), m%kext_Cabs(j+1), t)
            Csca(jw) = power_law_in_lambda(m%kext_Csca(j), m%kext_Csca(j+1), t)
            g_here   = (1.0_wp - t)*m%kext_gbar(j) + t*m%kext_gbar(j+1)
         end if
         if (present(gbar))   gbar(jw) = g_here
         if (present(albedo)) then
            albedo(jw) = 0.0_wp
            if (Cext(jw) > 0.0_wp) albedo(jw) = Csca(jw) / Cext(jw)
         end if
      end do
   end subroutine dust_extinction


   pure function power_law_in_lambda(y1, y2, t) result(y)
      ! Value at fraction t of a ln(lambda) interval whose ends are y1 and y2,
      ! taking the cross section to follow a local power law C ~ lambda^p --
      ! i.e. linear in ln C.  Where an end is not positive the power law is not
      ! defined, and the interpolation falls back to linear in ln(lambda).
      real(wp), intent(in) :: y1, y2, t
      real(wp) :: y
      if (y1 > 0.0_wp .and. y2 > 0.0_wp) then
         y = exp((1.0_wp - t)*log(y1) + t*log(y2))
      else
         y = (1.0_wp - t)*y1 + t*y2
      end if
   end function power_law_in_lambda


   ! Size-distribution-integrated extinction of the (active) model m, per H
   ! atom -- the calculation that produces the data/kext_*.dat curves which
   ! dust_extinction then serves to an RT host.
   !
   ! The size integral is the plain binned sum over each population, because
   ! dn(a) already carries the bin width:
   !   C_abs/H  = sum_pop sum_a dn(a) * Cabs(lambda, a)
   !   C_sca/H  = sum_pop sum_a dn(a) * Csca(lambda, a)
   !   C_ext/H  = C_abs/H + C_sca/H
   !   <cos>    = sum dn * Csca * g  /  sum dn * Csca      (scattering-weighted)
   ! A population whose Csca / gsca are unallocated contributes zero to those
   ! terms; in the astrodust model the PAHs are exactly that case, so they enter
   ! through absorption only and the astrodust grains carry all the scattering.
   !
   ! Units: all cross sections [cm^2/H]; gbar dimensionless.
   ! REQUIRES: m is the most recently built model (its grids == the globals).
   subroutine size_integrated_extinction(m, Cext, Cabs, Csca, gbar, albedo, status)
      type(dust_model_t), intent(in)  :: m
      real(wp),           intent(out) :: Cext(:), Cabs(:), Csca(:)   ! (NLAM) [cm^2/H]
      ! Scattering-weighted asymmetry; 0 where nothing scatters.
      real(wp), optional, intent(out) :: gbar(:)                     ! (NLAM)
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
      logical :: bad

      if (present(status)) status = 0

      bad = size(Cext) /= m%NLAM .or. size(Cabs) /= m%NLAM .or. size(Csca) /= m%NLAM
      if (present(gbar))   bad = bad .or. size(gbar)   /= m%NLAM
      if (present(albedo)) bad = bad .or. size(albedo) /= m%NLAM
      if (bad) then
         if (present(status)) then
            status = 1;  return
         else
            write(*,'(a,i0)') 'size_integrated_extinction: output arrays must be of size m%NLAM=', &
                              m%NLAM
            stop 1
         end if
      end if

      allocate(gnum(m%NLAM))
      Cabs = 0.0_wp;  Csca = 0.0_wp;  gnum = 0.0_wp

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

      Cext = Cabs + Csca
      if (present(gbar)) then
         gbar = 0.0_wp
         do jw = 1, m%NLAM
            if (Csca(jw) > 0.0_wp) gbar(jw) = gnum(jw) / Csca(jw)
         end do
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


   ! Load the precomputed size-integrated extinction curve that this model will
   ! serve through dust_extinction.  Called at the end of every builder, because
   ! the table path is relative to the directory the builder runs in and a host
   ! may change directory afterwards.
   !
   ! A path the CALLER named is a contract: failing to read it is a build
   ! failure, reported through ok = .false., because a host that asks for a file
   ! that is not there has a configuration error.  The model's DEFAULT path is
   ! only an offer: if it cannot be read, the model is left with kext_n = 0 and
   ! the build succeeds, so that a standalone driver which only computes
   ! emission still runs where no table has been generated yet -- and so that
   ! calc_kext.x can build the very model whose table it is about to write.
   ! An empty path names no file, so default_path = '' means the model has no
   ! default at all, and kext_path = '' asks for no table rather than for a file
   ! that cannot be opened.
   subroutine load_model_extinction_table(m, default_path, kext_path, ok, path)
      type(dust_model_t), intent(inout) :: m
      character(len=*),   intent(in)  :: default_path
      character(len=*), optional, intent(in) :: kext_path
      logical,            intent(out) :: ok
      ! The file that was tried, for the caller's error message.
      character(len=*),   intent(out) :: path
      real(wp), allocatable :: albedo_t(:)
      logical :: required, tok

      ok = .true.
      m%kext_n = 0;  m%kext_path = ''
      required = present(kext_path)
      if (required) then
         path = kext_path
      else
         path = default_path
      end if
      if (len_trim(path) == 0) return

      call load_kext_table(trim(path), m%kext_n, m%kext_lam, albedo_t, m%kext_gbar, &
                           m%kext_Cext, m%kext_Cabs, m%kext_Csca, tok)
      ! The albedo column is read for validation only; dust_extinction re-derives
      ! the albedo from C_sca/C_ext, which is written with more digits.
      if (allocated(albedo_t)) deallocate(albedo_t)
      if (tok) then
         m%kext_path = path
      else
         m%kext_n = 0
         ! A table named by the caller is required, and its failure reaches that
         ! caller as a build status.  The default one is optional -- the
         ! emission model is complete without it -- but dropping it in silence
         ! would take dust_extinction away from a host with nothing to show for
         ! it.  Announce the loss.  The usual causes are a table that has not
         ! been generated yet and one left behind by an earlier wavelength grid,
         ! neither of which the reader can distinguish from a corrupt file.
         if (.not. required) then
            write(*,'(a,a)') ' WARNING: could not read the default extinction table ', trim(path)
            write(*,'(a)')   '          The model is built, but dust_extinction has nothing to serve.'
         end if
         ok = .not. required
      end if
   end subroutine load_model_extinction_table


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
