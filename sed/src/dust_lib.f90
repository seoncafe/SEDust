module dust_lib
   ! RT-facing entry point for the model-agnostic dust thermal-emission
   ! library. A Fortran 3D radiative-transfer code links this module and:
   !
   !   use dust_lib
   !   type(dust_model_t) :: m
   !   call build_astrodust(m, qtab, sizedist, NT, T_lo, T_hi)   ! once
   !   ...
   !   do icell = 1, ncells
   !      ! ... assemble local mean intensity J_lam(:) on m%lam ...
   !      call dust_emission(m, J_lam, lamI_total [, lamI_chan])  ! per cell
   !   end do
   !
   ! EXTREME ULTRAVIOLET. The astrodust model's wavelength grid is the grid of
   ! the T-matrix Q table it was built from, and TWO tables ship side by side:
   !
   !   q_astrodust_P0.20_Fe0.00_1.400.dat      1129 wavelengths, 0.0912-39810 um
   !   q_astrodust_P0.20_Fe0.00_1.400_euv.dat  1762 wavelengths, 1.0e-4-39810 um
   !
   ! The second carries the same calculation below the Lyman limit, to 12398 eV;
   ! the first is the second with those wavelengths dropped, so over the 1129
   ! they share they are the same numbers cell for cell. WHICH ONE THE HOST
   ! PASSES AS qtab DECIDES WHAT ITS MODEL COVERS. Given the _euv table, the
   ! ionizing band is INSIDE the table, computed the same way as the rest of it,
   ! so a host that transports ionizing radiation reads it off the table and
   ! needs nothing further -- in particular it does not have to link
   ! libtmatrix.a.
   !
   ! Below 0.0912 um the _euv table's wavelengths are the DH21 dielectric
   ! function's own energy nodes rather than a resampling of them, so each
   ! absorption edge it carries stays an exact step between adjacent grid points
   ! instead of being averaged into a ramp. One added point is not a node,
   ! 0.0912*(1-1e-4): it resolves the Lyman-limit step of the radiation FIELD,
   ! not anything in the grain.
   !
   ! build_astrodust and build_dl07 take the optional lam_min [um]:
   !   call build_dl07(m, qtab, sizedist, NT, T_lo, T_hi, lam_min=6.21e-5_wp)
   ! Log-spaced points are then prepended from lam_min up to the table, no more
   ! coarsely than the table's own dln(lambda) at its short end (taken from the
   ! table rather than written down: 0.00794 on the _euv axis, 0.01156 on the
   ! other), and m%lam is that longer grid. OMITTING lam_min, or passing one the
   ! table already covers, prepends nothing.
   !
   ! WHAT lam_min DOES FOR ASTRODUST DEPENDS ON WHICH TABLE IT WAS GIVEN. On the
   ! 1129-wavelength table a lam_min below 0.0912 um builds a band, and the
   ! optics in it come from the injected route described next. On the _euv table
   ! there is nothing left to ask for: the DH21 dielectric function stops at
   ! 1.00003e-4 um, longward of that table's own first wavelength, so every legal
   ! lam_min either falls inside the table -- no band, no work -- or below the
   ! dielectric data, where it is refused rather than served with a frozen
   ! refractive index. DL07 can be extended below either table, to 6.205e-5 um,
   ! because the D03 optical constants reach further than both; its optics are
   ! Mie on those functions at every wavelength, so that extension needs no
   ! T-matrix either.
   !
   ! WHAT THE INJECTED ROUTE IS STILL FOR. When a band IS computed (n_euv > 0),
   ! the astrodust optics there come from the DH21 dielectric function by the
   ! T-MATRIX, on the same b/a = 1.400 oblate spheroid and the same
   ! random-orientation average as the Q table itself (spheroid_q, which takes
   ! the Rayleigh dipole limit below x = 0.1 and the geometric-optics limit
   ! above x = 50, exactly as the table's own driver does), so the band
   ! continues the table rather than gluing a different particle to it:
   ! measured at the seam, C_abs matched a log-log extrapolation of the table
   ! to 0.07% at worst over the whole size range, and <cos> agreed exactly
   ! where the table's own regime bounds put it at 0.
   !
   ! That spheroid is the ONLY thing in this library that needs the T-matrix,
   ! and it is not compiled in by default. The plain libsedust.a therefore has
   ! no T-matrix in it and links without one; asking it for the spheroid
   ! (euv_tmatrix = .true., the default whenever lam_min is given) is REFUSED
   ! with status 6 rather than answered with a sphere. A model built on the
   ! 1129-wavelength table with a lam_min below 0.0912 um reaches that refusal;
   ! one built on the _euv table never does, because it never forms a band. To
   ! get the spheroid:
   !
   !   cd tmatrix && make libtmatrix.a
   !   cd sed     && WITH_TMATRIX=1 make libsedust.a     # or ./build_lib.sh
   !
   ! and, once in the host, before the model is built:
   !
   !   use euv_astrodust_tmatrix, only: use_tmatrix_euv_band_optics
   !   call use_tmatrix_euv_band_optics()
   !
   ! (Or register an implementation of your own: sed_register_euv_band_optics
   ! takes any procedure matching the abstract interface euv_band_optics_i,
   ! both re-exported here, and sed_forget_euv_band_optics undoes it.)
   !
   ! That costs real time, and THAT COST IS WHY THE TABLE WAS EXTENDED INSTEAD:
   ! the band below is now precomputed once, not resolved per build. Measured
   ! here on 167 radii while it was still solved on the fly (gfortran 13.1,
   ! -O3), band only:
   !
   !   lam_min                 n_euv  band pts  T-matrix pts  threads    time
   !   0.0247968 um (50 eV)      113    18,871        12,204        1   527.5 s
   !                                                               2   264.2 s
   !                                                               4   136.7 s
   !                                                               8    69.2 s
   !   1.001e-4 um (12.4 keV)    590    98,530        41,900       72    73.5 s
   !
   ! so the band is built with OpenMP (one T-matrix workspace for each thread,
   ! 80 MiB apiece once its full-direct arrays are allocated; nothing is
   ! allocated when there is no EUV band). Peak resident memory grows with it:
   ! 64 / 110 / 197 / 376 MB at 1 / 2 / 4 / 8 threads, 3.0 GB at 72. The size
   ! loop carries no dependence, so the result does not depend on the thread
   ! count or the schedule -- the four runs above were byte-identical.
   !
   ! A host that computes a band and cannot afford this passes
   ! euv_tmatrix = .false.:
   !   call build_astrodust(m, qtab, sizedist, NT, T_lo, T_hi, &
   !                        lam_min=0.0124_wp, euv_tmatrix=.false.)
   ! which substitutes Bohren-Huffman Mie for the volume-equivalent SPHERE and
   ! runs in a tenth of a second. It is an APPROXIMATION: good to ~2% where the
   ! T-matrix regime holds, but converging on a fixed ~2% underestimate for
   ! grains in the geometric-optics limit rather than on zero, and leaving a
   ! visible step at the seam (median 1.2%, worst 16% in C_abs; <cos> jumps
   ! from 0.9 to 0 for a >~ 1 um). q_astrodust.f90 documents it size by size.
   ! On the _euv table the flag selects nothing, for the reason above -- there
   ! is no band to compute -- and neither shipped kext product goes through the
   ! sphere either, because calc_kext.x passes no lam_min for astrodust and
   ! reads the whole range off whichever table it was given.
   !
   ! The carbonaceous optics have no such seam: above the DL07 PAH cutoff
   ! (21.4 eV) they are Mie on the D03 graphite dielectric functions, which
   ! run to 6.2e-5 um.
   !
   ! The dielectric function must be the one the Q table was computed from --
   ! same porosity, iron fraction and axial ratio -- or the model changes
   ! material at the seam. build_astrodust takes it as the optional
   ! astrodust_index_path for that reason, and reports the file whenever the
   ! extension is active; omitted, it is the P = 0.20, f_Fe = 0.00, b/a = 1.400
   ! file both shipped Q tables were computed from. lam_min shorter
   ! than that file's own shortest wavelength (1.00003e-4 um = 12.4 keV) is
   ! refused rather than served with a frozen refractive index.
   !
   ! WHERE THE DATA IS. build_dust's data_dir is the data root for the whole
   ! library while the build runs: the dielectric functions, the default
   ! extinction curves and the stored cross-section tables all resolve inside
   ! it. It used to resolve the model directories only, the dielectric
   ! functions being compile-time paths relative to the WORKING directory, so a
   ! host that pointed data_dir elsewhere separated a model from the optical
   ! constants its optics were computed on -- and the open had no iostat, so
   ! the process died before any status could be set. Both are fixed: a host
   ! may call build_dust from any directory with an absolute data_dir, and a
   ! path that cannot be read arrives as a status, never as an abort. (A path
   ! the caller names itself is used as given and never composed with the root,
   ! so an absolute one stays absolute.)
   !
   ! ONE MODEL AT A TIME, FOR EMISSION. dust_emission solves through
   ! sed_grain_loop, which reads the module-level grids (lam, NLAM, NT,
   ! T_first) that the LAST build filled, so it can only answer for the most
   ! recently built model. Passing an older one returns status 3 (or stops,
   ! when status is omitted) rather than that other model's numbers, which is
   ! what used to happen in silence. dust_build_table inherits this, being
   ! built from dust_emission solves.
   !
   ! dust_extinction and size_integrated_extinction are NOT restricted: both
   ! read only the model argument, so a host may hold several models and query
   ! their optics in any order. (size_integrated_extinction's own header
   ! claimed the restriction for a while; it never applied to it.)
   !
   ! Two usage modes:
   !   (a) single-cell EXACT solve: dust_emission(m, J_lam, ...)  -- arbitrary J(lambda).
   !   (b) precomputed TABLE + interpolation, for when the cell field is a
   !       fixed reference SHAPE scaled by an intensity U:
   !         call dust_build_table(m, J_ref, U_grid, tab)        ! once
   !         call dust_emission_interp(tab, U, lamI_total [, lamI_chan]) ! per cell
   !
   ! NONLINEARITY CAVEAT: stochastic heating makes the emission a NONLINEAR
   ! functional of the full J(lambda), not of a scalar. The table is valid
   ! ONLY when the cell field is  U * (fixed J_ref shape).  For cells whose
   ! field SHAPE departs from J_ref (e.g. hardened spectra near hot stars),
   ! use the single-cell exact dust_emission instead. The table is built from
   ! exact solves, so it reproduces the U grid points to round-off, except that
   ! dust_emission_interp floors the stored emissivities at 1e-300 before the
   ! log, so a node whose true value is 0 comes back as 1e-300, not 0;
   ! interpolation between grid points is the usual smooth-in-U approximation.
   !
   ! dust_emission takes an optional final argument, status (integer):
   !   call dust_emission(m, J_lam, lamI_total [, lamI_chan] [, status])
   ! On return status = 0 means success, 1 an unknown m%stoch_method, and 2 a
   ! 'qm' model whose populations are missing their radii. When status is
   ! present a bad model is reported through it instead of stopping the
   ! process; when it is omitted such a model stops the run, as before.
   !
   ! EXTINCTION. dust_extinction is the extinction counterpart of
   ! dust_emission, so an RT host takes its opacity from the same model object
   ! and on the same wavelength grid (m%lam) as its emission:
   !   call dust_extinction(m, Cext, Cabs, Csca [, gbar] [, albedo] [, status])
   ! All three required outputs are (m%NLAM) cross sections per H atom
   ! [cm^2/H], with Cext = Cabs + Csca. Optional gbar is the
   ! scattering-weighted asymmetry <cos>; optional albedo is Csca/Cext, 0 where
   ! Cext underflows.
   !
   ! What it returns is the PRECOMPUTED size-integrated curve the model was
   ! built with -- /kext of data/<model>/sedust_<model>.h5, or a kext_*.dat
   ! table -- interpolated onto m%lam: log(lambda)-log(C) for the cross
   ! sections, linear in log(lambda) for <cos>, and the table value itself
   ! wherever a model wavelength coincides with a table wavelength.
   ! Nothing is integrated over grain size at call time. The size integral is
   ! size_integrated_extinction, which is what wrote those products and is what
   ! the standalone calc_kext.x calls; it is re-exported here for a host that
   ! wants to redo the integral rather than read the product.
   !
   ! The BUILDER loads the table, not dust_extinction, because a host is free
   ! to change directory once the model is built. build_dust and each builder
   ! take an optional kext_path:
   !   call build_dust(m, 'astrodust', '../data', &
   !                   kext_path='../data/astrodust/sedust_astrodust.h5')
   ! A file named that way is a contract: if it cannot be read the build FAILS.
   ! Omitting it falls to the model's default, which is the HDF5 product --
   !   astrodust   ../data/astrodust/sedust_astrodust.h5
   !   dl07        ../data/dl07/sedust_dl07.h5
   !   mrn         ../data/mrn/sedust_mrn.h5
   !   zubko       ../data/zubko/sedust_zubko.h5
   !   from_files  none
   ! -- with the model's EUV text table behind it for a tree built without
   ! HDF5, or one whose product carries no /kext. A default that cannot be read
   ! either way leaves the model with no extinction to serve (dust_extinction
   ! then reports status 2) but does NOT fail the build, so a driver that only
   ! wants emission still runs. That holds THROUGH build_dust as well: it
   ! forwards kext_path only when the caller gave one. It used to pass the
   ! product path unconditionally, which made the soft default hard -- a tree
   ! without the curve could then not build a model for emission alone -- and
   ! it applied a text fallback to two models out of three. One route now, for
   ! every model.
   !
   ! include_euv AND lam_min MEAN ONE THING EACH, FOR EVERY MODEL.
   ! include_euv selects the view: .false. cuts the model's own axis at
   ! i_lyman, the LAST node at or below 0.0912 um, so the grid a host gets
   ! COVERS the Lyman limit whatever model it named. lam_min states the
   ! shortest wavelength the model must COVER, and never truncates: astrodust,
   ! DL07 and MRN meet it by extending on the dielectric function their optics
   ! come from, zubko and from_files meet it when their own tables already reach and
   ! refuse it (status 4 out of build_dust) when they do not.
   !
   ! Both were previously model-dependent, and a host with one call and one
   ! floor met the difference: include_euv was an index cut for two models and
   ! a wavelength cut for the third, whose grid then began 1.2e-5 of itself
   ! INSIDE the Lyman limit -- past dust_extinction's edge allowance, so the
   ! same call served for two models was refused for the third. lam_min
   ! extended two models and truncated the other two.
   !
   ! The HDF5 product holds ONE axis and the index i_lyman at which it crosses
   ! the Lyman limit, so include_euv picks the view rather than the file. The
   ! text defaults are the EUV products for the same reason: their wavelength
   ! range CONTAINS the plain grid's, the two being the same size integral on
   ! the same optics with nodes that coincide over the overlap. One default
   ! therefore covers a host that transports ionizing radiation and one that
   ! does not. It is the NARROWER table that has to be named explicitly -- for
   ! astrodust 1129 nodes
   ! from 0.0912 um against 1762 from 1.0e-4 um, for DL07 the same 1129 against
   ! 1823 from 6.205e-5 um, and for zubko 865 from 0.0912 um against the 1201 of
   ! its own optics tables, from 1.0e-3 um. A model grid reaching outside the
   ! table it was given is refused (status 3), never extrapolated.
   !
   ! Every model carries its scattering optics, so every one writes a physical
   ! albedo into its table:
   !   astrodust  Csca and <cos> from the random-orientation T-matrix Q table,
   !              plus the graphite fraction xi_gra(a) of HD23 eq. 15 for the
   !              PAHs, which is the only part of a PAH that scatters
   !   dl07       Mie on the D03 astrosilicate and graphite dielectric
   !              functions (q_silicate_full / q_graphite_full); the
   !              carbonaceous grains scatter as the graphite sphere
   !   mrn        the same two calculations, on this model's own radius grid
   !   zubko      Q_sca and <cos> columns of the model's own ZDA optics tables
   !   from_files same, for whatever tables the descriptor names
   ! A population that genuinely does not scatter leaves its optics unallocated
   ! and contributes zero to the size integral. dust_extinction's status is 0
   ! on success, 1 if an output array is not of size m%NLAM, 2 if no table was
   ! loaded, and 3 if m%lam reaches outside the table; when status is omitted
   ! such a call stops the run. dust_emission adds status 3 = m is not the most
   ! recently built model.  size_integrated_extinction reads no table and is not
   ! restricted to the last-built model, so it uses status 1 alone.
   !
   ! DUST MASS. The cross sections above are per H nucleon. A host that carries
   ! a dust mass density instead converts with the model's own dust mass per H:
   !   Mdust_H = dust_mass_per_H(m)                ! [g/H]
   !   K_abs   = Cabs / Mdust_H                    ! [cm^2/g]
   ! It is one wavelength-independent constant, summed from the model's size
   ! distribution and the solid density of each population's material,
   ! M_dust/N_H = sum_pop rho_bulk * sum_a (4/3) pi a^3 dn(a). Every builder
   ! states those densities, so it is defined for every model in the tree:
   !   astrodust  rho_Ad = 2.74, rho_PAH = 2.0 g/cm^3 (HD23)
   !   dl07       astrosilicate 3.5, graphitic carbon 2.2 g/cm^3 (Draine & Li
   !              2007 sec. 2)
   !   mrn        the same two
   !   zubko      each component's density as its own ZDA optics file states it
   !   from_files the descriptor's rho, or the optics file's when that is 0
   ! It is the same number calc_kext.x normalizes the trailing K_abs column of
   ! every data/kext_*.dat by.
   !
   ! dust_build_table and dust_emission_interp take the same optional final
   ! status argument (0 = success); when present a bad argument is reported
   ! through it instead of stopping the process; when omitted such a call stops
   ! the run.
   !   dust_build_table:      1 size(J_ref) /= m%NLAM
   !                          2 size(U_grid) < 2
   !                          3 U_grid not positive-and-strictly-increasing
   !   dust_emission_interp:  1 U <= 0
   !                          2 size(lamI_total) /= tab%NLAM
   !                          3 lamI_chan present but not (tab%NLAM, tab%n_channel)
   !
   ! build_dust HAS ONE STATUS VOCABULARY, whatever the model, plus an optional
   ! `message` carrying the reason in words for a host that must print one line
   ! before a collective abort:
   !   0   built                      5   extinction table
   !   1   optics table               6   calorimetry
   !   2   size distribution          7   grid inconsistent between components
   !   3   dielectric function        8   EUV spheroid optics unavailable
   !   4   lam_min not coverable      9   model definition (config / descriptor)
   !   90  model name not one of those it codes
   !   91  from_files without config_path
   !   92  zubko_optics not one of 'zda' | 'mie_d03'
   ! The builders below it keep their OWN numbering, which is not the same
   ! one -- an unreadable extinction table is 5 for astrodust and DL07, 3 for
   ! MRN, 6 for zubko and 9 for from_files, and 6 means something else again for
   ! astrodust. A host calling them directly reads the list below; a host on
   ! build_dust reads the list above and does not branch on the model name.
   !
   ! The model builders take the same optional final argument, status (integer,
   ! 0 = success). When present, a missing or malformed input file is reported
   ! through it and the model is NOT built (so an RT host can recover); when
   ! omitted such a failure stops the process, which the CLI drivers rely on.
   ! The exact non-zero code only distinguishes the failing stage; the contract
   ! is simply "0 = built, non-zero = build failed". Codes for each builder:
   !   build_astrodust / build_dl07:  1 Q-table load failed
   !                                  2 size-distribution load failed
   !                                  5 an explicitly named kext_path failed to
   !                                    load
   !   build_dl07 only:               8 pah_xsec is not 'dl07' | 'ld01'.  Only
   !                                    reachable through build_dl07 itself;
   !                                    build_dust does not forward pah_xsec.
   !   build_astrodust only, and only when lam_min asks for the EUV band:
   !                                  3 astrodust dielectric function load failed
   !                                  4 lam_min below that dielectric function's
   !                                    shortest wavelength
   !                                  6 euv_tmatrix = .true. but the spheroid
   !                                    optics of that band are not available:
   !                                    no euv_band_optics_i is registered (a
   !                                    library built without the T-matrix), or
   !                                    the registered one could not compute
   !                                    the band
   !   build_mrn:     1 wavelength axis read failed
   !                  2 lam_min below the D03 dielectric functions' own
   !                    shortest wavelength (EUV band only)
   !                  3 an explicitly named kext_path failed to load
   !   build_zubko:   1 config read failed        2 fewer than 3 components
   !                  3 a component's optics read  4 grid inconsistent
   !                  5 a component's calorimetry read failed
   !                  6 an explicitly named kext_path failed to load
   !                  7 lam_min shorter than this model's own optics table
   !   build_from_files: 1 descriptor open   2 too many pop: lines
   !                     3 invalid channel   4 no pop: lines
   !                     5 optics read       6 grid inconsistent
   !                     7 size-dist read    8 calorimetry read failed
   !                     9 an explicitly named kext_path failed to load
   !                     10 lam_min shorter than the model's own optics tables
   !
   ! LINKING. One archive, and one .mod search path:
   !   cd sed && make libsedust.a          # or ./build_lib.sh -> sed/lib/
   !   gfortran -I<sed> my_rt.f90 -L<sed> -lsedust -fopenmp -o my_rt.x
   ! Nothing in this library reaches the T-matrix, so ../tmatrix need not even
   ! be built.
   !
   ! A host that wants the astrodust EUV band on the spheroid builds the
   ! T-matrix into the SAME archive, so the link line does not change:
   !   cd tmatrix && make libtmatrix.a
   !   cd sed     && WITH_TMATRIX=1 make libsedust.a
   !   gfortran -I<sed> my_rt.f90 -L<sed> -lsedust -fopenmp -o my_rt.x
   ! and calls use_tmatrix_euv_band_optics() once, as above. libtmatrix.a
   ! itself needs no OpenMP runtime; the threading is on the SEDust side.
   !
   ! The validated solver core (sed_grain_loop & helpers in sed_astrodust_mod)
   ! is untouched; this module only re-exports the model API and adds the
   ! table/interpolation layer.
   use constants,         only: wp
   use sed_paths,         only: sed_set_data_root, sed_get_data_root
   use sed_mathlib,           only: locate
   use sed_astrodust_mod, only: dust_model_t, build_dust, &
                                build_astrodust, build_dl07, build_mrn, &
                                build_zubko, build_from_files, &
                                dust_emission, dust_emission_single_teq, &
                                dust_extinction, size_integrated_extinction, &
                                dust_mass_per_H, euv_band_optics_i, &
                                sed_register_euv_band_optics, &
                                sed_forget_euv_band_optics
   implicit none
   private

   ! Re-exported model API.  build_dust is the one entry point: it takes the
   ! model by name and a data directory, and include_euv decides whether the
   ! grid carries the ionizing band.  The builders below it stay exported
   ! so that a host naming its own files keeps working unchanged.
   public :: dust_model_t, build_dust
   public :: build_astrodust, build_dl07, build_mrn, build_zubko, build_from_files
   ! Where the library looks for the data it was not handed a path to.
   ! build_dust sets this from its own data_dir and restores it, so a host on
   ! that entry point never calls these.  A host that calls the builders
   ! directly, and keeps its data somewhere other than <workdir>/../data, sets
   ! the root once before the first build.
   public :: sed_set_data_root, sed_get_data_root
   public :: dust_emission, dust_emission_single_teq, dust_extinction
   public :: size_integrated_extinction, dust_mass_per_H
   ! Optics of the astrodust EUV band; see the EXTREME ULTRAVIOLET note above.
   public :: euv_band_optics_i
   public :: sed_register_euv_band_optics, sed_forget_euv_band_optics
   ! Table API
   public :: dust_emis_table_t, dust_build_table, dust_emission_interp, dust_free_table
   ! Convenience accessors
   public :: dust_nlam, dust_lambda, dust_n_channel, dust_channel_name

   type :: dust_emis_table_t
      integer               :: NLAM = 0, n_channel = 0, NU = 0
      real(wp), allocatable :: U(:)            ! (NU)   intensity-scaling grid
      real(wp), allocatable :: logU(:)         ! (NU)   log(U), cached for interp
      real(wp), allocatable :: lam(:)          ! (NLAM) [um]
      real(wp), allocatable :: J_ref(:)        ! (NLAM) reference field shape (U=1)
      real(wp), allocatable :: total(:,:)      ! (NLAM, NU)
      real(wp), allocatable :: chan(:,:,:)     ! (NLAM, n_channel, NU)
   end type dust_emis_table_t

contains

   ! --- accessors -------------------------------------------------------
   pure integer function dust_nlam(m)
      type(dust_model_t), intent(in) :: m
      dust_nlam = m%NLAM
   end function dust_nlam

   pure integer function dust_n_channel(m)
      type(dust_model_t), intent(in) :: m
      dust_n_channel = m%n_channel
   end function dust_n_channel

   function dust_lambda(m) result(lam)
      type(dust_model_t), intent(in) :: m
      real(wp), allocatable :: lam(:)
      lam = m%lam
   end function dust_lambda

   function dust_channel_name(m, ic) result(name)
      type(dust_model_t), intent(in) :: m
      integer,            intent(in) :: ic
      character(len=16) :: name
      name = m%channel_name(ic)
   end function dust_channel_name

   ! --- emission table over an intensity grid ---------------------------
   subroutine dust_build_table(m, J_ref, U_grid, tab, status)
      ! Precompute lamI(lambda) for J = U*J_ref at each U in U_grid (which must
      ! be positive and strictly ascending). m must be the active model.
      type(dust_model_t),      intent(in)  :: m
      real(wp),                intent(in)  :: J_ref(:)    ! (NLAM) shape at U=1
      real(wp),                intent(in)  :: U_grid(:)   ! (NU)
      type(dust_emis_table_t), intent(out) :: tab
      ! Optional error report (0 = success); see the module header for codes.
      integer, optional,       intent(out) :: status
      real(wp), allocatable :: total(:), chan(:,:)
      integer :: iu, nu

      if (present(status)) status = 0
      nu = size(U_grid)

      if (size(J_ref) /= m%NLAM) then
         if (present(status)) then
            status = 1;  return
         else
            write(*,'(a,i0,a,i0)') 'dust_build_table: size(J_ref)=', size(J_ref), &
                                    ' /= m%NLAM=', m%NLAM
            stop 1
         end if
      end if
      if (nu < 2) then
         if (present(status)) then
            status = 2;  return
         else
            write(*,'(a,i0)') 'dust_build_table: need size(U_grid) >= 2, got ', nu
            stop 1
         end if
      end if
      if (U_grid(1) <= 0.0_wp .or. any(U_grid(2:nu) <= U_grid(1:nu-1))) then
         if (present(status)) then
            status = 3;  return
         else
            write(*,'(a)') 'dust_build_table: U_grid must be positive and strictly increasing'
            stop 1
         end if
      end if

      call dust_free_table(tab)
      tab%NLAM = m%NLAM;  tab%n_channel = m%n_channel;  tab%NU = size(U_grid)
      allocate(tab%U(tab%NU), tab%logU(tab%NU), tab%lam(tab%NLAM), tab%J_ref(tab%NLAM))
      allocate(tab%total(tab%NLAM, tab%NU))
      allocate(tab%chan(tab%NLAM, tab%n_channel, tab%NU))
      tab%U = U_grid;  tab%lam = m%lam;  tab%J_ref = J_ref
      tab%logU = log(tab%U)   ! cached once; the U bracket is grid-fixed at interp

      allocate(total(m%NLAM), chan(m%NLAM, m%n_channel))
      do iu = 1, tab%NU
         call dust_emission(m, U_grid(iu)*J_ref, total, chan)
         tab%total(:, iu)  = total
         tab%chan(:, :, iu) = chan
      end do
      deallocate(total, chan)
   end subroutine dust_build_table

   subroutine dust_emission_interp(tab, U, lamI_total, lamI_chan, status)
      ! Log-log interpolate the table at intensity U (per wavelength and
      ! channel). Reproduces the U grid points to round-off, except that the
      ! stored emissivities are floored at 1e-300 before the log, so a node
      ! whose true value is 0 comes back as 1e-300, not 0; clamps to the grid
      ! ends outside [U(1), U(NU)].
      type(dust_emis_table_t), intent(in)  :: tab
      real(wp),                intent(in)  :: U
      real(wp),                intent(out) :: lamI_total(:)      ! (NLAM)
      real(wp), optional,      intent(out) :: lamI_chan(:,:)     ! (NLAM, n_channel)
      ! Optional error report (0 = success); see the module header for codes.
      ! Validation is scalar/size-only to keep this on the hot path.
      integer, optional,       intent(out) :: status
      real(wp) :: ly(tab%NU)
      real(wp) :: lr, lUq
      integer  :: k, c, jlo

      if (present(status)) status = 0
      if (U <= 0.0_wp) then
         if (present(status)) then
            status = 1;  return
         else
            write(*,'(a,es12.4)') 'dust_emission_interp: need U > 0, got ', U
            stop 1
         end if
      end if
      if (size(lamI_total) /= tab%NLAM) then
         if (present(status)) then
            status = 2;  return
         else
            write(*,'(a,i0,a,i0)') 'dust_emission_interp: size(lamI_total)=', &
                                    size(lamI_total), ' /= tab%NLAM=', tab%NLAM
            stop 1
         end if
      end if
      if (present(lamI_chan)) then
         if (size(lamI_chan,1) /= tab%NLAM .or. size(lamI_chan,2) /= tab%n_channel) then
            if (present(status)) then
               status = 3;  return
            else
               write(*,'(a)') 'dust_emission_interp: lamI_chan must be (tab%NLAM, tab%n_channel)'
               stop 1
            end if
         end if
      end if

      ! tab%logU (= log(tab%U)) is cached at build time and is strictly
      ! ascending because U_grid is positive and strictly increasing. The U
      ! bracket is the same for every wavelength and channel, so locate it once
      ! and reuse jlo below. The interpolation arithmetic is the ascending
      ! branch of interp1, kept term-for-term (log-log, clamped to the ends).
      lUq = log(U)
      call locate(tab%logU, lUq, jlo)
      do k = 1, tab%NLAM
         ly = log(max(tab%total(k, :), 1.0e-300_wp))
         if (jlo == 0) then
            lr = ly(1)
         else if (jlo == tab%NU) then
            lr = ly(tab%NU)
         else
            lr = ly(jlo) + (ly(jlo+1)-ly(jlo))*(lUq-tab%logU(jlo))/(tab%logU(jlo+1)-tab%logU(jlo))
         end if
         lamI_total(k) = exp(lr)
      end do
      if (present(lamI_chan)) then
         do c = 1, tab%n_channel
            do k = 1, tab%NLAM
               ly = log(max(tab%chan(k, c, :), 1.0e-300_wp))
               if (jlo == 0) then
                  lr = ly(1)
               else if (jlo == tab%NU) then
                  lr = ly(tab%NU)
               else
                  lr = ly(jlo) + (ly(jlo+1)-ly(jlo))*(lUq-tab%logU(jlo))/(tab%logU(jlo+1)-tab%logU(jlo))
               end if
               lamI_chan(k, c) = exp(lr)
            end do
         end do
      end if
   end subroutine dust_emission_interp

   subroutine dust_free_table(tab)
      type(dust_emis_table_t), intent(inout) :: tab
      if (allocated(tab%U))     deallocate(tab%U)
      if (allocated(tab%logU))  deallocate(tab%logU)
      if (allocated(tab%lam))   deallocate(tab%lam)
      if (allocated(tab%J_ref)) deallocate(tab%J_ref)
      if (allocated(tab%total)) deallocate(tab%total)
      if (allocated(tab%chan))  deallocate(tab%chan)
      tab%NLAM = 0;  tab%n_channel = 0;  tab%NU = 0
   end subroutine dust_free_table

end module dust_lib
