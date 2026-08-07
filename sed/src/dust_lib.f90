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
   ! dust_emission also takes an optional lamI_pol (same shape and units as
   ! lamI_total):
   !   call dust_emission(m, J_lam, lamI_total [, lamI_chan] [, status] [, lamI_pol])
   ! It returns the INTRINSIC polarized emission -- what a population of
   ! perfectly aligned grains seen with their symmetry axes in the plane of
   ! the sky would radiate -- summed over the populations that carry both
   ! Cpol and falign. Populations without polarized optics (PAHs, and any
   ! model built without them) contribute zero. The geometric sin^2(gamma)
   ! projection and any turbulent depolarization are left to the radiative
   ! transfer and are NOT applied here. Being optional, it leaves every
   ! existing caller valid.
   !
   ! EXTINCTION. dust_extinction is the extinction counterpart of
   ! dust_emission, so an RT host takes its opacity from the same model object
   ! and on the same wavelength grid (m%lam) as its emission:
   !   call dust_extinction(m, Cext, Cabs, Csca [, gbar] [, Cpol_ext] &
   !                        [, Cbir_ext] [, albedo] [, status])
   ! All three required outputs are (m%NLAM) cross sections per H atom
   ! [cm^2/H], integrated over the size distribution of every population:
   ! Cext = Cabs + Csca. Optional gbar is the scattering asymmetry <cos>;
   ! optional albedo is Csca/Cext, 0 where Cext underflows -- derived here
   ! rather than by the caller so every host gets the same convention at those
   ! wavelengths. Optional Cpol_ext is the dichroic (polarized) extinction,
   ! sum dn*Cpol_ext*f_align [cm^2/H], and Cbir_ext the birefringent extinction
   ! from the same weighting; the size integral and the alignment weight are
   ! done here, but the sin^2(gamma) geometry factor and any turbulent
   ! depolarization are left to the radiative transfer, exactly as for
   ! lamI_pol. Populations without scattering or polarized optics (the PAHs)
   ! contribute zero to those terms and enter through absorption only.
   !
   ! WHERE THE NUMBERS COME FROM. The four SCALAR outputs -- Cext, Cabs, Csca,
   ! gbar -- are READ FROM A TABLE, one of the data/kext_*.dat products of
   ! calc_kext.x, attached to the model at build time and interpolated onto
   ! m%lam (log-log in the cross sections, linear in log(lambda) for the signed
   ! gbar). The transport optics of a dust model do not depend on the
   ! transport, so there is nothing to gain by re-running the size integral
   ! during an RT run: the table is that same integral, done once and recorded.
   ! A wavelength of m%lam that coincides with a table node takes the tabulated
   ! value unchanged, and nothing is extrapolated -- a grid running outside the
   ! table is refused (status 3), not served a frozen boundary value.
   !
   ! The two POLARIZED outputs -- Cpol_ext and Cbir_ext -- are NOT tabulated.
   ! They are computed from the model's own orientation-resolved optics on
   ! every call, because their f_align(a) weight is runtime state an RT host
   ! resets cell by cell through dust_set_alignment, which a table fixed at
   ! build time could not follow. This asymmetry is deliberate.
   !
   ! The first-principles size integral is still reachable under its own name,
   !   call size_integrated_extinction(m, Cext, Cabs, Csca [, gbar] &
   !                        [, Cpol_ext] [, Cbir_ext] [, albedo] [, status])
   ! with the identical argument list. It needs no table, computes all six
   ! outputs from the model's optics, and is what the standalone calculators
   ! (calc_kext.x, which writes the tables) use.
   !
   ! DUST MASS PER H. Both routines return cross sections PER H NUCLEON. A host
   ! that carries a dust mass density instead converts them with
   !   Mdust_H = dust_mass_per_H(m)          ! [g/H], one number per model
   !   kappa   = Cabs / Mdust_H              ! [cm^2 per gram of dust]
   ! It is the model's own size distribution weighted by the solid density of
   ! each population's material,
   !   M_dust/N_H = sum_pop rho_bulk * sum_a (4/3) pi a_cm^3 dn_pop(a),
   ! so the opacity it produces refers to the same grains the emission does.
   ! Being a property of the model it is wavelength- and temperature-
   ! independent; evaluate it once. The densities are the astrodust model's
   ! rho_Ad = 2.74 and rho_PAH = 2.0 g/cm^3, the DL07 model's astrosilicate 3.5
   ! and graphitic carbon 2.2 g/cm^3 (both Draine & Li 2007 sec. 2), and, for
   ! the Zubko and file-defined models, the density each optics file declares.
   ! The same number normalizes the K_abs column of every /kext product.
   !
   ! WHICH TABLE. build_dust and each builder take an optional kext_path naming
   ! the file. Omitting it takes that model's default, which is the model's
   ! HDF5 product:
   !   astrodust  ../data/astrodust/sedust_astrodust.h5
   !   dl07       ../data/dl07/sedust_dl07.h5
   !   zubko      ../data/zubko/sedust_zubko.h5
   !   from_files (none -- a file-defined model's product is named after the
   !               model, so a host wanting extinction must name the file)
   ! Behind each sits that model's EUV text table, for a tree built without
   ! HDF5 or a product carrying no /kext.  That order holds THROUGH build_dust
   ! as well: it forwards kext_path only when the caller gave one.  It used to
   ! pass the product path unconditionally, which made the soft default hard --
   ! a tree without the curve could then not build a model for emission alone --
   ! and it applied a text fallback to two models out of three.  One route now,
   ! for all four.
   !
   ! include_euv AND lam_min MEAN ONE THING EACH, FOR EVERY MODEL.
   ! include_euv selects the view: .false. cuts the model's own axis at
   ! i_lyman, the LAST node at or below 0.0912 um, so the grid a host gets
   ! COVERS the Lyman limit whatever model it named.  lam_min states the
   ! shortest wavelength the model must COVER, and never truncates: astrodust
   ! and DL07 meet it by extending on the dielectric function their optics come
   ! from, zubko and from_files meet it when their own tables already reach and
   ! refuse it (status 4 out of build_dust) when they do not.
   !
   ! Both were previously model-dependent, and a host with one call and one
   ! floor met the difference: include_euv was an index cut for two models and
   ! a wavelength cut for the third, whose grid then began 1.2e-5 of itself
   ! INSIDE the Lyman limit -- past dust_extinction's edge allowance, so the
   ! same call served for two models was refused for the third.  lam_min
   ! extended two models and truncated the other two.  The zubko non-ionizing
   ! view is therefore 866 nodes from 0.0899843 um, not 865 from 0.0912011.
   !
   ! The HDF5 product holds ONE wavelength axis and the index i_lyman at which
   ! it crosses the Lyman limit, so include_euv picks the view, not the file.
   ! The text tables behind them are the EUV ones for the same reason: each is
   ! the WIDEST grid its model has, so one file serves a host that transports
   ! ionizing radiation and a host that does not -- the latter's grid is a
   ! subset of the former's, on the same nodes. Every model has two: astrodust
   ! 1129 nodes from 0.0912 um against 1762 from 1.0e-4 um, DL07 the same 1129
   ! against 1823 from 6.205e-5 um, and zubko 865 from 0.0912 um against the
   ! 866 from 0.0899843 um against the 1201 of its own optics tables, from
   ! 1.0e-3 um. Both astrodust products
   ! carry the dichroic eighth column; the narrower text one also keeps the
   ! original narrower field widths. Naming a kext_path is how you point at a
   ! product of your own. A kext_path that cannot be read fails the build (see
   ! the codes below); a
   ! DEFAULT that cannot be read does not, since the emission-only drivers must
   ! still build and calc_kext.x builds a model precisely to write the table
   ! that is not there yet. dust_extinction then reports status 2.
   !
   ! The table is read at BUILD time, not at query time, because a host is free
   ! to change directory once the model is built.
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
   ! sed_grain_loop, which reads the module-level grids that the LAST build
   ! filled, so it can only answer for the most recently built model. Passing
   ! an older one returns status 3 (or stops, when status is omitted) rather
   ! than that other model's numbers, which is what used to happen in silence.
   ! dust_build_table inherits this, being built from dust_emission solves.
   ! dust_extinction and size_integrated_extinction are NOT restricted: both
   ! read only the model argument, so a host may hold several models and query
   ! their optics in any order. (size_integrated_extinction's own header
   ! claimed the restriction for a while; it never applied to it.)
   !
   ! status is 0 on success, 1 if an output array is not of size m%NLAM, 2 if
   ! no extinction table was loaded for this model, and 3 if m%lam runs outside
   ! the table; when it is omitted such a call stops the run.
   !
   ! All four builders carry scattering optics, so albedo and gbar are physical
   ! for every model: astrodust from the T-matrix Q table at every wavelength
   ! the model has, the ionizing band included, DL07 from Mie on the D03
   ! silicate and graphite functions, Zubko and file-defined models from the
   ! Q_sca and g columns of their own tables. Only astrodust carries polarized
   ! optics.
   !
   ! GRAIN ALIGNMENT. Both the polarized emission (lamI_pol) and the dichroic
   ! extinction (Cpol_ext) are weighted by an alignment efficiency
   ! f_align(a_eff), which build_astrodust initializes to the Hensley & Draine
   ! (2023) Table 1 fit. A host that wants a cell-dependent alignment state
   ! overrides it, either as that power law with its own parameters
   !     call dust_set_alignment(m, f_max, a_align, alpha_align [, status])
   !     f_align(a) = f_max / (1 + (a_align/a)**alpha_align)
   ! or as an arbitrary tabulated profile, interpolated in log(a) onto each
   ! population's radius grid and clamped at the ends of the table
   !     call dust_set_alignment_profile(m, aeff_in, falign_in [, status])
   ! for prescriptions that power law cannot express (a RAT-derived Rayleigh
   ! reduction factor, the GRADE-POL exponential). Both leave every population
   ! without polarized optics -- the PAHs, which HD23 take to be unaligned --
   ! untouched and contributing zero.
   !
   ! Calling neither leaves the HD23 alignment in place, so an existing host
   ! sees no change. The current state is readable off the model as
   ! m%align_fmax, m%align_a [um], m%align_alpha and m%align_tabulated; the
   ! three scalars describe the loaded efficiency only while align_tabulated
   ! is .false. They are reported, not applied -- assigning to them does not
   ! re-fill f_align; only the two setters do.
   !
   ! Alignment is a size WEIGHT and enters nowhere in the energy balance, so
   ! resetting it does not invalidate any P(T) solution and does not require a
   ! re-solve: it is one function evaluation on the size grid, and the total
   ! unpolarized emission is bit-for-bit unchanged. A single setter feeds both
   ! dust_emission's lamI_pol and dust_extinction's Cpol_ext, which therefore
   ! cannot fall out of step.
   !
   ! Division of labor, as for lamI_pol and Cpol_ext throughout: SEDust does
   ! the size-distribution integral and the alignment weight; the sin^2(gamma)
   ! projection onto the plane of the sky, any turbulent depolarization
   ! F_turb, and the position angle are the host's job.
   !
   ! status codes (0 = success; when omitted a bad argument stops the run):
   !   dust_set_alignment:          1 a_align <= 0
   !                                2 alpha_align <= 0
   !                                3 f_max outside [0, 1]
   !   dust_set_alignment_profile:  1 aeff_in/falign_in size mismatch or < 2 points
   !                                2 aeff_in not positive and strictly increasing
   !                                3 a falign_in value outside [-1, 1]
   ! A tabulated efficiency outside [-1, 1] is rejected rather than clamped:
   ! |f| <= 1 is the physical bound, so exceeding it is a caller error.
   ! NEGATIVE values are accepted on purpose -- a grain in the wrong internal
   ! alignment state has a negative Rayleigh reduction factor, flipping the
   ! polarization direction by 90 degrees, and clamping at zero would delete
   ! that effect.
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
   ! The model builders take the same optional final argument, status (integer,
   ! 0 = success). When present, a missing or malformed input file is reported
   ! through it and the model is NOT built (so an RT host can recover); when
   ! omitted such a failure stops the process, which the CLI drivers rely on.
   ! The exact non-zero code only distinguishes the failing stage; the contract
   ! is simply "0 = built, non-zero = build failed".
   !
   ! build_dust HAS ONE STATUS VOCABULARY, whatever the model, plus an optional
   ! `message` carrying the reason in words for a host that must print one line
   ! before a collective abort:
   !   0   built                      6   calorimetry
   !   1   optics table               7   grid inconsistent between components
   !   2   size distribution          8   EUV spheroid optics unavailable
   !   3   dielectric function        9   model definition (config/descriptor)
   !   4   lam_min not coverable      10  polarized optics
   !   5   extinction table
   !   90  model name not one of the four
   !   91  from_files without config_path
   ! The four builders below keep their OWN numbering, which is not the same
   ! one -- an unreadable extinction table is 10 for astrodust, 5 for DL07, 6
   ! for zubko and 9 for from_files. A host calling them directly reads the
   ! list below; a host on build_dust reads the list above and does not branch
   ! on the model name. Codes for each builder:
   !   build_astrodust / build_dl07:  1 Q-table load failed
   !                                  2 size-distribution load failed
   !                                  6 astrodust dielectric function load
   !                                    failed (EUV band only)
   !                                  7 lam_min below the dielectric function's
   !                                    own shortest wavelength (EUV band only)
   !                                  8 EUV polarized table load failed
   !                                    (build_astrodust only)
   !                                  9 the EUV band runs outside the
   !                                    wavelengths that table covers,
   !                                    0.0124-0.0912 um (build_astrodust only)
   !                                 11 euv_tmatrix = .true. but the spheroid
   !                                    optics of the EUV band are not
   !                                    available: no euv_band_optics_i is
   !                                    registered (a library built without
   !                                    the T-matrix), or the registered one
   !                                    could not compute the band
   !                                    (build_astrodust only)
   !                                  (build_astrodust also forwards sed_init's
   !                                   polarized-optics codes 3, 4, 5)
   !                                  build_dl07 only:
   !                                  5 an explicitly named kext_path failed to
   !                                    load
   !                                  build_astrodust only:
   !                                 10 an explicitly named kext_path failed to
   !                                    load. It is 10 here and 5 everywhere
   !                                    else because this builder's polarized
   !                                    codes already occupy 3, 4, 5, 8 and 9.
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
   ! WAVELENGTH GRID AND THE EUV. m%lam is the grid of the T-matrix Q table the
   ! model was built from, and TWO scalar tables ship side by side:
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
   ! needs nothing further. Below 0.0912 um that table's wavelengths are the
   ! DH21 dielectric function's own energy nodes, so each absorption edge stays
   ! an exact step between adjacent grid points; the one added point that is not
   ! a node, 0.0912*(1-1e-4), resolves the Lyman-limit step of the radiation
   ! FIELD rather than anything in the grain.
   !
   ! build_astrodust and build_dl07 take an optional lam_min [um] that prepends
   ! log-spaced points from lam_min up to the table, no more coarsely than the
   ! table's own spacing at its short end (taken from the table rather than
   ! written down: 0.00794 in dln(lambda) on the _euv axis, 0.01156 on the
   ! other). Omitting it, or passing one the table already covers, prepends
   ! nothing.
   !
   ! WHAT lam_min DOES FOR ASTRODUST DEPENDS ON WHICH TABLE IT WAS GIVEN. On the
   ! 1129-wavelength table a lam_min below 0.0912 um builds a band. On the _euv
   ! table there is nothing left to ask for: the DH21 dielectric function stops
   ! at 1.00003e-4 um, longward of that table's own first wavelength, so every
   ! legal lam_min either falls inside the table -- no band, no work -- or below
   ! the dielectric data, where it is refused rather than served with a frozen
   ! refractive index. DL07 can be extended below either table, to 6.205e-5 um,
   ! because the D03 optical constants reach further than both; its optics are
   ! Mie on those functions at every wavelength, so that extension needs no
   ! T-matrix either.
   !
   ! WHEN a band IS computed (n_euv > 0), the SCALAR astrodust optics there come
   ! from the DH21 dielectric function on the b/a = 1.400 oblate SPHEROID the Q
   ! table itself is made of, so the grain does not change shape at the seam;
   ! Cext, Cabs, Csca, gbar and albedo are complete over the whole grid.
   !
   ! That spheroid is the ONLY thing in this library that needs the T-matrix,
   ! and it is not compiled in by default. The plain libsedust.a therefore has
   ! no T-matrix in it and links without one; asking it for the spheroid
   ! (euv_tmatrix = .true., the default whenever lam_min is given) is REFUSED
   ! with status 11 rather than answered with a sphere. A model built on the
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
   ! A host that would rather have the far cheaper volume-equivalent sphere
   ! asks for it explicitly with euv_tmatrix = .false.
   !
   ! LINKING. One archive, and one .mod search path:
   !   cd sed && make libsedust.a          # or ./build_lib.sh -> sed/lib/
   !   gfortran -I<sed> my_rt.f90 -L<sed> -lsedust -fopenmp -o my_rt.x
   ! Nothing in this library reaches the T-matrix, so ../tmatrix need not even
   ! be built.
   !
   ! A host that DOES want the astrodust EUV band on the spheroid builds the
   ! T-matrix into the SAME archive, so the link line does not change:
   !   cd tmatrix && make libtmatrix.a
   !   cd sed     && WITH_TMATRIX=1 make libsedust.a
   !   gfortran -I<sed> my_rt.f90 -L<sed> -lsedust -fopenmp -o my_rt.x
   ! and calls use_tmatrix_euv_band_optics() once, as above. libtmatrix.a
   ! itself needs no OpenMP runtime; the threading is on the SEDust side.
   !
   ! The POLARIZED optics are NOT extended by default, and this is the one
   ! place where a caller can be misled by a zero. Cpol_ext and Cbir_ext below
   ! 0.0912 um require an EUV companion table computed for the spheroid itself,
   ! from the same dielectric function and the same first-principles core as
   ! the main table -- and that table is NOT distributed. Without it build_Cpol
   ! reports on error_unit that the band is "zero by omission, not by physics"
   ! and leaves it at zero. The true dichroism there is not small: it grows by
   ! a factor 2.9 below the Lyman limit, peaking at 3.6e-3 of C_ext near
   ! 20.6 eV with the sign opposite to the optical band. Generate the table
   ! with tmatrix/driver/run_q_jori.f90 (lam / lamfile mode) and pass the
   ! resulting .dat + .wave pair as qpol_euv_path / qpol_euv_wave_path. It
   ! stops at 0.0124 um (100 eV), and a polarized build with lam_min shortward
   ! of that is refused (status 9) rather than answered with an invalid
   ! large-x limit. See build_Cpol for what the table resolves and what it
   ! leaves at the geometric-optics zero.
   !
   ! The validated solver core (sed_grain_loop & helpers in sed_astrodust_mod)
   ! is untouched; this module only re-exports the model API and adds the
   ! table/interpolation layer.
   use constants,         only: wp
   use sed_paths,         only: sed_set_data_root, sed_get_data_root
   use sed_mathlib,           only: locate
   use sed_astrodust_mod, only: dust_model_t, build_dust, &
                                build_astrodust, build_dl07, build_zubko, build_from_files, &
                                dust_emission, dust_emission_single_teq, &
                                dust_extinction, size_integrated_extinction, &
                                dust_mass_per_H, &
                                dust_has_polarized_optics, &
                                euv_band_optics_i, &
                                sed_register_euv_band_optics, &
                                sed_forget_euv_band_optics, &
                                dust_set_alignment, dust_set_alignment_profile
   ! Aligned-grain polarized scattering optics for a polarized RT host. The
   ! init call (load_scatmat_aligned) runs once from serial code; the query
   ! calls are pure reads, safe from OpenMP threads. The loaded grids and arrays
   ! (scm_*) are exposed read-only so a host can inline its own lookups.
   use scatmat_aligned_mod, only: load_scatmat_aligned, free_scatmat_aligned, &
                                scatmat_band, extinction_matrix_aligned, &
                                mueller_matrix_aligned, mueller_matrix_random, &
                                mueller_matrix_total, scattering_cross_sections, &
                                scm_loaded, scm_nband, scm_nti, scm_nts, scm_nphi, &
                                scm_ntheta, scm_lambda, scm_theta_i, scm_theta_s, &
                                scm_phi, scm_theta_ran, &
                                scm_cos_theta_s, scm_cext_al, scm_cpol_al, scm_cbir_al, &
                                scm_csca_al, scm_csca_pol_al, scm_cext_tot, scm_csca_tot, &
                                scm_cext_ref, scm_csca_ref, scm_F_tot, scm_F_ref, scm_Z, &
                                scm_profile_name, scm_fmax, scm_a_align, scm_alpha, &
                                scm_profile_mismatch, scm_bytes
   implicit none
   private

   ! Re-exported model API
   ! One entry point for every coded model: it takes the model by name and a
   ! data directory, and include_euv decides whether the grid carries the
   ! ionizing band.  The four builders below it stay exported so that a host
   ! naming its own files keeps working unchanged.
   public :: dust_model_t, build_dust
   ! Where the library looks for the data it was not handed a path to.
   ! build_dust sets this from its own data_dir and restores it, so a host on
   ! that entry point never calls these.  A host that calls the four builders
   ! directly, and keeps its data somewhere other than <workdir>/../data, sets
   ! the root once before the first build.
   public :: sed_set_data_root, sed_get_data_root
   public :: build_astrodust, build_dl07, build_zubko, build_from_files
   public :: dust_emission, dust_emission_single_teq, dust_extinction
   public :: size_integrated_extinction
   public :: dust_mass_per_H
   ! Optics of the astrodust EUV band; see the WAVELENGTH GRID AND THE EUV note
   ! above.
   public :: euv_band_optics_i
   public :: sed_register_euv_band_optics, sed_forget_euv_band_optics
   public :: dust_has_polarized_optics
   public :: dust_set_alignment, dust_set_alignment_profile
   ! Re-exported aligned-scattering API (initialization + path queries)
   public :: load_scatmat_aligned, free_scatmat_aligned, scatmat_band, &
             extinction_matrix_aligned, mueller_matrix_aligned, &
             mueller_matrix_random, mueller_matrix_total, scattering_cross_sections
   ! Re-exported read-only storage (public, protected in scatmat_aligned_mod)
   public :: scm_loaded, scm_nband, scm_nti, scm_nts, scm_nphi, scm_ntheta, &
             scm_lambda, scm_theta_i, scm_theta_s, scm_phi, scm_theta_ran, &
             scm_cos_theta_s, scm_cext_al, scm_cpol_al, scm_cbir_al, scm_csca_al, &
             scm_csca_pol_al, scm_cext_tot, scm_csca_tot, scm_cext_ref, scm_csca_ref, &
             scm_F_tot, scm_F_ref, scm_Z, &
             scm_profile_name, scm_fmax, scm_a_align, scm_alpha, &
             scm_profile_mismatch, scm_bytes
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
