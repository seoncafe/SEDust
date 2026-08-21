# SEDust

Version 1.20.

A self-contained library for the optical properties and thermal emission of
interstellar dust: dielectric functions -> T-matrix / Mie cross sections ->
grain temperature distributions -> emergent infrared SED.

SEDust owns the dust physics. A radiative-transfer code that links it is
responsible for photon transport and geometry only: everything about the
grains, from cross sections through temperature distributions to emission and
polarization, is computed here and handed over as wavelength-resolved
quantities per H atom.

SEDust is **model-agnostic**. The HD23 astrodust+PAH model, the Draine & Li
(2007) carbonaceous+silicate model, and the Zubko et al. (2004) BARE-GR-S model
are handled as peers through one derived type (`dust_model_t`) and one emission
call (`dust_emission`). The solver core does not know which model it is running.

The whole package can also be linked into a Fortran 3D radiative-transfer code
as a static library, `libsedust.a`, with a two-step API: initialize once, then
solve one cell at a time. Every model builder and solver call takes an optional
`status` argument, so a missing input file or an invalid model is reported back
to the host instead of stopping the process.

## Which version

SEDust is split into two versions because most radiative-transfer codes either
do not carry polarization at all, or carry it only for spherical grains.

This is **version 1.20**, the polarized branch. Use it for polarized transfer
that accounts for the non-spherical (**spheroidal**) grain shape -- dichroic
extinction, birefringence, and scattering by aligned spheroidal grains -- for
which it provides the direction-dependent extinction matrix `K(theta_i)` and the
aligned-grain scattering matrix, in addition to the scalar optics and thermal
emission.

For transfer that does not carry polarization, use **version 1.00**, the scalar
branch: it computes the same unpolarized cross sections and emission through the
same `dust_emission` / `dust_extinction` API without the polarized optics or
their data tables. (Version 1.20 can also run scalar-only, by building a model
with `load_polarized_optics = .false.`, when a host wants both from one build.)

**1.20 is a superset of 1.00.** Every scalar capability of the scalar branch is
present here and computes the same numbers: the same builders, the same
`dust_emission` / `dust_extinction`, the same `lam_min` extension into the
ionizing band, the same `calc_kext.x` products. The differences are the added
polarized quantities and one consequence of adding them — the `status` codes of
`sed_init` / `build_astrodust` are numbered differently, because codes 3-5 are
taken by the polarized inputs, so what is `3` and `4` in v1.00 is `6` and `7`
here. Code that only checks `status /= 0` is unaffected.

The two versions are branches of one repository. `v1.00` is the default branch,
so this version has to be asked for by name:

```
git clone -b v1.20 git@github.com:seoncafe/SEDust.git # this version, 1.20
git clone git@github.com:seoncafe/SEDust.git          # the scalar version
```

An existing clone switches with `git checkout v1.20`.

## Layout

```
SEDust/
  sed/          the SED solver: cross sections, enthalpy, P(T), emission
    src/        library modules + drivers
    rt_example/ examples of linking libsedust.a into an RT code, one
                minimal and two showing the polarized quantities
  mc/           Draine & Anderson (1985) Monte Carlo solver (independent check)
  tmatrix/      Mishchenko T-matrix engine + driver; writes the Q table
  data/         one directory per dust model, plus what the models share
    astrodust/  everything that model owns: `sedust_astrodust.h5` (the
    dl07/       primary form -- its wavelength axis, cross-section tables
    zubko/      and extinction curve in one file), the same as text
                (`q_*.dat`, `kext_*.dat`), and, where the model IS a set of
                files, its definition (the ZDA config, optics, calorimetry)
    dielectric/ shared material data: the D03 / DH21 / D16 optical
                constants and the PAH cross sections.  A dielectric function
                is not one model's -- DL07 and Zubko read the same D03
                astrosilicate -- so it does not live inside a model directory
    release/    published reference tables (HD23, Draine) and the HD23 size
                distribution
  docs/         technical reports and the library user manual
  pyutil/       small Python helpers (radiation fields, SED from Cabs, and
                sedust_h5 -- the reader for data/<model>/sedust_<model>.h5)
```

Everything the code reads at run time ships with the package. There are no
paths outside this directory.

## Build and run

Requires `gfortran` (OpenMP for the parallel drivers). No autoconf, no
top-level configure; each subdirectory has its own `Makefile`.

```sh
# the T-matrix library.  Only the astrodust EUV band uses it, so it is needed
# by calc_sed.x astrodust, calc_kext.x, and a WITH_TMATRIX library build
cd tmatrix && make libtmatrix.a

# the SED solver.  Every executable is named for the quantity it computes --
# calc_sed.x the emission SED, calc_kext.x the size-integrated transport optics,
# calc_qtable.x the stored cross-section tables, calc_enthalpy.x the enthalpy
# tables -- so that none of them reads as the build tool.
cd ../sed
make                        # calc_enthalpy.x calc_sed.x calc_kext.x
./calc_sed.x astrodust      # astrodust+PAH SED at log U = 0.20 -> output/
./calc_sed.x dl07           # Draine & Li (2007) SED at U = 1   -> output/
./calc_sed.x zubko          # Zubko/ZDA SED at U = 1            -> output/
./calc_enthalpy.x           # astrodust enthalpy tables (add c2 for the
                            #   Stage-1 density-corrected prefactor)

# Every program in sed/ takes ONE vocabulary of run settings
# (sed/src/sed_run_options.f90) and declares which axes it has a referent for; a
# word of an axis a program does not have is refused by name rather than
# ignored.  calc_sed.x declares all of the following.  Any combination is a
# valid run; naming two values of one axis, or a setting the chosen solver does
# not read, is refused rather than silently dropped.  Every setting tags the
# output filename, in a fixed order whatever order the words are typed in, so a
# variant run cannot clobber the production one.
#
#   solver    heuristic (default) | draine | equil | qm | qm_dbcon | qm_stati
#   grid      euv          solve the ionizing band as well
#   field     mathis_orig  literal Mathis 1983 (4000 K dilution 1e-13, CMB 2.9 K)
#             logU=X       scale the Mathis field to U = 10^X
#             hardfield    fill the band below the Lyman limit with a diluted
#                          1e5 K blackbody, so the field carries photons above
#                          13.6 eV (implies euv)
#   emission  induced      multiply by the stimulated-emission factor (1 + J/B)
#             photcut      bound a bin's emission by its own enthalpy (heuristic)
#   qm sizes  nstate=N     enthalpy bins;   nisrf=N   field wavelengths
#   graphite  gra_d03_sphere | gra_d16_sphere | gra_d16_spheroid
#             the graphite of the PAH-to-graphite xi blend (astrodust and dl07,
#             which compute it; zubko reads the ZDA tables)
#   PAH xsec  dl07 (default) | ld01   which published carbonaceous absorption
#             the blend takes (dl07 model; calc_kext.x takes it too)
#   enthalpy  c2           astrodust Stage-1 density-corrected prefactor
#
./calc_sed.x astrodust qm_stati nstate=500  # -> sed_astrodust_qm_stati_ns500_*.dat
./calc_sed.x astrodust gra_d16_spheroid     # PAH xi blend on the D16 b/a=1.4 spheroid
./calc_sed.x astrodust mathis_orig          # -> sed_astrodust_morig_*.dat
./calc_sed.x dl07 lmc2_10 draine euv        # -> sed_dl07_lmc2_10_euv_draine.dat
./calc_sed.x zubko euv hardfield            # -> sed_zubko_euv_hardfield.dat

# the library, for embedding in an RT code.  No T-matrix in it: link with
#   -I. -L. -lsedust
make libsedust.a
# ... and with the astrodust EUV band on the spheroid of the Q table, in the
# same archive and on the same link line (the host then calls
# use_tmatrix_euv_band_optics() once before build_astrodust)
WITH_TMATRIX=1 make libsedust.a
make use_dustlib_scatmat.x  # reference consumer of the aligned-scattering API

# size-integrated transport optics: lambda, albedo, <cos>, C_ext/C_abs/C_sca per H
./calc_kext.x astrodust     # -> ../data/astrodust/kext_astrodust_MW.dat (+ dichroic column)
./calc_kext.x astrodust euv # the same model carried into the ionizing band
./calc_kext.x dl07 euv      # -> ../data/dl07/kext_dl07_MW_euv.dat
./calc_kext.x zubko         # -> ../data/zubko/kext_zubko_BARE_GR_S.dat (cut at the Lyman limit)
./calc_kext.x zubko euv     # -> ../data/zubko/kext_zubko_BARE_GR_S_euv.dat (the whole ZDA range)
./calc_kext.x from_files ../data/zubko/zubko_descriptor.txt

# the optics products.  ORDER MATTERS: calc_qtable.x lays down the wavelength
# axis and replaces data/<model>/sedust_<model>.h5, so it runs first and calc_kext.x
# puts /kext back.  Both also write the text products beside them.
./calc_qtable.x             # all three models -> ../data/<model>/q_*.dat
                            #                 -> ../data/<model>/sedust_<model>.h5
./calc_kext.x astrodust euv # ... then /kext into each of the three files
./calc_polarized_optics.x          # -> /polarized in ../data/astrodust/sedust_astrodust.h5

# polarized extinction alone, checked against the HD23 release
make calc_polext.x         && ./calc_polext.x

# checks against the published reference tables, from the same programs
./calc_kext.x dl07 ld01     # -> ../data/dl07/kext_ld01_MW.dat
#   The DL07 model with the Li & Draine (2001) carbonaceous absorption instead
#   of the Draine & Li (2007) one -- everything else the same, so a ratio of
#   kext_ld01_MW.dat to kext_dl07_MW.dat isolates that one cross section.
#   Draine's 2003 table (../data/release/kext_albedo_WD_MW_3.1_60_D03.all_2003)
#   was computed with the LD01 cross sections, so the ld01 run is the one that
#   should reproduce it, and calc_kext.x prints that comparison decade by
#   decade at the end of every dl07 run.  Measured on the two curves: the
#   carbonaceous C_abs of ld01 is 0.915 of dl07 at 1.047 um, where DL07 added
#   the PAH cation near-infrared absorption of Mattioda et al. (2005), and
#   1.109 at the 6.2 um C-C feature.
./calc_sed.x dl07 mathis_orig      # -> output/sed_dl07_mw31_60_morig.dat
./calc_sed.x zubko euv mathis_orig # -> output/sed_zubko_euv_morig.dat
#   The same models solved in the literal Mathis 1983 field convention
#   (4000 K dilution w_4000 = 1e-13, CMB 2.9 K) beside the default corrected
#   one (1.65e-13, 2.725 K), so that the convention behind each published
#   reference SED is settled by direct comparison.  The convention is fixed
#   BEFORE the model is built, so the CMB temperature in the cooling term and
#   the one in the field always agree.

# the Monte Carlo cross-check
cd ../mc && make && ./main_mc_sed.x run_sed.nml

# regenerating the scalar T-matrix Q tables (optional; both tables ship with SEDust)
cd ../tmatrix && make && ./run_tmatrix.x test   # smoke test -> output/..._euv.test.dat
./run_tmatrix.x                                 # full sweep -> output/..._euv.dat
make lyman_cut                                  # -> the 1129-wavelength companion

# orientation-resolved (polarized) Q table from first principles
./run_q_jori.x test                             # sample + full-sweep time estimate
# full sweep is ~16 h on one core; parallelize over wavelength windows:
#   ./run_q_jori.x range 1 400   (etc.)  then  ./run_q_jori.x merge output/...jw*.dat
# individual wavelengths, for a band the shipped table does not cover:
#   ./run_q_jori.x lam 0.0602 [ja=1:3] [tag=NAME]   ./run_q_jori.x lamfile PATH
#   ./run_q_jori.x lammerge STEM FILE [FILE ...]    (assembles the ja= windows)

# scattering matrix of randomly oriented grains (five optical bands ship with SEDust)
./run_scatmat.x 0.55                            # one wavelength; ./run_scatmat.x all for the grid

# fixed-orientation scattering matrix of ALIGNED grains (also five optical bands)
./run_scatmat_aligned.x test                    # one band, reduced grid: timing + OpenMP check
./run_scatmat_aligned.x                         # default UBVRI bands -> output/ (~8 min, 32 threads)
#   profile=FILE regenerates under a different alignment profile in minutes
```

Outputs are plain ASCII `.dat` files written to each subdirectory's `output/`;
`calc_kext.x` writes into `data/` instead.

## The ionizing band

The HDF5 product holds ONE wavelength axis and the index `i_lyman` at which it
crosses the Lyman limit, so `include_euv` picks the view rather than the file:
`.false.` (the default) returns `lambda(i_lyman:)` and the same rows of every
wavelength-indexed array.  For astrodust that is 1129 wavelengths
(0.0912--39810 um) out of 1762 (1.0e-4--39810 um), for DL07 1129 out of 1823,
for Zubko 866 out of 1201.

The scalar T-matrix text pair behind it, `q_astrodust_P0.20_Fe0.00_1.400.dat`
and `..._euv.dat`, is the same split as two files: the EUV one is the whole
sweep `run_tmatrix.x` writes and the other is that file with every wavelength
shortward of the Lyman limit dropped (`make lyman_cut`), row selection only, so
they agree cell for cell over the 1129 they share.

Below 0.0912 um the `_euv` table's wavelengths are the DH21 dielectric
function's own energy nodes rather than a resampling of them, so every
absorption edge stays an exact step between adjacent wavelengths instead of
being averaged into a ramp by a uniform 200-per-decade axis. Counting the edges
takes a threshold and the answer moves with it: 23 places where k jumps by more
than 0.5% across a pair of nodes closer than 5e-4 in relative energy, 14 of them
closer than 1e-4. The largest jumps are +211% at Fe K (7124 eV), +131% at O K
(538), +67% at Fe L (724) and +30% at C K (291). One point of that band is not a
dielectric node: 0.0912*(1-1e-4), placed to resolve the Lyman-limit step of the
*radiation field* (see below).

`build_astrodust` and `build_dl07` take an optional `lam_min` [um] for a host
needing to reach below the table it was given. It prepends log-spaced points
below the table, and is refused when it asks for a wavelength the model's own
dielectric data does not cover. What that does for astrodust depends on the
table:

- on the default 1129-wavelength table, a `lam_min` below 0.0912 um **does**
  build a band. Its optics come from the T-matrix on the same spheroid when
  `libtmatrix.a` is linked in (`WITH_TMATRIX=1 make libsedust.a`), and are
  refused with status 11 when it is not — never silently replaced by a sphere
  unless `euv_tmatrix = .false.` asks for one;
- on the `_euv` table there is nothing left to ask for: the DH21 dielectric
  function stops at 1.000032e-4 um, longward of that table's first wavelength,
  so every legal `lam_min` falls inside the table and prepends nothing.

DL07 can be extended below either table, to 6.205e-5 um, because the D03 optical
constants reach further than both; its band is Mie on those functions and needs
no T-matrix.

| product | rows | lambda [um] | what sets the floor |
|---|---:|---|---|
| `data/astrodust/kext_astrodust_MW.dat` | 1129 | 0.0912 - 39810 | non-EUV Q-table grid |
| `data/astrodust/kext_astrodust_MW_euv.dat` | 1762 | 1e-4 - 39810 | EUV Q-table grid |
| `data/dl07/kext_dl07_MW.dat` | 1129 | 0.0912 - 39810 | non-EUV Q-table grid |
| `data/dl07/kext_dl07_MW_euv.dat` | 1823 | 6.205e-5 - 39810 | the D03 dielectric functions |
| `data/zubko/kext_zubko_BARE_GR_S.dat` | 866 | 0.08998 - 1e4 | the ZDA optics table, cut at the Lyman limit |
| `data/zubko/kext_zubko_BARE_GR_S_euv.dat` | 1201 | 1e-3 - 1e4 | the ZDA optics table itself |

### Stored Q tables

Every model's `(lambda, a_eff)` cross sections are in `/qtable` of its HDF5
product, which is what ships.  `calc_qtable.x` writes them as text beside it as
well, and those files are what the table below names; they are regenerable, so
they are not tracked:

| model | non-EUV | EUV | a_eff grid |
|---|---|---|---|
| astrodust | `data/astrodust/q_astrodust_P0.20_Fe0.00_1.400.dat` (1129) | `..._euv.dat` (1762) | 169, from `data/astrodust/DH21_aeff` |
| DL07 | `data/dl07/q_dl07_{sil,gra,pah_neu,pah_ion}.dat` (1129) | `..._euv.dat` (1823) | 84, Draine's analytic grid, 3.548e-4 - 5.012 um |
| Zubko | `data/zubko/q_zubko_{sil,gra,pah}.dat` (866) | `..._euv.dat` (1201) | 121 (sil, gra) / 28 (PAH), the ZDA tables' own |

`sed/calc_qtable.x` writes all of these, and the HDF5 product of each model
with them; the astrodust scalar pair is the exception, computed by
`tmatrix/run_tmatrix.x` plus `make lyman_cut` and installed into
`data/astrodust/`. Every non-EUV file is the row subset of its EUV counterpart
that starts at the Lyman limit — verified, not asserted. In the HDF5 product the
pair is one array and its `i_lyman`, so there is no second file to keep in
step (§ *The HDF5 optics products* in the manual).

**These are our own tables.** For DL07 they are what the model already computes
at build time, written down rather than recomputed. For Zubko they are a genuine
recomputation: the model is still built on the ZDA optics tables under
`data/zubko/` by default, while these are what this tree's Mie gives on the
Draine (2003) optical constants (`data/dielectric/index_silD03`,
`index_CpaD03`, `index_CpeD03`). The ZDA headers name Draine's older `eps_Sil`
and `eps_Gra`, but that label does not survive a check: over all 121 x 1201
cells our silicate table reproduces the distributed one to a mean relative
2.4e-6 in Q_abs and 8.9e-6 in Q_sca, with `<cos>` agreeing to 8.5e-3 at worst,
and the agreement holds over the ZDA tables' whole range, out to 1e4 um. The
shipped tables were computed on D03 as well, whatever their header says. Both sets are
in the HDF5 product: `/qtable/{sil,gra,pah}` is the distributed one, which the
model is built on unless `build_zubko`'s `optics` argument asks otherwise, and
`/qtable/{sil,gra,pah}_mie_d03` is this recomputation.

The one text table anything still opens is the astrodust scalar pair. It is the
INPUT `calc_qtable.x` reads to write the HDF5 product, so that program cannot
read it back out of its own output; it is also the route a tree built without
HDF5 falls to. Everything else — the drivers, `libsedust.a`, the `rt_example`
consumers, the Python reader — opens the `.h5`.

No floor is a free choice: each is the shortest wavelength the data the model is
made of actually covers, and a shorter request is refused rather than served
with a refractive index frozen at the table boundary.

The DL07 extension below the table uses Mie on the volume-equivalent sphere; it
agrees with Draine's published table for the same model to 0.056-0.86% in C_ext
per decade. **There is no external reference for the astrodust ionizing band**:
the HD23 release stops at 12.4 eV. See `docs/SEDust_user_manual.pdf` for the full
accounting.

**The polarized optics did not follow the table into that band, deliberately.**
The `kext_*.dat` products above are scalar transport optics, except that the two
astrodust ones carry the dichroic extinction as an eighth column -- no other
model has polarized optics. On `kext_astrodust_MW.dat` every row of it is a
measured value; on `kext_astrodust_MW_euv.dat` the 633 rows below 0.0912 um are
exactly zero, which the header states as a deficit rather than a measurement.
The orientation-resolved tables the polarized quantities are read from stay
on the 1129-wavelength DH21 axis, 0.0912-39810 um, because polarized transfer is
run at the few wavelengths
an observation was made at rather than swept — see [Polarization in the ionizing
band](#polarization-in-the-ionizing-band) for what fills the gap when a run
needs it.

### A step in the field needs a point, like a step in the material

An ISRF is exactly zero below the Lyman limit and finite above it
(`radfield.f90`: `if (lambda(i) < 0.0912d0) J(i) = 0`). The dielectric energy
axis has no node between 13.595 and 14.000 eV, so a grid built from those nodes
alone leaves one cell 2.98% wide straddling that step; a consumer's trapezoid
integral then ramps up from zero across it and absorbs photons the field does
not carry. Measured on the shipped optics, that inflates the absorbed power of
the smallest grains by 1.74% and moves the emergent SED by up to 13%. The extra
point at 0.0912*(1-1e-4) brings it to 0.006%. It is the same move the extension
makes for every absorption edge — resolve a step with a close pair of points —
applied to the one step that belongs to the radiation field rather than to the
grain.

### The single-photon ceiling comes from the field

The stochastic solvers cap the enthalpy window at the energy of the hardest
single photon a grain can absorb, `max(13.6 eV, hc/lambda_hard)`, where
`lambda_hard` is the shortest wavelength at which the incident `J_lambda`
exceeds 1e-12 of its own peak (`hardest_photon_energy` in
`sed/src/radfield.f90`). **That ceiling is a property of the radiation field,
not of the wavelength grid**, and every shipped model now depends on the
distinction. The Q table reaches 1e-4 um (12.4 keV) and the Zubko/ZDA optics
tables 1e-3 um (1.24 keV), while a Mathis field carries nothing below 0.0912 um:
reading the grid instead of the field would raise the single-photon bound by
factors of 912 and 91 respectively, with no photon to justify either. That
coarsens the fixed-count enthalpy bins — and the refinement loops, guarded and
capped at ten iterations on the way down, would not walk it back — shifting the
emergent SED with nothing behind it.
`docs/EUV_EXTENSION_HOST_REGRESSION.md` records the measurement.

Before the Q table was carried into the ionizing band, astrodust and DL07 could
not tell the two apart: their grid floor *was* the Lyman limit (13.595 eV <
13.6 eV, so the constant stayed selected either way), and Zubko was the only
model that exercised the distinction. Now all three do.
`./calc_sed.x zubko [euv]` still prints both candidates side by side, and is the
quickest place to see them disagree.

### The photon that was not there (fixed)

`calc_P` builds the upward transition rates of the stochastic solver. Its
highest-bin term — the rate into the top enthalpy bin, driven by every photon
more energetic than that bin's gap — used to integrate **from `lambda(1)`, the
short end of the optics grid**, and to read `J_lambda` at the gap wavelength
through an interpolator that *clamps* outside its range
(`interp1` in `sed/src/sed_mathlib.f90`). The top bin sits above the
single-photon bound by construction, so its gap wavelength was always shorter
than the grid, and the clamp answered every time with `J_lambda` at the grid
edge instead of zero. On a grid ending at the Lyman limit that is the full
Mathis intensity there, 1.359 in its own units: **every top-bin transition was
driven by a photon flux the field does not carry.** The same clamp made the
integral run backwards when the gap exceeded the grid's range, subtracting a
rate instead of adding none.

Both limits now come from the field — `hardest_photon_energy`, the same routine
the enthalpy ceiling uses — and the band is skipped outright when the gap is
harder than the hardest photon. The sample count follows the band's width
instead of being a flat 51 points.

**This changed all three shipped models**, and for Zubko it has nothing to do
with the astrodust table: that model's grid has always started at 1e-3 um, so it
carried the same defect independently. Measured against the previous products:
astrodust up to 12.1% locally and −1.1% to +0.9% in integrated power, DL07 9.7%
and −2.8%, Zubko 3.0% (6.0% with `euv`). Energy conservation improved, from
+0.236% to −0.057% (astrodust S1, emitted/absorbed per H).

**Any SED produced before this fix includes photons that were not in the field**,
in the top enthalpy bin of every stochastically heated grain. Do not compare
against those products without accounting for it. The polarized emission
`lamI_pol` moves with the total, since both are the same P(T).

## Stochastic heating

Small grains do not reach a steady temperature: each absorbed ultraviolet
photon drives a large transient excursion, so the emission is set by the
temperature probability distribution `P(T)` rather than by a single `T_eq`.
SEDust provides three independent solvers, all reading the same cross sections
and enthalpy, and agreeing with one another to a few percent in every band from
the near-infrared through the sub-millimeter.

| `stoch_method` | solver | cooling treatment |
|---|---|---|
| `'heuristic'` | Guhathakurta & Draine matrix, look-ahead grid narrowing | continuous (**default**) |
| `'draine'`    | Guhathakurta & Draine matrix, iterative refinement | continuous |
| `'qm'`        | energy-space transition matrix, BiCG sparse solve | thermal-discrete (`dbdis`), thermal-continuous (`dbcon`), or exact-statistical (`stati`) |
| `'equil'`     | equilibrium temperature, no stochastic solve | n/a |

The three `'qm'` cooling treatments differ in how the downward rates are
built. `dbdis` (the reference) and `dbcon` weight the emission with a Planck
factor at each bin's temperature, treating the bin as thermal; `dbcon` further
collapses the cooling into a continuous term. `stati` is the exact-statistical
(microcanonical) treatment: the downward rate carries the vibrational
density-of-states ratio `g_f/g_i`, counted exactly by the Beyer-Swinehart
recursion over the grain's mode spectrum, and the emission kernel is
degeneracy-weighted rather than a bin blackbody -- no grain temperature enters
at all. The exact count is feasible only for the smallest carbonaceous grains:
the path is capped at `a <= 25 A`, and in the production configuration it runs
on PAHs with `N_C <= 116` (`a <= 6.3 A`) while every other grain, silicates
included, reverts to `dbdis` within the same run. The total SED it produces
agrees with the thermal-discrete one to 0.02% band-integrated, with the
largest pointwise differences confined below 1.5 um where the SED is orders of
magnitude under its peak. Select it with
`./calc_sed.x astrodust qm_stati [nstate=N]` -- or with the same word on
`./calc_sed.x dl07` and `./calc_sed.x zubko`, whose carbonaceous populations reach the
same solver.

The Monte Carlo solver in `mc/` follows Draine & Anderson (1985) and tracks
`T(t)` through individual absorption events. It shares no code with the matrix
solvers, which is what makes it a useful check on them.

Grains whose equilibrium enthalpy exceeds 150 eV are placed at `T_eq` and skip
the stochastic solve; the gate propagates forward in grain size.

## Polarization

For the astrodust model SEDust also computes the polarized cross sections of
aligned spheroidal grains, from an orientation-resolved spheroid table. By
default this is the Draine & Hensley (2021) table that ships in
`data/dielectric/`. Two quantities are available:

| Quantity | Where |
|---|---|
| polarized emission | optional `lamI_pol` argument of `dust_emission` |
| polarized extinction | optional `Cpol_ext` argument of `dust_extinction`, or the eighth column of `data/astrodust/kext_astrodust_MW.dat` |

Both are intrinsic values: the size integral and the alignment efficiency
`f_align(a)` are already applied, while the viewing geometry (the angle
between the field and the line of sight, and any turbulent depolarization)
is left to the caller. Codes that read only the first seven columns of the
extinction table are unaffected. `dust_extinction` and the table agree to the
precision the file is written with, so a code that links the library can take
its opacity from the call, on the model's own wavelength grid, and skip the
file entirely.

The polarized *extinction* covers the whole model grid: `Cpol_ext` (dichroism),
`Cbir_ext` (birefringence) and `Cpol` (polarized absorption) come from an
orientation-resolved table spanning all 1129 wavelengths from 0.0912 to
39810 um, so nothing has to be regenerated to use them. The default is the HD23
release file `data/astrodust/q_DH21Ad_P0.20_Fe0.00_1.400.dat.gz`; SEDust's own
regenerated table (`tmatrix/output/q_astrodust_jori_*.dat.gz`) is opt-in through
`qpol_path` and is the one carrying the birefringence block, so `Cbir_ext` is
zero with the default and nonzero with it.

The computed polarized extinction reproduces the released
`polarized_extinction.dat` to a median of 0.03%. The polarized emission
fraction reaches 17.2% at 154 um and 19.2% at 850 um.

The orientation-resolved table itself can also be regenerated from the astrodust
dielectric function with SEDust's own T-matrix engine, so the polarized optics
need not be taken from the release file. Each of the three size-parameter regimes
is computed from first principles — the Rayleigh polarizability, the
fixed-orientation amplitude matrix (optical theorem for extinction, phase-matrix
integral for scattering), and a projected-area-plus-Fresnel geometric-optics
limit — and matches the release to a few parts in 10^4 wherever grains carry
weight. Fed through `calc_polext`, the regenerated table reproduces
`polarized_extinction.dat` to a median of 0.06%, computed with no recourse to the
release optics. `tmatrix/run_q_jori.x` writes the table (a drop-in for the
release format); `oriented_cross_sections` is the same computation for a single
grain and wavelength, for an arbitrary point or a shape the table does not cover.
To use the regenerated table in a run, pass its path as `qpol_path` to
`build_astrodust` or `sed_init`; the default stays the release table.

The alignment efficiency `f_align(a)` can be replaced on an existing model
with `dust_set_alignment` (the HD23 power law) or `dust_set_alignment_profile`
(an arbitrary tabulated profile, for a RAT-derived reduction factor). Both are
size weights applied outside the temperature solution, so neither re-solves
`P(T)` and neither changes `lamI_total`.

For randomly oriented grains the full scattering (Mueller) matrix is also
computed, by `tmatrix/run_scatmat.x`, and stored for five optical bands
(approximately UBVRI) as 181 scattering angles by the six independent elements
`F11 F22 F33 F44 F12 F34`. Run `./run_scatmat.x all` for the full wavelength
grid if more bands are needed.

### Polarization in the ionizing band

The polarized products are computed **per band, for the wavelengths a run
actually needs** — the aligned scattering matrix ships five optical bands
(approximately UBVRI), and that is the pattern. The ionizing band is outside
what these products are used for, so shortward of 0.0912 um the polarized
optics are **not** part of the shipped model and no EUV companion is generated.

If `lam_min` widens the grid into that band anyway, `build_Cpol` looks for the
EUV companion table `q_astrodust_jori_euv_P0.20_Fe0.00_1.400.dat.gz`, does not
find it, and says so on stderr:

```
build_Cpol: no EUV polarized table (...)
            dichroic extinction and birefringence are zero below 9.120E-02 um
            (zero by omission, not by physics).
```

**Those zeros are a documented deficit, not a physical result.** Measured from
first principles (`tmatrix/driver/euv_polarized_optics.f90`), the true
alignment-weighted `|C_pol,ext| / C_ext` is 1.26e-3 at the 0.0912 um join, rises
to 3.64e-3 near 20.6 eV — a factor 2.9 — and decays to 1.6e-4 at 100 eV, **with
the sign opposite to the optical band**; the reversal falls near 0.106 um. A
sphere has exactly zero dichroism and zero birefringence, so the scalar EUV
optics (Mie on the volume-equivalent sphere) could not have supplied these
either.

### Computing the wavelengths you need

There is little reason to run polarized transfer at every wavelength; the five
bands above cover the usual case. What matters is that any other wavelength
*can* be computed. `tmatrix/run_q_jori.x` takes wavelengths directly:

```sh
./run_q_jori.x lam L1 [L2 ...] [ja=JA1:JA2] [tag=NAME]   # wavelengths [um] on the command line
./run_q_jori.x lamfile PATH    [ja=JA1:JA2] [tag=NAME]   # one per line, '#' comments allowed
./run_q_jori.x lammerge STEM FILE [FILE ...]             # reassemble the ja= windows
```

It writes a pair, `q_astrodust_jori_P0.20_Fe0.00_1.400.TAG.dat` and
`.TAG.wave`. Hand that pair to `sed_init` / `build_astrodust` as
`qpol_euv_path` and `qpol_euv_wave_path`, and `build_Cpol` reads it with a
double interpolation in log(lambda) and log(a). Because it interpolates, the
file needs **at least two wavelengths**.

`ja=JA1:JA2` splits the 169 radii across **processes**, not threads:
Mishchenko's solver hands the converged T-matrix to `AMPL` through
`COMMON /TMAT/` and keeps further working storage in COMMON blocks, two of them
blank, so the core is not re-entrant. Measured on this machine: 9m22s wall for
34m48s of processor time, 57 processes over three wavelengths.

**Convergence is certified, not assumed.** In `lam` / `lamfile` / `euv` mode
every node with size parameter x in (50, 60] is solved twice — `(DDELT, NDGS) =
(1e-3, 2)` and `(3e-3, 3)` — and rejected if the two disagree by more than 5%,
in which case the geometric-optics limit is used instead. Above x = 60 the
geometric-optics limit is used outright. This matters because **`IERR = 0` is
not a sufficient condition for convergence above x ~ 55**: at lam = 0.0912 um,
a = 0.9441 um (x = 65.04) the solver returns `IERR = 0` and a value ~50% away
from its neighbors. The older `range` and full-sweep modes keep the previous
x > 50 rule, so they still reproduce the shipped table byte for byte.

The certification has a price. Because x > 60 is left at the geometric-optics
zero for the dichroic extinction, about 46% of the band's dichroism is missing
at 100 eV (nothing at the join, ~8% at 0.031 um, ~20% at 0.022 um). The true
value lies between the certified 1.06e-4 and the 1/x extrapolation 3.14e-4 —
already ~1% of the optical-band signal, so the residual is a small fraction of
an already small number. `C_pol`, which drives polarized *emission*, has no such
gap: the geometric-optics limit gets it from the opaque-grain Fresnel surface
integral.

**Scattering matrices at other wavelengths are more awkward.**
`run_scatmat_aligned.x` already accepts wavelengths on the command line
(`./run_scatmat_aligned.x 0.44 0.55 0.79`), but three things are missing that
`run_q_jori.x` has: the output stem is fixed, so a custom run **overwrites** the
shipped five-band table; there is no merge mode, so bands cannot be *added* to
an existing table (a full recomputation costs ~4.5 min per band); and there is
neither a radius split nor the x > 50 certification above. Adding a band today
therefore means regenerating the whole table and accepting the uncertified
large-x treatment.

**The library reads tables; it never runs the T-matrix.** That is not a
code-size decision — it is that convergence cannot be certified at run time
above x ~ 55, and a silent wrong value is worse than an absent one. (A revision
of the T-matrix core is planned, so read this as the current policy rather than
a permanent one.)

The birefringence that converts linear into circular polarization on propagation
is the real part of the same forward-amplitude difference whose imaginary part
gives the dichroism, and it comes for free from the fixed-orientation amplitude
already computed. The regenerated table stores it as an optional 4th block, from
which the birefringence cross section follows; since no astrodust reference for it
exists, it is certified internally by Kramers-Kronig against the dichroism to a
median of about 0.1%. `dust_extinction` returns it through an optional `Cbir_ext`
argument, which is zero when the loaded table has no 4th block, as the release
table does not; consuming it in the transfer is the RT code's task.

Scattering by *aligned* grains is now computed from first principles for the
astrodust spheroid. `tmatrix/run_scatmat_aligned.x` builds the fixed-orientation
Mueller matrix `Z(theta_i; theta_s, phi)` of the DH21 oblate spheroid — from
Mishchenko's fixed-orientation amplitude in the T-matrix regime and the analytic
dipole below it — and size-integrates it over the astrodust distribution with the
alignment weight `f_align(a)`, writing five optical bands (approximately UBVRI) to
`output/scatmat_aligned_astrodust_P0.20_Fe0.00_1.400.dat`. The same file carries
the 4x4 extinction matrix `K(theta_i)` — the total `Cext`, the dichroic `Cpol`,
and the birefringent `Cbir` on the incidence-angle grid, from the forward
amplitudes and therefore exact at every propagation angle rather than interpolated
as `sin^2` — together with the two random-orientation matrices of the aligned and
of the full population, so the unaligned remainder follows by subtraction. A
cell's alignment enters only through a scalar `eta` (its local alignment scale)
and `theta_i = acos(k-hat . B-hat)`: the aligned optics scale as `eta` and the
unaligned extinction adds `Cext_tot - eta*Cext_ref`, exact by the linearity of the
size integral in `f_align`. The table ships gzipped (the reader opens `.gz`
directly) and `run_scatmat_aligned.x` regenerates it in minutes (`profile=FILE`
swaps the alignment profile).

The library reads it with the same lifecycle as the rest of the API: a
`scatmat_path` argument on `sed_init` / `build_astrodust` loads and integrates it
once (serial), and six pure-read query calls — `extinction_matrix_aligned`,
`mueller_matrix_aligned`, `mueller_matrix_random`, `mueller_matrix_total`,
`scattering_cross_sections`, and the band selector `scatmat_band` — serve a photon
path concurrently from OpenMP threads. `mueller_matrix_total` is the recommended
one: it returns the absolute combined phase matrix — the aligned part plus the
random-orientation remainder, correctly `1/(4 pi)`-normalized and rotated into the
grain frame — in a single call. `sed/rt_example/use_dustlib_scatmat.f90` is a
minimal two-cell reference consumer. This completes the material side of the Peest-formalism
contract for aligned-grain polarized transfer: SEDust returns every quantity in
matrix form, and the MoCafe-side consumption (frame rotations, direction sampling,
peel-off, the `exp(-K tau)` step) remains future work. The random-orientation
scatmat file above stays for unaligned use, and it does not limit far-infrared or
submillimeter polarized emission, where scattering is negligible. The PAH
component is treated as unaligned, and the DL07 and Zubko models have no polarized
optics.

A build can also skip the polarized optics entirely: `load_polarized_optics=.false.`
on `build_astrodust` / `sed_init` never opens the orientation-resolved table, leaves
`Cpol`/`Cpol_ext`/`Cbir_ext`/`falign` zero, and returns the scalar cross sections
and total SED bit-identical to a polarized build at zero alignment.
`dust_has_polarized_optics(m)` reports whether a built model carries polarized
optics at all.

## Documentation

- [`docs/astrodust_sed_report.pdf`](docs/astrodust_sed_report.pdf) — the astrodust+PAH pipeline, its validation
  against the HD23 release, and the resolution of the far-infrared offset.
- [`docs/SEDust_user_manual.pdf`](docs/SEDust_user_manual.pdf) — the `libsedust.a` API: model builders,
  channels, solver options, and how to link it into an RT code.
- [`docs/mc_pT_report.pdf`](docs/mc_pT_report.pdf) — the Monte Carlo algorithm, its adaptive-grid
  engines, and its validation against the matrix solvers.
- [`docs/sedust_polarization_implementation.pdf`](docs/sedust_polarization_implementation.pdf) — how the polarized optics are
  built: the orientation-resolved table, the derived cross sections, the
  aligned-grain scattering matrix and its two-layer API, the implementation
  decisions and their reasons, and the verification.
- [`docs/aligned_grain_polarization.pdf`](docs/aligned_grain_polarization.pdf) — background on grain alignment and
  polarized radiative transfer, and what a radiative-transfer code would need
  in order to use the polarized optics.

## References

- [Draine & Anderson 1985, ApJ, 292, 494](https://ui.adsabs.harvard.edu/abs/1985ApJ...292..494D/abstract) — the Monte Carlo temperature-history method the `mc/` solver follows.
- [Guhathakurta & Draine 1989, ApJ, 345, 230](https://ui.adsabs.harvard.edu/abs/1989ApJ...345..230G/abstract) — the transition-matrix solution for `P(T)` of stochastically heated grains.
- [Draine & Li 2001, ApJ, 551, 807](https://ui.adsabs.harvard.edu/abs/2001ApJ...551..807D/abstract) — the grain enthalpy used here, and the LD01 PAH cross-section vintage.
- [Weingartner & Draine 2001, ApJ, 548, 296](https://ui.adsabs.harvard.edu/abs/2001ApJ...548..296W/abstract) — the Milky Way size distributions the DL07 model is built on.
- [Zubko, Dwek & Arendt 2004, ApJS, 152, 211](https://ui.adsabs.harvard.edu/abs/2004ApJS..152..211Z/abstract) — the ZDA BARE-GR-S composition, size distributions, and calorimetry.
- [Draine & Li 2007, ApJ, 657, 810](https://ui.adsabs.harvard.edu/abs/2007ApJ...657..810D/abstract) — the DL07 PAH cross sections and emission model.
- [Mishchenko & Travis 1998, JQSRT, 60, 309](https://ui.adsabs.harvard.edu/abs/1998JQSRT..60..309M/abstract) — the T-matrix method used for the spheroidal astrodust optics.
- [Draine & Hensley 2021, ApJ, 909, 94](https://ui.adsabs.harvard.edu/abs/2021ApJ...909...94D/abstract) — the astrodust dielectric function.
- [Hensley & Draine 2023, ApJ, 948, 55](https://ui.adsabs.harvard.edu/abs/2023ApJ...948...55H/abstract) — the astrodust+PAH model this package reproduces.

## Author

Kwang-il Seon (KASI/UST)

---

Last updated: 2026-08-21 14:11 KST
