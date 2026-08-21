# SEDust — build, run, and the details

Companion to [README.md](README.md): how to build the package, every run
setting the drivers take, and the two parts of the physics that decide what a
run actually computes — how far into the ionizing band the grid reaches, and
which stochastic-heating solver resolves `P(T)`.

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
./calc_sed.x mrn            # MRN (1977) SED at U = 1           -> output/
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
./calc_sed.x mrn draine euv                 # -> sed_mrn_euv_draine.dat
./calc_sed.x zubko euv hardfield            # -> sed_zubko_euv_hardfield.dat

# the library, for embedding in an RT code.  No T-matrix in it: link with
#   -I. -L. -lsedust
make libsedust.a
# ... and with the astrodust EUV band on the spheroid of the Q table, in the
# same archive and on the same link line (the host then calls
# use_tmatrix_euv_band_optics() once before build_astrodust)
WITH_TMATRIX=1 make libsedust.a

# size-integrated transport optics: lambda, albedo, <cos>, C_ext/C_abs/C_sca per H
./calc_kext.x astrodust     # -> ../data/astrodust/kext_astrodust_MW.dat
./calc_kext.x astrodust euv # the same model carried into the ionizing band
./calc_kext.x dl07 euv      # -> ../data/dl07/kext_dl07_MW_euv.dat
./calc_kext.x mrn           # -> ../data/mrn/kext_mrn.dat
./calc_kext.x mrn euv       # -> ../data/mrn/kext_mrn_euv.dat
./calc_kext.x zubko         # -> ../data/zubko/kext_zubko_BARE_GR_S.dat (cut at the Lyman limit)
./calc_kext.x zubko euv     # -> ../data/zubko/kext_zubko_BARE_GR_S_euv.dat (the whole ZDA range)
./calc_kext.x from_files ../data/zubko/zubko_descriptor.txt

# the optics products.  ORDER MATTERS: calc_qtable.x lays down the wavelength
# axis and replaces data/<model>/sedust_<model>.h5, so it runs first and calc_kext.x
# puts /kext back.  Both also write the text products beside them.
./calc_qtable.x             # every model -> ../data/<model>/q_*.dat
                            #             -> ../data/<model>/sedust_<model>.h5
./calc_kext.x astrodust euv # ... then /kext into each of those files
./check_build_dust.x        # build_dust on HDF5 vs the builders on text

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

# regenerating the T-matrix Q tables (optional; both tables ship with SEDust)
cd ../tmatrix && make && ./run_tmatrix.x test   # smoke test -> output/..._euv.test.dat
./run_tmatrix.x                                 # full sweep -> output/..._euv.dat
make lyman_cut                                  # -> the 1129-wavelength companion
```

Outputs are plain ASCII `.dat` files written to each subdirectory's `output/`;
`calc_kext.x` writes into `data/` instead.

## The ionizing band

The HDF5 product holds ONE wavelength axis and the index `i_lyman` at which it
crosses the Lyman limit, so `include_euv` picks the view rather than the file:
`.false.` (the default) returns `lambda(i_lyman:)` and the same rows of every
wavelength-indexed array.  For astrodust that is 1129 wavelengths
(0.0912--39810 um) out of 1762 (1.0e-4--39810 um), for DL07 and MRN 1129 out of
1823, for Zubko 866 out of 1201.

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

`build_astrodust`, `build_dl07` and `build_mrn` take an optional `lam_min` [um] for a host
needing to reach below the table it was given. It prepends log-spaced points
below the table, and is refused when it asks for a wavelength the model's own
dielectric data does not cover. What that does for astrodust depends on the
table:

- on the default 1129-wavelength table, a `lam_min` below 0.0912 um **does**
  build a band. Its optics come from the T-matrix on the same spheroid when
  `libtmatrix.a` is linked in (`WITH_TMATRIX=1 make libsedust.a`), and are
  refused with status 6 when it is not — never silently replaced by a sphere
  unless `euv_tmatrix = .false.` asks for one;
- on the `_euv` table there is nothing left to ask for: the DH21 dielectric
  function stops at 1.000032e-4 um, longward of that table's first wavelength,
  so every legal `lam_min` falls inside the table and prepends nothing.

DL07 and MRN can be extended below either table, to 6.205e-5 um, because the D03
optical constants reach further than both; their band is Mie on those functions
and needs no T-matrix.

| product | rows | lambda [um] | what sets the floor |
|---|---:|---|---|
| `data/astrodust/kext_astrodust_MW.dat` | 1129 | 0.0912 - 39810 | non-EUV Q-table grid |
| `data/astrodust/kext_astrodust_MW_euv.dat` | 1762 | 1e-4 - 39810 | EUV Q-table grid |
| `data/dl07/kext_dl07_MW.dat` | 1129 | 0.0912 - 39810 | non-EUV Q-table grid |
| `data/dl07/kext_dl07_MW_euv.dat` | 1823 | 6.205e-5 - 39810 | the D03 dielectric functions |
| `data/mrn/kext_mrn.dat` | 1129 | 0.0912 - 39810 | non-EUV Q-table grid |
| `data/mrn/kext_mrn_euv.dat` | 1823 | 6.205e-5 - 39810 | the D03 dielectric functions |
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
| MRN | `data/mrn/q_mrn_{sil,gra}.dat` (1129) | `..._euv.dat` (1823) | 70, log-spaced 0.005 - 0.25 um, the power law's own cutoffs on its ends |
| Zubko | `data/zubko/q_zubko_{sil,gra,pah}.dat` (866) | `..._euv.dat` (1201) | 121 (sil, gra) / 28 (PAH), the ZDA tables' own |

`sed/calc_qtable.x` writes all of these, and the HDF5 product of each model
with them; the astrodust scalar pair is the exception, computed by
`tmatrix/run_tmatrix.x` plus `make lyman_cut` and installed into
`data/astrodust/`. Every non-EUV file is the row subset of its EUV counterpart
that starts at the Lyman limit — verified, not asserted. In the HDF5 product the
pair is one array and its `i_lyman`, so there is no second file to keep in
step (§ *The HDF5 optics products* in the manual).

**These are our own tables.** For DL07 and MRN they are what those models
already compute at build time, written down rather than recomputed -- the two
share the calculation and differ in the radius grid it runs on. For Zubko they are a genuine
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

**This changed all three models shipped at the time**, and for Zubko it has nothing to do
with the astrodust table: that model's grid has always started at 1e-3 um, so it
carried the same defect independently. Measured against the previous products:
astrodust up to 12.1% locally and −1.1% to +0.9% in integrated power, DL07 9.7%
and −2.8%, Zubko 3.0% (6.0% with `euv`). Energy conservation improved, from
+0.236% to −0.057% (astrodust S1, emitted/absorbed per H).

**Any SED produced before this fix includes photons that were not in the field**,
in the top enthalpy bin of every stochastically heated grain. Do not compare
against those products without accounting for it.

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
