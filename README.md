# SEDust

A self-contained library for the optical properties and thermal emission of
interstellar dust: dielectric functions -> T-matrix / Mie cross sections ->
grain temperature distributions -> emergent infrared SED.

SEDust is **model-agnostic**. The HD23 astrodust+PAH model, the Draine & Li
(2007) carbonaceous+silicate model, and the Zubko et al. (2004) BARE-GR-S model
are handled as peers through one derived type (`dust_model_t`) and one emission
call (`dust_emission`). The solver core does not know which model it is running.

The whole package can also be linked into a Fortran 3D radiative-transfer code
as a static library, `libsedust.a`, with a two-step API: initialize once, then
solve one cell at a time. Alongside the emission call, `dust_extinction` returns
the model's size-integrated extinction, scattering, and asymmetry cross sections
per H on the same wavelength grid, so a transfer code takes its opacity and its
emission from one model object. Those cross sections come from the precomputed
`data/<model>/kext_*.dat` table the builder loaded, which is the same product
`calc_kext.x` writes; the size integral behind it is `size_integrated_extinction`.
Every model builder and solver call takes an
optional `status` argument, so a missing input file or an invalid model is
reported back to the host instead of stopping the process.

## Which version

SEDust is split into two versions because most radiative-transfer codes either
do not carry polarization at all, or carry it only for spherical grains.

This is **version 1.00**, the scalar branch: unpolarized cross sections, thermal
emission, and the `dust_emission` / `dust_extinction` API for a radiative-transfer
host. It covers that common case.

Use **version 1.20** instead for polarized transfer that accounts for the
non-spherical (**spheroidal**) grain shape -- dichroic extinction, birefringence,
and scattering by aligned spheroidal grains. It adds the direction-dependent
extinction matrix and the aligned-grain scattering matrix on top of everything
here.

## Layout

```
SEDust/
  sed/          the SED solver: cross sections, enthalpy, P(T), emission
    src/        library modules + drivers
    rt_example/ minimal example of linking libsedust.a into an RT code
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
# by main_astrodust.x, calc_kext.x, and a WITH_TMATRIX library build
cd tmatrix && make libtmatrix.a

# the SED solver
cd ../sed
make                        # make_enthalpy.x main_astrodust.x main_dl07.x
                            # main_zubko.x calc_kext.x
./main_astrodust.x          # astrodust+PAH SED at log U = 0.20 -> output/
./main_dl07.x               # Draine & Li (2007) SED at U = 1   -> output/
./main_astrodust.x euv      # ... either driver on the _euv Q table's 1762 grid
./main_zubko.x              # Zubko/ZDA SED at U = 1            -> output/
./main_zubko.x euv          # ... on the whole ZDA range, not cut at the Lyman limit
./main_zubko.x euv hardfield  # ... and with a hard component added to the FIELD

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
./calc_kext.x zubko         # -> ../data/zubko/kext_zubko_BARE_GR_S.dat (cut at the Lyman limit)
./calc_kext.x zubko euv     # -> ../data/zubko/kext_zubko_BARE_GR_S_euv.dat (the whole ZDA range)
./calc_kext.x from_files ../data/zubko/zubko_descriptor.txt

# the optics products.  ORDER MATTERS: make_qtable.x lays down the wavelength
# axis and replaces data/<model>/sedust_<model>.h5, so it runs first and calc_kext.x
# puts /kext back.  Both also write the text products beside them.
./make_qtable.x             # all three models -> ../data/<model>/q_*.dat
                            #                 -> ../data/<model>/sedust_<model>.h5
./calc_kext.x astrodust euv # ... then /kext into each of the three files
./check_build_dust.x        # build_dust on HDF5 vs the builders on text

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
(0.0912--39810 um) out of 1762 (1.0e-4--39810 um), for DL07 1129 out of 1823,
for Zubko 865 out of 1201.

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
  refused with status 6 when it is not — never silently replaced by a sphere
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
| `data/zubko/kext_zubko_BARE_GR_S.dat` | 865 | 0.0912 - 1e4 | the ZDA optics table, cut at the Lyman limit |
| `data/zubko/kext_zubko_BARE_GR_S_euv.dat` | 1201 | 1e-3 - 1e4 | the ZDA optics table itself |

### Stored Q tables

Every model's `(lambda, a_eff)` cross sections are in `/qtable` of its HDF5
product, which is what ships.  `make_qtable.x` writes them as text beside it as
well, and those files are what the table below names; they are regenerable, so
they are not tracked:

| model | non-EUV | EUV | a_eff grid |
|---|---|---|---|
| astrodust | `data/astrodust/q_astrodust_P0.20_Fe0.00_1.400.dat` (1129) | `..._euv.dat` (1762) | 169, from `data/astrodust/DH21_aeff` |
| DL07 | `data/dl07/q_dl07_{sil,gra,pah_neu,pah_ion}.dat` (1129) | `..._euv.dat` (1823) | 84, Draine's analytic grid, 3.548e-4 - 5.012 um |
| Zubko | `data/zubko/q_zubko_{sil,gra,pah}.dat` (865) | `..._euv.dat` (1201) | 121 (sil, gra) / 28 (PAH), the ZDA tables' own |

`sed/make_qtable.x` writes all of these, and the HDF5 product of each model
with them; the astrodust scalar pair is the exception, computed by
`tmatrix/run_tmatrix.x` plus `make lyman_cut` and installed into
`data/astrodust/`. Every non-EUV file is the row subset of its EUV counterpart
that starts at the Lyman limit — verified, not asserted. In the HDF5 product the
pair is one array and its `i_lyman`, so there is no second file to keep in
step (§ *The HDF5 optics products* in the manual).

**These are our own tables.** For DL07 they are what the model already computes
at build time, written down rather than recomputed. For Zubko they are a genuine
recomputation: the model itself still reads the ZDA optics tables under
`data/zubko/`, while these are what this tree's Mie gives on the dielectric
functions those tables name (`data/dielectric/eps_Sil`, `eps_Gra`, added with
their provenance in `where_draine_eps.txt`). Over the wavelengths the dielectric
files cover, our silicate table agrees with the shipped one to a mean 3.5% in
Q_abs and 4.7% in Q_sca, with `<cos>` to 0.05; longward of 1e3 um the
dielectric functions run out and the two extension laws differ, which the file
headers state.

The one text table anything still opens is the astrodust scalar pair. It is the
INPUT `make_qtable.x` reads to write the HDF5 product, so that program cannot
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
`./main_zubko.x [euv]` still prints both candidates side by side, and is the
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
| `'qm'`        | energy-space transition matrix, BiCG sparse solve | thermal-discrete (`dbdis`) or thermal-continuous (`dbcon`) |
| `'equil'`     | equilibrium temperature, no stochastic solve | n/a |

The Monte Carlo solver in `mc/` follows Draine & Anderson (1985) and tracks
`T(t)` through individual absorption events. It shares no code with the matrix
solvers, which is what makes it a useful check on them.

Grains whose equilibrium enthalpy exceeds 150 eV are placed at `T_eq` and skip
the stochastic solve; the gate propagates forward in grain size.

## Documentation

- `docs/astrodust_sed_report.pdf` — the astrodust+PAH pipeline, its validation
  against the HD23 release, and the resolution of the far-infrared offset.
- `docs/SEDust_user_manual.pdf` — the `libsedust.a` API: model builders,
  channels, solver options, and how to link it into an RT code.
- `docs/mc_pT_report.pdf` — the Monte Carlo algorithm, its adaptive-grid
  engines, and its validation against the matrix solvers.

Rebuild any of them with `pdflatex <name>.tex` (run twice for cross-references).

## References

- Draine, B. T., & Anderson, N. 1985, ApJ, 292, 494
- Guhathakurta, P., & Draine, B. T. 1989, ApJ, 345, 230
- Draine, B. T., & Li, A. 2001, ApJ, 551, 807
- Weingartner, J. C., & Draine, B. T. 2001, ApJ, 548, 296
- Zubko, V., Dwek, E., & Arendt, R. G. 2004, ApJS, 152, 211
- Draine, B. T., & Li, A. 2007, ApJ, 657, 810
- Mishchenko, M. I., & Travis, L. D. 1998, JQSRT, 60, 309
- Draine, B. T., & Hensley, B. S. 2021, ApJ, 909, 94
- Hensley, B. S., & Draine, B. T. 2023, ApJ, 948, 55

## Author

Kwang-il Seon (KASI/UST)

---

Last updated: 2026-08-07 17:54 KST
