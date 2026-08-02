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

## Layout

```
SEDust/
  sed/          the SED solver: cross sections, enthalpy, P(T), emission
    src/        library modules + drivers
    rt_example/ examples of linking libsedust.a into an RT code, one
                minimal and two showing the polarized quantities
  mc/           Draine & Anderson (1985) Monte Carlo solver (independent check)
  tmatrix/      Mishchenko T-matrix engine + driver; writes the Q table
  data/         dielectric functions, the HD23 public release tables,
                the Zubko (ZDA BARE-GR-S) optical constants, and the
                size-integrated `kext_*.dat` transport-optics curves
  docs/         technical reports and the library user manual
  pyutil/       small Python helpers (radiation fields, SED from Cabs)
```

Everything the code reads at run time ships with the package. There are no
paths outside this directory.

## Build and run

Requires `gfortran` (OpenMP for the parallel drivers). No autoconf, no
top-level configure; each subdirectory has its own `Makefile`.

```sh
# the SED solver
cd sed
make                        # make_enthalpy.x main_astrodust.x main_dl07.x
                            # main_zubko.x calc_kext.x
./main_astrodust.x          # astrodust+PAH SED at log U = 0.20 -> output/
./main_dl07.x               # Draine & Li (2007) SED at U = 1   -> output/
./main_zubko.x              # Zubko/ZDA SED at U = 1            -> output/
./main_zubko.x euv          # ... with the field carried into the EUV

# the library, for embedding in an RT code
make libsedust.a            # link with:  -L. -lsedust -I.
make use_dustlib_scatmat.x  # reference consumer of the aligned-scattering API

# size-integrated transport optics: lambda, albedo, <cos>, C_ext/C_abs/C_sca per H
./calc_kext.x astrodust     # -> ../data/kext_astrodust_MW.dat (+ dichroic column)
./calc_kext.x astrodust euv # the same model carried into the ionizing band
./calc_kext.x dl07 euv      # -> ../data/kext_dl07_MW_euv.dat
./calc_kext.x zubko         # -> ../data/kext_zubko_BARE_GR_S.dat
./calc_kext.x from_files ../data/zubko/zubko_descriptor.txt

# polarized extinction alone, checked against the HD23 release
make calc_polext.x         && ./calc_polext.x

# the Monte Carlo cross-check
cd ../mc && make && ./main_mc_sed.x run_sed.nml

# regenerating the T-matrix Q table (optional; the table ships with SEDust)
cd ../tmatrix && make && ./run_tmatrix.x test   # then ./run_tmatrix.x for the full sweep

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

The astrodust wavelength grid is the T-matrix Q table's, and that table stops at
0.0912 um — 13.6 eV, the Lyman limit. A host that transports ionizing radiation
passes the shortest wavelength it needs as the optional `lam_min` to
`build_astrodust` or `build_dl07`; the grid is then carried down to it. **Omit
`lam_min` and nothing changes** — the unextended model is bit-identical to
before.

| product | rows | lambda [um] | what sets the floor |
|---|---:|---|---|
| `data/kext_astrodust_MW.dat` | 1129 | 0.0912 - 39810 | the Q-table grid (no extension) |
| `data/kext_astrodust_MW_euv.dat` | 1719 | 1.001e-4 - 39810 | the DH21 astrodust dielectric function |
| `data/kext_dl07_MW_euv.dat` | 1761 | 6.205e-5 - 39810 | the D03 dielectric functions |
| `data/kext_zubko_BARE_GR_S.dat` | 1201 | 1e-3 - 1e4 | the ZDA optics table itself |

No floor is a free choice: each is the shortest wavelength the data the model is
made of actually covers, and a shorter request is refused rather than served
with a refractive index frozen at the table boundary.

In the extended band the astrodust optics are Mie for the volume-equivalent
sphere on the same DH21 dielectric function the Q table was computed from. That
shape approximation is bounded by about 2% and converges to a constant -2.08% in
the large-particle limit, not to zero. **There is no external reference for this
band**: the HD23 release stops at 12.4 eV, exactly where the extension begins.
The DL07 model *can* be checked, and agrees with Draine's published table for
the same model to 0.056-0.86% in C_ext per decade. See
`docs/SEDust_user_manual.pdf` for the full accounting.

**The polarized optics do not extend with it.** The four `kext_*.dat` products
above are scalar transport optics; only the unextended
`data/kext_astrodust_MW.dat` carries the dichroic extinction, as its eighth
column. Below 0.0912 um the dichroism and the birefringence need a separate EUV
companion table, which does not ship — see [Polarization in the ionizing
band](#polarization-in-the-ionizing-band).

### The single-photon ceiling comes from the field

The stochastic solvers cap the enthalpy window at the energy of the hardest
single photon a grain can absorb, `max(13.6 eV, hc/lambda_hard)`, where
`lambda_hard` is the shortest wavelength at which the incident `J_lambda`
exceeds 1e-12 of its own peak (`hardest_photon_energy` in
`sed/src/radfield.f90`). **That ceiling is a property of the radiation field,
not of the wavelength grid.** The two coincide for astrodust and DL07, whose Q
table floor is the Lyman limit itself (13.595 eV < 13.6 eV, so the constant
stays selected); they differ by a factor of 91 for Zubko, whose DustEM tables
start at 1e-3 um (1239.8 eV) while a Mathis field carries nothing below
0.0912 um. Reading the grid there would coarsen the fixed-count enthalpy bins —
and the refinement loops, guarded and capped at ten iterations on the way down,
would not walk it back — shifting the emergent SED by 1-2% with no photon behind
it. `docs/EUV_EXTENSION_HOST_REGRESSION.md` records the measurement.

`./main_zubko.x [euv]` is the driver that keeps this path covered: Zubko is the
one shipped model whose optics grid reaches far past the band a transported field
occupies, so anything that reads the grid where it should read the field shows up
there and nowhere else. It prints both candidates side by side. **Nothing changes
for a model built without `lam_min`** — astrodust and DL07 keep the 13.6 eV
constant either way.

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

## Polarization

For the astrodust model SEDust also computes the polarized cross sections of
aligned spheroidal grains, from an orientation-resolved spheroid table. By
default this is the Draine & Hensley (2021) table that ships in
`data/dielectric/`. Two quantities are available:

| Quantity | Where |
|---|---|
| polarized emission | optional `lamI_pol` argument of `dust_emission` |
| polarized extinction | optional `Cpol_ext` argument of `dust_extinction`, or the eighth column of `data/kext_astrodust_MW.dat` |

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
release file `data/dielectric/q_DH21Ad_P0.20_Fe0.00_1.400.dat.gz`; SEDust's own
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

Shortward of 0.0912 um the polarized optics are **not** part of the shipped
model. If `lam_min` widens the grid into that band, `build_Cpol` looks for the
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

- `docs/astrodust_sed_report.pdf` — the astrodust+PAH pipeline, its validation
  against the HD23 release, and the resolution of the far-infrared offset.
- `docs/SEDust_user_manual.pdf` — the `libsedust.a` API: model builders,
  channels, solver options, and how to link it into an RT code.
- `docs/mc_pT_report.pdf` — the Monte Carlo algorithm, its adaptive-grid
  engines, and its validation against the matrix solvers.
- `docs/sedust_polarization_implementation.pdf` — how the polarized optics are
  built: the orientation-resolved table, the derived cross sections, the
  aligned-grain scattering matrix and its two-layer API, the implementation
  decisions and their reasons, and the verification.
- `docs/aligned_grain_polarization.pdf` — background on grain alignment and
  polarized radiative transfer, and what a radiative-transfer code would need
  in order to use the polarized optics.

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

Last updated: 2026-08-02 10:09 KST
