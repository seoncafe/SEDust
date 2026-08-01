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
emission from one model object. Every model builder and solver call takes an
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

# size-integrated transport optics: lambda, albedo, <cos>, C_ext/C_abs/C_sca per H
./calc_kext.x astrodust     # -> ../data/kext_astrodust_MW.dat
./calc_kext.x astrodust euv # the same model carried into the ionizing band
./calc_kext.x dl07 euv      # -> ../data/kext_dl07_MW_euv.dat
./calc_kext.x zubko         # -> ../data/kext_zubko_BARE_GR_S.dat
./calc_kext.x from_files ../data/zubko/zubko_descriptor.txt

# the Monte Carlo cross-check
cd ../mc && make && ./main_mc_sed.x run_sed.nml

# regenerating the T-matrix Q table (optional; the table ships with SEDust)
cd ../tmatrix && make && ./run_tmatrix.x test   # then ./run_tmatrix.x for the full sweep
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

Last updated: 2026-08-02 00:30 KST
