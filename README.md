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

## Layout

```
SEDust/
  sed/          the SED solver: cross sections, enthalpy, P(T), emission
    src/        library modules + drivers
    rt_example/ examples of linking libsedust.a into an RT code, one
                minimal and one showing the polarized quantities
  mc/           Draine & Anderson (1985) Monte Carlo solver (independent check)
  tmatrix/      Mishchenko T-matrix engine + driver; writes the Q table
  data/         dielectric functions, the HD23 public release tables,
                and the Zubko (ZDA BARE-GR-S) optical constants
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
make                        # make_enthalpy.x  main_astrodust.x  main_dl07.x
./main_astrodust.x             # astrodust+PAH SED at log U = 0.20 -> output/
./main_dl07.x               # Draine & Li (2007) SED at U = 1   -> output/

# the library, for embedding in an RT code
make libsedust.a            # link with:  -L. -lsedust -I.
make use_dustlib_scatmat.x  # reference consumer of the aligned-scattering API

# optical-property tables (extinction, albedo, <cos>, K_abs, polarized ext.)
make calc_kext_astrodust.x && ./calc_kext_astrodust.x
make calc_kext_dl07.x      && ./calc_kext_dl07.x

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

# scattering matrix of randomly oriented grains (five optical bands ship with SEDust)
./run_scatmat.x 0.55                            # one wavelength; ./run_scatmat.x all for the grid

# fixed-orientation scattering matrix of ALIGNED grains (also five optical bands)
./run_scatmat_aligned.x test                    # one band, reduced grid: timing + OpenMP check
./run_scatmat_aligned.x                         # default UBVRI bands -> output/ (~8 min, 32 threads)
#   profile=FILE regenerates under a different alignment profile in minutes
```

Outputs are plain ASCII `.dat` files written to each subdirectory's `output/`.

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

---

Last updated: 2026-07-21 09:48 KST
