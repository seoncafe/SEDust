# Fifteen-field SEDust emission benchmark

`validation_driver.f90` evaluates the astrodust model on the three radiation shapes and five absorbed-power targets used by the original timing driver. It builds one model with `NT=200`, evaluates each selected emission method five times for each field, and writes every total and channel spectrum. The code normalizes each radiation shape with the original full-grid logarithmic trapezoid, so source changes do not alter the incident fields.

The shapes are `S0` (4.78–100 eV, a 40,000 K Planck spectrum), `S1` (the same band with attenuation above 13.6 eV), and `S2` (4.78 <= E < 13.6 eV). Target absorbed powers are `3e-23`, `1e-22`, `1e-21`, `1e-20`, and `1e-19` erg/s/H. Model arguments are `sd_index=7`, `u_isrf=1`, `include_euv=.true.`, and `T_first=2.7–5000 K` with 200 points.

Build the library with the ifort precise settings used by the MoCHII link from the source version being measured. From the repository root:

```sh
FC=ifort sed/build_lib.sh
sed/tests/build_validation.sh sed/lib corrected /tmp/validation_corrected.x
```

Use `exact` for a library without the new absorbed-power and emission-bin diagnostic symbols. Use `corrected` for a library exporting `radiation_absorption` and the `emission_bin_power_*` symbols. `FC` defaults to `ifort`, with `-O3 -qopenmp -recursive -fp-model precise -fpp -w` and HDF5 in `/opt/hdf5_intel`. Set `FC=gfortran` explicitly for a reference build; its HDF5 default is `/opt/hdf5_gfortran`. Use the same compiler and numerical flags for every speed or bit-identity comparison.

Create an empty output directory and run with one OpenMP thread:

```sh
mkdir -p /tmp/validation_corrected
OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE /tmp/validation_corrected.x 200 \
  'heuristic,heuristic#1e-6' /tmp/validation_corrected /absolute/path/to/SEDust/data
```

The `#1e-6` suffix sets `emission_bin_power_cutoff=1e-6` for the second method. Both settings run on the same initialized model. For a source version without this option, use `heuristic` alone.

`lambda.bin` starts with two 32-bit integers (`NLAM`, `NCHANNEL`), followed by two 64-bit arrays (`lambda`, `Cabs`). Each `field_S*_p*.bin` is one 64-bit `J_lambda` array. Each `spec_*.bin` holds the total `lambda*I_lambda` array followed by `NCHANNEL` arrays in Fortran column order. The `r1` through `r5` suffixes identify repeated evaluations. `timing.tsv` reports the original absorbed power, the absorbed power over the represented field support, emitted power, call times, repeat identity, and the total and omitted `P_i*kappB_i` sums with and without grain-number weighting.

To compare source variants, build separate libraries and executables from the pinned source versions, run each against the same data directory, and compare `lambda.bin` and all `field_*.bin` files byte for byte before comparing spectra. Keep the `exact` setting at zero cutoff when testing arithmetic identity. The optional emission-bin threshold should be judged from its omitted power and wavelength-band errors, as well as its timing.

The small Planck-grid test exercises 240,000 ordinary values, the all-zero Wien tail, zero temperature, an unsorted grid, nonpositive wavelengths, and points on both sides of `x=700` and `x=1e-4`. After building a library with `planck_grid_t`, run:

```sh
ifort -O3 -qopenmp -recursive -fp-model precise -I sed/lib \
  sed/tests/test_planck_grid.f90 sed/lib/libsedust.a -o /tmp/test_planck_grid.x
/tmp/test_planck_grid.x
```

The expected result is `BIT_DIFFERENCES=0`. The same test also passed with Intel Fortran `-O3 -fp-model precise`. Intel Fortran's default fast floating-point model changed 28,625 of the 240,000 ordinary bit patterns in this comparison through compiler arithmetic transformations; that result should not be described as a physics or algorithm difference. Use the project's precise floating-point setting for Intel builds when bit identity is required.

## Absorption boundary check

Run `FC=ifort sed/tests/check_radiation_support.sh` (ifort is the default). It compiles a small checked build
in a new temporary directory and verifies absorption boundaries, the cached
wavelength grid, equilibrium temperatures, the top-bin rate of a two-bin
grain (every absorbed photon raises it to the top bin, Draine & Li 2001
eq. 16, so 4 pi P(2) kappB(2) = P(1) times the support absorption, checked
as P(2)/P(1) to 1e-14), and
normalized probabilities for a dark field. It also checks that an isolated significant sample has zero continuum bandwidth in both absorption and the photon bound, reports `RADIATION_SUPPORT_ISOLATED` through the shared pure support query, and increments the public `isolated_radiation_support_count` through absorption/transition checks. Pure support queries do not increment the counter. The test does not build a dust model.

The same script then runs `test_absorption_balance.f90`, the energy closure of
`calc_P` for one grain: the steady state of the recursion emits exactly what its
upward rates absorb, and with the bin-overlap rates every bin below the top
absorbs the support integral `4 pi radiation_absorption`, so
`4 pi sum_(u>=2) P(u) kappB(u)` must equal `(1 - P(nT))` times that integral.
It covers 72 cases (six model grains, three field shapes, NT = 40 and 200, the
cached grid and another) with bin widths from 5e-10 of the softest field photon
to 700 times the hardest; the tolerance is 1e-14 and the worst case measured
2.2e-15 with ifort and 2.0e-15 with gfortran (2026-10-30, after the
physical constants became the exact SI values).

The benchmark is deliberately narrow: it measures the heuristic astrodust
solver on the stated fields. It does not establish convergence of the full
model with enthalpy resolution or the accuracy of every dust model.

## Emission threshold errors

After the ifort library build, the following test checks negative, too-large,
and NaN thresholds before model initialization. It also verifies that zero
cutoff reaches the existing model-validity check.

```sh
ifort -O3 -qopenmp -recursive -fp-model precise -I sed/lib sed/tests/test_emission_cutoff.f90 \
  sed/lib/libsedust.a /opt/hdf5_intel/lib/libhdf5_fortran.a \
  /opt/hdf5_intel/lib/libhdf5.a -lsz -lz -ldl \
  -o /tmp/test_emission_cutoff.x
/tmp/test_emission_cutoff.x
```

## Calling from OpenMP applications

The mutable dust model, `emission_bin_power_cutoff`, and emission diagnostics
are shared module state. Applications must serialize calls to `dust_emission`
from host OpenMP threads and change the option only between solves. The
library's own OpenMP loops remain supported. The isolated-support counter
uses an atomic update so internal parallel support checks cannot lose counts.
Read that cumulative counter before and after a solve; it counts checks,
not distinct fields or grains, and must not be reset during active calls.


## CMB accounting

`test_cmb_field_balance.f90` checks the full-field convention: a zero field
leaves the lowest enthalpy bin populated, and a pure CMB field or a field
without CMB obeys the finite-bin GD energy-closure relation. It does not test
exact Planck detailed balance, which this continuous-cooling approximation
does not enforce.

`test_cmb_population_invariance.f90` exercises `build_astrodust` and
`dust_emission` with both `heuristic` and `draine`. It perturbs the legacy
`kappCMB` arrays in memory and requires the emitted spectrum to remain exactly
unchanged. The complete supplied field is the only heating source. It
compiles the dust-model sources and `dust_lib` and links HDF5, and takes the
absolute data-directory path as its first argument; it does not write model
data.

```sh
cd sed
make test_cmb_field_balance.x && ./test_cmb_field_balance.x
make test_cmb_population_invariance.x
./test_cmb_population_invariance.x "$PWD/../data"
```

## HD23 field and current CMB checks

`test_hd23_field.f90` checks the UV boundaries of the final journal HD23
field, the three diluted blackbodies, `U` scaling of stellar light, the
unscaled 2.725 K CMB, and both Mathis conventions. The HD23 field follows
Hensley & Draine (2023), equations (9) and (10); `hd23` is an explicit
`calc_sed.x` option while corrected Mathis remains the default.

`test_qm_invalid_bins.f90` checks that `qm_solve_grain` refuses fewer than
five enthalpy bins, fewer than two transition-matrix wavelengths, and fewer
than two input wavelengths, returning `solved = .false.` and a zero spectrum.
`test_qm_rate_closure.f90` builds the Zubko BARE-GR-S model with its stored
optics and the U=1 Mathis field plus the hard EUV component. For one selected
PAH size, cooling method, and transition-grid count, it requires `QM_SOLVED`,
finite nonnegative emission, zero raw negative transition rates, and emitted
power within 1% of absorbed power. It verifies the reported actual cooling
method: `dbcon` and `dbdis` must remain selected, while a requested `stati`
may use `dbdis` when its mode count cannot be resolved. It also checks that a
negative absorption cross section returns `QM_INVALID_INPUT`, zero emission,
and zero raw negative rate count. The 1% bound covers finite energy-bin and
integration errors while rejecting the previously measured 3-4% small-PAH
energy discrepancy.
`mc/tests/test_mc_supplied_field.f90` checks that the Monte Carlo grain setup
uses a supplied field unscaled, independent of `U_isrf`, and that its default
equals an explicit `J_Mathis` field.

Each test has a make target in its own directory's `Makefile`; it compiles
only the modules the test uses.

```sh
cd sed
make test_hd23_field.x      && ./test_hd23_field.x
make test_qm_invalid_bins.x && ./test_qm_invalid_bins.x
make test_qm_rate_closure.x
./test_qm_rate_closure.x "$PWD/../data/zubko" 1 dbcon 150
./test_qm_rate_closure.x "$PWD/../data/zubko" 1 dbdis 150
./test_qm_rate_closure.x "$PWD/../data/zubko" 1 stati 150
cd ../mc
make test_mc_supplied_field.x && ./test_mc_supplied_field.x
```

On 2026-10-30, with gfortran 16.2.0 and the default `HDF5_PREFIX`, both v1.00
and v1.20 gave: `test_cmb_field_balance.x` passed the zero-field population
check with relative energy closures `-7.7716e-16` (pure CMB) and `-3.3307e-16`
(field without CMB) after the physical constants became the exact SI values
(`8.8818e-16` and `-9.9920e-16` before); `test_cmb_population_invariance.x` left the `heuristic`
and `draine` spectra unchanged; `test_qm_invalid_bins.x` and
`test_mc_supplied_field.x` passed; `test_hd23_field.x` passed once its
`uv_intensity` ended the UV segment at `lambda <= 0.246 um`, where `J_HD23`
ends it (a first version of the test stopped at that boundary). These focused
checks do not validate full SED or MC figure products.

For `test_qm_rate_closure.x`, argument 2 is the one-based PAH size index;
argument 3 is `dbcon`, `dbdis`, or `stati`; argument 4 is the transition-grid
count (`0` selects the full optical grid). Sizes 10, 14, and 20 can be used
for additional rate checks. The test reads the existing Zubko data and does
not write model products.

## Resolution and Rayleigh-limit diagnostics

A pre-repair single-grain diagnostic of 2026-10-30 (its driver and record are kept with
the manuscript calculations, outside this repository) used the v1.00
astrodust HDF5 input, a neutral 35-atom PAH of radius
`0.000421697 um`, HD23 `U=10^0.2`, `NT=200` from 2.7 to 5000 K, and
thermal-discrete `qm` with `nstate=500`. Emitted-to-absorbed power ratios were
0.8439877849 (solution rejected), 0.9683433756, and 0.9962086841 when the
radiation field used 200, 500, and 1129 wavelengths. The matching ionized PAH
at 200 wavelengths gave 1.0068300802. This single-grain diagnostic requires
a full-spectrum convergence comparison before drawing model-wide conclusions.
Those ratios describe the earlier transition-rate implementation and are not
current validation values; the corresponding astrodust calculation must be
remeasured with the repaired rates before numerical claims are made about it.

`tmatrix/tests/test_rayleigh_zero_scattering.f90` checks the vacuum,
lossless, and absorbing Rayleigh limits and optional scattering moments. It
needs only the driver object, not `libtmatrix.a`:

```sh
cd tmatrix
make test_rayleigh_zero_scattering.x && ./test_rayleigh_zero_scattering.x
make clean && make BUILD=debug test_rayleigh_zero_scattering.x
```

Both versions passed in the release build and in the `BUILD=debug` build
(GNU `-fcheck=all` and invalid, zero, and overflow floating-point traps). Before the correction, a vacuum-contrast test produced
NaN albedo and phase-matrix coefficients from 0/0. The corrected convention
returns finite zeros when no scattering occurs; no optical tables were
regenerated.
