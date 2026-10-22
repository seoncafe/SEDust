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
2.0e-15 with both ifort and gfortran.

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
