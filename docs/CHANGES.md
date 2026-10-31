# CHANGES

A dated record of code and documentation changes, newest first. Each entry names
the files it touches and, where a change moves numbers, states the measurement
that was made on the tree itself.

## 2026-10-31

### Photon support and MC diagnostics

Construct the photon CDF on the input wavelength nodes and invert the
nonnegative linear number density within each positive-mass interval.
Zero-probability intervals are skipped. Exact plateau probabilities use an
explicit endpoint limit without zero divided by zero. For the default Mathis
field, photon integration begins at 0.0912 micrometers; the cooling grid is
retained. Supplied radiation fields retain their own positive support.

At radius 0.001 micrometers and U = 1, the original sampler drew 1,881 of
2,000,000 photons outside the field support, carrying 2.281094% of sampled
absorbed energy. The repaired sampler drew zero unsupported photons in
2,000,000 samples at each of five graphite radii. The sample mean agrees
with the exact integral of the represented photon density within 0.596
standard errors in all five cases. Independent tests with actual astrodust
and PAH cross sections, Mathis and HD23 fields, and three radii found zero
unsupported samples for both source versions; those paths did not have the
standalone zero-UV defect.

The observed maximum temperature in the 500,000-event, seed-101 graphite
calculation at 10 angstroms is 420.738401 K, compared with 7808.555649 K from the
original sampler. These are trajectory extrema, not histogram edges or
universal temperature limits. The full graphite size set and all dependent
MC emission and temperature-comparison products were recalculated in both
branches. All 43 tracked MC tables are byte-identical between branches.
The 2,913 emission-grain diagnostics have zero omitted and outside-grid time;
the raw radiative closure residual including stored heat spans -1.083501%
to +0.075948%, before spectral rescaling. Finite cooling steps and the
0.5 K temperature floor remain documented numerical approximations.

Use constants-module values for photon energies, the 150 eV equilibrium
gate, pi, and the optional D85 Boltzmann constant; derive the MC
Stefan--Boltzmann coefficient from those fundamental constants. No grain
changed equilibrium classification in the measured production cases.
`diag_Teq_vs_Trad.dat` has a header and a sixth population column, is
truncated at the start of each total solve, and has deterministic population
and radius ordering. Its plotting script reads the population column.

A partial `a_um_list` runs only the specified positive radii. Both drivers
reject malformed or nonphysical inputs; U = 0 remains valid for a CMB-only
field. Large emission/event-absorption ratios use exponent notation.
`make test_photon_support.x` and `make test_driver_inputs` provide targeted
support, mean-energy, boundary, and driver-input checks. The existing
supplied-field and cooling-interval tests passed in both branches.

Current measurements, preserved products, source and input hashes, and
installation receipts are retained in
`SEDust_JKAS/results/20261031_mc_photon_support` beside the branches.

### Astrodust report: removed passages restored with remeasured values

`docs/astrodust_sed_report.tex` regains the passages the final rewrite of
2026-10-31 removed: the measurement-status and energy-budget paragraphs of
the introduction, the execution-cost timings, the original-Mathis band
statistics, the scaling test, the jori=1 substitution test (avenue 3), the
induced-emission analysis with its two tables, the history of the
`tab:mathis_c` values, the first D03-sphere graphite-continuum table, the
mode-by-mode PAH Drude tuning record (`tab:pah_tune`), the library structure
note, and the conclusions with the diagnostic record and next steps. Every
restored number was remeasured on 2026-10-31 from `sed/output` and from
untracked runs of the settled source (`mathis_orig`, `induced`,
`mathis_orig induced`, `c2 induced`, `c2 mathis_orig induced`, `draine`,
`draine hd23`), which reproduce the tracked astrodust products byte for byte.
The jori=1 substitution table and its -15.1% to -11.3% result, the DL07 eq. 2
attribution, and the PAH tuning fit cannot be remeasured from this tree and
are marked with their dates. Corrected on the way: the gross/net emission
derivation (B J/B_env equals n_T J, not J), the independent Teq check
(remeasured: rms 0.0016% over 161 astrodust grains, not 0.14%), the
single-grain ratio statements (mean 0.99992; plateau 1/1.1866), the
T-window algorithm switch (run time, not compile time), the HD23-field
shoulder maximum (+27% at 35 um), the feature-band range (2-21%), the
`morig` file-name token, and the timings (Xeon Gold 6154, 2026-10-31).

### mc_pT_report figures and numbers remeasured from the final MC products

`docs/mc_pT_report.tex` now describes the current engine and quotes numbers
measured from the tracked `mc/output` and `sed/output` products and the run
logs of `SEDust_JKAS/results/20261031_final_recalculation`. The PAH and
large-grain single-grain comparisons, the cap-grain cold-end tally, and the
graphite wall times were measured with a separate build whose sources match
the settled manifest; their inputs and outputs are in
`SEDust_JKAS/results/20261031_mc_report_remeasure`. All six MC figures
(`fig24_5`, `fig24_6`, `sed_mc_vs_ref`, `pT_compare_astrodust`,
`pT_compare_pah`, `pT_compare_large`) were regenerated. Measured on
2026-10-31: default-mode MC/reference band ratios 0.961 (1-5 um), 0.998,
1.001, 1.001 (5-30, 30-300, 300-3000 um), and the all-MC mode within 0.0003
of them; T_rad/Teq of the 265 equilibrium-gate grains between 0.99871 and
1.00008; grain-by-grain closure between -1.084% and +0.076%, unchanged to
2.1e-5 between 5e4 and 5e5 events; no waiting interval needed a second step
chunk.

The report also records a defect of the standalone graphite setup, not
corrected because the numerical sources are settled: `invert_cdf` returns the
first grid wavelength (1e-3 um) for the zero-probability CDF node, so
`grain_setup_graphite` samples 0.097% of its photons below 0.0912 um, where
`J_Mathis` is zero. These photons (13.6 eV to 1.24 keV) raise the mean sampled
photon energy by 1.8% and produce the high-temperature tails of the
Figure 24.6 graphite distributions (10 A reaches 7800 K). The
`grain_setup_from_cabs` path of the SED builder and `main_pT_compare.x`
starts at 0.0912 um and is not affected.

`mc/scripts`: `run_fig24_5.sh` and `run_fig24_6.sh` use `set -eo pipefail`,
keep the full output of each run in `<tag>.log`, and pass the current
diagnostic lines (`emit/event_abs`, `mc_interval`, `mc_energy_single`,
`mc_warmup`, errors); the `run_fig24_6.sh` comment now states its 500000
events. The plot scripts take `--data`/`--out` arguments, render with
LaTeX, and no longer carry the 2 K grid, `N_burn`, or single-bin claims;
`plot_fig24_5.py` counts photon absorptions without the end-of-interval
record, `plot_fig24_6.py` shows the full temperature range, and
`plot_sed_compare.py` labels the reference as `calc_sed.x astrodust` (GD
matrix, `heuristic` narrowing).

### Settled-source data, emission, and manuscript recalculation

Rebuilt both branches outside the working trees from the 112-file settled
source manifest and regenerated the optical data affected by constant changes.
All six dust models passed `check_build_dust` in each branch. Astrodust, DL07,
MRN, and the Zubko D03 Mie options have changed numerical HDF5 datasets;
THEMIS and G18D remain numerically identical. The astrodust PAH maximum
relative absorption change is 6.03821748e-4 at 0.004259018 micrometers.
An isolated calculation changing only `hc_evum` reproduces that extremum.

Regenerated branch SED and MC products with the current sources. The 68 common
SED spectra are byte-identical between branches. Of 43 tracked MC tables,
42 are byte-identical; the remaining equilibrium diagnostic has identical
numeric rows and differs only in row order. Across the nine MC emission runs,
all 2,913 recorded grain diagnostics have zero omitted and outside-grid time.
The signed radiative closure residual including stored heat spans
-1.08352% to +0.075952% before spectral rescaling. The 0.5 K floor and
finite cooling-step error remain numerical limitations. The fifteen distinct
graphite calculations have a maximum absolute energy residual of 1.71728%.

Corrected the final-run validator to distinguish neutral and ionized PAH
population blocks when pairing warm-up and energy records. The six earlier
successful runs were validated and reused; no numerical sources were changed
for this correction. Regenerated the astrodust report figures and reconciled
its current input conditions and QM cooling-method counts.

The complete measurements, preserved products, installed-file hashes, and
manuscript checks are recorded in
`SEDust_JKAS/results/20261031_final_recalculation` beside the branches.
The MC temperature report is maintained separately.

## 2026-10-30

### QM radiative rates: positivity, cooling opacity, and actual-method counts

`sed/src/stoch_qm.f90` now interpolates opacity and radiation between two
neighboring wavelengths with nonnegative convex weights. The three-point
`PARAB_INTERPL` algorithm, from Aigen Li and B. T. Draine, could overshoot
below zero between positive opacity samples. A Zubko PAH heating rate of
-5.016826684e-11 s^-1 and cooling rate of -0.002970407157 s^-1 were traced to
specific negative interpolated cross sections, not floating-point round-off.
The original algorithm remains named `PARAB_INTERPL` for state-density
interpolation; its source is Draine's `BDraine_src/src/parab_interpl.f`,
which is not distributed with SEDust.

Cooling always uses the original opacity and field grid shared by the output
spectrum. `nisrf=N` reduces only the photon-absorption grid. Its field scale
preserves absorbed power using the exact integral of the product of the two
linear interpolants, including cross terms. Coarse cooling opacity missed
narrow PAH infrared features and caused a 3-4% emission excess in the smallest
grain; preserving absorbed power alone did not preserve its cooling power.

Rates are checked before the continuous-cooling collapse. Negative or
nonfinite physical rates return `QM_INVALID_RATE`; GTH no longer clips them.
The optional `raw_negative_count_out` and verbose diagnostic retain the count
across the attempted state grids. The optional `actual_cooling_out` returns
the successful method. Callers in `sed_astrodust.f90` count actual `dbcon`,
`dbdis`, and `stati` solves in the OpenMP summaries; `calc_sed.f90` records
those counts in each spectral header. Rejected QM attempts and equilibrium
calculations are not counted as successful QM methods. Numeric output columns
are unchanged. QM retains its own mode-counting enthalpy and accepts no
enthalpy table.

Measured with gfortran, 200 states, Zubko BARE-GR-S neutral PAH radius 3.5 A,
Mathis U=1 plus the hard EUV field, and `nisrf=150`:

| Cooling method | Emitted/absorbed before | Emitted/absorbed after |
| --- | ---: | ---: |
| `dbcon` | 1.03011961 | 0.99658405 |
| `dbdis` | 1.03075130 | 0.99876866 |
| `stati` | 1.04049800 | 1.00013548 |

The new `test_qm_rate_closure.x` Makefile target checks actual Zubko optical
inputs, raw rate signs, successful-method reporting, emitted/absorbed power,
and rejection of negative input opacity. Fifteen parameterized checks at
radii 3.5, 9.09832, 13.8236, and 25.8886 A passed with zero negative raw rates
and power discrepancies below 1%. The existing invalid-bin checks passed in
both branches. Builds and numerical comparisons ran outside the repositories;
tracked spectra and manuscript figures were not regenerated.

The long-wavelength PAH difference between `qm_dbcon` and discrete cooling
is a limit of the continuous-cooling approximation and is documented in
README_HOWTO, not removed by these repairs. The final comparison report and
source hashes are recorded with the read-only trace evidence.


Validation of the complete model spectra used identical before/after inputs.
With `nisrf=150`, the Zubko PAH 1-5 um after/before ratios were 0.963428
(`dbcon`), 0.962629 (`dbdis`), and 0.958863 (`stati`); the 300-3000 um ratios
were 1.002228, 1.002439, and 1.002411. Astrodust PAH and DL07 carbonaceous
5-30 um ratios were 1.003400 and 1.003373. Full-field runs and all component
bands are recorded in `SEDust_JKAS/results/20261030_qm_rate_repair/REPORT.md`
and its tables. All ten after runs recorded zero negative raw rates and
zero invalid-rate refusals. The fourteen complete spectra were plotted and
visually inspected; only diagnostic outputs were written.

### Energy-space (qm) solver: cooling quadrature, field grid, direct solve, refusals

All in `sed/src/stoch_qm.f90` (identical in v1.00 and v1.20), with the calls
in `sed/src/sed_astrodust.f90`.

1. Cooling-rate quadrature. `simpson_therm_cooling` and
   `simpson_intrabin_therm_cooling` returned zero whenever their first three
   samples were zero. That happened whenever the photon-energy interval was
   hundreds of kT wide (enthalpy bins wider than the photons that cool them),
   because those samples fall far outside the Planck peak at a few kT. The
   integral now stops at E1 + 60 kT and its first step is at most kT, so the
   peak is always sampled. The omitted tail is below 1e-17 for an integrand
   falling as x^m e^-x with m <= 7 (stated at the routine). A zero sum is
   accepted only on a grid finer than kT. The heating and exact-statistical
   integrators had the same early exit; they now accept a zero sum only once
   the step resolves the field grid (`photon_grid_spacing`).
2. Transition-matrix wavelength grid. The matrix now uses the full field grid
   by default (`qm_nisrf_max = 0`); before, it used 200 points picked evenly
   in grid index. An explicit `nisrf=N` coarsens the grid to N points equally
   spaced in ln(lambda), keeping the absorbed power (`coarsen_absorbed_field`;
   README_HOWTO describes the approximation).
3. Steady state. The BiCG solver and its sparse-storage routines are removed.
   The steady state is now found by GTH state reduction
   (`stationary_distribution_gth`), which adds only nonnegative terms. The
   back substitution is rescaled so that probabilities spanning more than the
   double range do not overflow. Cost, measured with gfortran -O3 (one
   thread): 0.002, 0.03, 0.22, 2.1 and 9.7 s, and 0.3, 2, 8, 32 and 72 MB of
   work array, at nstate = 200, 500, 1000, 2000 and 3000.
4. No unconverged solution is used. A grain whose UMIN/UMAX window is still
   moving after the last pass is refused; the window loop now allows 30 passes
   (it was 10, after which the last state was accepted as is). The NaN
   replacement by `pstate_gd`, which was always zero, is removed.
   `qm_solve_grain` returns a status (`atom_count_cap`, `no_steady_state`,
   `window_not_converged`, `energy_not_conserved`). Every outcome prints a
   diagnostic line, the atom-count cap included. The callers print
   `[QM refusals: ...]` by reason; the refused grains are still solved by
   the GD matrix.
5. `qm_solve_grain` no longer takes `nt_wide, t_wide, h_wide` or `teq`. None
   was read: the bin temperatures come from the solver's own mode-counting
   enthalpy (`T_from_U_bisect`, `enthalpy_DL01`), by design, and comments at
   that routine and at the bin construction now say so. With these arguments
   removed and nothing else changed, `calc_sed.x astrodust qm` and
   `calc_sed.x zubko qm` (one thread) gave byte-identical outputs in both
   branches.

Measured with gfortran 16.2.0 (OpenMP, 24 threads before and 36 after; the counts do not depend on it). Counts are
grains returned to GD out of the grains not in equilibrium.

| run | refused before | refused after (reason) | P_em/P_abs of solved grains |
|---|---|---|---|
| `astrodust qm` | 33 | 2 (atom_count_cap) | 0.9856-1.0003 |
| `astrodust qm nstate=500` | 10 | 2 (atom_count_cap) | 0.9950-1.0003 |
| `astrodust qm_stati nstate=500` | 10 | 2 (atom_count_cap) | 0.9950-1.0003 |
| `astrodust qm_dbcon nstate=500` | 10 | 2 (atom_count_cap) | 0.9919-1.0002 |
| `dl07 mw31_60 qm` | 14 | 1 (atom_count_cap) | 0.9863-1.0003 |
| `dl07 mw31_60 qm_stati mathis_orig nstate=150` | 20 | 2 (atom_count_cap) | 0.9806-1.0022 |
| `zubko qm` | 10 | 2 (atom_count_cap) | 0.9915-1.0027 |
| `zubko euv hardfield qm_dbcon nisrf=150` | 49 | 1 (atom_count_cap) | 0.9930-1.0313 |

Band powers qm/heuristic after the change, `astrodust qm`, for the bands
1-5 / 5-30 / 30-300 / 300-3000 um:
- S1: 0.9988 / 1.0009 / 0.9996 / 1.0002;
- S2: 0.7380 / 0.9015 / 0.9994 / 1.0097;
- PAH: 0.9928 / 0.9989 / 0.9945 / 1.0006.

The S2 difference was 0.7392 / 0.9026 / 0.9998 / 1.0097 before. It is the
intended difference between the solver's own mode enthalpy and the model's
S2 H(T), not a solver error. The tracked `sed/output` qm products are not
regenerated by this entry. v1.00 and v1.20 give byte-identical outputs for
all 23 files of these runs.

### Exact SI physical constants in one module; QM fallback counts

`sed/src/constants.f90` now holds h, c, k_B and the electron-volt as the exact
SI defining values (h = 6.62607015e-27 erg s, c = 2.99792458e10 cm/s,
k_B = 1.380649e-16 erg/K, 1 eV = 1.602176634e-12 erg), their SI-unit forms
(`h_si`, `c_si`, `k_si`), and the derived `hc`, `hc_ergum`, `hc_evum`
(1.2398419843 eV um) and `hc_k` (1.4387768775 cm K). The modules that carried
their own rounded copies now `use` these: `p_sub.f90` (h of CODATA 2010),
`radfield.f90` (h, k_B of CODATA 2010), `sed_astrodust.f90` (h of CODATA 2010;
1 eV = 1.60218e-12 erg; hc = 1.98645e-16 erg cm), `stoch_qm.f90` (h = 6.62607e-27,
c = 2.99792e10, k_B = 1.38065e-16, pi = 3.141593, 1 eV = 1.60218e-12 erg,
hc/k = 1.43877), `enthalpy_v2.f90` (k_B of CODATA 2010 and hc/k = 1.4387687,
a wrong value; the correct one is 1.4387769), `calc_sed.f90` (1 eV of CODATA
2010, hc = 1.23984193 eV um), and the dielectric-table readers
`q_astrodust.f90`, `q_silicate.f90`, `q_graphite.f90` (hc = 1.23984 eV um).
The tests `test_planck_grid`, `test_radiation_support`,
`test_absorption_balance`, `test_hd23_field` and `validation_driver` take the
same constants, and `pyutil/sed_from_cabs.py` mirrors them. `pah_ioniz.f90`
(a port of Draine's grain-charging code) and `tmatrix/` are unchanged.

Measured with gfortran 16.2.0, OpenMP 6 threads, against the tracked
`sed/output` files (which the unchanged code reproduces byte for byte):
`calc_sed.x astrodust` changes by at most 1.4e-5 at a point (PAH, 0.40 um;
points above 1e-6 of the peak) and 2e-7 in the total; `calc_sed.x dl07
mw31_60` by at most 6.9e-6 at a point and 2e-7 in the total;
`calc_sed.x astrodust qm` by at most 3.2e-3 at a point (S2, 48 um), with
band changes of at most 3.8e-4 (1-5 um). The tracked `sed/output` files are
not regenerated by this entry. v1.00 and v1.20 give byte-identical outputs.

The dl07 and mrn EUV optics are computed through `q_silicate`/`q_graphite`,
whose wavelength axis moved by 1.6e-6; `check_build_dust.x` now finds the
stored `data/dl07` and `data/mrn` EUV products off by up to 1.1e-3 (C_abs) and
5e-3 (<cos>), and fails those 8 checks until the products are regenerated
(`calc_qtable.x`, then `calc_kext.x <model> euv`). Their non-EUV views and all
other models pass.

`check_radiation_support.sh` passes; its worst closure is now 2.2e-15 with
ifort (2.0e-15 before) and 2.0e-15 with gfortran. `test_planck_grid` gives
`BIT_DIFFERENCES=0`; `test_cmb_field_balance` closures are -7.8e-16 and
-3.3e-16 (8.9e-16 and -1.0e-15 before). `sed/tests/README.md` is updated.

`stoch_qm.f90`'s single-grain diagnostic line printed `natom` with `i5`,
which gave `*****` for every grain above 1e5 atoms; the integers are now
written with `i0`. The QM summary lines in `sed_astrodust.f90` now report what
solved each grain: `[QM batch: stoch= eeq_gate= qm_fail_to_gd= qm_fail_to_eq=]`.
The batch line did not count QM refusals at all, and the `qm` line of `sed_grain_loop`
counted a refusal solved by the GD narrowing as `fail_to_eq`. For
`calc_sed.x astrodust qm` the batch now reports 33 refusals solved by GD
(31 with the old constants, from the old line's stoch and equil counts); the
stderr diagnostic prints 2 fewer `ok=F` lines, because the atom-count cap
(`NATOM_QM_MAX`) returns before it (two carbonaceous grains of 0.0133 and
0.0141 um).

### P(T) comparison energy diagnostic label

`mc/src/main_pT_compare.f90` now labels its printed ratio
`emit/event_abs`. The denominator counts discrete photon absorption only;
continuous heating is excluded. This ratio is not a total energy closure
test. The arithmetic and saved P(T) products are unchanged.

### SED products regenerated without the stored CMB cooling term

The tracked `sed/output/sed_*.dat` products were made before `calc_P` stopped
subtracting the stored `kappCMB` integral from the cooling rate (the CMB is
already in `J_Mathis` and `J_HD23`). They were regenerated with the current
`calc_sed.x` (gfortran, OpenMP 8 threads), one run for each tracked file, with
the arguments read back from the filename tag. Every regenerated file keeps its
name, wavelength grid and columns; the column header now reads
`erg s^-1 sr^-1 H^-1`. All new spectra are identical, data line for data line,
to the committed code with only the `kappCMB` subtraction removed, so the
changes below are this correction alone. All 71 files are installed. Largest
point change (points above 1e-12 of the column peak): 7.5e-4
(`sed_dl07_mw31_60_morig.dat`, silicate column, 63.8 um). Largest change of
the total column integrated over 1-5, 5-30, 30-300 and 300-3000 um:
-4.1e-5, +4.1e-5, -8.6e-6 and -1.0e-4. Equilibrium runs are unchanged.
The new spectra are byte-identical in v1.00 and v1.20 wherever both
branches track the file. `sed_astrodust_hd23_{S1,S2,PAH}.dat` were already
current and are unchanged.

Four files change more than the 1e-4 level of the other products: `sed_g18d.dat` and `sed_g18d_morig.dat`
(PAH column -2.7% and -3.6% near 2 mm; total over 300-3000 um -3.8e-5 and
-5.2e-5), `sed_zubko_euv_hardfield_qm_dbcon_nisrf150.dat` (PAH hot tail
-13.7% at 0.084 um, where the PAH column is 1.3e-12 of its peak; the
smallest PAHs fail the QM solve under the hard field and are solved by the
GD matrix), and `sed_zubko_draine_morig_logU-.50_induced.dat` (PAH column
-1.3e-3 at 4.7 mm). These changes are also the `kappCMB` correction alone. The large relative
changes occur where the column is a small fraction of its peak (about 1e-5
of the PAH peak near 2 mm for G18D, 1.3e-12 at 0.084 um for the hard field).

`docs/astrodust_sed_report.tex` numbers and figures were recomputed from the
regenerated products (`docs/make_figs.py`, all eight figures redrawn); no
value printed by `make_figs.py` moved at the precision the report quotes, the
largest change being the maximum S1_C1 30-300 um residual, 23.33% -> 23.34%
(quoted as 23.3%), and the energy-budget and `hd23` tables are unchanged. The
induced-emission comparisons and the x1.186 scaling test, which the report had
kept from older code states, were re-measured with the current `calc_sed.x`
(FIR baseline of the induced-vs-HD23 table -15.83% -> -14.18%, scaled FIR
residual -0.2% -> +1.8%), and its QM-solver paragraph now gives the measured
band ratios (Stage S2 26% below the heuristic solver in 1-5 um) and the 29
grains that fall back to the GD matrix.

### Stationary-start comparison under prior constants

The preserved set `SEDust_JKAS/results/20261030_mc_stationary` completed
12 SED/P(T) calculations and 15 graphite calculations using the source
snapshot preceding the concurrent physical-constant conversion. It is not
a set of current-tree products. The 14 tracked MC tables and manuscript
selection remain unchanged pending the explicit numerical-sources-settled
notice and a final rebuild.

All 2913 SED grain trajectories have zero omitted and out-of-grid time.
The 2940 warm records contain actual durations and event counts; realized
counts range from 76 to 50583. Stored/absorbed fractions span
-6.02975e-5 to 1.29309e-4, and signed raw closure residuals span
-1.08352% to +0.0759518%. Across the saved component/band comparisons,
integrated-power ratios range from 0.998900412 to 1.002057569. These changes
include finite Monte Carlo fluctuations because preparation changes the
recorded random sequence.

The 200/100/50 angstrom graphite first occupied bin edges are
16.2743/9.69634/4.28594 K; these are sample-dependent bin edges, not
physical minimum temperatures. The 200 angstrom probability below
15 K changes from 2.01842e-5 to zero in this finite sample. For
1000 independent 10 angstrom starts, the bin-edge CDF difference from
a 500000-event reference is 0.0385893, within the start-only 99% DKW
half-width of 0.0514700. The reference has separate finite sampling
uncertainty. Individual warm times, full profiles, spectral bands, and
diagnostic limits are in the preserved set's `MEASUREMENTS.md`.

### Deterministic stationary Monte Carlo recording

`mc/src/mc_engine.f90` now prepares each simulated grain from the radiative
equilibrium temperature of the configured grain and supplied radiation field.
For a positive discrete absorption rate, the preparation duration is the fixed
`T_w = max(10 * tau_relax, 100 / rate_event)`; with zero event rate, it is
`10 * tau_relax`. The second term gives 100 expected events, not a required
number of realized events. The logged event count is the actual number before
`T_w`. Because the sampled event rate is independent of temperature, photon
arrivals form a homogeneous Poisson process. The preparation integrates the
current waiting interval only to `T_w`, without forcing an event at that time.
At positive rate, the recorded engine samples a fresh exponential residual
waiting time, consistent with Poisson memorylessness. The fixed duration
avoids selecting a recording start solely because a photon arrived.

`mc/src/main_mc.f90`, `mc/src/main_mc_sed.f90`,
`mc/src/main_pT_compare.f90`, and `mc/src/mc_sed.f90` start recording
from the prepared temperature and random state. Preparation contributes no
histogram residence, elapsed time, event absorption, emission,
or stored-energy change to recorded outputs. `N_events` and `t_max` constrain
the recorded interval only. An explicit positive `T_init` in `main_mc.x`
requests a transient experiment and bypasses preparation; the default is the
stationary procedure. The comparison driver starts all three engines from the
same prepared state. The `[mc_warmup]` record reports initial equilibrium
temperature, recording temperature, preparation time, and actual event count;
`mc_sed` also reports composition.

The buffered option measures its complete sample count with a two-pass replay
from that same state, then allocates the measured capacity. The limit is
1,000,000,000 samples, approximately 16 GB for the two real sample arrays.
The OpenMP SED builder serializes capacity measurement and buffered execution
to bound their aggregate buffer allocation. The drivers reject insufficient
capacity and replay-count mismatches. Input file open and malformed-namelist
failures now report an error and exit with nonzero status; `main_mc.x` also
rejects an unknown engine name. `docs/mc_pT_report.tex` and `README_HOWTO.md`
describe the current preparation and recording boundaries.

The fresh numerical production under this deterministic start is in progress.
No prior numerical product is attributed to this procedure. The entries below
retain the methods and measured values of their earlier checkpoints, including
CMB starts and runs without an excluded preparation interval where stated.

### Earlier completed MC products and manuscript verification

The earlier completed calculation set was
`SEDust_JKAS/results/20261030_mc_energy_consistent`. Its numbers were
measured before the deterministic stationary start above and are historical
for that method; they are not current stationary results.
Twelve SED/distribution runs and fifteen graphite runs completed in external
GNU Fortran builds. All accepted MC intervals have zero omitted and
out-of-grid time. The two version trees produce identical matched products.
The then-refreshed default tracked SED had final/preserved integrated powers
1.000007013, 1.000010176, 1.000158469, and 0.998882880 in 1-5, 5-30,
30-300, and 300-3000 um. Preserved tracked files lack full provenance;
these are product differences, not an isolated estimate of one correction.
The paper default ratios from that calculation were 1.000006503,
1.000012414, 1.000147607, and 0.998794121. Full temperature profiles and CDFs are retained. Their
largest changes include missing cold probability and different bin edges:
0.890341078 (fixed) and 0.608986554 (two-pass/buffered) for the paper.

Every simulated grain has a stored-energy diagnostic. Across paper runs,
maximum absolute stored/absorbed fractions are 0.000262176 (50,000 events)
and 0.0000308162 (500,000 events). The largest astrodust grain, initialized
at equilibrium, has stored/absorbed = -0.00000728550 in the paper all-MC
run, compared with 0.0903643 in the controlled CMB-start experiment.
The signed raw closure residual after subtracting stored heat spans
-0.0107978 to 0.000686893 in paper runs; spectral rescaling does not
remove this integration error. A matched 5000-event test gives residuals
-0.696281%, -0.221798%, and -0.056514% for maximum fractional temperature
steps 0.05, 0.025, and 0.01. The fast default remains 0.05.

Trajectory recording now returns an explicit capacity failure instead of
silently truncating the history. The single-grain driver reserves space
for event records and time samples; recording-disabled trajectories are
unchanged by this check. Focused tests verify overflow, reset, and identical
RNG, elapsed time, and energy when recording is disabled. All fifteen
refreshed graphite calculations retain raw energy accounting; the largest
absolute closure residual is 0.563324% of total absorption. The manuscript's
10 angstrom graphite example uses the measured 468.208 carbon atoms and
367.879 K after a 10 eV photon from 2.725 K. Manuscript text describes the
current methods, while these change records preserve diagnostic comparisons.

### Complete Monte Carlo waiting intervals and temperature histograms

The adaptive wavelength split uses the same `hc_erg_um` as the photon-event
energies. Its value is derived from the exact SI Planck constant and the speed
of light, converted to erg um: 1.9864458571489286e-12. The local constants in
`mc_sed` and `main_pT_compare` were 1.98644582e-13, a factor of ten too small
for the stated 1% temperature-jump threshold. They are removed. The engine's
Planck constant is also updated from 6.62606957e-27 to 6.62607015e-27 erg s;
this changes event energies by 8.7533e-8 in relative terms. A unit test checks
that a 1 um photon has energy 1.239841984332 eV. The default temperature-jump
threshold remains 0.01; its dimensional conversion is consistent with the
sampled photons. This correction changes the physical division between event
and continuous heating and requires affected MC products to be recalculated.

The Monte Carlo engines complete each sampled waiting interval before applying
its photon event. The fixed-grid and adaptive engines use chunks of 4000 and
300 steps, respectively, with a finite budget of 1,000,000 steps for one
interval. A nonfinite or non-advancing step, exhausted budget, temperature
outside the histogram grid, or exhausted sample capacity returns an explicit
status. Drivers stop before accepting a spectrum or grain comparison.
`cool_segment` reports the time actually integrated. No environment variables
or new run-input format are introduced.

All completed intervals, including the initial trajectory, contribute to the
histogram, elapsed time, absorbed event energy, and emitted energy. The
obsolete `N_burn` argument is removed from `mc_run_engine_2pass` and
`mc_run_engine_buffered`; source callers must remove that positional argument.
The two-pass probe covers all temperatures that the replay bins. The fixed
histogram covers 0.5–5000 K. Its column structure is unchanged, but its bin
edges and therefore the interpretation of saved MC probabilities change.
Residence outside the grid is counted and fails explicitly; it is not moved
to an edge bin. Buffered calculations require storage for the complete
trajectory and fail explicitly on capacity exhaustion. `main_pT_compare`
measures the required capacity from its identical-seed two-pass calculation
and checks the buffered count against it, with a limit of 1,000,000,000
samples (16 GB for the two real sample arrays).

The optional engine outputs `floor_time_out`, `out_of_grid_time_out`,
`omitted_time_out`, and `sample_count_out` separate numerical boundaries from
elapsed time. `status_out` reports `MC_OK`, `MC_STEP_BUDGET`,
`MC_INVALID_STEP`, `MC_INVALID_WAIT`, `MC_HISTOGRAM_RANGE`, or
`MC_SAMPLE_CAPACITY`. Omitting `status_out` causes an unsuccessful engine
call to stop. `mc_sed` reports the independently integrated prescribed-field
absorption, sampled event absorption, continuous heating, emitted power,
elapsed time, omitted time, floor residence, and histogram mass for each
simulated grain. Its bolometric rescale is not an acceptance test for either
the trajectory or the wavelength distribution.

At the imposed 0.5 K floor, a derivative directed toward lower temperature
leaves the grain at the floor for the rest of the waiting interval, including
its full residence and emission. This is an unresolved-cold-tail approximation,
not a continuous-heating equilibrium or a CMB temperature boundary. In the
controlled interval-only diagnostic at U=1.585, retaining the original
wavelength split to isolate time integration, the HD23 S2 grains at
0.000354813389233575 and 0.000375837404288444 um have mathematical
continuous-heating roots of 0.199355 and 0.258852 K under the implemented
QPL extrapolation. More than 99.998% of their CMB absorption power belongs
to discrete events. QPL is held at its 1 K table value below that temperature;
the cold probability density is not physically validated by this approximation.

Direct matched-grain calculations used 50,000 events, seeds 13354 and 14363,
and the same S2 optics, radiation field, and preserved original cutoff.
These controlled values isolate interval integration; they do not describe
the final corrected-cutoff production settings. The omitted time in the replay
was 751746640.014091 and 43451931.8374593 s. Completing the intervals
restored those times while keeping the absorbed event energy and final RNG
state identical. With the complete-histogram implementation, omitted and
out-of-grid time are zero, probability mass differs from unity by less than
1.3e-14, and the exact 0.5 K residence fractions are 0.4843276265 and
0.4015012335. The required sample counts are 61,208,568 and 59,317,056.
The full CDF changes relative to the truncated calculations are 0.00032106026
and 0.000025638214. The measured emission from bins below 1 K is less than
8e-9 of total emission for these two cases; this does not validate their cold
probability distributions or establish the approximation for other inputs.
Independent trajectory emission still exceeds absorbed energy by about 0.7%
in these cases, so exact energy conservation is not claimed.

With the physically consistent cutoff and exact Planck constant, both
matched grains use lambda_c=39,360 um. Their continuous-heating terms are
5.4944010e-19 and 5.8199754e-19 erg cm^-2 s^-1; the formal roots of the
constant-QPL continuation are about 0.0297602 K, outside its physical
validation range. Exact floor residence fractions in the final 50,000-event
runs are 0.4868630592 and 0.4091529688. Omitted and out-of-grid time are
zero, probability mass is within 3e-14 of unity, and sample counts are
61,210,768 and 59,324,376. Full CDF changes from the original truncated,
incorrect-cutoff calculations are 0.0028159573 and 0.0075563804. Separate
runs with the corrected cutoff but the truncating engine omit
8.9003996997e8 and 5.73672694165e7 s in replay; completing the same
intervals restores those times. Final temperature bins below 1 K contribute
7.17293e-9 and 7.79503e-9 of histogram-estimated emitted power. These
fractions describe the direct cold emission, not a bound on its influence
on later heating excursions. Full inputs and profiles are recorded in
`SEDust_JKAS/results/20261030_mc_interval_diagnostics/correct_hc/`.

Focused tests in both versions cover chunk-size independence, all three
engines, pure cooling, continuous heating, floor residence, long waits,
finite observation time, invalid steps and waiting rates, budget exhaustion,
truncated-probe range errors, and buffered capacity errors. Tests and builds
run outside the source repositories. The matched-grain measurements and full
profiles are in `SEDust_JKAS/results/20261030_mc_interval_diagnostics/`.

### Recalculated MC spectra and complete temperature distributions

The following measurements describe the preserved
`20261030_mc_physical_cutoff` diagnostic set. A subsequent independent energy
check found that interpolation in the photon enthalpy inverse loses energy
for small photon increments. These numbers are retained as a diagnostic
record, not accepted final results. The analytic-inverse calculations are
recorded separately below when complete.

External GNU builds completed the three SED calculations and four-grain
P(T) comparison for each version and for the paper inputs (12 driver runs).
All accepted grain trajectories have zero omitted and out-of-grid time.
The two version outputs agree for matching inputs. The tracked comparisons
use the unchanged Mathis-field namelists; the paper uses piecewise HD23
at U=10^0.2. These input sets must not be conflated.

For the tracked default, high-count, and all-MC spectra, respectively,
the 30–300 um power changes against the preserved tables are +0.015896%,
+0.009812%, and +0.014894%; the 300–3000 um changes are -0.111246%,
-0.110643%, and -0.097555%. For the matched paper comparison, the
high-count spectrum changes by +0.007532% and -0.119996% in those bands.
The paper's two-pass P(T) CDF changes by as much as 0.608986 across the
four grains, at 10 Angstrom; the tracked set's maximum is 0.649580.
This is distinct from the smaller matched cap-grain CDF changes above.
The fixed-grid cumulative-mass changes reach 0.890341 and 0.906744,
respectively, and include recovery of previously unrecorded cold residence.
The preserved fixed-grid distributions are not silently renormalized.
The saved column structure stays unchanged while the fixed bin edges change.

In that earlier calculation, the two-pass and buffered distributions agreed for the same seeds and
complete sample counts. Full density/CDF comparisons, input conditions,
independent raw absorption/emission diagnostics, and component band changes
are recorded in
`SEDust_JKAS/results/20261030_mc_physical_cutoff/MEASUREMENTS.md`.
Previously tracked tables lack complete input/source metadata, so their
changes are reported as a product refresh rather than an isolated estimate
of one correction. Earlier calculation directories are preserved.

### Analytic photon enthalpy inversion and initial thermal energy

Photon absorption in the MC engines now inverts the same analytic U(T) used
before absorption. The tabulated enthalpy locates a bracket; safeguarded
Newton iteration with bisection satisfies the analytic energy residual.
Below the table the bracket extends to the imposed 0.5 K floor, and above
the table it expands instead of silently capping temperature. Invalid or
unresolved targets return MC_INVALID_STEP before the event is accepted.
When the bracket endpoints are adjacent representable temperatures, the
endpoint with the smaller enthalpy residual is accepted only if the relative
residual is at most 1e-12. All other targets retain the stricter residual
criterion. A positive-energy photon cannot cool the grain; roundoff-sized
inverse differences are handled before updating the state.

The defect was measured with the preserved HD23 all-MC 5.011872336 um
astrodust calculation. The cumulative photon-jump loss was
2.65530331314e-4 erg. Stored heat increased by 2.69460436767e-4 erg, while
absorbed minus emitted energy was 5.34990014640e-4 erg. Thus finite initial
heating alone did not explain the deficit. The prior inverse interpolated
log(T) linearly in U even though each event began from analytic U(T), and
small photon increments repeatedly lost energy through that mismatch.
The matched Mathis calculation had a cumulative loss of 2.54406670312e-4 erg.
The separate signed closure residual after accounting for these terms was
below 8e-10 erg in both diagnostic runs.

At this earlier checkpoint, all-MC grains above the equilibrium-energy
threshold started at radiative equilibrium, while smaller grains started at
the CMB temperature. No preparation interval was excluded at that stage:
residence, time, absorption, and emission used the same intervals. The optional
stored_energy_out diagnostic returns analytic U(T_final)-U(T_initial);
mc_sed appends stored_energy_erg to its independent diagnostic record.
The spectrum and P(T) data-column structure is unchanged.

Focused tests cover the analytic inverse over 0.5–12000 K, multiple materials
and radii, zero and small photon increments, the observed adjacent-temperature
case, invalid targets, and stored-energy consistency between engines and
chunk sizes. Full measured results for the refreshed calculations are recorded
in SEDust_JKAS/results/20261030_mc_energy_consistent/ when validated.

### Explicit resolution controls and the zero-scattering Rayleigh limit

The three-dimensional radiative-transfer setting retains the fast `heuristic`
solver and a 200-temperature-point model. The energy-space `qm` method retains
200 enthalpy bins and a 200-wavelength radiation-field cap by default. Existing
`nstate=N` and `nisrf=N` CLI settings permit a deliberate resolution study;
for example, `astrodust qm nstate=500 nisrf=1129 hd23 logU=0.2` uses all
1129 source wavelengths when that grid is available. The effective count is
the lesser of the source-grid count and `qm_nisrf_max`. Library callers may set
the public `qm_nstate_default` and `qm_nisrf_max` module variables before
parallel solves and must not change them concurrently. A larger setting does
not establish spectral convergence by itself; compare the resulting spectra
and energy balance at successive resolutions.

A focused diagnostic used one neutral PAH at `a=0.000421697 um`, 35 atoms,
`NT=200` over 2.7–5000 K, HD23 `U=10^0.2`, and thermal-discrete QM with
`nstate=500`. Its emitted-to-absorbed power ratios for 200, 500, and 1129
radiation wavelengths were 0.8439877849 (rejected solution), 0.9683433756,
and 0.9962086841. The corresponding ionized PAH at 200 wavelengths gave
1.0068300802. This isolates a sampling problem for one grain and charge
state; it is not an accuracy estimate for a complete dust spectrum. The fast
default and all existing reference products remain unchanged by the explicit
resolution option.

The Rayleigh dipole limit now returns zero albedo and zero phase-matrix
moments when the refractive index has no contrast with vacuum and both
extinction and scattering vanish. Previously the normalization divided by
zero and produced NaNs. Focused checked GNU tests in both versions passed
vacuum sphere, oblate, and prolate cases, lossless and absorbing scattering,
the spherical Rayleigh formula, and optional moments. No optical tables were
regenerated.

### Explicit HD23 radiation field and complete CMB accounting

`sed/src/radfield.f90` implements the final journal version of Hensley &
Draine (2023, ApJ 948, 55), equations (9) and (10). The three diluted
blackbodies have `(T, W)` values `(3000 K, 7.0e-13)`,
`(4000 K, 1.65e-13)`, and `(7500 K, 1.0e-14)`. With wavelength in um,
the UV `lambda*u_lambda` terms in erg cm^-3 are `1.287e-9*lambda^4.4172`
for 0.0912–0.110 um, `6.825e-13*lambda` for 0.110–0.134 um, and
`2.373e-14*lambda^-0.6678` for 0.134–0.246 um. The code converts
`u_lambda` to `J_lambda` in W m^-3 sr^-1, scales stellar light by `U`,
and adds the unscaled 2.725 K CMB. The field is piecewise, as in Draine's
`mmp_radiation_field.f` and in `J_Mathis`: the UV power laws up to and
including 0.246 um, the three blackbodies above it. (The first version of
this routine, written the same day, added the blackbodies at every
wavelength because eq. 9 is printed as a sum; the review of 2026-10-30
replaced it, since the MMP83 UV fits already hold the stellar Wien tails,
and the release energy budget closes only with the piecewise form: emitted
over absorbed power of the HD23 release at log U = 0.2 is 1.000 (astrodust)
and 1.015 (PAH) with it, against 0.974 and 0.917 with the sum.)
The `J_HD23(U, lambda, J)` routine is public.

`calc_sed.x` accepts `hd23` as an explicit field setting for every dust
model and rejects its combination with `mathis_orig`. A default run still
uses the corrected Mathis field; `mathis_orig` retains the published 1983
convention. `hd23` adds a deterministic filename tag. Each SED run also
writes the spectrum supplied to the solver as
`output/field_<model>[_<tag>].dat`, including any `hardfield` replacement.
The emission headers identify the field and `U`; units of `lambda*I_lambda/N_H`
are corrected to erg s^-1 sr^-1 H^-1.

The Monte Carlo input files now accept `radiation_field = 'mathis'`,
`'mathis_orig'`, or `'hd23'`. The supplied field is used by both the
stochastic grains and the equilibrium grains, fixing the path that rebuilt
a Mathis field after the caller selected another spectrum. The MC drivers
write the field they supplied to `*_field.dat`. The complete field includes
CMB heating, so the GD cooling calculation does not subtract the stored
partial CMB integral a second time.

Fresh external GNU builds of both versions used `WITH_TMATRIX=1`,
`FC=/opt/miniconda3/bin/gfortran`, and
`HDF5_PREFIX=/data/opt/hdf5_gfortran`. The focused HD23 test passed its
UV boundaries, unscaled CMB, `U` scaling, and both Mathis conventions.
The zero-field population check passed. The measured relative finite-bin
energy-closure residuals were `8.8818e-16` for pure CMB and `-9.9920e-16`
for a field without CMB in each version. The subsequent manuscript calculations regenerated the affected SED and MC
figures under `../SEDust_JKAS/results/20261030_hd23`. That directory records
actual fields, run settings, frozen sources, checks, and reference limitations.
Older benchmark outputs elsewhere remain historical.

### Energy budget conclusion withdrawn

`docs/astrodust_sed_report.tex` (and its PDF) no longer states that the HD23
release fails to conserve energy. That conclusion, reached on 2026-06-12,
used our default `mathis` field, which follows Draine's 2008
`mmp_radiation_field.f`: it corrects the 4000 K dilution factor to 1.65e-13
but keeps the Mathis (1983) `W = 4e-13` at 3000 K. HD23 eq. (9) uses the
Draine (2011) value `7e-13` (see "Explicit HD23 radiation field" above for
the `J_HD23` routine and the `hd23` setting). The 2026-06-12 argument that a
field difference would shift astrodust and PAH equally holds only for a
uniform change of `U`; the 3000 K factor changes the field shape, and raises
the absorbed power of the release `C_abs` by 9.3% (astrodust) and 1.8% (PAH).

Measured on 2026-10-30 from the release tables (`C_abs` = extinction -
scattering, power-law extrapolation over 0.0912-0.1 um, emission
`4 pi int lambda I_lambda dln lambda` of the `irem` column, log U = 0.2),
emitted over absorbed power:

| field | astrodust | PAH |
|---|---|---|
| `mathis` (W3000 = 4e-13) | 1.093 | 1.033 |
| `hd23`, piecewise (W3000 = 7e-13) | 1.000 | 1.015 |
| eq. (9) as a literal sum | 0.974 | 0.917 |

Our SED with the HD23 field (`./calc_sed.x astrodust hd23`, gfortran,
NT = 200, heuristic; outputs `sed/output/sed_astrodust_hd23_{S1,S2,PAH}.dat`)
against the release, median of (ours - HD23)/HD23 with the band power ratio
in parentheses (`docs/make_figs.py` conventions):

| band (um) | astrodust S1 | astrodust S2 | PAH |
|---|---|---|---|
| 1-5 | +71.24% (1.387) | +164.35% (1.879) | +14.75% (1.022) |
| 5-30 | +3.70% (1.043) | +10.97% (1.156) | +4.34% (0.969) |
| 30-300 | -0.13% (0.999) | -0.37% (0.999) | +0.32% (1.001) |
| 300-3000 | -0.42% (0.994) | -1.06% (0.986) | +3.57% (1.021) |

The astrodust S1 residual at the 133.4 um peak is -0.69% (-9.80% with the
default field), and the 0.1-3000 um power ratio is 0.999 (astrodust) and
0.982 (PAH), against 0.914 and 0.965 with the default field. The 1-5 um
astrodust excess is the same with either field and remains open. The report
now carries both fields in Table `tab:ad_energy`, adds Table
`tab:hd23_bands`, retitles the section `sec:phaseE-resolution` ("the 3000 K
dilution factor of the heating field"), and rewrites the introduction
summary, the conclusions table and the closing discussion accordingly. The
default `mathis` numbers elsewhere in the report are kept and labeled as such.

### Count the illuminating CMB once

The GD callers in `sed_astrodust.f90` and `main_pT_compare.f90` no longer
subtract the stored CMB integral from cooling. `J_lam` already contains any
background radiation, so the previous subtraction duplicated that heating.
This corrects heuristic, iterative, and unsuccessful QM solves that use GD.
The public partial `kappCMB` arrays remain available, but production emission
does not consume them. The optional low-level `calc_P` argument is documented
as a separate continuous-heating term only when excluded from `Jfield`.

Added focused full-field closure and population-invariance tests. In fresh
GNU builds of both versions, the pure-CMB closure residual was 8.8818e-16 and the
CMB-free-field residual was -9.9920e-16. Both emission methods
were unchanged when the stored partial CMB arrays were perturbed. Existing
spectra and quoted benchmark values must be remeasured before serving as
validation of the corrected solver.

## 2026-10-22

### Intel Fortran builds

`sed/Makefile`, `mc/Makefile`: compiler family detection (`ifort`, `ifx`,
`mpiifort`, `mpiifx`) with the Intel flags of `sed/build_lib.sh`.
`mc/src/main_mc_sed.f90`, `main_pT_compare.f90`: `use omp_lib` (IFX provides no
intrinsic `omp_lib`). `mc/src/mc_sed.f90`: the BLOCK-local variables of the
parallel region are declared at routine scope and listed in `private`.
`sed/src/check_build_dust.f90`: the ZDA round-trip check allows two ulps,
because IFX 2026.1 reads 2 of the 145321 rows of `suvSil_121_1201.dat` one ulp
from the correctly rounded value. With IFX 2026.1.0 all `sed/` targets and `mc/`
build and the read-only checks pass.

### Other corrections

`sed/src/p_sub.f90` (comments): the bin representative energy H_u versus DL01
Sec. 4.1 and its validity range. `pyutil/sed_from_cabs.py`: `calc_P` and the
absorption integral follow `p_sub.f90` (192 cases: upward rates bit-identical,
P to 2.5e-14). `docs/SEDust_user_manual.tex`: the 1e-12 support-floor paragraph
no longer places 1e-10 below 1e-13. `sed/tests/*.sh` are executable.

### Port the v1.00 emission changes to v1.20

Transfer `5f68a42..710503b` and the current v1.00 working-tree updates,
including the bin-representative-energy qualification in `p_sub.f90` and
executable modes for the two new test scripts. Preserve the v1.20 polarized
optics, `Cpol` and `lamI_pol`, the aligned-scattering interfaces, and their
build targets.

`sed_astrodust.f90` was merged by hand where the v1.00 patch did not apply:
the emission-bin cutoff skip sits before both the total and the polarized
accumulation in `sed_grain_loop`, the prepared Planck grid is in the `draine`
and `qm` OpenMP shared lists beside the polarized variables, an unreadable
D03 file in `sed_init_dl07` returns status 7 (v1.20's code for the D03
refusal; 4 in v1.00), and `error_unit` stays imported. `sed/build_lib.sh`
now lists `dustem_io`, which `sed_astrodust.f90` uses and which the script
lacked here as it did in v1.00.

Validation on v1.20, gfortran, outside the repository: all 18 `sed/`
Makefile targets and `mc/` build, with no compiler warnings in the merged
modules; `test_scalar_mode`, `test_euv_extension`, `test_scatmat_aligned`,
`test_dustem_extinction`, `test_dustem_polarized_extinction`,
`check_build_dust`, `test_h5_storage_policy`, `test_size_distribution` and
`check_data_root` pass. Astrodust, NT = 200, `include_euv = .false.`, Mathis
field at U = 1.585 and 100, `heuristic` and `draine`: emitted over absorbed
power (both from the same wavelength grid, absorption with the size-integrated
C_abs of `dust_extinction`) was 0.99643-1.00031 before the transfer and is
1.00048-1.00050 after. The total power moves by +0.019% to +0.409%; the
polarized emission moves differently, e.g. heuristic at U = 100, 1-3 um:
total +0.58%, polarized -1.34%; the emission polarization fraction at 153 um
and 851 um changes by at most 3.5e-4 (absolute). With
`emission_bin_power_cutoff = 1e-6` (heuristic, U = 100) the omitted P*kappB
fraction is 1.4e-7 and the largest change relative to the peak is 3.8e-8 in
the total and 1.1e-9 in the polarized emission; a nonzero cutoff with
`draine` returns status 5 and a negative one status 4.

The absorption solver uses the finite-bin overlap construction of Draine
& Li (2001), equations (15), (16), and (28), on the represented radiation
support. Shared support integration, prepared Planck coefficients,
carbonaceous-optics reuse, and charge-specific PAH cooling are transferred
with the new focused tests in `sed/tests/`. The optional
`emission_bin_power_cutoff` defaults to zero and applies to the heuristic
emission sum, not to heating rates or the probability distribution.
Applications must serialize calls that use the shared mutable model state.

The GD enthalpy grid and cooling remain in use. Geometric bin edges do not
generally make `H(T)` the arithmetic bin centers of DL01 Sec. 4.1; this is
not the complete DL01 statistical treatment. Energy-closure tests do not
establish spectral convergence or bound relative polarized-emission errors.

`SRC_DUST_MODEL` and `SED_DUST_MODEL_SRC` replace the old build-list names.
The MC source list omits the redundant graphite-sphere declaration and
keeps the v1.20 polarized modules in dependency order. The existing
`SRC_POLH5` dependency correction is retained. Historical v1.20 entries
below are preserved. Documents use the updated solver terminology.
The manual and astrodust report also correct an existing claim that the
relative radiation intensity bounds the temperature probability. The support
floor is a numerical convention; excluded absorbed power must be assessed
with the actual cross sections and wavelength intervals.

No v1.00 timings, energy budgets, or spectral changes are adopted as v1.20
measurements. The 71 `sed/output/sed_*.dat` products were regenerated with the
merged code (gfortran; 66 change, the five `equil` products are byte-identical)
and the astrodust report re-measured from them: S1_C1 against HD23, 30-300 um
median -7.57% -> -7.37%, peak -10.06% -> -9.80%; Table mathis_c, the energy
budget (PAH ours/HD23 0.963 -> 0.962; astrodust 0.914, which the committed
products had not matched at 0.911), the Ad + PAH three-way and feature tables,
the conclusions, and the timings (Xeon Gold 6254, one thread: one
`dust_emission` 0.98 s, `calc_sed.x astrodust` 2.56 s). Two claims that were
already wrong before the transfer were corrected from the same products: the
60-300 um residual after the 1.093 amplitude is flat to +-3.3%, not +-2%, and
in the MC report the MC-to-matrix band ratios are NIR 0.959, MIR 0.998, FIR
1.001, sub-mm 1.002 (band-integrated), not "0.97/1.01/1.00/0.99, 1-3% in every
band". The tracked MC products (`mc/output`: three SEDs, `pT_compare_ad`,
the T_rad/T_eq diagnostic) and the four MC figures were regenerated. The
untracked text `kext_*_euv.dat` of astrodust, dl07, mrn and zubko were
regenerated as well: they predated the shipped HDF5 `/kext`, which `calc_kext.x`
reproduces bit for bit, and `pyutil/tests/test_sedust_h5.py` now passes.
Optical tables and the HDF5 products are unchanged.

## 2026-09-08

### G18 Model D gets a scattering asymmetry, computed from first principles

`tmatrix/driver/spheroid_asymmetry_table.f90` (new),
`tmatrix/driver/effective_medium.f90` (new),
`tmatrix/driver/dustem_optical_tables.f90` (new),
`tmatrix/driver/compare_mie_spheres.f90` (new),
`tmatrix/driver/read_index.f90`, `tmatrix/Makefile`,
`data/g18d/oprop/G_amCBE_0.3333x.DAT` (new),
`data/g18d/oprop/G_aSil2001BE6pctG_0.4x.DAT` (new),
`data/dielectric/index_amcBE_ZMCB96` (new),
`data/dielectric/eps_suvSil` (new),
`data/g18d/where.txt`, `data/g18d/sedust_g18d.h5`,
`data/g18d/kext_g18d.dat`, `data/g18d/kext_g18d_euv.dat`,
`sed/src/calc_qtable.f90`, `sed/src/calc_kext.f90`,
`sed/src/sed_astrodust.f90`, `sed/src/sedust_product.f90`,
`sed/src/dust_model_mod.f90`, `sed/src/dust_lib.f90`,
`sed/rt_example/use_dustlib.f90`, `docs/make_g18d_asymmetry_figs.py` (new),
`docs/SEDust_user_manual.tex`, `README_HOWTO.md`.

`g18d` used to reach a radiative-transfer host with `<cos theta> = 0` at all
800 wavelengths while its albedo reached 0.38 — the host would scatter
isotropically. The zero was honest (`m%gsca_complete` was `.false.` and every
driver said so) but it was still a zero, and for a model whose largest grains
reach `x = 2 pi a / lambda ~ 100` in the ultraviolet it is wrong. The cause:
the DustEM distribution ships no `G_amCBE_0.3333x.DAT` and no
`G_aSil2001BE6pctG_0.4x.DAT` — Guillet et al. (2018) computed no asymmetry
parameter, and DustEM reads a `G_` file only under its `pdr` run keyword,
which this model does not set — and those two populations hold 90% of the dust
mass.

**SEDust now computes those two tables**, in DustEM's own `G_` layout, on each
population's own `Q`-table radius grid and on `oprop/LAMBDA.DAT`, from the
same physics the distributed `Q_` tables are built on: prolate spheroids of
axis ratio 1/3 (BE amorphous carbon of Zubko et al. 1996) and 0.4 (the WD01
"smoothed UV" astrosilicate carrying 6% by volume of the same a-C as Maxwell
Garnett inclusions, their Eq. 28), in random orientation, with `Q` and `g`
defined against the volume-equivalent radius. Nothing else about the model
changed: `Q_abs` and `Q_sca` are still the distributed tables, read unchanged.
All three populations now carry a `g`, so `m%gsca_complete` is `.true.` and the
`<cos>` column of `kext_g18d[_euv].dat` is a measurement.

Regime policy, and the fraction of each table it produces, over the radii the
size distribution reaches:

| | Rayleigh, `g = 0` | T-matrix | sphere stand-in |
|---|---|---|---|
| `amCBE_0.3333x` (a <= 1.540 um) | 65.4% | 30.0% | 4.6% |
| `aSil2001BE6pctG_0.4x` (a <= 0.292 um) | 65.1% | 33.6% | 1.3% |

* `x < 0.1`: `g = 0` exactly. Rayleigh dipole scattering is forward-backward
  symmetric; the correction is O(x^2). The refractive index is not read.
* `0.1 <= x <= 50` and the solver converges: the random-orientation T-matrix,
  which is the analytic form of the 768-direction HEALPix orientation average
  of Guillet et al. Eq. (19). Tolerance 1e-3, retried at 1e-2 and 3e-2; a
  retried solve is kept only where it also reproduces the published `Q_abs`
  and `Q_sca` of that cell to 1%, because a converged solution of an
  ill-conditioned system that misses the extended-precision reference by more
  is a false convergence. 24130 / 369 / 145 cells (a-C) and 24381 / 150 / 60
  (silicate) were accepted at the three tolerances.
* otherwise: the volume-equivalent sphere's Mie `g` at the same `m`. This is
  the one approximation, and the same class of stand-in Guillet et al. used in
  the ultraviolet (their Eqs. 11-14, after Min et al. 2003).

The double-precision convergence boundary, measured: `a_V/lambda = 2.65` (a-C)
and `2.74` (silicate), `x_max = 16.6` and `16.8`. Guillet et al. ran the
extended-precision `amplq.lp.f` and quote `a_V/lambda ~ 3` for `b/a = 1/3`
(their footnote 4), so the stand-in covers only what lies past a boundary
close to theirs. `IERR = 3` is the multipole order running to `NPN1 = 100`
without the convergence test being met — the ill-conditioning of the
extended-boundary-condition method for elongated particles, not a storage
limit; a failing solve costs 0.85 s whatever `x` is.

Measured, on this tree:

* **The stand-in bound.** It is entered only at `x >= 5.35` (a-C) and
  `x >= 10.2` (silicate). Over that band, and over the model's radii, the
  sphere and the spheroid differ in `g` by at most **0.0146** (a-C) and
  **0.0094** (silicate). Anywhere in the T-matrix region the gap reaches 0.170
  and 0.183, but that is at `x ~ 1` and `x ~ 4`, where the T-matrix answers.
* **What the stand-in actually decides.** Cell counts overstate it: the size
  distributions are steep power laws (alpha = -4.67 and -3.08) and the cells
  regime 3 supplies are the largest radii. Weighted by each cell's
  contribution to the population's `C_sca` (dn/dln a times a^2 Qsca over the
  model's own radius range) it carries 4.8% of the a-C scattering at the Lyman
  limit, 2.1% at 0.55 um and 0% beyond 2 um, and 16% of the silicate's at the
  Lyman limit and 0% longward of 0.2 um. With the bounds above it can move
  each population's `<cos theta>` by at most 7e-4 (a-C) and 1.5e-3 (silicate,
  at the Lyman limit only).
* **Check 1**, our Mie with `index_amcBE_ZMCB96` against DustEM's own BE a-C
  sphere tables `Q_amCBEx` / `G_amCBEx`, 33500 cells over 0.0912-800 um:
  median 0.17% in `Q_abs`, 0.051% in `Q_sca`, 3.5e-6 in `g`; max 3.2%, 4.2%,
  6.2e-3. 541 of those cells carry an exact 0 in `G_amCBEx` where Mie gives
  `g < 0`, which a strongly absorbing sphere does around `x ~ 0.5`: that
  reference is floored, so those cells are counted and reported apart.
* **Check 2**, our Mie with `eps_suvSil` against Draine's `suvSil_81`, 14060
  cells: median 0.016%, 0.032%, 1.7e-6. The outliers are a defect of the
  reference, shown without a second code: at `x < 0.01`, `Q_abs/a` cannot
  depend on the radius, and `suvSil_81` puts 83 of 4205 such cells more than
  1% off the median over radii, worst 23.6%, in the 6-10 um silicate feature;
  ours obeys the scaling to 0.1%.
* **Check 3**, our T-matrix `Q` against the published spheroid tables where the
  solver converged: a-C median 0.15% (`Q_abs`) and 0.12% (`Q_sca`), max 3.0%
  and 4.1% over 20651 in-model cells. Silicate median 0.39% and 1.6%, but with
  a max of 38% confined to 0.09 < lambda < 0.2 um and largest at the smallest
  radii, where the calculation is pure Rayleigh — so it is `m(lambda)`, not the
  shape, the mixing rule or the orientation average. Outside that band the
  median is 0.35%. Inverting the published `Q_abs` for the implied
  `Im[(m^2-1)/(m^2+2)]` gives a smooth curve that none of the three Draine
  silicate tables in this tree reproduces, so **the astrosilicate index Guillet
  et al. used in the far ultraviolet is not one this tree holds**. Its effect on
  `g` is bounded directly: swapping the matrix for Draine's *unsmoothed*
  astrosilicate, which changes `Q_abs` by up to 36% at 0.145 um — the same size
  as the disagreement — moves `g` by at most 0.016.
* **Maxwell Garnett**, checked against its own derivation: `f = 0` gives
  `eps_m`, `f = 1` gives `eps_i`, identical components give `eps_m` for every
  `f`, all exactly, and the Clausius-Mossotti closure
  `(eps-eps_m)/(eps+2eps_m) = f (eps_i-eps_m)/(eps_i+2eps_m)` has a residual of
  1.7e-17.
* **Emission is unchanged.** `sed/output/sed_g18d.dat` and `sed_g18d_morig.dat`
  are byte-identical to their previous copies: the emission side never reads
  `g`.
* **Wall clock** 636 s for both populations on 36 threads (401 s + 235 s),
  ~37 MB of T-matrix workspace per thread.

`read_index.f90` was rewritten around a `refractive_index_t` type so that two
materials can be held at once, which the composite silicate needs, and it now
reads Draine's wavelength-ordered `eps_*` layout as well as his energy-ordered
`index_*` one. Header lines are recognized rather than counted. `load_index`
and `interp_m`, which every other driver in that directory uses, act on a
single table owned by the module and are **bit-identical** to the old reader:
compared over 200001 wavelengths on `index_DH21Ad_P0.20_0.00_1.400`,
`max|dn| = max|dk| = 0`. Its misleading `expected >= 6000` diagnostic, which
never matched its own `ndata < 2` test, is gone.

`data/dielectric/` gained copies of `index_amcBE_ZMCB96` and `eps_suvSil`
(130 KB together, byte-identical to `Grain/opt/`) so that the two `G_` files
can be regenerated from the release alone.

Figures: `docs/figs/g18d_mie_amCBEx.pdf`, `g18d_mie_suvSil.pdf`,
`g18d_spheroid_vs_dustem.pdf`, `g18d_sphere_vs_spheroid.pdf` and
`g18d_gbar.pdf`, from `docs/make_g18d_asymmetry_figs.py`.

### One `check_build_dust.x` in both trees; `dust_extinction` status 4 checked

`sed/src/check_build_dust.f90` (new -- the same file as v1.00's),
the narrower DustEM-only duplicate it replaces (deleted),
`sed/src/sed_astrodust.f90`, `sed/src/calc_kext.f90`, `sed/Makefile`,
`README_HOWTO.md`, `docs/SEDust_user_manual.tex`, `docs/CHANGES.md` (the
2026-09-05 and 2026-09-06 entries name the surviving program).

v1.20 carried a DustEM-only duplicate of this check, written when THEMIS and
G18D were added on 2026-09-05. It compared those two models' HDF5 products
against their DustEM text tables and nothing else, so v1.20 had no
product-vs-text check at all on astrodust, DL07, MRN or Zubko. v1.00 carried
`check_build_dust.x`, which makes that comparison for all six models on both
wavelength grids and adds the exact ZDA round trip. The two are now **one
program**: `sed/src/check_build_dust.f90`, byte-identical in the two trees and
built as `check_build_dust.x` in both. The narrower duplicate is deleted, and
v1.20 gains the four models it was not checking.

No library interface was changed to make the port fit. Every routine the
program calls -- `build_dust`, `build_astrodust`, `build_dl07`, `build_mrn`,
`build_zubko`, `build_dustem`, `size_integrated_extinction`,
`dust_mass_per_H`, `dust_extinction`, `read_zubko_optics`,
`read_sedust_qtable` -- takes the same arguments in the same order in both
trees for the arguments this program passes, so the one file compiles in both
against each tree's own `sed_astrodust.f90`.

Folded in from the deleted program: the two routes must name the **same
populations** in `gsca_missing`, not merely reach the same `gsca_complete`.
`build_dust` reaches the product through `build_dustem` with `qtable_path` but
does not pass that string back, so `check_dustem` asks the product route for it
directly, by the same call `build_dust` makes.

New check, `check_undefined_asymmetry`: `dust_extinction` status 4 is a
**warning**, not an error -- `gbar` was asked for from a model that carries no
scattering asymmetry, so the `gbar` returned is 0 because the asymmetry is
undefined and not because the scattering is isotropic. The check takes the
astrodust model, which has an asymmetry, calls `dust_extinction` as built and
again on a copy with `gsca_complete` forced `.false.`, and requires status 0
then 4, `gbar` identically 0 on the second call, and `C_ext`, `C_abs`, `C_sca`
bitwise identical between the two. `compare` no longer counts status 4 as a
failure either; it requires both routes to reach the same status, whichever it
is, which is what a model with no asymmetry needs.

`dust_extinction` now **writes that zero itself** (`sed/src/sed_astrodust.f90`)
instead of passing on whatever the extinction table's `g` column holds. The
routine's own header already promised the zero; it was true only because
`calc_kext.x` writes a zero `g` column for such a model, which is a property of
the file and not of the routine. A curve from another vintage or another tool
could hold the asymmetry of the populations that do carry one over the
scattering of all of them -- the asymmetry of nothing, which reads as a small
`g` rather than as a missing one. `size_integrated_extinction` refuses that
ratio for the same reason. **No shipped number moves**: the `g` column of
`kext_g18d.dat` and `kext_g18d_euv.dat` is 0 at every wavelength already, and
those are the only shipped curves whose model has `gsca_complete = .false.`

Measured, `./check_build_dust.x all`:

* astrodust, DL07, MRN and Zubko, both grids, against `TOL = 1e-5`: the size
  integrals agree to 1.2e-7 at worst (DL07 `C_sca`), $\langle\cos\theta\rangle$
  to 5.1e-6 (DL07), the served curves to 1.8e-8 in `C_ext`/`C_abs` and 6.4e-7
  in $\langle\cos\theta\rangle$, and `M_dust/N_H` exactly, at 0. The wavelength
  grids agree exactly except for the DL07 and MRN EUV grids, at 2.2e-16.
* THEMIS and G18D, both grids, against `TOL_STORED = 1e-10`: grid, size
  integrals, $\langle\cos\theta\rangle$ and `M_dust/N_H` all exactly 0; the
  served curve 4.1e-13 to 4.8e-13, the printed precision of the text table.
  The G18D half of that line was measured at 14:19 today and no longer
  reproduces; see the note below.
* the ZDA round trip: `/qtable/{sil,gra,pah}` equals the distributed
  `suvSil_121_1201.dat`, `Gra_121_1201.dat` and `PAH_28_1201_neu.dat` outright,
  max$|{\rm d}Q_{\rm abs}|$ = max$|{\rm d}Q_{\rm sca}|$ = max$|{\rm d}g|$ = 0.
* the status-4 check: astrodust as built gives status 0 with
  max$|\bar g|$ = 0.782; with `gsca_complete` forced `.false.` it gives status
  4, $\bar g$ identically 0, and max$|{\rm d}C|$ = 0 across `C_ext`, `C_abs`
  and `C_sca`.

The stale claim in `docs/SEDust_user_manual.tex` that four of the THEMIS/G18D
checks stand over tolerance at 1.9e-8 and 2.6e-8 is corrected: both sides of
that comparison were recomputed on 2026-09-06 and now agree to 4.5e-13. The
figure is kept in the text as history, stated as history.

**G18D is not verifiable in either tree at the moment, and the reason is not
this program.** `data/g18d/oprop/G_amCBE_0.3333x.DAT` and
`G_aSil2001BE6pctG_0.4x.DAT` -- the two asymmetry tables the DustEM
distribution does not publish and SEDust computes -- were written on
2026-09-08 (v1.20 at 14:23, v1.00 at 14:27), and neither tree's
`data/g18d/sedust_g18d.h5` has been rewritten since 2026-09-06. The text route
therefore finds a `g` for every scattering population
(`gsca_complete = .true.`) while the product route does not
(`gsca_complete = .false.`), and ten G18D assertions report that disagreement
in each tree: $\langle\cos\theta\rangle$, `gsca_complete`, both
`g`-population comparisons and the differing `dust_extinction` status, once per
grid. `C_ext`, `C_abs`, `C_sca`, the wavelength grid and `M_dust/N_H` still
agree exactly. Re-run `./check_build_dust.x g18d` once the product has been
rewritten from the new tables. The five other models, the ZDA round trip and
the status-4 check pass in both trees, with the numbers listed above.

## 2026-09-06

### HDF5 products back to 64-bit; `regime` alone stays an 8-bit code

`sed/src/sedust_h5.f90`, `sed/src/calc_qtable.f90`, `sed/src/calc_kext.f90`,
`sed/src/calc_polarized_optics.f90`, `sed/src/check_build_dust.f90`,
`sed/src/test_h5_storage_policy.f90`, `pyutil/migrate_h5_regime_int8.py` (new,
replacing `pyutil/migrate_h5_float32.py`), `sed/Makefile`, `README_HOWTO.md`,
`docs/SEDust_user_manual.tex`, every `data/<model>/sedust_<model>.h5`, and every
file under `sed/output/`.

The 32-bit conversion of the previous entry is undone. Two storage classes now,
not three:

| class | file type | datasets |
|---|---|---|
| values | 64-bit float | every computed quantity -- `Q_ext`, `Q_abs`, `Q_sca`, `g`, `albedo`, the `C_*` cross sections, `K_abs`, `Z`, `F_tot`, `F_ref`, `C_polext`, the aligned `C_*_al` -- and the coordinate axes `lambda`, `a_eff`, `theta`, `theta_i`, `theta_s`, `phi`, `band_lambda` they are tabulated on |
| codes | 8-bit integer | `regime` |

A number read back out of a product is again the number that was computed, so
nothing lands on disk that a later comparison has to be given slack for, and an
axis stays strictly ascending however finely it is spaced. `regime` keeps its
byte: three distinct values, which eight bytes would carry while saying they are
continuous.

Mechanism: the writers `h5_write_1d/2d/3d/6d` and the `write_real` they share no
longer take a `single` argument, and `h5dcreate_f` is back to
`H5T_NATIVE_DOUBLE`. `h5_write_2d_i8` and `h5_read_2d_int` stay -- those are the
`regime` path. The memory type on write and read was double throughout and still
is, so no reader changed here either and a caller passes and receives
`real(real64)`.

The products were **restored from the pre-conversion commit**, not widened back
from 32-bit. Widening would have returned 64-bit numbers that merely equal the
32-bit values, and the digits the conversion destroyed do not come back that
way. Every value dataset is therefore bit-identical to the pre-conversion
product.

`sed/output/` was **regenerated**, not restored: the header lines of several
files were corrected in the conversion commit and had to survive. Every
tracked file was rerun from the restored products, and the result compared
against the pre-conversion commit file by file.

Measured on this tree:

* every dataset of the six products: 174 of 175 are 64-bit float, the one
  `regime` is 8-bit integer; no product lost a dataset,
  `/polarized/scatmat_aligned/Z` included.
* the six products come to 122.4 MB on disk (astrodust 86.76, dl07 9.77,
  mrn 6.47, zubko 16.31, themis 2.16, g18d 0.94 MB).
* `./check_build_dust.x` passes at 1e-10 again, every reported deviation
  exactly 0.00E+00.
* `sed/output/`, all 74 tracked files rerun and compared against the
  pre-conversion commit: 68 byte-identical, 3 differing only in the header
  lines that commit corrected (the two THEMIS `CM20` channels now named
  `CM20_plaw-ed` and `CM20_logn`, and the G18D `aSil2001BE6pctG_0.4x` no
  longer truncated). Of the remaining three, `sed_g18d_morig.dat` has no
  pre-conversion counterpart -- the conversion commit created it -- and
  `sed_dl07_mw31_60_qm.dat` and `sed_dl07_mw31_60_qm_stati_morig_ns150.dat`
  differ in their numbers by 1.7e-4 and 2.3e-4 relative, because they were
  last written at commit f213897 and never rerun afterwards.
  Built from the pre-conversion sources against the pre-conversion products,
  those two runs give exactly what this tree now gives, so the regenerated
  files are the first that match the code and data beside them. The runs are
  deterministic: three consecutive runs give one md5.

Tolerance back at the value it had before the conversion, the reason for
loosening it having gone:

* `check_build_dust.f90` `TOL_STORED`, the DustEM product-vs-text gate:
  1e-6 -> **1e-10**.

`./test_h5_storage_policy.x` still guards the policy, now asking of every
dataset that a value be 64-bit float and `regime` 8-bit integer.
`pyutil/migrate_h5_regime_int8.py` replaces `pyutil/migrate_h5_float32.py`: it
stores the `regime` codes of an existing product as one byte each and requires
every other dataset to come out bit-for-bit identical.

## 2026-09-05

### HDF5 products store computed quantities as 32-bit; axes stay 64-bit

`sed/src/sedust_h5.f90`, `sed/src/calc_qtable.f90`, `sed/src/calc_kext.f90`,
`sed/src/calc_polarized_optics.f90`, `sed/src/check_build_dust.f90`,
`sed/src/test_h5_storage_policy.f90` (new),
`pyutil/migrate_h5_float32.py` (new), `sed/Makefile`, `README_HOWTO.md`,
`docs/SEDust_user_manual.tex`, and every `data/<model>/sedust_<model>.h5`.

Three storage classes, decided by what a dataset is:

| class | file type | datasets |
|---|---|---|
| computed quantities | 32-bit float | `Q_ext`, `Q_abs`, `Q_sca`, `g`, `albedo`, the `C_*` cross sections, `K_abs`, `Q_re`, `Z`, `F_tot`, `F_ref`, `C_polext`, the aligned `C_*_al` |
| coordinate axes | 64-bit float | `lambda`, `a_eff`, `theta`, `theta_i`, `theta_s`, `phi`, `band_lambda` |
| codes | 8-bit integer | `regime` |

The quantities carry their own uncertainty far above float32's 1.2e-7 relative
resolution, so the other four bytes buy nothing. The axes are the exception
because their problem is distinctness, not accuracy: the astrodust wavelength
grid resolves each X-ray absorption edge with a pair of points 6.7e-7 apart in
relative wavelength, six to eleven representable float32 values, and they are
0.06% of the payload.

Mechanism: the writers `h5_write_1d/2d/3d/6d` and the `write_real` they share
take an optional `single` argument, `.false.` by default so that an untagged
call site costs disk space and never accuracy. Only the FILE type moves; the
memory type on write and read stays `H5T_NATIVE_DOUBLE`, so no reader changed
and a caller still passes and receives `real(real64)`.

Products already in the tree were converted in place by
`pyutil/migrate_h5_float32.py` rather than regenerated: nothing about the
storage type of a number requires recomputing the number, and rerunning
`calc_qtable.x` would replace the file and take `/polarized` with it. The script
rewrites dataset by dataset (HDF5 leaves a hole where a dataset is replaced, so
only a rewrite reclaims the space), preserving values, shape, chunking,
compression, shuffle and every attribute.

Measured on this tree:

* the six products fall from 122.4 MB to 54.1 MB on disk
  (astrodust 86.76 -> 39.16, dl07 9.77 -> 4.06, mrn 6.47 -> 2.68,
  zubko 16.31 -> 6.77, themis 2.16 -> 0.99, g18d 0.94 -> 0.45 MB). The astrodust
  product keeps all 41 datasets, `/polarized/scatmat_aligned/Z` included, at its
  own chunking; the drop also clears the 50 MB warning GitHub prints on push.
* migrated values against the originals: worst relative change **5.96e-8**,
  which is one unit in the last place of a float32 significand. Two datasets
  reach below float32's normal range and lose their smallest entries to zero:
  the long-wavelength end of `kext/C_sca` (down to 6.2e-40 cm^2/H against a peak
  of 1.3e-21), and `/polarized/scatmat_aligned/Z`, 82064
  of 10179440 entries below 1.2e-38 against a peak of 3.9e-14 -- thirty orders
  of magnitude down, in a table only ever summed or interpolated linearly. Both are checked against
  the peak of their own array as well as elementwise.
* `kext_*.dat` text products regenerated from the migrated files: the E-format
  columns (`C_ext`, `C_abs`, `C_sca`, `K_abs`, `C_polext`) move by at most
  **9.97e-8** relative; the fixed-point `albedo` and `<cos>` columns of
  `kext_astrodust_MW.dat` move by one unit in the last printed digit (1.0e-6
  absolute), the underlying values by 1.1e-8 as the 13-digit `_euv` companion
  shows. `kext_ld01_MW*.dat` are unchanged, that route reading no stored table.
* `test_dustem_extinction.x` still agrees with DustEM's own `EXT_*.RES` to
  5.0e-7, unchanged: that comparison reads the DustEM text tables on both sides.

Tolerances re-tuned, because a 32-bit product and a text table of 7 to 13
significant digits cannot agree better than float32 resolution:

* `check_build_dust.f90` `TOL_STORED`, the DustEM product-vs-text gate:
  1e-10 -> **1e-6**. Measured after the change: 3.4e-8 (THEMIS and G18D, both
  views); `M_dust/N_H` still agrees exactly.

New: `./test_h5_storage_policy.x` reads the file types of every product in the
tree and fails unless each is 64-bit float (axes), 8-bit integer (`regime`) or
32-bit float (everything else). That is what keeps an axis from silently
becoming single later.

## 2026-08-20

### Programs renamed `calc_*`, the three SED drivers merged into `calc_sed.x`

`sed/src/calc_sed.f90` (new), `sed/src/calc_enthalpy.f90`,
`sed/src/calc_qtable.f90`, `sed/Makefile`, and in v1.20
`sed/src/calc_polarized_optics.f90`.  Every executable now carries the quantity
it computes in its name, `calc_<quantity>`, so that none of them reads as the
build tool:

| was | is |
|---|---|
| `main_astrodust.x` | `calc_sed.x astrodust` |
| `main_dl07.x` | `calc_sed.x dl07 [submodel]` |
| `main_zubko.x` | `calc_sed.x zubko` |
| `make_enthalpy.x` | `calc_enthalpy.x` |
| `make_qtable.x` | `calc_qtable.x` |
| `make_polarized.x` (v1.20) | `calc_polarized_optics.x` |

`calc_polarized_optics.x` is named for both of its products -- the
orientation-resolved cross sections AND the random and aligned scattering
matrices; `calc_qpol.x` would have named only the first.

The three model drivers were one program in three copies: same solver, same
radiation field, same emission term, same writer, differing in which model is
built and which channels are written.  `calc_sed.x <model> [settings ...]` is
that program, with the model as a required first argument.  The outputs are
named the same way for all three, `output/sed_<model>[_<submodel>][_<tag>]
[_<stage>].dat`:

| was | is |
|---|---|
| `astrodust_irem_ours_<tag><stage>.dat` | `sed_astrodust[_<tag>]_<stage>.dat` |
| `dl07_sed_ours_<submodel>[_<tag>].dat` | `sed_dl07_<submodel>[_<tag>].dat` |
| `zubko_sed_ours[_<tag>].dat` | `sed_zubko[_<tag>].dat` |

The existing outputs were RENAMED, not recomputed: a change of name changes no
number.  `calc_sed.x astrodust`, `calc_sed.x dl07` and `calc_sed.x zubko euv`
reproduce the renamed files byte for byte in both trees.

### One command-line system for every program in `sed/`

`sed/src/sed_run_options.f90` (rewritten), `sed/src/sed_apply_options.f90`
(new).  The settings are grouped into independent axes, and each program
DECLARES which axes it has a referent for before reading a single argument.  A
word of an undeclared axis is refused, naming the axis, instead of being
ignored:

    $ ./calc_kext.x dl07 qm
     calc_kext: 'qm' selects the stochastic-solver axis, which calc_kext has no
     referent for.

and when the axis exists in the program but not in the model it was asked
about, the model is named: `calc_sed: 'c2' selects the Stage-1-enthalpy axis,
which calc_sed's zubko has no referent for.`

| program | axes |
|---|---|
| `calc_sed.x` | subject, solver, grid, field, emission, matrix size, graphite (astrodust, dl07), PAH vintage (dl07), Stage-1 enthalpy (astrodust) |
| `calc_kext.x` | subject, grid, PAH vintage (dl07), Zubko size distribution and optics set (zubko) |
| `calc_qtable.x` | subject |
| `calc_enthalpy.x` | Stage-1 enthalpy (`c2`) |
| `calc_polarized_optics.x` (v1.20) | subject |
| `calc_polext.x` (v1.20) | none |

Where a setting LANDS is a second module, `sed_apply_options.f90`: a program
that computes an enthalpy table or a cross-section table links no solver and no
radiation field, and would otherwise have to drag both in to get the shared
command line.

`calc_enthalpy.x c2` writes `output/enthalpy_S1_c2.dat` and no Stage-2 table --
`enthalpy_S2` does not read that setting, so a tagged copy would be a second
name for the same numbers.

### `ld01`: the LD01 carbonaceous vintage as an axis of the DL07 model

`sed/src/qpah.f90`, `sed/src/sed_astrodust.f90`, `sed/src/calc_kext.f90`,
`sed/src/calc_sed.f90`.  `qpah.f90` gains `qpah_xsec_vintage` (`'dl07'` |
`'ld01'`) and the dispatcher `qpah_abs`, in the same shape as
`qpah_graphite_source`; `sed_init_dl07` calls the dispatcher, and `build_dl07`
takes a `pah_xsec` argument that also selects the vintage's own extinction
curve (`kext_ld01_MW*.dat`, `/kext_ld01`) so that what `dust_extinction` serves
is the size integral of the very cross sections the model was built on.  The
stored DL07 cross-section tables are the DL07 vintage, so an `ld01` build
solves every optic from the dielectric functions instead of reading a table
that cannot answer the request.

    ./calc_kext.x dl07 ld01 [euv]     -> ../data/dl07/kext_ld01_MW[_euv].dat
    ./calc_sed.x  dl07 ld01           -> output/sed_dl07_mw31_60_ld01.dat

Measured on `mw31_60`: the two curves share their `C_sca` and `<cos>` columns
bit for bit -- the vintages differ in the carbonaceous ABSORPTION only -- and
the carbonaceous `C_abs` of LD01 against DL07 runs from 0.915 at 1.047 um,
where DL07 adds the PAH cation near-infrared resonance after Mattioda et al.
(2005), to 1.109 at the 6.2 um C-C feature, reaching 1.000 by 100 um.

### `calc_kext_dl07.x` retired; the 2003-table comparison moved to the figure

The program computed the WD01/D03 size integral on the wavelengths of Draine's
2003 `kext_albedo` table, in both carbonaceous vintages.  With the `ld01` axis
above, `calc_kext.x dl07` and `calc_kext.x dl07 ld01` compute the same two
curves on the model's own grid, and `calc_kext.x` already prints the
point-by-point comparison against that table.  What is gone is its private
400-point size grid over 3.5e-4 - 3.0 um; the products use the model's own
84-point grid over 3.548e-4 - 5.012 um, which is the model definition.

Measured before removing it, at the seven wavelengths the two grids share (no
interpolation): `C_ext` agrees to 0.06-0.16%, worst at 0.0912 um where the
smallest grains weigh most.  Interpolating the product onto all 553 overlapping
reference wavelengths widens the spread to 0.3% in the continuum and 0.9% at the
9.7 and 11.3 um resonances, which is the cost of resampling a narrow feature and
not a difference between the two calculations.  Figure 1 of the JKAS paper
(`SEDust_JKAS/python/plot_optics_combined.py`) now reads the products; the
rendered panels are the same figure, with the curves reaching further at both
ends because the model grid does.

Removed: `sed/src/calc_kext_dl07.f90`, its make target,
`sed/output/kext_albedo_{dl07,ld01}_ours.dat`.

### `cmp_refs.x` retired -- and it was mixing two CMB conventions

The program solved the DL07 and Zubko SEDs at U = 1 in both Mathis-field
conventions.  With `mathis_orig` and `euv` reachable from the driver, its four
files are two pairs of driver runs, so it is redundant -- but it was also
wrong, and the driver is not.

`cmp_refs.f90` set `use_mathis_corrected` AFTER calling the builder.  The
builders bake the CMB temperature into `kappCMB` through `cmb_temperature()`,
and `calc_P` subtracts that term as the net cooling rate, so each model's
SECOND convention was evaluated against a cooling term built for the first.
Measured: `cmp_dl07_corr.dat` (built and run corrected) reproduces
`calc_sed.x dl07` to 1.000000000 in every point, and `cmp_zubko_orig.dat`
(built and run uncorrected) reproduces `calc_sed.x zubko euv mathis_orig` the
same way; the two mismatched files differ from their driver counterparts by up
to 9.3e-6 (DL07) and 6.4e-5 (Zubko, whose grid reaches 1e4 um where the CMB
term dominates).  `calc_sed.x` applies the convention BEFORE the build, so
field and cooling term always agree.

Removed: `sed/src/cmp_refs.f90`, its make target,
`sed/output/cmp_{dl07,zubko}_{corr,orig}.dat`.  Figure 10 of the JKAS paper now
reads `sed_dl07_mw31_60[_morig].dat` and `sed_zubko_euv[_morig].dat`; the
rendered figure differs from the old one in 2 pixels out of 226300.

Earlier entries in this file name the programs by their current names, so that
what they describe can still be run.

### One vocabulary of run settings for the three model drivers

`sed/src/sed_run_options.f90` (new), `sed/src/calc_sed.f90`,
`sed/src/calc_sed.f90`, `sed/src/calc_sed.f90`, `sed/Makefile`.

`calc_sed.x dl07` read `[model] [euv] [qm|draine]` and `calc_sed.x zubko` read
`[heuristic|draine|qm|equil] [euv] [hardfield]`, while `calc_sed.x astrodust` read
thirteen settings. The solver, the radiation field, the emission term and the
transition-matrix sizes are all shared code, so most of that difference was in
the drivers alone: the same setting was reachable for one model and not for
another that would have read it. The three drivers now take ONE vocabulary,
parsed and tagged by the new module, grouped into independent axes:

    solver    heuristic (default) | draine | equil | qm | qm_dbcon | qm_stati
    grid      euv
    field     mathis_orig, logU=X, hardfield
    emission  induced, photcut
    qm sizes  nstate=N, nisrf=N
    graphite  gra_d03_sphere | gra_d16_sphere | gra_d16_spheroid
    enthalpy  c2

Any combination across axes is a valid run. Two are refused instead of being
resolved silently: two values of one axis, and a setting the chosen solver does
not read (`nstate=`/`nisrf=` without a `qm` solver, `photcut` with any solver
but `heuristic`, which is the only one whose bin sum reads it —
`sed/src/sed_astrodust.f90` line 1463). Newly reachable: `equil` and
`heuristic` by name on `calc_sed.x astrodust`; `qm_dbcon`, `qm_stati`, `equil`,
`mathis_orig`, `logU=`, `induced`, `photcut`, `nstate=`, `nisrf=`, `hardfield`
and the graphite axis on `calc_sed.x dl07`; the same minus the graphite axis on
`calc_sed.x zubko`.

Two settings are NOT offered on every driver, because the model has no referent
for them. `c2` scales the astrodust Stage-1 enthalpy prefactor
(`sed/src/enthalpy_astrodust.f90` line 74, reached only from `sed_init`), and
DL07 and Zubko carry no astrodust component — DL07 uses `enthalpy_DL01`, Zubko
the ZDA calorimetry tables. The `gra_*` axis names the graphite of the
PAH-to-graphite xi blend, which `build_zubko` never computes: that model's PAH
absorption is read from its own optics tables (`sed/src/sed_astrodust.f90`
lines 2882–2900), and `qpah` appears in that builder only in a comment.

On `calc_sed.x dl07` a named graphite also switches the optics off the stored
tables: those were computed once with this model's own `'d03_sphere'` graphite,
so reading them back would leave the request with no effect at all. The driver
then passes `stored_q_dir=''` and every optic is solved from the dielectric
functions, which is also what puts two graphite variants on the same route.
Measured on `mw31_60`: `gra_d03_sphere` reproduces the stored-table default in
every printed digit of all four columns, so the route change costs nothing on
its own; against it, the carbonaceous channel of `gra_d16_spheroid` spans
0.61–16.7 and of `gra_d16_sphere` 0.33–15.4 in ratio over the grid.

Every setting tags the output filename, in one fixed order whatever order the
words are typed in, and a default run carries no tag — so the production
filenames are unchanged (`sed_dl07_<model>.dat`, `sed_zubko.dat`,
`sed_astrodust_<stage>.dat`). Regression, both trees: the default runs of
all three drivers, plus `calc_sed.x dl07` on `lmc2_10`, `smc`, `euv` and `qm`,
`calc_sed.x zubko` on `draine`, `equil`, `euv`, `euv hardfield` and `qm`, and
`calc_sed.x astrodust` on `mathis_orig` and `qm`, are byte-identical to the same
runs of the previous build.

### PAH cation near-infrared continuum: `exp[-(0.1x)^2]`

`sed/src/qpah.f90`, `qpah_dl07`. The near-infrared continuum that DL07 eq. 2
adds to a PAH cation (after Mattioda et al. 2005a) is

    C_abs/N_C = 3.5e-19 * 10^(-1.45/x) * exp[-(0.1x)^2] cm^2,   x = 1/lambda[um],

with `0.1x` inside the square. The routine carried `exp(-0.1 x^2)`, which is a
different, much narrower Gaussian. With the corrected factor the astrodust PAH
SED at log U = 0.20 rises by 0.99% bolometrically (0.1--3000 um integral of
`sed/output/sed_astrodust_PAH.dat`), and the DL07 model at U = 1 rises by
0.41% in its carbonaceous channel and 0.28% in the total
(`sed/output/sed_dl07_mw31_60.dat`); the DL07 silicate channel is
unchanged to the printed digits. The astrodust S1 and S2 outputs do not move at
all — the term lives in the PAH population only.

The correction also shifts the PAH row of the report's Mathis-convention table
(`docs/astrodust_sed_report.tex`, `tab:mathis_c`) by +0.76 pp on its own.

### D16 spheroid graphite table: random-orientation average, not jori = 1

`sed/src/q_graphite_d16.f90`. Draine's `qlib_gra_D16MGemt_1.400` is
orientation-resolved, and on the convention its own header states — the same one
the HD23 astrodust table uses — `jori = 1` is **k** parallel to **a**, a single
orientation, not the random-orientation average. The module read `jori = 1` and
labeled it `Q_rand`. It now forms

    <Q>_random = (1/3) Q(jori=2) + (2/3) Q(jori=3)

on load, which is what an unaligned carbonaceous population needs. Reading
`jori = 1` for it overstates `Q_abs` by ~22% in the Rayleigh limit, because
`jori = 1` carries no **E** parallel to **a** component at all. The
300--3000 um PAH-only band ratio of the D16-spheroid variant against HD23
`PAH_irem.dat` at log U = 0.20 moves from 1.21 to 1.05 (the 1.05 re-measured
here from `sed/output/sed_astrodust_sphdgra_PAH.dat`; the 1.21 is the
value the earlier, `jori = 1`, build gave). The production D16-sphere variant is
unaffected — it is a sphere table with no orientation axis — and remains the
closest to HD23 in that band, at 0.99.

### The graphite of the xi blend is a named source, not a boolean

`sed/src/qpah.f90`. The logical `qpah_use_d03_graphite` is replaced by the
character setting

    qpah_graphite_source = 'd16_sphere'    ! default; D16 turbostratic, MG EMT, sphere
                         | 'd03_sphere'    ! Draine 2003 graphite, Mie, (1/3)|| + (2/3)perp
                         | 'd16_spheroid'  ! the D16 material on the b/a = 1.4 oblate spheroid

with an unrecognized value stopping the run rather than being taken as one of
the two the boolean could express. A new routine `q_graphite_xi_blend_abs`
dispatches on it. Above the DL07 PAH cutoff (x >= 17.25) all three settings
still take graphite from the D03 dielectric functions, because there the blend
needs radii below the smallest the D16 tables carry; that branch is a property
of the tables, not a preference, and it was previously entangled with the
boolean.

Callers updated: `sed/src/calc_sed.f90`, `sed/src/calc_qtable.f90`,
`sed/src/sed_astrodust.f90` (`build_astrodust` and `build_dl07`).

`sed/src/calc_sed.f90` exposes the two non-default settings as the
arguments `gra_d03_sphere` and `gra_d16_spheroid`, which tag their output files
`sed_astrodust_d03gra_*.dat` and `sed_astrodust_sphdgra_*.dat`.
Naming both is refused. Both touch the PAH population only, so the S1 and S2
files they write are byte-identical to the production ones.

### Exact-statistical cooling kernel in the energy-space solver

`sed/src/stoch_qm.f90`. A third cooling treatment joins `dbdis`
(thermal-discrete, the reference) and `dbcon` (thermal-continuous):
`stati`, the microcanonical kernel, whose downward rate carries the degeneracy
ratio `g_f/g_i` instead of a Planck factor, with `g(U)` counted exactly by the
Beyer--Swinehart recursion over the vibrational mode spectrum. Emission is the
matching degeneracy-weighted kernel rather than a bin blackbody. New routines:
`calc_glu`, `stat_cooling_funct`, `intrabin_stat_cooling_funct`,
`simpson_stat_cooling`, `simpson_intrabin_stat_cooling`,
`calc_stat_cooling_afi`, `beyer_swinehart`, `T_from_U_modes`,
`thermal_energy_modes`, `compute_debye_lngu`. `build_enthalpy_bins_qm` now also
returns `lngu` and a flag saying whether it came from the exact count or from
the Debye limit; `build_transition_matrix` takes `lngu` and selects the downward
kernel on `method`. The 20 kT cutoff of the thermal kernel is not applied to
`stati`, which has no Planck suppression to justify it.

The count is feasible only for small carbonaceous grains: the solver caps the
`stati` path at a <= 25 A (`A_STATI_MAX_CM`), and within that cap the
feasibility guards hand grains with N_C >= 137 (a >= 6.7 A) back to `dbdis`,
so in the production configuration `stati` actually runs on PAH grains with
N_C <= 116 (a <= 6.3 A) and every other grain reverts to `dbdis` within the
same run. Silicate uses the thermal kernel throughout: the count that survives
the feasibility guards still leaves a residual non-monotonicity in ln g that the
exponential of the difference amplifies into a spurious mid-infrared spike.

Selected by `./calc_sed.x astrodust qm_stati [nstate=N]`, which tags its output
`sed_astrodust_qm_stati_*`. Measured against the thermal-discrete run at
`nstate=500`: the PAH SED integrates to 1.00033 times the `qm` one over
0.1--3000 um, and the band ratios stay within 0.25% over 1--3000 um (5--15 um
0.9990, 15--60 um 1.0025, 60--300 um 1.0003, 300--3000 um 0.9989); the
sub-1 um band, whose absolute level is negligible, differs by 11%. The S1 and S2
files of the two runs are byte-identical, since only PAH grains take the new
kernel. The `dbdis` and `dbcon` paths are unchanged apart from the added `lngu`
argument.

`beyer_swinehart` continues the counted density of states above the overflow
point with the thermodynamic asymptotic `d ln g = dE/kT + d ln(dE)`. The
coefficient there is hc/k = 1.43877 cm K. Draine's `dens_states.f` carries
1.48377 at this one place, a digit transposition — every other use of hc/k in
that source is 1.43877. The branch is not reached for the grains the
exact-statistical path runs on, so the correction changes no output of the
released solver.

### Two reference-comparison programs

- `sed/src/calc_kext_dl07.f90` (`make calc_kext_dl07.x`). Size-integrated
  C_ext/H, albedo, `<cos>` and C_abs/H for the WD01/D03 Milky Way model on the
  wavelength grid of Draine's `kext_albedo_WD_MW_3.1_60_D03.all_2003`, computed
  twice on one and the same size grid, size distribution, silicate optics,
  graphite scattering and charge mixing — once with the LD01 and once with the
  DL07 carbonaceous absorption cross sections, so the only difference between
  the two outputs is that cross section. Writes
  `sed/output/kext_albedo_{ld01,dl07}_ours.dat`. The 2003 table was built with
  the LD01 cross sections, so the LD01 run is the one that reproduces it. The
  program needs no enthalpy, no radiation field and no stored Q table, and
  therefore links without HDF5.
- `sed/src/cmp_refs.f90` (`make cmp_refs.x`). The DL07 (MW R_V = 3.1,
  b_C = 6e-5) and Zubko ZDA BARE-GR-S emission SEDs at U = 1, each solved in
  both Mathis-field conventions — the corrected 4000 K dilution
  w_4000 = 1.65e-13 and the w_4000 = 1e-13 printed in Mathis et al. (1983) — so
  that the convention behind each published reference SED can be settled by
  direct comparison. Writes `sed/output/cmp_{dl07,zubko}_{corr,orig}.dat`. Both
  models are built through `dust_lib` on this tree's own released optics
  products, the same ones `calc_sed.x dl07` and `calc_sed.x zubko` use.

Both targets are additions to `sed/Makefile`; the default `make` target is
unchanged, and `make cleanall` now removes them.

### `docs/make_figs.py`

New. Regenerates all eight figures of `astrodust_sed_report.tex` from the tree
it sits in and nothing else:

    figs/q_vs_lambda.pdf             Q_ext and Q_abs at four grain sizes
    figs/tau_residual.pdf            our tau_Ad/N_H vs the release extinction.dat
    figs/d03_comparison.pdf          HD23 Ad+PAH vs WD01/D03 (2003, 2009)
    figs/sed_stages.pdf              the three enthalpy stages vs astrodust_irem
    figs/cabs_per_grain_vs_hd23.pdf  single-grain C_abs vs the release Q table
    figs/sed_total_with_pah.pdf      Ad + PAH vs HD23 model_irem
    figs/sed_total_d16_threeway.pdf  the three graphite sources of the xi blend
    figs/polarization.pdf            p_max/N_H vs the release polarized_extinction

Two figures that used to need drivers of their own are now computed in the
script from products already in the tree: the optical-depth residual comes from
the Q table and the size distribution directly, and the single-grain C_abs comes
from the same random-orientation Q table `sed_init` loads, on that table's own
wavelengths, so no interpolation enters the comparison. The polarized figure
forms

    C_pol^ext = 0.5 pi a^2 [ Q_ext(E perp a) - Q_ext(E || a) ]   at k perp a,

which is the sign the release file and the polarized branch's own C_polext
product carry, so p_max comes out positive over the optical and the residual is
measured without an overall sign flip.

Each figure prints the numbers the report quotes for it, so a run of the script
is also the measurement the text is checked against. All text goes through
LaTeX, so no label carries a non-ASCII character. `sed_total_d16_threeway` is
skipped, with a message, until the two non-production graphite runs of
`calc_sed.x astrodust` have been made.

### `mc/` build and namelists

`mc/Makefile`. The Monte Carlo drivers that link the dust-model modules
(`main_mc_sed.x`, `main_pT_compare.x`) failed at link: `SED_PIPELINE` was
missing `sed_paths.f90`, `sedust_h5.f90`, `sedust_product.f90`,
`kext_table.f90` and `q_component.f90`, and the link line carried no HDF5
libraries, so `sed_data_path`, `read_sedust_grid` and `read_sedust_qtable` came
out undefined. The module list is now the same set and order as `SRC_PIPE` in
`sed/Makefile`, and the same `HDF5` / `HDF5_PREFIX` switch is honored, with
`-cpp` added so `-DSEDUST_HDF5` reaches the sources. `make HDF5=0` compiles the
HDF5 paths out, as in `sed/`.

`mc/*.nml`. `qtable_path` pointed at `../tmatrix/output/`, which the data
reorganization emptied; all eight namelists now name
`../data/astrodust/q_astrodust_P0.20_Fe0.00_1.400.dat`, where the table is.

### `pyutil`

`pyutil/tests/test_sedust_h5.py`. Two checks were asserting the wrong contract.
The Zubko file's `q_zubko_*` text products are the D03 Mie recomputation, which
the HDF5 file stores as the `*_mie_d03` groups, so the text comparison is mapped
onto those; the `sil`/`gra`/`pah` groups hold the ZDA tables as distributed,
whose text form is a different format this test does not parse. The expected
component list is now stated per model instead of being waived for astrodust.
The `i_lyman` check now tests the covering rule the Fortran writer's
`lyman_index` actually implements — `i_lyman` is the **last** node at or below
the Lyman limit, so the non-ionizing view starts at or just below 0.0912 um and
a host whose transport floor is the limit interpolates inside the table — rather
than the bracketing rule it had. The suite passes in both trees (65 checks in
v1.00, 71 in v1.20, the extra six covering the polarized branch).

`pyutil/sedust_h5.py`. The module docstring stated the opposite `i_lyman`
definition; corrected to the covering rule above.

### `docs/astrodust_sed_report.tex`

Every residual table and every quoted number was re-measured against the current
code, with `docs/make_figs.py` as the measurement, and the text now says which
of the older numbers are historical records of a superseded state rather than
present-day measurements:

- `tab:mathis_c` (both columns) re-measured; the original-Mathis variant was
  regenerated with `./calc_sed.x astrodust mathis_orig`. The FIR-peak row moves from
  -15.83%/-8.82% to -14.33%/-7.57% and the PAH row from -4.40%/+0.02% to
  -8.16%/-3.67%; the entry explains that the astrodust rows moved with the
  enthalpy correction and the PAH row with that plus the DL07 eq. 2 correction
  above.
- The three-way graphite table of the xi blend is now a single run of the
  present code for all three columns, so they can be read digit for digit
  against each other; the older D03-sphere table in the PAH section is kept and
  labeled as the record of the state in which it was taken.
- The PAH feature-strength comparison is stated feature by feature with the
  band each ratio is taken over, defined as one Drude FWHM to either side of the
  DL07 Table 1 mode center.
- The `jori` convention of the D16 spheroid table, and its effect on the sub-mm
  band ratio, are written up in the xi-blend section.
- References to drivers that no longer exist (`sed/tau_check.x`,
  `sed/src/dump_per_grain.f90`, `tmatrix/driver/fallback.f90`) are replaced by
  what the tree does carry (`docs/make_figs.py`, `rayleigh_limit` in
  `tmatrix/driver/asymptotic_optics.f90`).

### `docs/SEDust_user_manual.tex`

`gra_d03_sphere` and `gra_d16_spheroid` added to the `calc_sed.x astrodust` argument
list; `calc_kext_dl07.x` and `cmp_refs.x` documented among the standalone
programs; the `kext` product paths corrected to the directory each model now
keeps its products in; the
Zubko non-ionizing product's point count and first wavelength corrected to 866
and 0.08998 um, which is the covering rule stated above.
