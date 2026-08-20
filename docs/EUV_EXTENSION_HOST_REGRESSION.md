# Host-code regression of the EUV extension

**Date:** 2026-08-01
**SEDust revision measured:** `536759d` ("EUV extension"), plus the working-tree
changes to `sed/src/sed_astrodust.f90` and `sed/src/zubko_io.f90`
**Baseline:** the SEDust copy inside MoCafe v2.00 as it stood before the update
(MoCafe git HEAD)
**Host:** MoCafe v2.00, `examples/dustemis/model_compare_{astrodust,dl07,zubko}.in`
**Status: SUPERSEDED, 2026-08-03.** The 2026-08-01 fix was incomplete. It
repaired three of the four places that confuse the wavelength grid with the
radiation field; the fourth, `calc_P` in `sed/src/p_sub.f90`, was found on
2026-08-03 and is described in section 10. Every host measurement in sections
4-6 was taken with that fourth site still defective, so those numbers no longer
describe the library and have to be remeasured. Sections 1-9 are left as
written, as the record of the first defect and of a fix that looked complete.

**Layout note, 2026-08-07.** The paths in this document predate the HDF5
migration and the data-layout rewrite that followed it. The optics are now one
HDF5 file per model under `data/<model>/`, with the text products beside them;
`data/kext_*.dat` and `tmatrix/output/q_astrodust_*.dat` referred to here now
live under `data/astrodust/`. None of that changes the measurements above, which
are about how the wavelength grid and the radiation field are confused in the
solver, not about where files sit. See `HDF5_MIGRATION_STATUS.md`.

What the 2026-08-01 entry claimed, and what is no longer true: that with the
fix in place a Mathis field leaves the zubko output md5 identical to the
baseline, and astrodust and dl07 byte-identical. With the fourth site repaired
as well, a Mathis field moves zubko by 3.05 per cent (6.00 per cent for its EUV
variant) even though its grid length never changes, and all three models move.
Section 10 gives the numbers.

At the time of measurement: no SEDust source or data file was modified. The
reverted-bound library of section 5 was built in a scratch directory and
deleted afterwards.

## 1. Executive summary

The EUV extension stays inactive unless the host passes `lam_min`, and for
`astrodust` and `dl07` the update is a no-op down to roundoff, exactly as
designed. For `zubko` it is not. Without any host code change and without
`lam_min`, the dust emission shifts by 1-2% across the bright part of the SED,
by up to 5.4% in the emission image, and the dust temperature rises by 0.07 K
(median, everywhere positive).

The cause is confirmed rather than inferred (section 5): reverting only the three
single-photon-bound expressions and nothing else reproduces the baseline to
8e-13, the same roundoff level the other two models sit at.

Those three expressions intend to bound the energy of *the hardest photon the
radiation field carries*. What they actually measure is *the short-wavelength
end of the model's own optics grid*. The two coincide for astrodust and dl07,
whose grid stops at the Lyman limit; they differ by a factor of 91 for zubko,
whose DustEM tables start at 1.0e-3 um (1.24 keV) while the field being
transported carries nothing below 0.0912 um.

## 2. What was synchronized

Nine sources brought into the MoCafe copy, after which they are byte-identical
to this tree:

| file | change |
|---|---|
| `sed_astrodust.f90` | EUV grid extension, `planck_integration_grid`, single-photon bounds |
| `stoch_qm.f90` | single-photon bound in `qm_solve_grain` |
| `q_astrodust.f90` | new: astrodust Mie optics from the DH21 dielectric function |
| `qpah.f90` | D03 graphite above the PAH cutoff |
| `q_silicate.f90` | `silicate_index_lambda_range` |
| `q_graphite.f90` | `graphite_index_lambda_range` |
| `q_graphite_d16_sphere.f90` | `D16_LAM_MIN_UM` |
| `zubko_io.f90` | module header only |
| `dust_lib.f90` | module header only |

One data file was added, `data/dielectric/index_DH21Ad_P0.20_0.00_1.400`
(355 kB), read only when the EUV band is active.

## 3. How the comparison was made

Both MoCafe binaries were built with `-no-ipo`, differing only in which
`libsedust.a` they link: one from MoCafe git HEAD (baseline), one from the
sources above. Runs used `par%iseed = 1234`, `OMP_NUM_THREADS=1`, 8 MPI ranks
and `par%use_master_slave = .false.`, the four conditions under which MoCafe is
reproducible run to run.

Six inputs were run on each side, twelve runs in all: the three
`model_compare_*.in` as shipped, where the transport optics come from
`kext_file`, and the same three with `kext_file` blanked so that the transport
optics come from the SEDust model object instead. Every dataset in every output
file was compared by value with h5py; HDF5 container metadata makes a byte
comparison meaningless.

## 4. Results

Size-integrated extinction is unchanged, to every digit printed, from the
values recorded before the update:

| model | nlam | Cext(0.55 um) | albedo | gbar |
|---|---|---|---|---|
| astrodust | 1129 | 2.978E-22 | 0.6978 | 0.6693 |
| dl07 | 1129 | 4.874E-22 | 0.6724 | 0.5373 |
| zubko | 1201 | 4.306E-22 | 0.5435 | 0.4313 |
| from_files | 1201 | 4.245E-22 | 0.5433 | 0.4305 |

Dust emission, largest relative difference over all datasets of a file:

| model | observed image | dust SED / Tdust | J_lambda |
|---|---|---|---|
| astrodust | identical | 1.7e-13 | 6.0e-16 |
| dl07 | identical | 2.8e-12 | 5.6e-16 |
| zubko | identical | **5.9e-02** | 5.9e-16 |

The observed image is identical because these inputs place the dust emission in
its own output file; the scattered and direct light never touch the affected
code. The `J_lambda` agreement at 6e-16 shows the transport itself is untouched:
the difference enters only where the tallied field is turned into emission.

For zubko, in detail:

| quantity | difference |
|---|---|
| Tdust | +0.0705 K median, +0.025 to +0.138 K, max 1.4% relative (T range 6.3-25.5 K) |
| SED_emergent, ten brightest bins | -0.9% to -1.9% |
| DustEmis_image, above 1e-3 of peak | 1.9% median, 5.4% max |
| total emitted energy | -0.08% |

The change is a systematic redistribution, not noise: the temperature rises in
every cell, the bright SED bins all fall, and the total energy is nearly
conserved.

## 5. Confirming the cause

Three expressions introduced by the update raise a single-photon energy bound
from a fixed 13.6 eV to the harder of 13.6 eV and hc/lambda at the short end of
the grid:

| site | expression |
|---|---|
| `stoch_qm.f90:1866-1868`, `qm_solve_grain` | `u_photon_max = HC_CGS / minval(isrf_wl_full)` |
| `sed_astrodust.f90:982`, `sed_grain_loop` | `max(U_UV1_ERG, HC_ERG_UM/lam(1))` |
| `sed_astrodust.f90:1410`, `narrow_iterative` | `max(U_UV1_ERG, HC_ERG_UM/lam(1))` |

A library was built from the updated sources with only these three reverted to
their previous form, everything else left as is, and the two zubko cases rerun.
Against the baseline:

```
gm_zubko_dustsed.h5    maxrel = 2.0e-13
sed_zubko_dustsed.h5   maxrel = 8.2e-13
gm_zubko.h5            maxrel = 0
sed_zubko.h5           maxrel = 0
```

That is the roundoff level astrodust and dl07 already sat at, so these three
expressions account for the entire zubko change and nothing else in the update
contributes to it.

## 6. Why the two definitions diverge for zubko

Both `minval(isrf_wl_full)` and `lam(1)` are the model's wavelength grid, which
is the grid of its optics table. They are not where the radiation field actually
has photons -- a field that is zero over the first half of the grid moves the
bound just as much as one that is not.

| model | grid floor | implied photon energy |
|---|---|---|
| astrodust | 0.0912 um (T-matrix Q table) | 13.595 eV |
| dl07 | 0.0912 um (same Q table) | 13.595 eV |
| zubko | 1.000e-3 um (DustEM tables) | 1239.8 eV |

For astrodust and dl07 the grid floor is the Lyman limit itself, and
13.595 eV < 13.6 eV, so `max()` returns the old constant and nothing moves --
which is precisely the invariance the code comments claim. For zubko the same
`max()` returns 1.24 keV, 91 times the old bound, although in these runs the
transported field stops at 0.0912 um and no photon above 13.6 eV exists.

The consequence is numerical rather than a change of physics. `umax` sets the
upper edge of the grain-enthalpy bin set at fixed bin count, so a 91-fold wider
starting window means coarser bins. The refinement loop in `qm_solve_grain` does
contract it -- `umax = 0.8*ub(jcut) + 0.2*umax` whenever the top-bin probability
falls below `PMIN_UP_QM` -- but it stops at the guards `umax > 1.02*umaxlo` and
`umax > 1.01*ub(jcut)`, and it is capped at `MAX_ITER = 10`. Starting 91 times
too high therefore does not return to the same endpoint, and the residual
coarseness is what the 1-2% shift measures.

## 7. Verdict and options

Removing the hard-wired 13.6 eV is right: when the grid is carried into the EUV
the bound must follow, or the temperature excursions driven by the hardest
photons are clipped. The defect is in what stands in for "hardest photon."

The bound should come from the field, i.e. the shortest wavelength at which
J_lambda is nonzero, rather than from the extent of the optics grid. On an
EUV-extended astrodust grid illuminated to its short end the two agree, so the
intended behavior is preserved; on a grid that reaches far past the illumination
-- zubko today, and any model whose optics table is wider than the transported
band -- the field-based bound is both the physically meaningful one and the one
that keeps the enthalpy bins as fine as before.

Three sites need the same treatment. `qm_solve_grain` has the field in hand
(`isrf_full`) and can take the minimum over the wavelengths where it is nonzero.
`sed_grain_loop` and `narrow_iterative` read `lam(1)` from module state and would
need the same quantity passed in or stored alongside the grid.

Until that is done, `zubko` results computed with this revision and with the
previous one differ at the level quoted in section 4, and the difference is a bin-grid
artifact rather than an improvement.

## 8. The fix

*Added 2026-08-01, after the measurement above.*

`sed/src/radfield.f90` gained

```fortran
pure function hardest_photon_energy(lam_um, J_lam) result(u_photon)
```

which returns `hc / lam_hard`, `lam_hard` being the shortest wavelength at which
`J_lambda` exceeds `J_REL_FLOOR = 1e-12` of its own peak. `J_lam` may be in any
units, since only ratios within the array are used. The three sites of section 5
now read:

| site | expression |
|---|---|
| `stoch_qm.f90`, `qm_solve_grain` | `u_photon_max = hardest_photon_energy(lam_um, j_lam_si)` |
| `sed_astrodust.f90`, `sed_grain_loop` | `max(U_UV1_ERG, hardest_photon_energy(lam, J_lam))` |
| `sed_astrodust.f90`, `narrow_iterative` | `max(U_UV1_ERG, hardest_photon_energy(lam, J_lam))` |

**This table was incomplete.** A fourth site, `calc_P` in `sed/src/p_sub.f90`,
reads the same quantity from the grid and was not changed here. Section 10.

**Why the threshold is 1e-12 of the peak and not an exact `> 0`.** Exact zeros
must be excluded — `J_Mathis` returns exactly zero below 0.0912 um — but an
exact positivity test would accept denormal or roundoff-level residue in the
unilluminated bins of a field handed over by a host, and push the bound back up
by orders of magnitude. A component a fraction 1e-12 below the peak of
`J_lambda` contributes at most that fraction of the photon absorption rate;
allowing two decades for the run of `C_abs` across the illuminated band, the
enthalpy states it could populate carry probability <~ 1e-10 of the peak, at or
below the 1e-13 tail level (`PMIN_UP`, `PMIN_UP_QM`) at which both solvers
already truncate the window. Such a component cannot change a resolved
excursion, and 1e-12 still sits eighteen decades above the residue it is meant
to reject.

The error is one-sided on purpose. Both refinement loops *expand* the window
while the top bin is still populated, so an underestimated bound is corrected as
the loop runs; their contraction is guarded (`umax > 1.02*umaxlo`,
`umax > 1.01*ub(jcut)`) and can stop short of where it began, so an overestimate
is not corrected at all. `MAX_ITER = 10` bounds the correction in either
direction. Overestimating is the dangerous direction, which is why the threshold
is set high enough to reject junk rather than as low as floating point allows.

Reproducing the section 5 comparison, now with the fix in place:

| condition | fixed 13.6 eV vs. field-based |
|---|---|
| Mathis field (stops at 0.0912 um), zubko | output md5 identical — baseline fully recovered |
| field carrying an EUV component, zubko | max relative difference 0.46 |
| astrodust, dl07 | byte-identical |

The first row is exact rather than the 8e-13 of section 5: for a field that
stops at the Lyman limit the `max()` selects `U_UV1_ERG` itself, so the code
path is the same one the hard-coded constant took. `./calc_sed.x zubko` prints this
directly,

```
$ ./calc_sed.x zubko euv
 optics grid short end :  1.0000E-03 um  ->   1.23984E+03 eV
 hardest photon in the field:    1.35946E+01 eV
 single-photon bound in use :    1.36000E+01 eV
```

and `./calc_sed.x zubko euv hardfield`, which puts a hard component in the band
below the Lyman limit, gives `3.70141E+02 eV` for each of the last two lines.

## 9. Keeping the path covered

The defect escaped `calc_sed.x astrodust` and `calc_sed.x dl07` because for those two
models the grid floor and the field's short end are the same number. No driver
in the tree ran the emission solvers on a model whose optics grid reaches past
the illuminated band, which is the reason a grid-based bound survived review.
`sed/src/calc_sed.f90` closes that gap and is built by the default `make`:

```
./calc_sed.x zubko [settings ...]
```

It takes the run settings the three model drivers share (solver, grid, field,
emission term, transition-matrix sizes); `euv` and `hardfield` are the two that
matter here. It is the emission counterpart of `./calc_kext.x zubko`, which
exercises only the extinction size integral. `euv` builds the model on the whole
ZDA optics range instead of cutting it at the Lyman limit; the Mathis field
still stops there, so the bound must stay at 13.6 eV however far the optics grid
extends. `hardfield` then fills the band below the Lyman limit with a diluted
1e5 K blackbody, and the bound must follow the field up. Running `euv` and
`euv hardfield` is the regression.

## 10. The fourth site

*Added 2026-08-03.*

Section 8 listed three places that took the single-photon bound from the
wavelength grid. There were four. `calc_P` in `sed/src/p_sub.f90` integrates
the transitions into the highest enthalpy bin -- every photon harder than that
bin -- and took the lower limit of that integral from the grid as well:

```fortran
dlnwav = log(wavl/lambda(1))/(nwav-1)     ! lambda(1): where the grid ends
```

Two consequences, both measured rather than inferred:

**A photon that was not there.** `interp1` in `sed_mathlib.f90` clamps outside
its range: below `x(1)` it returns `y(1)`. `wavl = hc/ener` for the top bin is
always shorter than the grid floor, so on the old 0.0912 um grid every top-bin
transition was fed `J(0.0912) = 1.359` -- the value at the Lyman limit -- for a
band in which `J_Mathis` is exactly zero. The emission solvers were absorbing
photons the field does not carry.

**An integral that ran backwards.** When `ener > u_hard` the span is negative,
so `dlnwav < 0` and the loop subtracted transition rate instead of adding it.

Both are repaired. The lower limit is now `hardest_photon_energy(lambda,
Jfield)`, the same function section 8 introduced; `Jwavl` is set to zero when
`ener > u_hard`; and the loop is skipped entirely in that case. The fixed count
of 51 sample points became `min(301, max(51, ceiling(span/0.02) + 1))`, so the
sampling follows the width of the illuminated band instead of the width of the
grid.

### What changed, decomposed

Four builds were run to separate the causes: old grid with the old `calc_P`
(the committed results), old grid with only the ghost photon removed, old grid
with the whole fix, and the 12.4 keV grid with the whole fix. Largest relative
difference over each file:

| step | S1 | S2 | PAH |
|---|---:|---:|---:|
| ghost photon removed | **9.98%** | **12.08%** | **3.93%** |
| sampling rule | 0.0000% | 0.0000% | 0.0000% |
| grid extension itself | 0.0088% | 0.0085% | 1.011% |

The decomposition is complete: the steps sum to the total. Nearly all of the
change is the ghost photon. The sampling rule contributing nothing is a result,
not a null: once the lower limit is right, 51 points were already converged.
Extending the table to 12.4 keV moves the astrodust SED by 0.009 per cent.

### Energy conservation

Emitted power over absorbed power, the same field and the same size
distribution on both sides:

| | S1+PAH | S2+PAH |
|---|---:|---:|
| before (ghost photon present) | +0.236% | -0.125% |
| after | **-0.057%** | -0.108% |

The old code is the worse of the two, which is the direction a spurious
absorption term predicts.

### Which models are affected

All three, and zubko independently of anything to do with the astrodust table.
Its grid starts at 1.0e-3 um and its length does not change here, yet its SED
moves by 3.05 per cent (6.00 per cent for the EUV variant) purely from this
fix -- for the reason section 6 already gave, now applied to a fourth site.
DL07 moves most in integrated power, -2.81 per cent.

**Any SED produced before 2026-08-03 contains photons the radiation field did
not carry**, in the transitions into the top enthalpy bin. That includes the
host measurements in sections 4-6 of this document.

### What still has to be remeasured

Sections 4-6 record MoCafe runs (8 MPI ranks, `iseed = 1234`,
`OMP_NUM_THREADS = 1`, `examples/dustemis/model_compare_*.in`). They cannot be
redone inside SEDust, and the SEDust copy inside MoCafe is older than this fix,
so the host has to be updated first. Specifically:

- the dust-emission table of section 4 (astrodust 1.7e-13, dl07 2.8e-12, zubko
  5.9e-02) and the zubko detail table below it;
- the section 8 verification table, whose "md5 identical" and "byte-identical"
  rows no longer hold;
- section 9's coverage claim, which should now say that a driver exercising the
  top-bin integral on a model whose grid reaches past the illuminated band is
  what was missing.

Size-integrated extinction (the first table of section 4) does **not** need
remeasuring: it is read from the `kext_*.dat` tables and never passes through
`calc_P`. Its `nlam` entries do change, to 1762 for astrodust and dl07, because
the Q table now carries the grid to 12.4 keV.

### How it escaped the 2026-08-01 review

Section 9 asked why no driver caught the first defect and answered: for
astrodust and dl07 the grid floor and the field's short end are the same
number. That answer was right and also incomplete. The search that produced the
section 8 table looked for expressions computing a *bound* -- `umax`, a window
edge. `calc_P` does not compute a bound; it computes an integral, and its lower
limit is the same physical quantity wearing a different name. Grepping for the
concept rather than for the idiom would have found it.

