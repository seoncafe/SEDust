# Plan: orientation-resolved Q from the T-matrix engine, and its comparison against the HD23 astrodust table

**Date:** 2026 July 20

This is the planning document that guided the work. It states the goal, the physics to be
computed, the engine changes required, the output format, the comparison strategy, the
staging, and the open questions.

**Status: implemented (2026 July 20).** All four stages are done, verified against HD23, and
committed; the geometric-optics regime was carried through the projected-area plus Fresnel
model, and the full-grid table is generated and shipped. The as-built account, with the
validation numbers, is in `sedust_polarization_implementation.pdf` §7. This document is kept
for the reasoning and the staging it records.

**Engine status (2026 August 3).** The T-matrix engine this plan calls has since been
rebuilt as a modern-Fortran library. The numbers are unchanged — every product reproduces
byte for byte — but the file layout, the entry-point names, and the ownership of state are
no longer what Sections 3 and 9 describe. The library lives in `tmatrix/src/`:
`tmatrix_kinds.f90`, `tmatrix_status.f90`, `tmatrix_types.f90`, `tmatrix_core.f90`,
`tmatrix_full_direct_bridge.f90`, `fixed_orientation_amplitude.f90`,
`scattering_matrix_expansion.f90`, `tmatrix_api.f90`, and `lpd.f` (the only fixed-form
file, LAPACK/BLAS). Mishchenko's unmodified sources are kept as non-compiled provenance in
`tmatrix/reference/upstream/`; `tmd_one.f` and `ampl_oriented.f` are gone. The astrodust
policy stays outside the library, in `driver/read_index.f90`, `driver/asymptotic_optics.f90`
and the regime selection in `driver/spheroid_optics.f90`. Entry points map as

| this plan | now |
|---|---|
| `TMD_ONE` | `tmatrix_eval(work, a_um, lambda_um, nr, ki, options, result)` |
| `TMD_ONE_SCATMAT` | `tmatrix_eval_scattering_matrix(work, ..., result, scatmat)`, `scatmat` of type `tmatrix_scatmat_t` (`al1..al4, be1, be2, lmax, nmax_tm`) |
| `AMPL(NMAX, DLAM, ...)` | `ampl(work, scatmat, tl, tl1, pl, pl1, alpha, beta, vv, vh, hv, hh, ierr)` — `nmax` and the wavelength come from `scatmat%nmax_tm` and the workspace, not from the argument list |
| `MATR` | `scatmat_from_moments` (module `scattering_matrix_expansion`; returns arrays instead of printing) |
| `HOVENR` | `vdm_hovenier_test` (same module, returns a status) |

No `COMMON`, `EQUIVALENCE` or mutable `SAVE` remains in the library. The converged
T-matrix lives in `work%full_tstore(NPN6,NPN4,NPN4,8)`, fourth index 1-4 for
`Re T11,T12,T21,T22` and 5-8 for the imaginary parts, element for element the old
`COMMON /TMAT/ RT11..IT22`; `GSP` writes its own arrays (`full_gsp_*`) and so no longer
overwrites it, which retires the save-and-restore described in Section 3.2. Evaluation is
one active caller for each workspace, and different workspaces may be used at the same
time; `ampl` takes its workspace `intent(in)` and only reads the stored T-matrix, so
several threads may read one workspace concurrently. A workspace is 83,926,096 byte
(80.04 MiB), which is what a thread costs. Misuse is caught rather than silently wrong: the
workspace carries a generation counter `tm_tag`, incremented on every evaluation whether it
succeeds or fails, `tmatrix_eval_scattering_matrix` copies it into `scatmat%state_tag`, and
`ampl` returns `ierr = 4` when the two disagree (`ierr`: 0 success, 1 angle out of range,
2 no valid T-matrix, 3 `nmax_tm` out of range, 4 generation mismatch).

Three capabilities of Mishchenko's original are absent from the library — polydisperse size
integration (`DISTRB`/`POWER`), the equal-surface-area radius conversion `RAT /= 1`
(`SAREA`/`SAREAC`/`SURFCH`), and the droplet shape `NP = -3` (`RSP4`/`DROP`). None of them
was lost in this conversion; all three were already missing from the earlier f77 wrapper, so
nothing in this plan used them.

What follows is the plan as it stood in July 2026. It is corrected only where it asserts
something about the code that is no longer true; the reasoning, the staging, and the
anticipated file list are left as written.

**Scope note.** The orientation-resolved Q table covered here is one stage toward the
actual goal: computing the polarized radiative transfer of Seon (2018, ApJ 862, 87) —
dust scattering plus dichroic extinction by aligned grains — with properly derived
grain optics. The scattering side, the fixed-orientation Mueller matrix of aligned
grains, is the subject of `tmatrix_aligned_scattering_plan.md`.

**Added after this plan (see §9).** Three things the plan does not describe were built
later and are recorded at the end of this document so it does not read as the whole
story: named-wavelength generation (`lam` / `lamfile` / `lammerge`), the
extreme-ultraviolet axis below the Lyman limit, and a convergence certification for
size parameters above x = 50.

---

## 1. Motivation

The astrodust optics in SEDust split into two channels with two different sources:

- **Total intensity (Cabs, Csca, g)** already come from SEDust's own T-matrix run,
  `tmatrix/output/q_astrodust_P0.20_Fe0.00_1.400.dat`, read through
  `q_table_mod::load_q_table`. This is a random-orientation average and is computed from
  first principles here.
- **Polarization (Cpol, Cpol_ext)** still come from the HD23-shipped orientation-resolved
  table, `data/dielectric/q_DH21Ad_P0.20_Fe0.00_1.400.dat.gz`, read through
  `q_table_jori_mod::load_q_table_jori`. SEDust does not yet produce this table itself.

The reason the polarization channel still depends on the HD23 file is narrow: the current
driver `tmatrix/driver/run_tmatrix.f90` writes only the random-orientation-averaged
`Q_ext, Q_abs, Q_sca` and does not emit the orientation-resolved (jori = 1, 2, 3) blocks
that polarization needs. The T-matrix itself is orientation-independent and is already
computed and stored inside every engine call; only the fixed-orientation cross-section
layer on top of it is missing.

The goal of this plan is to add that layer, produce an orientation-resolved table in the
same format the HD23 file uses, and compare the two so that the polarization channel can
also stand on a first-principles calculation.

---

## 2. Physics: what must be computed

For each wavelength lambda and effective radius a_eff, three orientations of the incident
wave relative to the spheroid symmetry axis `a` are required. The convention is the one
already documented in `sed/src/q_table_jori.f90`:

| jori | geometry        | Rayleigh-limit polarizability seen |
|------|-----------------|------------------------------------|
| 1    | k parallel to a | alpha_b (E transverse to the axis) |
| 2    | k perp a, E parallel to a | alpha_a (E along the axis) |
| 3    | k perp a, E perp a        | alpha_b                    |

Each orientation needs `C_ext`, `C_abs`, and `C_sca`, because the HD23 table stores all
three quantities for all three orientations. The derived combinations follow the existing
convention in `q_table_jori.f90`:

- `Q_pol = 0.5 * (Q(jori=3) - Q(jori=2))`  (polarization cross section)
- `Q_ran = (Q(1) + Q(2) + Q(3)) / 3`  (three-point orientation average)

Two cross sections come from the fixed-orientation amplitude matrix; the third is a
difference:

- **C_ext(jori)** follows from the optical theorem,
  `C_ext = (4*pi/k) * Im[S_forward]`, where `S_forward` is the forward-scattering amplitude
  for the chosen orientation and incident polarization.
- **C_sca(jori)** is the integral of the fixed-orientation phase matrix over scattering
  angle. This is the one genuinely new piece of physics.
- **C_abs(jori) = C_ext(jori) - C_sca(jori)**.

An important consequence for staging: polarized *emission*, which is what a radiative-
transfer host consumes in the far-infrared and sub-millimeter, is governed by C_abs(jori)
(Kirchhoff). In that band the grains sit deep in the Rayleigh regime (x = 2*pi*a/lambda is
of order 1e-3 for a ~ 0.1 um at lambda ~ 100-850 um), where C_abs is analytic per
orientation and C_sca is negligible. The emission-critical band is therefore the easy one.
Polarized *extinction* (dichroism) needs C_ext(jori) across the optical, which is the
T-matrix regime.

---

## 3. Engine work, by size-parameter regime

`run_tmatrix.f90` already branches on x into three regimes; each branch needs an
orientation-resolved counterpart that returns Q(jori=1,2,3) for ext, abs, and sca.

### 3.1 Rayleigh regime, x < 0.1

`tmatrix/driver/asymptotic_optics.f90` already forms the spheroid polarizability components
`alpha_a` (E along the axis) and `alpha_b` (E perpendicular) and then averages them as
`(1/3) alpha_a + (2/3) alpha_b`. The orientation-resolved values are the pre-average
quantities: emit them per jori instead of collapsing them. This is the cheapest branch and
covers the far-infrared and sub-millimeter emission that matters first.

Cross-check built into the physics: at small x, jori = 1 and jori = 3 both see alpha_b, so
`Q(jori=1)` and `Q(jori=3)` must be equal in this limit and diverge only as retardation
grows with x.

### 3.2 T-matrix regime, 0.1 < x < 50

The T-matrix is solved inside every scattering-matrix call. When this was written that call
was `TMD_ONE_SCATMAT` in `tmatrix/src/tmd_one.f`, and Mishchenko's fixed-orientation
amplitude routine `AMPL` with its helper `VIGAMPL` sat unbuilt in `tmatrix/src/ampld.lp.f`,
reading the converged T-matrix out of `COMMON /TMAT/`. One caution learned while
implementing this (Stage B): the random-orientation expansion `GSP`, called at the end of
`TMD_ONE_SCATMAT`, reused the `/TMAT/` storage as scratch through an `EQUIVALENCE` and so
destroyed the T-matrix it read. The block therefore did not hold the converged T-matrix on
return unless it was saved before `GSP` and restored after, and the implementation kept an
intact copy in a second common block for exactly this reason. Neither the overlay nor the
restore exists any more: the library keeps the T-matrix in `work%full_tstore` and gives
`GSP` its own arrays (see the engine status above). The work is:

1. Add `AMPL` and `VIGAMPL` to the build (`tmatrix/Makefile:46`, `SRC_TM`), or fold them
   into `tmd_one.f`. Keep the Mishchenko routine names unchanged so they stay checkable
   against the original.
2. Add an F90 wrapper that (i) runs `TMD_ONE_SCATMAT` to fill `COMMON /TMAT/`, then (ii)
   calls `AMPL` at beta = 0 (jori = 1) and beta = 90 degrees (jori = 2, 3) in the forward
   direction, then (iii) applies the optical theorem for C_ext(jori) and integrates the
   fixed-orientation phase matrix for C_sca(jori), giving C_abs(jori) by difference.

No T-matrix is recomputed; the added cost is the amplitude evaluation and the angular
integral, both small next to the solve that already ran.

### 3.3 Geometric-optics regime, x > 50

`geometric_optics_limit` needs an orientation-resolved counterpart, or this branch is left
as a documented approximation with its domain of validity stated at the code site. The
polarized contribution here is small, so this branch is a candidate for deferral; the
decision is recorded as an open question in section 7.

New routines are to be named for the physics they compute, for example
`rayleigh_oriented`, `tmatrix_oriented`, and `cext_from_forward_amplitude`.

---

## 4. Output format

`q_table_jori_mod::load_q_table_jori` takes the table path as an argument, so if the new
output is written in the same format as the HD23 file, the whole polarization side reads it
with a single path change. That format is 12 header lines followed by the free-format stream

```
((Q(jw, jr, jori), jw = 1..1129), jr = 1..169), jori = 1..3
```

written once for Q_ext, once for Q_abs, once for Q_sca. The grid axes stay in the companion
files `DH21_wave` and `DH21_aeff`.

The existing column-format `q_astrodust_...dat` is left as it is; the orientation-resolved
table is a second, parallel output.

---

## 5. Comparison strategy

Four layers, ordered so that the ones that do not rely on trusting the HD23 file come first.

1. **Internal consistency, and a measurement of the three-point average error.** Compare
   the new `(Q(1) + Q(2) + Q(3)) / 3` against the existing `q_astrodust` table, whose
   random average is the exact continuous orientation average the engine computes. The
   difference is the error of the three-point average that the HD23 table uses; the comment
   at `sed/src/sed_astrodust.f90:1587` refers to exactly this approximation. In the Rayleigh
   regime the two must agree to rounding.

2. **Rayleigh-to-T-matrix continuity.** At the x -> 0.1 boundary the oriented Rayleigh
   closed form and the oriented T-matrix result must agree. This checks the two independent
   code paths against each other with no external reference.

3. **Direct comparison against HD23.** Compare the new table's Q(jori) against
   `q_DH21Ad_...` entry by entry over (lambda, a_eff, jori). Secondary check: at small x,
   Q(jori=1) approx Q(jori=3), diverging as x grows.

4. **End to end.** Run `sed/src/calc_polext.f90` against the new table (change `F_Q` only)
   and reproduce `data/release/polarized_extinction.dat`, reusing the existing deviation
   statistics. Then feed the new table into `sed_astrodust` for the polarized emission the
   host consumes.

---

## 6. Staging

Ordered to bring the radiative-transfer host's need forward.

- **A.** Rayleigh oriented C_abs, giving polarized emission in the far-infrared and
  sub-millimeter. Validated by layers 1 and 4. This alone moves the polarized emission that
  the MoCafe host uses onto a first-principles footing.
- **B.** T-matrix oriented C_ext, giving polarized extinction across the optical. Validated
  by layers 2 and 3.
- **C.** C_sca(jori) from the phase-matrix integral, completing C_abs over all wavelengths;
  settle the geometric-optics branch.
- **D.** Full sweep over the 1129 x 169 grid with three orientations, write the table,
  update the documentation. The amplitude calls are cheap next to the T-matrix solve, so the
  run time is close to the current `q_astrodust` sweep.

---

## 7. Risks and open questions

- **C_sca(jori) integral** is the only genuinely new physics. It should be checked against a
  Mishchenko published test case before it is trusted.
- **Geometric-optics branch:** implement the oriented variant, or keep it as a documented
  approximation. Its polarized contribution is small, but the extinction table is not
  complete without it.
- **HD23 jori = 1 convention.** The HD23 release defines jori = 1 as k parallel to a, a
  single orientation rather than a random-orientation average. This matches the SEDust
  convention, so layer 3 compares like with like.

---

## 8. Files touched (anticipated)

- `tmatrix/Makefile` — add `AMPL`, `VIGAMPL` to the build.
- `tmatrix/src/ampld.lp.f` — Mishchenko fixed-orientation amplitude routine, names unchanged.
- `tmatrix/driver/run_tmatrix.f90` — orientation branches, new output writer.
- `tmatrix/driver/asymptotic_optics.f90` — emit the pre-average Rayleigh components per jori.
- `sed/src/calc_polext.f90` — comparison run against the new table (path change only).
- New comparison utility for layers 1-3 (table difference and the three-point-average error).

---

## 9. Built after this plan

The plan above ends at the full DH21-grid sweep. Three additions came later; the
implementation account is in `sedust_polarization_implementation.tex` (§7.6-§7.8).

### 9.1 Named wavelengths

`range JW1 JW2` takes *indices* into a fixed wavelength file, so it can only reach
wavelengths that file already carries. Polarized transfer needs named wavelengths
instead, so `run_q_jori.x` gained three modes:

```
./run_q_jori.x lam L1 [L2 ...] [ja=JA1:JA2] [tag=NAME]   # wavelengths [um]
./run_q_jori.x lamfile PATH    [ja=JA1:JA2] [tag=NAME]   # one per line, '#' comments
./run_q_jori.x lammerge STEM FILE [FILE ...]             # reassemble ja= windows
```

Each writes `q_astrodust_jori_P0.20_Fe0.00_1.400.TAG.dat` and its companion
`.TAG.wave` axis; the default TAG is `lamL` for one wavelength, `lamLMIN_LMAX_nN`
for several. The pair is handed to `sed_init` / `build_astrodust` through
`qpol_euv_path` and `qpol_euv_wave_path`. The reader interpolates in log(lambda),
so at least two wavelengths are needed to fill a band.

`ja=JA1:JA2` splits the 169 radii over separate **processes**, not threads. The
reason at the time was the engine: Mishchenko's solver passed the converged T-matrix
to `AMPL` through `COMMON /TMAT/` and kept further working storage in COMMON (two
blocks blank), so the core was not re-entrant. That obstacle is gone — the library
holds all of it in a caller-owned workspace, one active caller for each workspace and
different workspaces concurrently — but a workspace costs 80.04 MiB, so memory, not
re-entrancy, now sets how many may run at once; the driver still splits over
processes. Measured before the library conversion: three wavelengths (0.0124, 0.0602,
0.0912 um) over 57 processes took 9m22s wall for 34m48s of processor time — a speedup
of 3.7, not 57, because the few radii just below the x = 60 ceiling carry nearly all
the cost.

`euv` selects a wavelength *axis file* and therefore cannot be combined with
`lam` / `lamfile` / `lammerge`, which carry their own axis; the driver rejects the
combination.

### 9.2 The extreme-ultraviolet band, and what is currently zero

`Cpol`, `Cpol_ext` and `Cbir_ext` cover the whole DH21 axis (1129 nodes,
0.0912-39810 um) — from the release table
`data/dielectric/q_DH21Ad_P0.20_Fe0.00_1.400.dat.gz` by default, or from the
regenerated `tmatrix/output/q_astrodust_jori_P0.20_Fe0.00_1.400.dat.gz` through
`qpol_path`, which carries the 4th block and so is what makes `Cbir_ext` nonzero.
Below the Lyman limit `build_Cpol` looks for the companion table
`q_astrodust_jori_euv_P0.20_Fe0.00_1.400.dat.gz`, **which is not shipped**. When it
is absent the extreme-ultraviolet dichroism and birefringence stay zero and the
reader says so on stderr ("zero by omission, not by physics").

That zero is a documented deficit, not the physics. First-principles size integrals
(`tmatrix/driver/euv_polarized_optics.f90`) give the alignment-weighted
`|C_pol,ext| / C_ext` as 1.26e-3 at the 0.0912 um seam, *rising* to 3.64e-3 near
20.6 eV (a factor 2.9) before decaying to 1.6e-4 at 100 eV, with the sign opposite
to the optical band; the reversal is near 0.106 um. `|C_bir,ext| / C_ext` falls
monotonically from 7.2e-3 at the seam and changes sign near 51 eV. A sphere has
exactly zero dichroic extinction and exactly zero birefringence, so scalar EUV
optics on the volume-equivalent sphere could not have supplied these entries at
all — the shape has to be retained.

With the companion table present, `build_Cpol` interpolates it in log(lambda) and
log(a); a grid outside the table's 0.0124-0.0912 um coverage is rejected with
`status = 9` rather than extrapolated. The axis stops at 0.0124 um (100 eV) because
above x = 60 only the opaque geometric-optics limit remains, and it needs the chord
optical depth 4 Im(m) x to be large: that quantity is 3.7 at x = 50 there and falls
below 1 by 200 eV.

### 9.3 Certification above x = 50

Section 3.3 left the x > 50 branch as a documented approximation. The
extreme ultraviolet pushes size parameters up, so `euv`, `lam` and `lamfile` now
attempt the T-matrix up to x = 60 (`X_TM_MAX`) under a two-setting certification:
each node in (50, 60] is solved at the production `(DDELT, NDGS) = (1e-3, 2)` and
again at `(3e-3, 3)`, and is kept only if the two agree to `TOL_CHK = 5%` on every
channel; a node that fails falls to geometric optics. The plain and `range` sweeps
keep the older x > 50 rule, so they still reproduce the shipped table byte for byte.

- `NDGS = 4` is unusable as the check setting: `NGAUSS = NDGS*NMAX` then exceeds
  `NPNG1 = 300`, so it would fail for want of storage, not of convergence.
- **`IERR = 0` is not sufficient at x >~ 55.** At lambda = 0.0912 um and
  a = 0.9441 um (x = 65.04) the solver reports `IERR = 0` and returns a value 50%
  away from its neighbors. At the seam, x = 51.7 and 54.7 pass the two-setting
  check while 57.9, 61.4, 65.0 and 68.9 fail.
- The polarized channels are differences of order-unity numbers, so a relative
  comparison alone is meaningless near their zeros; `Q_POL_FLOOR = 1e-3` sets the
  absolute scale below which they need not agree (about `|C_pol|/C_ext ~ 5e-4`
  after the size integral, ~1% of the V-band dichroism).

**What it costs.** Leaving x > 60 at the geometric-optics value (which carries
`Q_pol,ext = Q_bir,ext = 0` exactly) loses part of the dichroic extinction: nothing
at the seam, ~8% of the band total at 0.031 um, ~20% at 0.022 um, ~46% at
0.0124 um. The true value lies between the certified figure (lower bound 1.06e-4)
and a 1/x continuation (upper bound 3.14e-4). `C_pol` — the absorption dichroism
that drives polarized emission — has no such gap, because the geometric-optics
limit obtains it from the opaque-grain Fresnel surface integral.

### 9.4 Why the library still reads tables

No T-matrix is solved at run time. The reason is not code size but that convergence
cannot be certified at x >~ 55 from inside a transport loop (§9.3), and a silently
wrong optic is worse than a missing one. Tables are generated offline, where a node
can be rejected and the rejection counted in the file header. The core has since been
rebuilt as a library (see the engine status above); that changed the code structure,
not this reason, so the policy stands.
