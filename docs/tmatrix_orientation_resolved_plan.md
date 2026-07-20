# Plan: orientation-resolved Q from the T-matrix engine, and its comparison against the HD23 astrodust table

**Date:** 2026 July 20

This is a planning document, not a specification or a report of completed work. It states
the goal, the physics to be computed, the engine changes required, the output format, the
comparison strategy, the staging, and the open questions. Nothing here has been implemented
yet.

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

The T-matrix is solved inside every `TMD_ONE_SCATMAT` call (`tmatrix/src/tmd_one.f`).
Mishchenko's fixed-orientation amplitude routine `AMPL` (`tmatrix/src/ampld.lp.f:535`), with
its helper `VIGAMPL` (`:822`), reads the converged T-matrix from `COMMON /TMAT/` but is not
currently in the build. One caution learned while implementing this (Stage B): the
random-orientation expansion `GSP`, called at the end of `TMD_ONE_SCATMAT`, reuses the
`/TMAT/` storage as scratch through an EQUIVALENCE and so destroys the T-matrix it reads. The
block therefore does not hold the converged T-matrix on return unless it is saved before
`GSP` and restored after; the implementation keeps an intact copy in a second common block
for exactly this reason. The work is:

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
   random average is the exact continuous orientation average from `tmd_one`. The
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
