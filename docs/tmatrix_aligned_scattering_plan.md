# Plan: fixed-orientation scattering matrix of aligned astrodust grains

**Date:** 2026 July 20

**Status: implemented (2026 July 20).** All three stages of Section 6 are built
and verified. As-built summary: the engine
(`tmatrix/driver/scattering_matrix_oriented.f90`), the shared physics module and
production driver (`aligned_population_optics.f90`, `run_scatmat_aligned.f90` ->
`run_scatmat_aligned.x`), and the library API (`sed/src/scatmat_aligned.f90`,
module `scatmat_aligned_mod`) with a reference consumer
(`sed/rt_example/use_dustlib_scatmat.f90`). Anchors (Section 5): closure
`INT Z11 dOmega` vs `C_sca(jori)` `2.7e-7`; optical theorem vs `C_ext(jori)`
`2e-14`; orientation average vs GSP `F(Theta)` `1.4e-6` including the `F12`/`F34`
signs; analytic dipole `2e-16`; dipole-vs-T-matrix continuity at `x ~ 0.1`
`0.5-1.1%`; `phi` mirror `3.5e-7`; equatorial mapping `8.8e-16`; size-integrated
`K` vs independent jori-table integrals `Cext 1.6e-6`, `Cpol 2.0e-5`,
`Cbir 1.3e-7`; Rayleigh `K(theta_i)` sin^2 law exact to `3.3e-16`. Library checks
on the production table: reader round-trip and symmetry reconstruction bitwise 0,
eta algebra exact, `K(90)` vs jori integral vs `dust_extinction` `Cpol 1.6e-5` /
`Cbir 8.0e-5`, `rt_example` closures `0.997-1.001`. Production run: five UBVRI
bands (0.36, 0.44, 0.55, 0.64, 0.79 um), ~8 min at 32 threads, file 138 MB
(gzip 38 MB); the `x > 50` omitted scattering fraction measured at most `3.1e-15`.
The full account is in `docs/sedust_polarization_implementation.tex`
(Section "Aligned-grain scattering matrix and the extinction matrix").

## 1. Motivation

Seon (2018, ApJ 862, 87) modeled the optical polarization of edge-on galaxies with
radiative transfer that includes both dust scattering and dichroic extinction. Doing
that calculation properly requires the angle-resolved scattering (Mueller) matrix of
*aligned* spheroidal grains, not only the random-orientation matrix: an aligned grain
polarizes scattered light differently depending on how the incident ray is inclined
to its symmetry axis, and this is the scattering-side counterpart of the dichroic
extinction SEDust already provides.

SEDust v1.20 supplies, for the astrodust spheroid,

- orientation-resolved cross sections C_ext, C_abs, C_sca (jori = 1, 2, 3), computed
  from first principles and validated against HD23;
- the birefringence optic (4th table block);
- the *random-orientation* scattering matrix (`run_scatmat.x`).

What is missing is the fixed-orientation Mueller matrix Z as a function of both the
incidence direction relative to the grain axis and the scattering direction. This
plan adds it. The building blocks already exist: `AMPL` (src/ampl_oriented.f) returns
the 2x2 complex amplitude matrix for arbitrary incidence angles, scattering angles,
and grain Euler angles, reading the converged T-matrix from COMMON /TMAT/, which
`TMD_ONE_SCATMAT` leaves valid (the /TMATK/ restore).

## 2. Physics

### Geometry and Stokes convention

Grain symmetry axis a-hat along z (Euler ALPHA = BETA = 0). Incident direction at
polar angle theta_i from the axis, azimuth 0. Scattered direction (theta_s, phi).
Because the grain is axisymmetric, only the azimuth difference matters, so

    Z = Z(theta_i; theta_s, phi),    16 elements, in um^2 sr^-1.

Stokes vectors are defined in Mishchenko's convention: the (v, h) = (theta-hat,
phi-hat) meridional basis of each propagation direction in the grain frame, with
Q = I_v - I_h. The radiative-transfer code rotates Stokes vectors into and out of
this frame at each scattering event (the grain frame's z is the local alignment
axis, i.e. the magnetic field direction).

The phase matrix follows from the amplitude matrix S = [[S11,S12],[S21,S22]] =
[[VV,VH],[HV,HH]] by the standard bilinear combinations (Mishchenko, Travis & Lacis
2002, Eqs. 2.106-2.121). AMPL returns amplitudes carrying the dimension of length,
so |S|^2 is an area, and Z is a differential scattering cross section.

### Symmetries (storage reduction and checks)

- The plane containing the axis and the incident direction is a symmetry plane of
  the scattering problem: phi -> 360 - phi maps Z onto itself with sign flips of the
  two off-diagonal 2x2 blocks (elements 13, 14, 23, 24, 31, 32, 41, 42). Store
  phi in [0, 180] only.
- The oblate spheroid is symmetric under z -> -z: theta_i -> 180 - theta_i maps onto
  the stored range with theta_s -> 180 - theta_s and phi unchanged (the numerical
  verification rejected the phi -> 180 - phi guess written before implementation;
  the off-diagonal 2x2 blocks flip sign, the exact pattern is in the table header).
  Store theta_i in [0, 90] only. Physically this also encodes that the alignment
  axis is headless.
- At theta_i = 0 the problem is axisymmetric: Z must become independent of phi and
  take the six-element block structure. This is a check, not an assumption.

### Size-parameter regimes

- x = 2 pi a_eff / lambda < 0.1: analytic electric-dipole Mueller matrix from the
  spheroid polarizability tensor (alpha_a, alpha_b of asymptotic_optics.f90).
  The scattered amplitude is S proportional to k^2 (n-hat x (n-hat x alpha E)); all
  16 elements follow in closed form.
- 0.1 <= x <= 50: T-matrix via TMD_ONE_SCATMAT + AMPL.
- x > 50: not implemented. For the optical bands shipped (0.36-0.79 um) this regime
  is reached only by grains carrying a fraction ~1e-14 of the scattering weight;
  the driver skips such bins and reports the skipped weight so the omission is
  provably negligible. It stops if the skipped weight is not negligible.

### Size integration and alignment weighting

The physically correct decomposition for partial alignment (the same picture used
for polarized extinction) is a perfectly aligned fraction f_align(a) plus a randomly
oriented remainder:

    Z_total = INT da n(a) [ f_align(a) Z_oriented(a) + (1 - f_align(a)) F_random(a) ]

The driver therefore writes two products per band:

1. **Aligned part** Z_al(theta_i; theta_s, phi) = INT da n(a) f_align(a) Z(a; ...),
   absolute units um^2 sr^-1 per H.
2. **Unaligned remainder** F_unal(Theta) = the run_scatmat random-orientation matrix
   accumulated with weight n(a) (1 - f_align(a)) instead of n(a), same 6-element
   format and normalization as the existing scatmat file.

One TMD_ONE_SCATMAT call per (band, size) serves both: the GSP expansion
coefficients feed the unaligned accumulation, and the /TMAT/ left valid by the
/TMATK/ restore feeds the AMPL loop.

f_align defaults to the HD23 power-law fit (`falign_hd23` in sed/src/q_table_jori.f90,
reused directly). A regenerated run can substitute a different profile; the run
used is recorded in the header.

### Alignment dependence at run time

The RT host must be able to vary the alignment degree — between runs (a different
f_align profile) and across space (cell-to-cell). Both follow from the linearity of
the size integral in f_align, without approximating the physics:

- **Cell-to-cell scaling.** For a local alignment degree written as
  f_cell(a) = eta * f_ref(a) with one scalar eta in [0, 1] per cell, the aligned
  part scales exactly: Z_al,cell = eta * Z_al,ref, and the unaligned remainder is
  F_unal,cell(Theta) = F_tot(Theta) - eta * F_ref(Theta), where F_tot is the
  n(a)-weighted random-orientation matrix (the existing scatmat product) and
  F_ref is the same matrix weighted by n(a) f_ref(a). So three stored integrals
  give the exact scattering optics of any eta.
- **Profile changes.** A different f_ref(a) (a_align, alpha_align, f_max, or a
  tabulated RAT profile) requires re-integrating over size. The generator therefore
  also writes a size-resolved intermediate (single-size Z and F on the size grid,
  not tracked in git); re-integrating it under a new profile takes seconds and
  needs no new T-matrix computation. If the intermediate file is available, the
  library performs this integration itself at load time with the model's current
  f_align, so a profile set through `dust_set_alignment` /
  `dust_set_alignment_profile` is honored exactly.

## 3. Products

- `tmatrix/driver/scattering_matrix_oriented.f90` — module with
  `mueller_matrix_fixed_orientation` (T-matrix regime; AMPL + bilinears) and
  `rayleigh_mueller_matrix_oriented` (analytic dipole limit).
- `tmatrix/driver/run_scatmat_aligned.f90` -> `run_scatmat_aligned.x` — size-
  integrated production driver; CLI like run_scatmat (wavelengths, `test`, default
  UBVRI bands 0.36, 0.44, 0.55, 0.64, 0.79 um). Writes
  `output/scatmat_aligned_astrodust_P0.20_Fe0.00_1.400.dat` (aligned part) and
  `output/scatmat_unaligned_astrodust_P0.20_Fe0.00_1.400.dat` (remainder).
- `tmatrix/driver/compare_scatmat_aligned.f90` — verification program (Section 5).
- Makefile targets for both; docs updated after verification.

### Angular grid of the aligned table

theta_i = 0(5)90 [19 values] x theta_s = 0(1)180 [181] x phi = 0(5)180 [37],
16 elements per node. About 127k rows per band. The 1-degree theta_s step keeps
the forward peak resolved at x ~ 5; theta_i and phi vary smoothly and tolerate
5-degree steps. Header records, for every theta_i, the closure integral
C_sca_al(theta_i) = INT Z11 dOmega (unpolarized incidence) for normalization and
consistency. Final file size and whether to track it in git are decided after
generation.

Alongside Z the products carry the extinction-matrix elements on the same
theta_i grid: C_ext(theta_i) for the two polarizations and the corresponding
birefringent real parts, taken from the forward amplitudes. In the Rayleigh
limit these must reduce to the analytic sin^2 dependence between the jori
values (anchor); at theta_i = 0 and 90 degrees they must reproduce the jori
table exactly.

## 4. Consumption by radiative transfer

### Library contract: alignment input -> matrices for the Peest formalism

The division of labor with the radiative-transfer code (MoCafe) is: **MoCafe
supplies the alignment assumptions — the alignment profile and the local
magnetic-field geometry — and SEDust returns the matrices; MoCafe then carries
out the polarized transfer** following the framework of Peest et al. (2017,
A&A 601, A92) as extended to aligned spheroids, dichroism, and birefringence by
Peest et al. (2023, A&A 673, A112, MCPOL). Frame rotations, azimuth sampling,
peel-off, and the exp(-K tau) transfer step are RT-side; every material
quantity in matrix form is SEDust-side.

The MoCafe side of this contract does not exist yet; MoCafe will be updated
separately to consume this interface. Until then, the deliverables that define
the contract for that future development are the library API itself, the
handoff document (mocafe_polarization_handoff.md, which becomes the
development specification), and an rt_example program acting as a minimal
reference consumer of every call.

The alignment state is specified by the three parameters of the RAT-motivated
power law already adopted in SEDust (following the Hoang et al. picture of
radiative-torque alignment: an efficiency ceiling, a critical size, and a
rolloff sharpness),

    f_align(a) = f_max / (1 + (a_align / a)^alpha_align),

set through the existing `dust_set_alignment(m, f_max, a_align, alpha_align)`
(or an arbitrary tabulated profile through `dust_set_alignment_profile`),
optionally scaled cell-to-cell by the scalar eta. Given these inputs the
library returns, per wavelength, with theta_i the angle between the propagation
direction and the local alignment axis (the magnetic field):

    K(theta_i)                  4x4 extinction matrix, um^2 per H: C_ext(theta_i)
                                on the diagonal, the dichroic C_pol(theta_i) in
                                the IQ block, the birefringent C_bir(theta_i) in
                                the UV block
    Z_al(theta_i; theta_s, phi) phase matrix of the aligned population,
                                16 elements, um^2 sr^-1 per H
    F_unal(Theta)               randomly oriented remainder, 6 elements
    C_sca_al(theta_i), C_sca_unal    (albedo and polarization-dependent optical
                                depth follow from K and these)

computed with the f_align(a) currently set on the model — not with a profile
frozen into a shipped file. Before this plan, the three alignment inputs
reached only the polarized extinction and polarized emission channels, and no
scattering-matrix API existed.

K(theta_i) is stored on the theta_i grid from the same forward amplitudes the
generator computes anyway, which makes the extinction matrix exact at every
propagation angle. This replaces the sin^2(psi) interpolation between the
psi = 0 and 90 degree values that the three-orientation jori table forces, and
reduces to that law analytically in the Rayleigh limit — which serves as a
verification anchor.

### API structure: initialize once, query along the path

A Monte Carlo photon crosses many cells, and the polarized-optics calls sit in
the innermost loop. The existing library already follows the required shape
for unpolarized transport — `sed_init` / `build_astrodust` load and integrate
once, and `dust_extinction` / `dust_emission` then serve each cell's queries from
memory. The polarized-scattering layer adopts the same lifecycle, split into
two strictly separated layers:

- **Initialization (once per run, before the photon loop).** File reads, the
  size integrals under the model's current f_align, unit conversions, grid
  bookkeeping (cos theta abscissae, the phi-mirror extension), and any sampling
  aids (e.g. cumulative distributions of Z11 over theta_s per theta_i) are all
  done here and kept in memory. This layer is not thread-safe and is called
  from serial code.
- **Path queries (per step / per scattering event).** Pure reads of the
  in-memory structures: no file access, no allocation, no size integration, no
  state mutation — safe to call concurrently from OpenMP photon threads. The
  cell dependence enters only through scalars the RT already has in hand: the
  local eta (linear scaling, applied by the caller or passed as an argument
  that multiplies a preloaded array element) and theta_i, which the RT obtains
  from one dot product k-hat . B-hat. Nothing about a cell requires
  precomputation inside SEDust.

In addition to the interpolating query routines, the loaded grids themselves
are exposed read-only (the pattern already used for the extinction table), so
a host that wants to inline its own interpolation or build its own sampling
tables can take the arrays once at startup and never call back into the
library during transport.

At each scattering event off the aligned population: express the photon direction
in the local grain frame (z = B-hat) to get theta_i, sample (theta_s, phi) from the
Z11 sky brightness, transform the Stokes vector with the full 4x4 Z (meridional-
basis rotations on both sides), and renormalize by C_sca_al(theta_i). The unaligned
population scatters with F_unal exactly as with the existing scatmat file. Extinction
along the ray already uses C_ext +/- C_pol (dichroism) and C_bir (birefringence).
The handoff recipe in docs/mocafe_polarization_handoff.md is updated with this.

## 5. Verification anchors

A. **Closure vs the jori table (single size).** INT Z11 dOmega at theta_i = 0 must
   equal C_sca(jori=1); at theta_i = 90, INT (Z11 + Z12) dOmega = C_sca(jori=2) and
   INT (Z11 - Z12) dOmega = C_sca(jori=3) (incident E parallel / perpendicular to
   the axis are pure +Q / -Q states in the meridional basis). Reference values from
   the 4-block table already validated against HD23.
B. **Optical theorem wiring.** (4 pi / k) Im S11(forward) and Im S22(forward) at
   theta_i = 0, 90 must reproduce C_ext(jori) exactly as tmatrix_oriented.f90
   computes them — validates the Euler-angle/geometry parameterization.
C. **Random-orientation average.** Averaging Z over grain orientations (Gauss
   quadrature in cos beta, uniform alpha) at fixed scattering geometry must
   reproduce the GSP-expansion F(Theta) of TMD_ONE_SCATMAT for the same size.
   This is the strongest end-to-end check of the bilinear construction.
D. **Rayleigh limit.** The analytic dipole matrix against the T-matrix matrix at
   x near 0.1 (continuity), and against closed-form Rayleigh expressions at
   theta_i = 0 (exactness of the dipole implementation).
E. **Symmetries.** Numerical verification of the phi mirror relation, the
   equatorial mapping, and the theta_i = 0 block structure.
F. **Size-integrated closure (production table).** Header C_sca_al(theta_i = 0, 90)
   against the same f_align-weighted size integrals computed independently from the
   4-block jori table (the calc_polext-style quadrature).

## 6. Staging

- **Stage 1** — engine module + compare program, anchors A-E at single sizes.
  Gate: all anchors pass.
- **Stage 2** — production driver, test mode, anchor F, then the 5-band production
  run (process-parallel per band if needed).
- **Stage 3** — documentation: README, sedust_polarization_implementation.tex,
  aligned_grain_polarization.tex, SEDust_user_manual.tex, and the MoCafe handoff.

## 7. Cost estimate

Per (band, size) in the T-matrix regime: one TMD_ONE_SCATMAT (seconds at most) plus
~127k AMPL evaluations (tens of microseconds each) — minutes per size at the largest
x. With ~100 populated sizes per band and 5 bands, a few hours serial; the driver
prints progress and the bands can run as separate processes.
