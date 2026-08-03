# libtmatrix

`libtmatrix.a` provides a SEDust-independent Fortran API for one randomly
oriented particle.  It accepts size, wavelength, real/imaginary refractive
index, and Mishchenko shape options; it does not read a dielectric table or
apply astrodust-specific asymptotic limits.

Build and validate it from this directory:

```sh
make libtmatrix.a example test smoke test-full-direct-parallel test-api-parallel
```

Source layout is separated by responsibility:

- `src/` contains the production library; `lpd.f` is its only remaining
  fixed-form source.
- `driver/` contains the size-parameter regime selection (`spheroid_optics.f90`,
  `asymptotic_optics.f90`) and the standalone astrodust table generator.
- `reference/upstream/` contains the unmodified Mishchenko standalone sources
  and their parameter include files. They are provenance only and are never
  compiled.

An external consumer compiles with the module directory and links the archive:

```sh
gfortran -I/path/to/tmatrix consumer.f90 \
  -L/path/to/tmatrix -ltmatrix -o consumer.x
```

Use `tmatrix_api` and follow this lifecycle:

```fortran
type(tmatrix_workspace_t) :: work
type(tmatrix_options_t) :: options
type(tmatrix_result_t) :: result

call tmatrix_workspace_init(work)
call tmatrix_eval(work, a_um, lambda_um, n_real, n_imag, options, result)
if (result%status /= 0) then
   ! inspect result%message, result%legacy_ierr, result%lapack_info
end if
call tmatrix_workspace_finalize(work)
```

`shape=-1` is a spheroid and `aspect_ratio=b/a`; values above one are oblate
in the convention used by the Mishchenko implementation.  The raw library does
not decide which size-parameter regime a particle belongs to.  That physics
lives in `driver/spheroid_optics.f90`, whose `spheroid_q` covers the whole
size-parameter range — Rayleigh dipole limit, T-matrix, geometric optics limit —
and returns the regime flag.  `driver/run_tmatrix.f90` is then only the
standalone astrodust-table consumer that supplies the grid and the refractive
index.  Any other code that wants the same optics should call `spheroid_q`
rather than reimplement the regime selection.

## Random-orientation scattering matrix

`tmatrix_eval` returns cross sections only. `tmatrix_eval_scattering_matrix`
runs the identical calculation and additionally fills a `tmatrix_scatmat_t`:

```fortran
type(tmatrix_scatmat_t) :: scatmat

call tmatrix_eval_scattering_matrix(work, a_um, lambda_um, n_real, n_imag, &
                                    options, result, scatmat)
```

- `scatmat%al1 … al4`, `be1`, `be2` are the generalized-spherical-function
  expansion coefficients of the random-orientation scattering matrix in the
  Mishchenko / Hovenier convention, each `tmatrix_expansion_size` (201)
  elements long. Element 1 is the `l = 0` term, the normalization is
  `al1(1) = 1`, and entries beyond `lmax+1` are zero.
- `scatmat%lmax` is the truncation order; the highest populated index is
  `lmax+1`.
- `scatmat%nmax_tm` is the multipole truncation order of the converged
  T-matrix, that is the highest `N` and the highest azimuthal order `M` for
  which the workspace holds T-matrix elements.

`result` carries the same cross sections and the same status as
`tmatrix_eval`. A nonzero `result%status` leaves `scatmat` all zero with
`lmax = 0` and `nmax_tm = 0`.

`scattering_matrix_expansion` turns those coefficients into angles and back:

- `scatmat_from_moments(a1, a2, a3, a4, b1, b2, lmax, npna, theta, f11, f22, f33, f44, f12, f34)`
  evaluates the six scattering-matrix elements on `npna` angles equally
  spaced from 0 to 180 degrees. `f11` carries the normalization of `a1`;
  `f22`, `f33`, `f44`, `f12`, and `f34` are not divided by `f11`, so the
  degree of linear polarization for unpolarized incident light is
  `-f12/f11`.
- `vdm_hovenier_test(l1, a1, a2, a3, a4, b1, b2, kontr, lviol)` applies the
  van der Mee & Hovenier necessary conditions to `l1 = lmax+1` coefficients
  and returns `kontr = 1` when they hold, `kontr = 2` when they do not, with
  `lviol` the lowest offending `l` (or `-1`). It is meant for coefficients
  accumulated over a size distribution; for a single size the conditions are
  those the expansion routine already satisfies.

Both are stateless array arithmetic: no workspace, no shared state, safe to
call concurrently.

## Fixed-orientation amplitude matrix

`ampl` returns the 2x2 amplitude matrix `(vv, vh, hv, hh)` of the same
particle in a fixed orientation. Angles are in degrees: `tl`, `tl1` the
incidence and scattering polar angles, `pl`, `pl1` the corresponding
azimuths, `alpha`, `beta` the Euler angles orienting the particle symmetry
axis. Each element carries the dimension of length.

```fortran
call ampl(work, scatmat, tl, tl1, pl, pl1, alpha, beta, vv, vh, hv, hh, ierr)
```

**`ampl` computes nothing from particle parameters.** It reads the converged
T-matrix that the numerical core left in `work%full_tstore`. The truncation
order and the wavelength of the multipole sum are therefore not arguments:
they come from `scatmat%nmax_tm` and from the workspace, both written by the
evaluation that produced the stored T-matrix. A caller cannot hand them a
value belonging to some other particle.

**The one obligation left to the caller** is that `work` and `scatmat` belong
together: the same `work` must have completed a successful
`tmatrix_eval_scattering_matrix`, and `scatmat` must be what that call
returned. Nothing else is required, and a broken pairing is reported rather
than silently computed:

| `ierr` | meaning |
| --- | --- |
| 0 | success |
| 1 | an angular parameter is outside its allowable range (`alpha` and the azimuths outside 0–360, `beta` and the polar angles outside 0–180). The original reported this through `WRITE` and `STOP`; a library does neither |
| 2 | the workspace holds no converged T-matrix — nothing has been evaluated on it yet, the last evaluation failed, or it has been re-initialized |
| 3 | `scatmat%nmax_tm` outside `1 … 80`, i.e. outside the stored T-matrix block |
| 4 | `scatmat` belongs to a T-matrix that a later evaluation on this workspace replaced |

Every nonzero `ierr` returns a zero amplitude matrix.

Code 4 is what makes interleaving safe. Each evaluation stamps the workspace
with a new generation tag and copies it into the `scatmat` it returns, so any
`scatmat` obtained earlier stops matching. `tmatrix_eval` returns no `scatmat`
and so can never be the evaluation an `ampl` call refers to, but it advances
the tag all the same, and so does a failed evaluation, because a failure may
already have overwritten part of the storage. A caller that sweeps many
scattering directions for one particle therefore needs no bookkeeping of its
own; a caller that interleaves particles either re-evaluates before each
`ampl` call or keeps one workspace for each particle, and learns immediately
if it does neither.

`ampl` takes both `work` and `scatmat` as `intent(in)` and keeps all its
scratch in automatic local storage, so several threads may call it
concurrently on one workspace — useful when sweeping many scattering
directions for one particle. Evaluating a particle, by contrast, writes the
workspace and admits only one caller.

## Fixed-orientation phase matrix

`phase_matrix_from_amplitude` turns that amplitude matrix into the 4x4 phase
(Mueller) matrix of the Stokes vector `(I, Q, U, V)` written in the meridional
basis `e_1 = theta-hat` (V), `e_2 = phi-hat` (H), which is the basis `ampl`
rotates its amplitudes into:

```fortran
call phase_matrix_from_amplitude(vv, vh, hv, hh, z)   ! z(4,4)
```

These are Mishchenko's `Z11 … Z44` bilinears, Eqs. (13)–(29) of Appl. Opt. 39,
1026 (2000), which `reference/upstream/ampld.lp.f` writes inline in its main
program. `z` carries the square of the units of the amplitude matrix, so with
the wavelength in microns it is a differential scattering cross section in
`um^2 sr^-1`.

The routine reads no workspace and holds no state, and the orientation reaches
it only through the amplitudes it is handed, so an orientation average may
accumulate it at any `alpha`, `beta` and any number of threads may call it at
once.

## Concurrency contract

`tmatrix_eval` is the default full-direct, reentrant backend. It allocates the
remaining T-matrix and GSP components on first use, then evaluates entirely
against the caller-owned `tmatrix_workspace_t`. Different workspaces may be
used concurrently by OpenMP or another threading runtime; a single workspace
must not be used by more than one active call.

No default evaluation enters a legacy bridge or an OpenMP critical region. The
standalone table driver and `examples/evaluate_one.f90` both call this default
API. `make test-full-direct-parallel` compares serial and concurrent default
results with 1, 2, 4, and 8 independent workspaces; `make test-api-parallel`
does the same through the public API regression harness; and
`make test-scattering-matrix-parallel` does the same for
`tmatrix_eval_scattering_matrix` followed by `ampl`, comparing IEEE bit
patterns rather than a tolerance.

## Modern full-direct core

The default full-direct backend is implemented in
`src/tmatrix_core.f90`. It is a free-form module whose public
`tmd_one_full_direct` procedure and every numerical helper are module
procedures, so their internal calls have explicit interfaces. The unmodified
Mishchenko standalone sources it was translated from are kept, uncompiled, in
`reference/upstream/`.

The module owns the former `tmd.par.f` limits as private parameters and uses
`wp` (`real64`) and `real32` for the original `REAL*8` and `REAL*4` storage.
Formulae, loop order, storage ranks, and caller-owned workspace layout are
unchanged. The original implicit scalar typing is intentionally retained in
this first mechanical migration; making all local declarations explicit is a
separate, regression-gated cleanup rather than a numerical change hidden in
the structural conversion.

`src/fixed_orientation_amplitude.f90` and
`src/scattering_matrix_expansion.f90` are migrated the same way. They keep
Mishchenko's routine names `AMPL`, `VIGAMPL`, and the recurrences of `MATR`
and `HOVENR` so that they stay checkable against the original sources, while
returning status instead of writing or stopping.

`src/fixed_orientation_phase_matrix.f90` carries the `Z11 … Z44` bilinears of
the same `ampld.lp.f`. They are inline statements of its main program rather
than a routine of its own, so there is no original name to keep; the
expressions and their operation order are unchanged, and the declarations are
explicit.

The Makefile compiles this core with `-frecursive -fno-inline`.

`-fno-inline` used to be a crash workaround, and is not one any more. The
scratch arrays that made an inlined `-O2` build overflow the OpenMP worker
stack now belong to the workspace: the deepest full-direct call chain needs
about 190 KB of stack inlined and 175 KB not inlined, measured with
`-fstack-usage`, against about 3.1 MB and 1.7 MB before the move. The default
worker stack on the reference platform sat between those two figures, which is
what the flag was really hiding. An inlined build now passes the 1/2/4/8-thread
suites down to `OMP_STACKSIZE=256K`.

The flag is retained for bit reproducibility. Inlining changes where the
compiler contracts multiply-add pairs, which moves the last bits of the
convergence test in `TMD_ONE_FULL_DIRECT`; 2 of the 6258 rows of a shipped
astrodust table slice then shift by about 1e-6 in `Q_abs`, `Q_sca`, albedo and
`g`, with no flag changes. That is not a physics difference, but it does mean
an inlined build no longer reproduces the shipped table byte for byte.
`-ffp-contract=off` does not restore it. Dropping `-fno-inline` buys about 3%
of wall time; keeping it is the deliberate trade.

After its first successful default evaluation, one workspace accounts for
83,926,096 bytes (80.04 MiB) of explicit storage:

- 14,400 bytes for fixed `FAC` and `SS` tables;
- 38,938,880 bytes for CT/CTT, Bessel, and GSP `B1/B2` components;
- 16,588,800 bytes for eight T-matrix single-precision fields;
- 1,280,000 bytes for 16 TMATR double-precision work fields;
- 24,729,600 bytes for six GSP D-array fields; and
- 2,374,416 bytes for the `TT`, `TMATR0`/`TMATR`, and `GSP` scratch arrays
  (640,000 + 960,000 + 774,416) that were automatic arrays until they moved
  here.

Every thread needs its own workspace, so that figure multiplies by the thread
count. `tmatrix_eval_scattering_matrix` and `ampl` add nothing to it: the
expansion coefficients travel in the caller's `tmatrix_scatmat_t`, and the
amplitude evaluation keeps its own scratch in automatic local storage.

`tmatrix_eval` is the only evaluation entry point. The fixed-form migration
backends and their COMMON-state bridges have been removed now that the
free-form core reproduces the 120-case golden reference and the shipped
astrodust table byte for byte; `tests/golden/legacy_cases.dat` remains the
recorded reference those backends produced.

`libtmatrix.a` contains only the modern public types/API, full-direct core and
bridge, and the LAPACK/BLAS routines required by that core. It requires no
OpenMP runtime. `make test-production-package` checks these archive boundaries
and also builds the external-consumer example without `-fopenmp`.

Those LAPACK/BLAS routines keep their standard names (`zgetrf_`, `zgemm_`,
`xerbla_`, ...) and are visible in the archive. A host that links its own
LAPACK alongside `libtmatrix.a` should control link order: this `xerbla_`
returns to its caller instead of printing and stopping, so binding a host's
LAPACK to it would silence that host's argument-error diagnostics.

