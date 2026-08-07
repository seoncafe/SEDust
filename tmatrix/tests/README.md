# T-matrix regression tests

`golden/legacy_cases.dat` is the direct-call reference for the 120 fixed
single-particle cases in `test_cases.f90`: spheroid, cylinder, and Chebyshev
shapes across five radii, four wavelengths, and weak/strong absorption. It was
recorded from the fixed-form Mishchenko implementation this library was
translated from, built with `gfortran` and the default `tmatrix/Makefile` flags
on the reference platform. That implementation is no longer in the tree, so
this file is the only surviving statement of what it computed. Do not
regenerate it to make a changed implementation pass: explain and approve any
intentional numerical change.
Reference SHA-256 values are `1a127149875971f2232f35f46072b32491d7a75efdb889be75d2b70da11935b1`
for `legacy_cases.dat` and `6ce24b9519d5779ff7972c923c66d79503036d15757a9f4ecff26830da6bfe8c`
for `standalone_test.dat`.

`make test` runs the public API comparison against that reference, the
full-direct workspace comparison, and the API invalid-input/status tests. The
standalone `run_tmatrix.x test` command is kept as a separate Q-table smoke
test because it exercises the size-parameter regime selection in
`driver/spheroid_optics.f90`, which sits outside the raw library API.

`test_parallel_api.x` and `test_parallel_full_direct.x` exercise the same
default `tmatrix_eval` entry point with one workspace for each OpenMP caller.
`make test-full-direct-parallel` executes the independent-workspace check with
1, 2, 4, and 8 threads; `make test-api-parallel` runs the corresponding public
API harness at the same thread counts. Only those parallel-test executables
compile and link with `-fopenmp`; the production library and external example
do not.

`test_parallel_scattering_matrix.x` covers the same 120 cases through
`tmatrix_eval_scattering_matrix` followed by four `ampl` calls on each case,
each passing back the `scatmat` that call returned, and
`make test-scattering-matrix-parallel` runs it at 1, 2, 4, and 8 threads.
It compares IEEE bit patterns, not a tolerance: both runs execute the same
code on the same inputs, so any difference would be state leaking between
workspaces rather than rounding. The `ampl` status of every call is compared
too, so a workspace that lost or replaced its T-matrix would show up as a
status mismatch rather than as amplitudes.

`test_parallel_spheroid_optics.x` is the only harness that enters the driver
layer. It runs `spheroid_q` on 48 cases spanning all three size-parameter
regimes — 16 Rayleigh, 24 T-matrix, 8 geometric optics, asserted from the
returned flags so a case list that stopped covering a regime is reported rather
than passed — and calls `rayleigh_limit` and `geometric_optics_limit` directly
with their orientation-resolved outputs. Those optional outputs are the point:
the cross-section-only path through `spheroid_q` never enters the surface
quadrature of `fresnel_opaque_absorption`, so a harness that stopped at
`spheroid_q` passes even when that routine's locals are shared between threads.
`make test-spheroid-optics-parallel` runs it at 1, 2, 4, and 8 threads, again on
IEEE bit patterns, and compares the regime flag and the library status as well
as the efficiencies.

`golden/standalone_test.dat` is the 49-row output of `run_tmatrix.x test`,
including the historical table layout and flags. `make smoke` writes its
candidate only to `/tmp`, compares it byte-for-byte to this table, and removes
the temporary file. Regenerate this reference only after an explicitly reviewed
numerical/output-format change.

It was last regenerated when the astrodust wavelength axis was carried from
0.0912 um (13.6 eV) down to 1.0e-4 um (12398 eV), taking the table from 1129 to
1762 wavelengths (`data/astrodust/DH21_wave_to_12keV`). The subset takes seven
samples along each axis, so a longer axis picks different wavelengths — the
stride moved from 188 to 293 and none of the sampled wavelengths survived. That
is a change of which points are sampled, not of what the code computes at them:
the full sweep on the extended axis reproduces the full sweep on the old one
byte for byte over all 190801 rows at lambda >= 0.0912 um, checked separately by
building the driver against the old axis and running it. The new subset spans
1.0e-4 to 3.85e4 um and covers all three regimes more evenly than the old one
(14 T-matrix / 23 Rayleigh / 12 geometric-optics rows, against 13 / 34 / 2).
