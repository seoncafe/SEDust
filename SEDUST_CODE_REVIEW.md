# SEDust Code Review

Date: 2026-07-15

Reviewed revision: `8e8a39d` (`main`)

Scope: `sed/`, `mc/`, `tmatrix/`, `pyutil/`, build files, and the public documentation

## Verification and disposition (2026-07-15)

Every finding below was independently re-verified against revision `8e8a39d`
by reading the cited code. The review is factually accurate almost throughout:
the cited line numbers were correct, and 16 of the 18 findings were confirmed
as described. Two corrections: the failure H5 predicts is unreachable (see the
H5 row), and H4's severity is overstated because nothing in SEDust imports
`pyutil`, so no validation path is affected. M4(d) is partially correct: on a
step-cap hit the fixed engine keeps the full elapsed time, so the
"photon applied early" mechanism is specific to the adaptive engines.

Fixes were applied in two passes on 2026-07-15: commit `4b82459`
(H2, H3, M2, M3, M6, P3, and the two solve-time stops of P4), followed by a
cleanup pass (M1, H4, H5, the M4 comment, M5 validation, M7 guards).

| Finding | Verdict | Disposition |
|---|---|---|
| H1 global model state | Confirmed. Real for multi-model or externally threaded RT use; harmless for the shipped CLI, which builds one model once. | Deferred (design refactor). The one-active-model constraint is stated in the code and in this document. |
| H2 qm radius grid | Confirmed: a latent crash/garbage path for `'qm'` with Zubko or file-defined models, with no guard. | Fixed (`4b82459`): `grain_pop_t%aeff` added, every builder fills it, the qm branch reads the population radii, and `dust_emission` guards (status=2). Verified end to end with `build_zubko` + `'qm'` (three populations with different size grids, status=0). |
| H3 T-matrix convergence flag | Confirmed; inherited from the original `tmd.lp.f`. Never triggered in the shipped table: 47,574 calls, flags 0/10/20 only. | Fixed (`4b82459`): `IERR = 5` on refinement-loop exhaustion; header and driver legends updated. Test-mode flag distribution unchanged. |
| H4 Python Mathis field | Numbers confirmed (1.0e-13 / 2.9 K vs the Fortran default 1.65e-13 / 2.725 K); severity overstated, since `pyutil` has no callers inside SEDust. | Fixed (cleanup pass): `corrected=True` default matching the Fortran `use_mathis_corrected`; `corrected=False` reproduces the CLI `mathis_orig`. Unit docstrings in `sed_from_cabs.py` made consistent (SI). |
| H5 uninitialized `Pmax` | Declaration confirmed, but the claimed failure is unreachable: the first-iteration `TMAX <= TMIN` exit cannot occur given how `UMIN`/`UMAX` are constructed, so `Pmax`/`P_out` are always assigned before any read. | Defensive `Pmax`/`P_out` initialization added (cleanup pass). |
| M1 Planck overflow | Confirmed. Under default IEEE execution the overflow yields the correct 0 limit, so shipped results were right; the routine failed only under floating-point traps. | Fixed (cleanup pass): three-branch `bbody` (overflow guard, `sinh` identity for small `x`); the production range is bit-identical by construction. |
| M2 MC straddling bin | Confirmed exactly as described: the interval crossing `lam_c` was dropped from both integrals. | Fixed (`4b82459`): the crossing interval is split at `lam_c` for the heating integral, the event integral, and the event CDF (normalization matches by construction; sampled wavelengths stay at or below `lam_c`). Measured effect, same seed and single thread: +0.036% integrated MC SED, at most 0.3% at any wavelength. |
| M3 zero event rate | Confirmed as latent; unreachable with the shipped inputs (the Q table's blue edge is the Lyman limit, and grains whose cutoff clamps that low take the equilibrium gate). | Fixed (`4b82459`): setup guard with a clear stop and message. |
| M4 step control | (a)-(c) confirmed; (d) partially correct (see above). The energy effect is removed by the absorbed-power rescale in `mc_sed`, and the step caps are not reached in practice (cooling converges in tens of steps against caps of 300/4000). | The misleading "same logic" comment corrected and the dead `DLNT_MAX` in `cool_segment` removed (cleanup pass). Unifying the two tolerances is deferred. |
| M5 table API | Confirmed, both halves. | Validation added and the exactness comments corrected (cleanup pass): optional `status` on `dust_build_table` and `dust_emission_interp`, same pattern as `dust_emission`. The repeated log/bracket work (the performance half) is deferred. |
| M6 unconditional diagnostics | Confirmed, including one line to stderr for every grain in the qm path. | Fixed (`4b82459`): `sed_verbose`/`qm_verbose` guards driven by `dust_model_t%verbose`; the library path is silent by default, the CLI drivers are unchanged. |
| M7 loader trust | Confirmed, all four hazards. | Guards added (cleanup pass): `MAXP` checked before writing, channel >= 1, wavelength-grid consistency across populations, at least 2 radii. Negative test (channel 0) stops cleanly; the shipped descriptor does not false-trip. |
| M8 table consistency checks | Confirmed as missing; the shipped table itself is clean (all 190,801 rows scanned: no NaN/Inf, no negative cross sections, albedo and asymmetry in range). | Deferred (writer/loader validation). |
| P1 solver workspaces | Confirmed: two `NT x NT` matrices (about 640 KB at `NT=200`) allocated and zeroed for every stochastic solve, inside the OpenMP region for the draine/qm paths. | Deferred. |
| P2 Planck precompute | Confirmed: `bbody` is re-evaluated inside the size loop of `build_kappB`. | Deferred (startup cost only). |
| P3 build flags / archive | Confirmed, all three points. | Fixed (`4b82459`): portable default flags with `-march=native` kept as a commented option, and the archive is removed before `ar` so a dropped source cannot leave a stale member. `-w` is retained for now. |
| P4 `stop` in library code | Confirmed: 29 stops in the library sources, of which 2 are reachable at solve time from `dust_emission`. | The two solve-time stops fixed (`4b82459`): `dust_emission` validates up front with an optional `status`, and the qm sparse overflow now reports the grain unsolved (existing GD fallback) instead of stopping. The reader/builder stops, which fire only at model-construction time on bad input files, are deferred. |
| P5 numeric kinds | Not separately verified. | Deferred. |

Regression evidence for the fixes: the astrodust CLI outputs (default and qm,
6 files) are byte-identical to the pre-fix reference, and remain so after the
portable-flag rebuild and the cleanup pass; `main_dl07.x`, the `rt_example`
link test, and `run_tmatrix.x test` all pass; the shipped Q table is untouched.

## Executive summary

SEDust has a clear physical pipeline, useful independent solver paths, and unusually detailed in-source numerical notes. The default astrodust heuristic run also completes successfully with bounds checking enabled. However, the public library abstraction is not yet as self-contained or model-agnostic as its API suggests. Several routines still depend on mutable module-global state, and this creates correctness and thread-safety risks for an RT embedding.

The highest-priority issues are:

1. `dust_model_t` is not a self-contained model: emission uses global wavelength, temperature, radius, and solver state belonging to the most recently built model.
2. The QM solver cannot safely process Zubko or generic file-defined populations because their per-population radius grids are not stored, while the QM path indexes the global `aeff` array.
3. The T-matrix wrapper can exhaust its Gaussian-quadrature refinement loop and still return `IERR=0`, marking an unconverged result as successful.
4. The Python Mathis field has drifted from the production Fortran defaults, so Python and Fortran validation paths no longer construct the same radiation field.
5. A rare early-exit path in iterative temperature-window narrowing reads an uninitialized `Pmax` and may accept an undefined probability distribution.

Before using the exact solver concurrently over a large 3D grid, I recommend making all grid/state dependencies explicit, adding input and shape validation at the public API boundary, silencing per-grain diagnostics by default, and adding small regression tests that compare energy balance and stored baselines.

## Review method

The review included:

- source inspection of the solver, Monte Carlo, T-matrix, I/O, Python, and public API paths;
- strict GNU Fortran syntax compilation with `-std=f2018 -Wall -Wextra -Wimplicit-interface -Wsurprising` for modern Fortran code, and a legacy-mode warning build for the bundled T-matrix code;
- a debug build with `-fcheck=all -fbacktrace`;
- a default astrodust heuristic run with `OMP_NUM_THREADS=1`, which completed all S1, S2, and PAH stages;
- a floating-point-trap build with `-ffpe-trap=invalid,zero,overflow`, which reproduced an overflow in the Planck function during initialization;
- Python bytecode compilation of `pyutil/*.py`.

No automated unit or regression test target was found. The repository contains comparison drivers, plots, reports, and baseline data, but they are not wired into a pass/fail test suite.

## Findings

### High priority

#### H1. `dust_model_t` still depends on the most recently built global model

**Evidence**

- `sed/src/sed_astrodust.f90:1483-1490` explicitly states that module-global grids are the active model and that only one model may be active.
- `dust_emission` accepts a model object but immediately writes the global `stoch_method` at line 1911.
- `sed_grain_loop` uses global `T_first`, `lam`, `aeff`, `NLAM`, and `NT` rather than fields passed from `m` (for example lines 681-695 and 840-870).
- `calc_Teq` relies on integration weights cached globally by `p_sub_setup`; a later builder replaces those weights.
- `apply_induced_factor` uses global `NLAM` and `lam` at lines 1476-1478 rather than the model wavelength array.

**Impact**

- Building model B invalidates subsequent emission calls on a still-live model A.
- A shape mismatch can produce a bounds error; equal-sized but different grids can silently produce physically wrong output.
- Concurrent calls are not safely reentrant. Even calls using the same model write shared state (`stoch_method`), which is an OpenMP data race under the language memory model.
- The API and README describe a model object, but correctness depends on hidden call order.

**Recommendation**

Pass a model or solver-context object all the way into `sed_grain_loop`, `calc_Teq`, `calc_P`, narrowing helpers, and induced-emission handling. Store integration weights in that context. Remove writes to global solver configuration during emission. If a full refactor is deferred, add a model generation/token check and a documented serial-only guard so stale or concurrent calls fail clearly instead of returning wrong results.

#### H2. QM emission uses a radius grid that Zubko and file-defined models do not own

**Evidence**

- `grain_pop_t` has `dn`, cross sections, and enthalpy, but no radius array (`sed/src/dust_model_mod.f90:29-37`).
- `build_zubko` and `build_from_files` explicitly allocate `m%aeff(0)` because grids are said to be held per population (`sed/src/sed_astrodust.f90:1662-1665` and 1833-1835), but no per-population radius field exists.
- In the QM branch, `sed_grain_loop` unconditionally evaluates `a_cm_qm = aeff(ir) * UM2CM` (`sed/src/sed_astrodust.f90:864-870`). This is the module-global `aeff`, not population data.
- The same global `aeff` is listed as shared in the OpenMP region at line 843.

**Impact**

For `m%stoch_method='qm'`, Zubko and generic file-defined models can index an unallocated, zero-length, or stale radius array. With checking enabled this should fail; without checking it may pass a garbage radius to the enthalpy-mode builder and return invalid emission.

**Recommendation**

Add `aeff(:)` to `grain_pop_t`, populate it in every builder, and pass `aeff_pop` into `sed_grain_loop`. Do not use a model-level radius grid for models whose populations have different grids. Add a regression test that executes `build_zubko` and `build_from_files` with all four solver choices.

#### H3. T-matrix quadrature non-convergence can be reported as success

**Evidence**

In `tmatrix/src/tmd_one.f:170-199`, the wrapper increases `NGAUSS` through `NPNG1`. It jumps to label 155 when `DSCA` and `DEXT` converge, but ordinary exhaustion of the `DO 150` loop falls through to the same label. No flag records whether the convergence condition was reached, and `IERR` remains zero.

The driver treats `IERR=0` as a converged T-matrix point and writes `flag=0` (`tmatrix/driver/run_tmatrix.f90:165-178`).

**Impact**

An unconverged optical property can be indistinguishable from a valid one in the shipped Q table. It will not trigger either fallback, and downstream SED calculations have no way to identify it.

**Recommendation**

Track a logical `CONVERGED_GAUSS`, set it only on the convergence branch, and return a new nonzero `IERR` when the loop exhausts. Extend the driver flag legend and report counts of every fallback/error code. Also reject non-finite or non-physical results before writing a nominally successful row.

#### H4. Python and Fortran construct different default Mathis radiation fields

**Evidence**

- Production Fortran defaults to the corrected 4000 K dilution factor `1.65e-13` and a 2.725 K CMB (`sed/src/radfield.f90:13-18`, 46-49).
- Python uses `1.0e-13` and 2.9 K (`pyutil/radiation_fields.py:16`, 33-38), and its module docstring still claims it mirrors the Fortran implementation.
- `pyutil/sed_from_cabs.py:17-20` describes `J_lam` as CGS, while its function documentation and calculations use the same SI Planck convention as Fortran.

**Impact**

Python-versus-Fortran comparisons include a real input-physics difference, especially in optical heating and long-wavelength CMB behavior. This can be mistaken for a solver discrepancy.

**Recommendation**

Expose a shared, explicit choice such as `corrected=True` in Python, make the corrected values the default, and add a small cross-language fixture containing wavelengths and expected `J_Mathis` values. Correct the unit documentation so the public convention is unambiguous.

#### H5. Iterative narrowing can read uninitialized `Pmax`

**Evidence**

`Pmax` is declared but not initialized in `narrow_iterative` (`sed/src/sed_astrodust.f90:1182-1187`). The iteration may exit at line 1207 when `TMAX <= TMIN`, before `Pmax` is assigned at line 1214. After the loop, lines 1278-1283 test `Pmax` and may set `converged=.true.`.

**Impact**

For a collapsed/clamped temperature window, behavior is undefined. The caller may treat an uninitialized `P_out` as a valid stochastic distribution. The path is likely rare for the shipped default model, but it is reachable with extreme fields, narrow temperature grids, or custom enthalpy tables.

**Recommendation**

Initialize `Pmax=0`, initialize `P_out=0`, and track whether `calc_P` has completed at least once. Only accept a best-effort result when a valid distribution was actually computed and normalized.

### Medium priority

#### M1. The Planck implementation overflows instead of taking the analytic limit

**Evidence**

`bbody` directly evaluates `exp(hc_kB/(T*lambda_m))` (`sed/src/radfield.f90:99-104`). A debug run with overflow trapping failed during `build_kappB` at this line. The Python version already uses `expm1`, showing the intended stable approach.

**Impact**

Normal IEEE execution often converts the overflowed denominator to infinity and returns zero, so the default run succeeds. Nevertheless, the routine fails under floating-point traps, generates avoidable exceptions, and loses precision for very small exponent arguments through `exp(x)-1` cancellation.

**Recommendation**

Validate `T>0` and `lambda>0`; return zero above a safe exponent threshold; use `expm1(x)` in the remaining range; and use a Rayleigh-Jeans expansion for very small `x` if the wavelength range is extended further.

#### M2. Monte Carlo cutoff integration drops the interval crossing `lam_c`

**Evidence**

In `grain_setup_from_cabs`, continuous heating includes an interval only when its left endpoint is at or above the cutoff, while the event rate includes it only when its right endpoint is at or below the cutoff (`mc/src/mc_engine.f90:182-190`). The one interval that straddles `lam_c` is included in neither term. The event CDF uses a similar endpoint rule at lines 199-206.

**Impact**

Absorbed power and photon rate are systematically underestimated by the omitted bin. The adaptive cutoff usually lies between tabulated wavelengths, so this is the common case rather than an exceptional one. The size of the bias depends on local grid spacing and spectral structure.

**Recommendation**

Interpolate the integrand at `lam_c` and split the crossing interval between continuous and stochastic contributions. Build the event CDF from the same split integral so `rate_event` and the CDF normalization are identical by construction.

#### M3. A zero stochastic-event rate is converted into a later division by zero

**Evidence**

When the event integral is zero, `Fmax` is replaced by 1 only for CDF construction (`mc/src/mc_engine.f90:208-214`); `g%rate_event` remains zero. Every engine later calls `rng_exp(rng, g%rate_event)`, which divides by the rate (`mc/src/mc_rng.f90:87-93`).

**Impact**

Custom cross sections with no absorption below the cutoff, a degenerate wavelength grid, or an extreme cutoff can cause infinity or a floating-point exception. The fabricated CDF does not make the physical event process valid.

**Recommendation**

Treat `rate_event<=0` as a continuous-heating/equilibrium-only case, or return a clear setup status. Validate positive, finite wavelength and cross-section inputs before constructing the CDF.

#### M4. Adaptive MC engines use a much coarser cooling step and silently truncate long segments

**Evidence**

- The fixed engine limits a step to a 5% fractional temperature change (`DT_FRAC_MAX=0.05`, `mc/src/mc_engine.f90:454-500`).
- The shared adaptive `step_advance` uses `DLNT_MAX=0.5`, allowing a 50% change (`mc/src/mc_engine.f90:625-665`). Its comment says it reuses the same trajectory math, but it does not use the tightened limit.
- Adaptive engines stop an interval after `MAX_STEPS=300` (`mc/src/mc_engine.f90:795-797`, 844-856, 899-910, and 968-970) without reporting that `tau_seg < dt_act`. They then continue to the photon event using the shortened elapsed time.
- The fixed engine has the same silent-cap pattern with a larger cap of 4000 (`mc/src/mc_engine.f90:461-472`).

**Impact**

The adaptive and fixed solvers integrate different trajectories. If a cap is hit, the next photon is applied too early and histogram/emission time is lost, biasing both `P(T)` and energy accounting.

**Recommendation**

Use one step-control implementation and tolerance for all engines. Return a status when the step cap is reached, dynamically continue with a safe fallback, and count/report cap hits. Add convergence tests over step tolerance and `MAX_STEPS`.

#### M5. The table interpolation API lacks validation and repeats identical work per wavelength

**Correctness evidence**

`dust_build_table` does not validate the lengths of `J_ref`, positivity/monotonicity of `U_grid`, or a minimum grid size (`sed/src/dust_lib.f90:84-108`). `dust_emission_interp` immediately computes `log(tab%U)` and `log(U)` (`sed/src/dust_lib.f90:110-136`). Zero/negative values therefore produce infinities or NaNs. Output array shapes and channel indices are also unchecked.

Replacing exact zeros with `1e-300` means a table lookup is not strictly exact at a zero-valued grid point, despite the API comment.

**Performance evidence**

Every cell call allocates/deallocates `lU`, recomputes `log(tab%U)`, and performs the same binary bracket search separately for every wavelength and every channel. For an RT loop, this is avoidable work in the hottest public API path.

**Recommendation**

Validate all inputs when building the table. Store `logU` and preferably logged spectra in `dust_emis_table_t`. Locate the two bracketing U indices once per cell, then interpolate complete wavelength/channel slices with array operations. Handle exact endpoints and true zeros explicitly.

#### M6. Library calls emit unconditional diagnostics from inner solver paths

**Evidence**

`sed_grain_loop` writes a diagnostic line after every population solve (`sed/src/sed_astrodust.f90:734-736`, 832-833, 909-911, and 928). The QM grain routine writes one line per grain to unit 0 (`sed/src/stoch_qm.f90:2280-2308`).

**Impact**

An exact per-cell RT loop can produce multiple lines per population per cell; QM can produce hundreds more. I/O can dominate runtime, interleave across threads, and make the library unsuitable for large production grids.

**Recommendation**

Add a verbosity/log callback in the solver context, default it to silent for library calls, and aggregate counters so a driver can request one summary after a batch.

#### M7. Generic loaders trust dimensions and descriptor counts too far

**Evidence**

- `build_from_files` increments `npop` and writes fixed arrays before checking the `MAXP=16` limit (`sed/src/sed_astrodust.f90:1774-1797`).
- Channel values are not checked to be positive before they are used as `Jchan(:,ic)` indices (`sed/src/sed_astrodust.f90:1806-1815`, 1922-1923).
- Zubko and generic builders assume all population wavelength grids match the first one, but do not compare `nwave` or wavelength values before looping over global `NLAM` (`sed/src/sed_astrodust.f90:1643-1674`, 1817-1844).
- Radius grids are assumed to contain at least two points when endpoint `dln(a)` values access element 2 (`sed/src/sed_astrodust.f90:1713-1718`, 1864-1869).

**Impact**

Malformed or simply heterogeneous custom inputs can cause out-of-bounds access or silent grid mixing. This weakens the advertised data-driven model path.

**Recommendation**

Validate descriptor counts before incrementing, require channels in `1:n_channel`, require at least two strictly increasing radii/wavelengths/temperatures, verify identical wavelength grids across populations (or resample explicitly), and return structured errors instead of calling `stop` from library code.

#### M8. Optical tables are not checked for physical consistency

**Evidence**

`run_tmatrix` computes `qabs=qext-qsca` and writes values without checking finiteness, non-negativity, `0<=albedo<=1`, or `-1<=g<=1` (`tmatrix/driver/run_tmatrix.f90:154-182`). `load_q_table` trusts row ordering and overwrites axis values without verifying that every row belongs to the expected grid (`sed/src/q_table.f90:74-113`).

**Impact**

A corrupt, unconverged, shuffled, or partially merged table can flow directly into negative cross sections and logarithmic interpolation. Row-count checks alone do not establish table integrity.

**Recommendation**

Validate every generated row, record failures explicitly, and make the loader verify axis consistency, monotonicity, finite values, and `Qext ~= Qabs + Qsca` within tolerance. Reject extra rows as well as missing rows.

### Performance and maintainability opportunities

#### P1. Reuse stochastic-solver workspaces

`calc_P` allocates and zeros two `NT x NT` matrices for every grain solve (`sed/src/p_sub.f90:104-117`, 186). `dust_emission` and `sed_grain_loop` also allocate temporary spectra and probability arrays on every call/population/thread. In a many-cell RT run, allocator traffic is significant even before matrix work begins.

Introduce a per-thread `solver_workspace_t` containing matrices, spectra, probability arrays, and narrowing buffers. Allocate it once per thread or solver instance and reuse it. Longer term, exploit the triangular transition structure to reduce storage and unnecessary zeroing.

#### P2. Precompute Planck values during model construction

`build_kappB` evaluates `bbody(T,w)` inside the size loop even though it depends only on temperature and wavelength (`sed/src/sed_astrodust.f90:1424-1435`). Precompute `B(NW_INT,NT)` once, then integrate each cross-section column against it. This mainly improves startup and repeated model construction.

#### P3. Make the build warning-clean and reproducible

The main Makefiles use `-w`, hiding all warnings, and the SED build uses `-march=native`, which can make a static library unusable on older compute nodes. Provide separate release and debug profiles, for example:

- release: optimized, portable architecture baseline, warnings enabled;
- native: opt-in host-specific optimization;
- debug/check: `-O0 -g -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow`;
- CI: strict standard/warning compilation.

The `libsedust.a` recipe also updates an existing archive with `ar rcs`; if a source/object is removed from `LIB_SRC`, its old member can remain in the archive. Recreate the archive from a clean object list or delete the archive immediately before archiving.

#### P4. Replace `stop` in library modules with status returns

Many readers and solver helpers terminate the entire process. This is reasonable for standalone drivers but not for an embedded RT library. Return an error code/message (or use a caller-provided error handler) at public/library boundaries, and let standalone programs decide whether to terminate.

#### P5. Standardize numeric kinds and constants

Strict compilation reports many default-real literals being converted to `real(wp)`, obsolescent statement functions in `enthalpy_v2.f90`, and real equality comparisons. Most are not immediate defects, but consistent `_wp` literals and modern internal/pure functions would reduce accidental single-precision evaluation and make future kind changes safer.

The bundled T-matrix code intentionally mixes `REAL*8` calculations with large `REAL*4` storage arrays. Given the requested convergence tolerance, this should be documented and covered by a precision-sensitivity comparison before tightening `DDELT`.

## Testing recommendations

The following tests would cover the most consequential risks with modest runtime:

1. **Public API shape/error tests**: wrong `J_lam` length, invalid U grid, invalid channel, one-point axes, non-monotonic grids, and zero event rate.
2. **Model-lifetime test**: build A, build B, call A again. The final design should either reproduce A exactly or reject the call explicitly.
3. **Solver/model matrix**: astrodust, DL07, Zubko, and file-defined models crossed with `heuristic`, `draine`, `qm`, and `equil`.
4. **Thread test**: compare serial and concurrent exact calls under ThreadSanitizer-capable tooling or an OpenMP stress test.
5. **Energy balance**: absorbed versus emitted bolometric power for representative small, transition, and equilibrium grains under multiple U values.
6. **MC convergence**: vary event count, step tolerance, histogram size, engine type, and thread count; require statistical agreement with confidence intervals.
7. **T-matrix status test**: deliberately force insufficient `NPNG1`/quadrature capacity and require a nonzero error flag.
8. **Cross-language radiation test**: compare Fortran and Python `J_Mathis` and Planck values at fixed wavelengths.
9. **Baseline spectra**: convert selected existing comparisons into numeric tolerances rather than plot-only checks.
10. **Debug CI run**: execute at least a small astrodust solve with bounds, uninitialized-value, and floating-point checks.

## Suggested implementation order

1. Fix H2, H3, and H5 because they can directly produce invalid results or undefined behavior.
2. Refactor the model/context boundary in H1 before scaling the library to multi-model or concurrent RT use.
3. Synchronize Python radiation fields and stabilize the Planck function.
4. Repair MC cutoff/zero-rate handling and unify cooling-step control.
5. Add validation and regression tests around every public builder and emission entry point.
6. Optimize table interpolation and workspace reuse after correctness tests are in place.

## Positive observations

- The repository includes independent GD, QM, and Monte Carlo approaches, which is a strong basis for cross-validation.
- The source comments document many numerical decisions and previously observed discrepancies in enough detail to make future regression tests practical.
- Most modern Fortran code uses explicit interfaces through modules and compiled successfully under strict syntax checks.
- The default astrodust heuristic path completed under `-fcheck=all`, providing a useful starting baseline.
- Data dependencies are packaged inside SEDust, and the standalone drivers use clear, reproducible relative paths.
