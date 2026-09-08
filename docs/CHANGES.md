# CHANGES

A dated record of code and documentation changes, newest first. Each entry names
the files it touches and, where a change moves numbers, states the measurement
that was made on the tree itself.

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

`sed/src/check_build_dust.f90`, `sed/src/sed_astrodust.f90`,
`docs/SEDust_user_manual.tex`.

v1.20 carried `test_dustem_product.x`, written when THEMIS and G18D were added
on 2026-09-05. It compared those two models' HDF5 products against their DustEM
text tables and nothing else, so v1.20 had no product-vs-text check at all on
astrodust, DL07, MRN or Zubko. v1.00 carried `check_build_dust.x`, which makes
that comparison for all six models on both wavelength grids and adds the exact
ZDA round trip. The two are now **one program**: `sed/src/check_build_dust.f90`,
byte-identical in the two trees and built as `check_build_dust.x` in both. The
narrower duplicate is deleted, and v1.20 gains the four models it was not
checking.

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
`sed/src/check_build_dust.f90`, `sed/src/test_h5_storage_policy.f90`,
`pyutil/migrate_h5_regime_int8.py` (new, replacing
`pyutil/migrate_h5_float32.py`), `sed/Makefile`, `README_HOWTO.md`,
`docs/SEDust_user_manual.tex`, every `data/<model>/sedust_<model>.h5`, and every
file under `sed/output/`.

The 32-bit conversion of the previous entry is undone. Two storage classes now,
not three:

| class | file type | datasets |
|---|---|---|
| values | 64-bit float | every computed quantity -- `Q_ext`, `Q_abs`, `Q_sca`, `g`, `albedo`, the `C_*` cross sections, `K_abs` -- and the coordinate axes `lambda`, `a_eff` they are tabulated on |
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

* every dataset of the six products: 150 of 151 are 64-bit float, the one
  `regime` is 8-bit integer; no product lost a dataset.
* the six products come to 45.0 MB on disk (astrodust 9.31, dl07 9.77,
  mrn 6.47, zubko 16.31, themis 2.16, g18d 0.94 MB).
* `sed/output/`, all 87 tracked files rerun and compared against the
  pre-conversion commit: 79 byte-identical, 6 differing only in the header
  lines that commit corrected (the two THEMIS `CM20` channels now named
  `CM20_plaw-ed` and `CM20_logn`, the G18D `aSil2001BE6pctG_0.4x` no longer
  truncated, and the `d16gra` / `sphdgra` DL07 headers naming the graphite
  actually used). The remaining two, `sed_dl07_mw31_60_qm.dat` and
  `sed_dl07_mw31_60_qm_stati_morig_ns150.dat`, differ in their numbers by
  1.7e-4 and 2.3e-4 relative -- because they were last written at commit
  78a7431 and never rerun afterwards. Built from the pre-conversion sources
  against the pre-conversion products, those two runs give exactly what this
  tree now gives, so the regenerated files are the first that match the code
  and data beside them. The runs are deterministic: three consecutive runs
  give one md5.

Tolerances back at the values they had before the conversion, the reason for
loosening them having gone:

* `check_build_dust.f90` `TOL_STORED`, the DustEM product-vs-text gate:
  1e-6 -> **1e-10**.
* `check_build_dust.f90` `check_zubko_zda`: back from a float32-ulp bound to
  the exact `da /= 0 .or. ds /= 0 .or. dg /= 0`, the stored ZDA tables having to
  equal the distributed files outright again.

At `TOL_STORED = 1e-10` the size integrals, `M_dust/N_H` and the wavelength
grid of THEMIS and G18D agree exactly, at 0, and the ZDA round trip is exact.
Four checks are reported over tolerance, and none of them is about storage:
they compare `/kext` in the product against the distributed
`kext_themis_euv.dat` and `kext_g18d_euv.dat`, which predate the `/kext`
group beside them. `/kext` recomputes bit-identically; the text curves do not,
and recomputing both from one run brings them to 4.5e-13 of each other against
the 1.9e-8 (THEMIS) and 2.6e-8 (G18D) that stand now. The text curves are
what wants regenerating; they are left alone here.

`./test_h5_storage_policy.x` still guards the policy, now asking of every
dataset that a value be 64-bit float and `regime` 8-bit integer.
`pyutil/migrate_h5_regime_int8.py` replaces `pyutil/migrate_h5_float32.py`: it
stores the `regime` codes of an existing product as one byte each and requires
every other dataset to come out bit-for-bit identical.

## 2026-09-05

### HDF5 products store computed quantities as 32-bit; axes stay 64-bit

`sed/src/sedust_h5.f90`, `sed/src/calc_qtable.f90`, `sed/src/calc_kext.f90`,
`sed/src/check_build_dust.f90`,
`sed/src/test_h5_storage_policy.f90` (new),
`pyutil/migrate_h5_float32.py` (new), `sed/Makefile`, `README_HOWTO.md`,
`docs/SEDust_user_manual.tex`, and every `data/<model>/sedust_<model>.h5`.

Three storage classes, decided by what a dataset is:

| class | file type | datasets |
|---|---|---|
| computed quantities | 32-bit float | `Q_ext`, `Q_abs`, `Q_sca`, `g`, `albedo`, the `C_*` cross sections, `K_abs` |
| coordinate axes | 64-bit float | `lambda`, `a_eff` |
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

* the six products fall from 45.0 MB to 18.7 MB on disk
  (astrodust 9.31 -> 3.76, dl07 9.77 -> 4.06, mrn 6.47 -> 2.68,
  zubko 16.31 -> 6.77, themis 2.16 -> 0.99, g18d 0.94 -> 0.45 MB);
  no product lost a dataset.
* migrated values against the originals: worst relative change **5.96e-8**,
  which is one unit in the last place of a float32 significand. Two datasets
  reach below float32's normal range and lose their smallest entries to zero:
  the long-wavelength end of `kext/C_sca` (down to 6.2e-40 cm^2/H against a peak
  of 1.3e-21). Both are checked against
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
  1e-10 -> **1e-6**. Measured after the change: 3.4e-8.
* `check_build_dust.f90` `check_zubko_zda`, which required the stored ZDA tables
  to equal the distributed files **exactly**: now within one float32 ulp of each
  array's own peak. Measured: max |dQ_abs| 1.2e-7, |dQ_sca| 2.4e-7, |dg| 3.0e-8,
  against peaks of about 2, 4 and 0.5.

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

`mc/Makefile`. The Monte Carlo drivers that pull in the SED pipeline
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
