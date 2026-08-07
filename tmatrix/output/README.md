# tmatrix/output -- what the T-matrix drivers write

The products here are what `../driver/run_*.f90` writes: the output of a sweep,
tracked so that regenerating one is optional.

**This is not where the solver reads them from.** Everything a dust model owns
lives in that model's own directory, so the installed copies the SED solver and
the polarized radiative-transfer API read are under `../../data/astrodust/`.
A sweep produces a file here; installing it means copying it there. The two can
therefore differ, and do: `q_astrodust_jori_P0.20_Fe0.00_1.400.dat.gz` here is
SEDust's own regenerated orientation-resolved table (opt-in through
`qpol_path`, and the only one with a birefringence block), while the default
`../../data/astrodust/q_DH21Ad_P0.20_Fe0.00_1.400.dat.gz` is the HD23 release
table.

## The two wavelength axes do not match, and that is deliberate

| product | wavelength axis | length |
|---|---|---|
| `q_astrodust_P0.20_Fe0.00_1.400.dat` | `../../data/astrodust/DH21_wave` | 1129 |
| `q_astrodust_P0.20_Fe0.00_1.400_euv.dat` | `../../data/astrodust/DH21_wave_to_12keV` | 1762 |
| `q_astrodust_jori_P0.20_Fe0.00_1.400.dat.gz` | `../../data/astrodust/DH21_wave` | 1129 |

The two SCALAR products are deliberately separate.  The plain file preserves
the 1129-wavelength non-ionizing grid.  The `_euv` file prepends the Draine &
Hensley (2021) dielectric function's own nodes below 0.0912 um, reaching
1.0e-4 um (12398 eV).  Its final 1129 wavelength blocks are value for value
the plain file; the extension changed no existing cell.

The ORIENTATION-RESOLVED (polarized) table was NOT extended with it. Polarized
transfer is run at a handful of wavelengths, not on a whole axis, so the
wavelengths a given calculation needs are generated on demand
(`run_q_jori.x lam L1 L2 ...`, or its `euv` mode) and supplied through
`sed_init`'s `qpol_euv_path` / `qpol_euv_wave_path`. Sweeping the 633 new
wavelengths for all 169 radii would cost hours of solver time for entries no
planned calculation reads.

What this means at run time: when an EUV run selects the 1762-point scalar
table, the model's polarized block covers 0.0912-39810 um and the 633
wavelengths below 0.0912 um have no
polarized entry. `build_Cpol` locates that block from the polarized table's own
coverage, reports the uncovered band on `error_unit` as "zero by omission, not
by physics", and leaves `C_pol_ext` and `C_bir_ext` at zero there. It is a
stated deficit, not a physical value: a b/a = 1.4 spheroid has non-zero
dichroic extinction throughout that band. Supply a companion table for the
wavelengths a calculation needs, and that band is filled from it.

## Scattering-matrix products are independent of the axis

`scatmat_astrodust_*.dat` and `scatmat_aligned_astrodust_*.dat(.gz)` are
computed at the five UBVRI wavelengths 0.36, 0.44, 0.55, 0.64 and 0.79 um
(`BANDS` in `../driver/run_scatmat_aligned.f90`; the command line for
`run_scatmat.x`). None of them is a node of either wavelength axis -- the
refractive index at each is interpolated from the dielectric function directly
-- so extending the axis leaves these products unchanged. Their size comes from
the angular grid (theta_i x theta_s x phi times the 4x4 Mueller matrix), not
from a wavelength count.

## Smoke-test products

`*.test.dat` are the `test`-mode outputs of the drivers. They take their
sampling stride from the length of the wavelength axis their driver sweeps, so
each one follows the axis of the product it samples: `q_astrodust_*.test.dat`
was regenerated on the 1762-point axis, while `q_astrodust_jori_*.test.dat`
stays on the 1129-point one and is unchanged. `scatmat_*.test.dat` use fixed
wavelength lists and are unaffected.
