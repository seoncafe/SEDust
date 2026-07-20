# MoCafe polarized RT: development specification

**Date:** 2026 July 20

**Role of this document.** This is the development specification for MoCafe's
future polarized radiative-transfer update. The SEDust side of the interface is
built and verified; the MoCafe side described here does not exist yet. The
document tells the MoCafe developer exactly what to call, when, in what order,
with what units, and what to do with each result. It is written against the
library API in `sed/src/scatmat_aligned.f90` (module `scatmat_aligned_mod`), the
extinction/emission calls in `sed/src/sed_astrodust.f90`, and the reference
consumer `sed/rt_example/use_dustlib_scatmat.f90`.

The physics behind the quantities is in `docs/aligned_grain_polarization.pdf`;
how they are built is in `docs/sedust_polarization_implementation.pdf` §8; the
calling interface is in `docs/SEDust_user_manual.pdf` §5–§6. See §11 below for
which document answers which question.

---

## 1. Division of labor

The contract implements the Peest et al. (2017, 2023) Monte Carlo formalism as
extended to aligned spheroids, dichroism, and birefringence, and is the route to
redoing the polarized transfer of Seon (2018).

- **MoCafe supplies the alignment assumptions and the field geometry.** The
  alignment profile `f_align(a)` (default HD23, or set through the SEDust
  setters), the local magnetic-field direction `B-hat` in every cell, and a
  scalar alignment scale `eta ∈ [0, 1]` for each cell.
- **SEDust supplies the matrices.** Every material quantity in matrix form: the
  4×4 extinction matrix `K(theta_i)`, the aligned phase matrix
  `Z_al(theta_i; theta_s, phi)`, the random-orientation remainder, the scattering
  cross sections, and the polarized emission vector.
- **MoCafe transports.** Frame rotations, direction sampling, peel-off, and the
  `exp(-K tau)` transfer step. None of these is a SEDust responsibility.

The cell dependence reaches SEDust through only two scalars MoCafe already holds:
`eta`, and `theta_i = acos(k-hat . B-hat)` (one dot product).

---

## 2. The calls: signatures and units

### Isotropic optics and the emission source (on `m%lam`, cm²/H)

```fortran
call dust_extinction(m, Cext, Cabs, Csca, gbar=gbar, &
                     Cpol_ext=Cpol_ext, Cbir_ext=Cbir_ext)   ! cm^2/H
call dust_emission(m, J_lam, lamI_total, lamI_pol=lamI_pol)  ! per cell, emission
```

`Cbir_ext` is the birefringent extinction (zero when the loaded table has no 4th
block); `lamI_pol` is the **intrinsic** polarized emission (size integral and
`f_align` done; `sin^2 gamma`, turbulent depolarization, and position angle left
to MoCafe — see §4 and the user manual §6.7).

### Aligned scattering optics (µm² per H; note the unit change)

```fortran
call scatmat_band(lambda_um, iband, exact)                    ! pick the band once
call extinction_matrix_aligned(iband, theta_i, eta, kmat)     ! kmat(4,4)  [um^2/H]
call mueller_matrix_aligned(iband, theta_i, theta_s, phi, z)  ! z(4,4) [um^2 sr^-1/H], eta=1
call mueller_matrix_random(iband, big_theta, f_tot, f_ref)    ! 6-elt random matrices
call scattering_cross_sections(iband, theta_i, eta, csca_aligned, csca_unaligned)
```

| call | in | out (units) |
|---|---|---|
| `scatmat_band` | `lambda_um` [µm] | `iband`; `exact` (`.true.` if matched to 1e-3) |
| `extinction_matrix_aligned` | `iband`, `theta_i` [deg, 0–180 folded to 0–90], `eta` | `kmat(4,4)` [µm²/H] |
| `mueller_matrix_aligned` | `iband`, `theta_i`, `theta_s` [deg, 0–180], `phi` [deg, 0–360] | `z(4,4)` [µm² sr⁻¹/H] at `eta=1` |
| `mueller_matrix_random` | `iband`, `big_theta` [deg] | `f_tot(6)`, `f_ref(6)`, α₁-normalized |
| `scattering_cross_sections` | `iband`, `theta_i`, `eta` | `csca_aligned`, `csca_unaligned` [µm²/H] |

`kmat` structure: diagonal `eta*Cext_al`; `K(1,2)=K(2,1)=eta*Cpol_al` (dichroism);
`K(3,4)=eta*Cbir_al`, `K(4,3)=-eta*Cbir_al` (birefringence), in the Mishchenko
meridional basis (`Q = Iv - Ih`) of the grain frame.

**Unit caveat.** The aligned-scattering quantities are in **µm²** per H, while
`dust_extinction` returns **cm²**/H. Multiply the `dust_extinction` outputs by
`1e8` before combining them with the aligned matrices (§4).

### Exposed storage (public, protected) for hosts that inline lookups

`scm_lambda`, `scm_theta_i`, `scm_theta_s`, `scm_phi`, `scm_theta_ran`,
`scm_cext_tot`, `scm_csca_tot`, `scm_cext_ref`, `scm_csca_ref`, `scm_F_tot`,
`scm_F_ref`, `scm_Z`, and the sizes `scm_nband`, `scm_nti`, `scm_nts`,
`scm_nphi`, `scm_ntheta`, plus the flag `scm_loaded`.

---

## 3. Initialization checklist (once per run, serial)

1. Build the model with the aligned table path:
   ```fortran
   call build_astrodust(m, qtab, sizedist, NT, T_lo, T_hi, status=st, &
                        qpol_path=QPOL, scatmat_path=SCATMAT)
   ```
   `scatmat_path` failing to load is an error (`st = 3`); a missing `qpol_path`
   only disables polarized emission. Check `st == 0` and `scm_loaded`.
2. Optionally set the alignment: `dust_set_alignment(m, f_max, a_align, alpha)` or
   `dust_set_alignment_profile(m, aeff_in, falign_in)`. **Heed status code 4**: it
   means the loaded scattering table's recorded profile differs from the one just
   set (non-fatal — the emission and `Cpol_ext` channels re-integrate live and
   honor the new profile, but the loaded scattering matrices no longer match it).
   The sanctioned runtime variation of the aligned scattering optics is `eta`, not
   a change of profile; a different profile means regenerating the table (§5).
3. Compute the isotropic optics once on `m%lam` with `dust_extinction`
   (`Cext, Cabs, Csca, gbar, Cpol_ext, Cbir_ext`). Keep them; convert to µm² (×1e8)
   the copies you will combine with `K_al`.
4. Select the working bands with `scatmat_band` and cache the `iband` values, so
   no wavelength search runs in the photon loop.

This layer is **not** thread-safe. Do it before the photon loop.

---

## 4. In each cell and at each scattering event

`eta` is the cell's alignment field; `theta_i = acos(k-hat . B-hat)`.

### Extinction matrix, with the double-counting subtraction

The isotropic `Cext` from `dust_extinction` already contains the aligned
population's reference-orientation average `Cext_ref`. To make extinction
direction-dependent without counting that population twice, **subtract the
reference piece and add the direction-resolved matrix**:

```
K_total = (Cext_iso - eta*Cext_ref) * I  +  eta * K_al(theta_i)
```

with `Cext_iso` the isotropic total at the band (µm²/H) and
`Cext_ref = scm_cext_ref(iband)`. At the incidence average `<K_al> = Cext_ref`,
so the two cancel and the diagonal returns `Cext_iso` — nothing is double-counted.
Implemented literally in `use_dustlib_scatmat.f90` (add the scalar
`Cext_iso - eta*Cext_ref` to the four diagonal elements of `eta*K_al`).

### Population selection

```fortran
call scattering_cross_sections(iband, theta_i, eta, csca_aligned, csca_unaligned)
```
`csca_aligned = eta*Csca_al(theta_i)`; `csca_unaligned = Csca_tot - eta*Csca_ref`.
Their sum is the total scattering cross section, and the ratio gives the albedo
and the probability that an interaction is with the aligned versus the unaligned
population.

### Scattering

- **Aligned population.** Rotate the photon Stokes vector into the grain frame
  (`z = B-hat`) using the meridional rotation, sample `(theta_s, phi)` from the
  `Z_al,11` sky brightness, apply the full 4×4 `Z_al` (scaled by `eta`), rotate
  back out. `mueller_matrix_aligned` returns `Z_al` at `eta=1`; scale it yourself.
- **Unaligned population.** Work in the scattering plane exactly as the existing
  random-orientation `scatmat` file is used: `mueller_matrix_random` returns
  `f_tot`, `f_ref`; the absolute remainder matrix is
  `F_unal = Csca_tot*f_tot - eta*Csca_ref*f_ref` (six elements
  `F11 F22 F33 F44 F12 F34`).

### Peel-off

Evaluate the same phase matrix (`Z_al` for the aligned case in the grain frame,
`F_unal` in the scattering plane for the unaligned case) at the forced scattering
geometry toward the observer, and weight by the albedo.

### Emission source term

`j_I = lamI_total`, and the polarized emissivities follow from the intrinsic
`lamI_pol` with the field geometry MoCafe applies:
```
(j_Q, j_U) = lamI_pol * sin^2(gamma) * F_turb * (cos 2*phi, sin 2*phi)
```
with `gamma` the angle between the field and the line of sight and `phi` the
position angle of the projected field (user manual §6.7). The same
`sin^2(gamma) F_turb` multiplies `Cpol_ext` when the extinction matrix is
assembled for the far-infrared-only path (no scattering).

---

## 5. Regenerating the table

The shipped table is `f_align = falign_hd23` (`f_max = 1`, `a_align = 0.0749 um`,
`alpha = 1.80`). To use a different profile:

```
cd tmatrix
make run_scatmat_aligned.x
./run_scatmat_aligned.x profile=my_falign.dat        # two columns: a_eff[um], f_align
```

Cost: the five UBVRI bands (0.36, 0.44, 0.55, 0.64, 0.79 µm) take about 8 minutes
at 32 threads (about 26 minutes at 8); the file is 138 MB, 38 MB gzipped. A single
band or a reduced-grid `test` run is seconds to a minute. Point `scatmat_path` at
the regenerated file. `run_scatmat_aligned.x` with no argument writes the default
UBVRI file; a list of wavelengths writes those bands.

---

## 6. Thread-safety rules

- **Initialization is serial, once.** `load_scatmat_aligned` (via `scatmat_path`),
  the size integrals, and `free_scatmat_aligned` mutate module state and must run
  from serial code before the photon loop.
- **Path queries are concurrent.** `scatmat_band`, `extinction_matrix_aligned`,
  `mueller_matrix_aligned`, `mueller_matrix_random`, `scattering_cross_sections`,
  and `dust_extinction` are pure reads — they write only their own output
  arguments, do no I/O and no allocation — so they are safe to call from OpenMP
  photon threads. The generator's OpenMP correctness (1-thread vs N-thread bitwise
  identical) is already verified; the library layer adds no shared state.

---

## 7. File format and symmetry relations

`tmatrix/output/scatmat_aligned_astrodust_P0.20_Fe0.00_1.400.dat` (ASCII), one
`# lambda =` block per band, each with three data blocks (all per H):

| block | rows | columns |
|---|---|---|
| K | `theta_i` grid (19: 0,5,…,90) | `theta_i` `Cext_al` `Cpol_al` `Cbir_al` `Csca_al` [µm²/H] |
| F | 181: `Theta` = 0,1,…,180 | `Theta` `F_tot(11,22,33,44,12,34)` `F_ref(…)`, α₁-normalized |
| Z | 19×181×37 nested `theta_i→theta_s→phi` | `theta_i theta_s phi  Z11 Z12 … Z44` [µm² sr⁻¹/H] |

Grids: `theta_i` 0(5)90 (19), `theta_s` 0(1)180 (181), `phi` 0(5)180 (37). Block
header scalars per band: `Cext_tot`, `Csca_tot`, `Cext_ref`, `Csca_ref` [µm²/H].

**Eta contract.** `Z_al,cell = eta*Z_al`; `K_al,cell = eta*K_al`; the unaligned
remainder is `F_unal = Csca_tot*F_tot - eta*Csca_ref*F_ref` (absolute units); the
unaligned population adds an isotropic `Cext_tot - eta*Cext_ref` (its `Cpol` and
`Cbir` are zero). Exact by the linearity of the size integral in `f_align`.

**Symmetries** (used to reconstruct the unstored ranges; the library does this):
- `phi → 360-phi`: the two off-diagonal 2×2 blocks (elements 13,14,23,24,31,32,41,42)
  flip sign; store `phi ∈ [0,180]`.
- `theta_i → 180-theta_i`: maps to `theta_s → 180-theta_s`, `phi` unchanged, with
  the same off-diagonal-block sign flip; store `theta_i ∈ [0,90]`.
- `theta_i = 0`: `Z(0; theta_s, phi) = Z(0; theta_s, 0) R(phi)`, the `phi=0`
  matrix being six-element block-diagonal.

---

## 8. What MoCafe must implement itself

SEDust returns the matrices in the grain frame; everything below is MoCafe's:

1. **Frame rotations.** The meridional rotation of the photon Stokes vector into
   the grain frame (`z = B-hat`) before an aligned scatter and back out after it;
   the position-angle rotation that splits polarization into `Q` and `U`.
2. **Direction sampling.** Drawing `(theta_s, phi)` from `Z_al,11` (aligned) or
   `Theta` from `F_unal,11` (unaligned), with the peel-off/forced-scattering
   bookkeeping.
3. **The `exp(-K tau)` transfer step.** Prefer the Peest et al. (2023) closed
   form for a segment of constant density and field: `I`/`Q` evolve through
   `cosh`/`sinh` of `Cpol n s`, `U`/`V` through `cos`/`sin` of `Cbir n s`, all
   scaled by `exp(-Cext n s)`.
4. **Cell inputs.** A magnetic-field direction and an `eta` in every cell (a
   turbulent component is required — Kim et al. 2023 show a smooth field cannot
   reproduce the observed low polarization fraction and narrow dispersion), and a
   vectorized radiation-field tally **only** if alignment is to be computed
   internally rather than assumed (that is the larger, later stage).
5. **Imagery.** Assembling Stokes maps from the peel-off contributions and the
   emission source term, evaluated per observer direction (the angles `gamma`,
   `phi` depend on viewing direction, so they cannot be cached per cell for an
   all-sky output).

---

## 9. Limits of the released data

These are limits of the optics as published, not of the implementation.

- **Birefringence is computed; the release has none.** The regenerated
  orientation-resolved table carries a 4th block (real part of the
  forward-amplitude difference), certified by Kramers-Kronig against the
  dichroism to ~0.1%. `dust_extinction` returns it as `Cbir_ext`, and the aligned
  `K(theta_i)` carries it at `theta_i = 90`. With the release-format default it is
  zero (circular polarization then vanishes).
- **PAHs contribute nothing** (`f_align = 0`, an HD23 assumption); `build_dl07`
  and `build_zubko` have no polarized optics at all, and free any polarized arrays
  left by a prior astrodust build so `lamI_pol` returns zeros rather than stale
  values.
- **Far-infrared and submillimeter polarized emission needs no scattering
  matrix.** There the albedo is essentially zero, so the transfer needs only the
  extinction matrix and the emission vector, both fully exported. The HAWC+ 154 µm
  milestone lies inside this regime and should be attempted first; the scattering
  optics are what the optical update of Seon (2018) additionally requires.

---

## 10. Verification the SEDust side already passed

| Check | Result |
|---|---|
| Engine closure `INT Z11 dOmega` vs `C_sca(jori)` | `2.7e-7` |
| Optical theorem vs `C_ext(jori)` | `2e-14` |
| Orientation average vs GSP `F(Theta)` (incl. F12/F34 signs) | `1.4e-6` |
| Dipole vs closed form; dipole vs T-matrix at `x~0.1` | `2e-16`; `0.5–1.1%` |
| Symmetries (phi mirror; equatorial mapping) | `3.5e-7`; `8.8e-16` |
| Size-integrated `K` vs jori integrals (Cext/Cpol/Cbir) | `1.6e-6` / `2.0e-5` / `1.3e-7` |
| Library: reader round-trip, symmetry reconstruction, eta algebra | bitwise 0 / bitwise 0 / residual 0 |
| Library `K(90)` vs jori vs `dust_extinction` (Cpol/Cbir) | `1.6e-5` / `8.0e-5` |
| `use_dustlib_scatmat.x` closures | `0.997–1.001` |
| `x > 50` omitted scattering fraction | `≤ 3.1e-15` |
| `calc_polext` regression (unchanged) | median `0.0296%` |

---

## 11. Which document to read

| Question | Document |
|---|---|
| Why do this; alignment physics; the transfer equation; the aligned scattering matrix physics and its place in the Peest formalism | `docs/aligned_grain_polarization.pdf` |
| How MoCafe calls the library: signatures, units, argument ranges, lifecycle | `docs/SEDust_user_manual.pdf` §5–§6 |
| How the aligned scattering optics were built and verified | `docs/sedust_polarization_implementation.pdf` §8 |
| The consumption recipe and this contract | this file |

---

## 12. References

- Draine, B. T., & Hensley, B. S. 2021b, ApJ, 919, 65 (DH21b). `2021ApJ...919...65D`
- Hensley, B. S., & Draine, B. T. 2023, ApJ, 948, 55 (HD23). `2023ApJ...948...55H`
- Jones, T. J., et al. 2020, AJ, 160, 167. `2020AJ....160..167J`
- Kim, J., et al. 2023, AJ, 165, 223. `2023AJ....165..223K`
- Peest, C., Camps, P., Stalevski, M., Baes, M., & Siebenmorgen, R. 2017, A&A, 601, A92. `2017A&A...601A..92P`
- Peest, C., et al. 2023, A&A, 673, A112. `2023A&A...673A.112P`
- Planck Collaboration 2020, A&A, 641, A12. `2020A&A...641A..12P`
- Reissl, S., Wolf, S., & Brauer, R. 2016, A&A, 593, A87 (POLARIS I). `2016A&A...593A..87R`
- Seon, K.-I. 2018, ApJ, 862, 87. `2018ApJ...862...87S`
