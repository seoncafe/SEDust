# SEDust polarization: session handoff note

**Date:** 2026 July 19

This is a short handoff note, not a specification. It records what was completed, what
was decided and why, how far the numbers can be trusted, and where to start next. The
detailed material lives elsewhere; see §6 for which document answers which question.

---

## 1. What is done, what is left

Completed in this session:

- **Orientation-resolved optics.** `sed/src/q_table_jori.f90` reads the three orientation
  blocks of the DH21 astrodust spheroid table and forms the polarized and
  random-orientation combinations, plus the HD23 alignment fit `falign_hd23(a)`.
- **Polarized emission.** `sed_init` builds `Cpol` (polarized absorption) and
  `falign_ad`; the grain loop accumulates the polarized channel at every site, so
  stochastic heating is treated exactly; `dust_emission` gained an optional `lamI_pol`
  output of the same shape and units as `lamI_total`.
- **Polarized extinction.** `sed_init` also builds `Cpol_ext`, and
  `calc_kext_astrodust.x` writes the size-integrated result as an eighth column,
  `C_polext/H` [cm^2/H], of `data/kext_astrodust_MW.dat`.
- **Random-orientation scattering matrix.** Computed and stored for five optical
  bands; see §4.

Two setters change the alignment efficiency on an existing model, both
re-exported by `dust_lib`:

```fortran
call dust_set_alignment(m, f_max, a_align, alpha_align [, status])
call dust_set_alignment_profile(m, aeff_in, falign_in [, status])
```

The first installs the HD23 power law `f_max / (1 + (a_align/a)**alpha_align)`,
with `a_align` in um. The second interpolates an arbitrary tabulated
`(aeff_in [um], falign_in)` profile onto each population's radius grid, linearly
in `log a` — the route for a RAT-derived `R(a)` or the GRADE-POL exponential.

**Neither requires re-solving `P(T)`.** The alignment efficiency is a size
weight applied *outside* the temperature solution, so it never enters `P(T)` and
`lamI_total` is bit-for-bit unchanged across a setter call. A scan over
alignment parameters costs one `dust_emission` call per point, not a model
rebuild.

`falign_in` must lie in `[-1, 1]` and is **rejected, not clamped**, outside it.
Negative values inside the range are deliberately allowed: slow internal
relaxation gives `Q_X = -0.1`, a negative reduction factor that flips the
polarization direction by 90 degrees, which is the accepted origin of
polarization parallel to B at millimeter wavelengths. Clamping at zero would
delete that physics silently.

Both exported quantities are **intrinsic**: the size-distribution integral and the
`f_align` weighting are done, but the `sin^2(gamma)` geometric factor, the turbulent
depolarization, and the field position angle are left to the consumer.

Not started: polarized transfer on the consumer side. That is the next session.

---

## 2. Decisions and their reasons

**`f_align` is the HD23 analytic fit, not the release column.** The alignment column of
`size_distribution.dat` differs from the fit by up to 0.32%. The fit is definitive; do
not mix the two sources.

**`Cabs` keeps our exact T-matrix orientation average while `Cpol` comes from the release
table's 1/3 trace (the modified picket fence approximation).** This mixture was chosen
deliberately, and v1.20 deliberately leaves it in place: it keeps the existing SEDs
byte-identical, and in the far infrared the two averages agree closely. Measured on the shared grid, `|trace-average / exact - 1|` for
`Q_abs` has a median of **0.022%** beyond 30 um, with a worst case of 2.0% among the
grains carrying 99.999% of the geometric cross section.

**That mixture is not acceptable in the ultraviolet**, where the median rises to **0.21%**
and the worst case to **9.5%**. The emission path is safe because polarized emission is a
far-infrared phenomenon, but the polarized extinction column spans the full wavelength
grid. Revisit this before using column 8 in the ultraviolet or the optical. The clean
resolution is to recompute `Cabs` on the same convention, which breaks byte-identity of
the existing SEDs and so needs to be a separate, deliberately validated step. Summary for
a consumer: **harmless in the optical and the infrared, unresolved in the ultraviolet.**

**Polarized extinction is exported as a table column, not a library accessor.** A
`dust_cpol_ext` accessor was considered and rejected: `dust_lib` is an emission-facing
API and exposes no extinction surface at all, so such an accessor would be the only
extinction entry point in the module. The extinction table is already how a host consumes
extinction. If `dust_lib` later gains a full extinction API, `Cpol_ext` should join it
then.

**The module array `Cpol` was not renamed to `Cpol_abs`.** It appears at roughly 40 sites
and its name is part of the public surface (`grain_pop_t%Cpol`, the optional
`set_pop(..., Cpol_in=...)`). The declaration carries an explicit note instead: `Cpol` is
the absorption one, `Cpol_ext` its extinction twin.

---

## 3. Verification

| Check | Result |
|---|---|
| Column 8 against the released `polarized_extinction.dat` | median 0.0296%, max 0.94%, over the 992 points with lambda >= 0.11 um (the reference turns negative below that) |
| Column 8 against `calc_polext.x` | agree to 7e-8 relative, the precision the ASCII format carries; two independent routes to the same quantity |
| Columns 1-7 of `kext_astrodust_MW.dat` | byte-identical to the pre-change file, 1129 lines |
| `astrodust_irem_ours_S1/S2/PAH.dat`, `dl07_sed_ours_mw31_60.dat` | byte-identical to the reference outputs |
| Build | zero warnings |

Intrinsic polarization fraction `lamI_pol / lamI_total`, Mathis ISRF at U = 1.585:

| lambda | fraction |
|---|---|
| 12 um | 0.00% |
| 100 um | 13.67% |
| 154 um | 17.15% |
| 350 um | 18.79% |
| 850 um | 19.15% |
| median over 300-3000 um | 19.18% |

The mid-infrared zero is correct: there the emission is dominated by stochastically
heated PAHs, which HD23 take to be unaligned.

**The 19.15% versus Planck's 22% is not a bug.** Planck measures
`p_max = 22.0 (+3.5, -1.4) %` at 353 GHz. Ours is the maximum before any geometric
factor, and every real line of sight reduces it, so the model cannot reach the observed
maximum. This is a known property of the HD23 model: the axial ratio b/a = 1.4 is the
minimum non-sphericity that DH21b found the starlight polarization efficiency to require,
and it was not tuned to maximize the submillimeter polarization fraction. Do not "fix" it
by adjusting `f_align` or the axial ratio without treating that as a change to the dust
model itself.

---

## 4. Limits of the released data

These are limits of the optics as published, not of the implementation.

**No circular polarization or birefringence.** The table header offers `Q_ext`, `Q_abs`
and `Q_sca` only. None of them is a phase retardation, so `C_circ` cannot be built and
must be set to zero. Obtaining it would mean computing the fixed-orientation amplitude
matrix ourselves; `tmatrix/src/ampld.lp.f` exists in the tree for this and is currently
excluded from the build.

**No scattering matrix for aligned grains.** The release gives the integrated `Q_sca`
with no angular information, so scattering by an aligned spheroid cannot be modeled from
it. Reaching optical polarized transfer would require computing wavelength-resolved
Mueller matrices (`ampld.lp.f`, DDSCAT, ADDA or CosTuuM) and a matrix-based scattering
treatment on the consumer side, with a data volume of `lambda x theta x orientation x 16`
components.

**The random-orientation scattering matrix, on the other hand, is available.** This is a
different object from the previous paragraph — do not read the two as contradicting each
other. It describes how an *unaligned* astrodust population scatters, and carries no
information about grain orientation with respect to B, but it is what a first optical
model needs for the scattering source term, and it already improves on Seon (2018), whose
phase matrix was spherical Mie.

| | |
|---|---|
| file | `tmatrix/output/scatmat_astrodust_P0.20_Fe0.00_1.400.dat` |
| produced by | `tmatrix/driver/run_scatmat.f90` → `run_scatmat.x` |
| wavelengths | 0.36, 0.44, 0.55, 0.64, 0.79 um (approximately UBVRI), one block each |
| angles | 181 rows, theta = 0, 1, ..., 180 deg, 0 = forward |
| columns | `F11 F22 F33 F44 F12 F34` — the six independent elements of a randomly oriented particle with a plane of symmetry |
| normalization | `(1/2) integral F11 dcos(theta) = 1`; multiply by the `Csca` in each block header for absolute units |

Degree of linear polarization for unpolarized incident light is `-F12/F11`. The five bands
were chosen to cover the planned scope (V band first, one or two more optical bands
possible; no FIR or UV use foreseen). The full 1129-point sweep costs about 1.7 hours and
was not run. Extending needs no code change: `./run_scatmat.x all`, or an explicit list of
wavelengths, at about 5.5 s per wavelength.

**Far-infrared polarized emission needs neither.** In that regime the albedo is
essentially zero, so the transfer needs only the extinction matrix and the emission
vector, both fully determined by what we export. The HAWC+ 154 um comparison lies inside
this regime, so the far-infrared milestone is complete without any scattering matrix and
should be attempted first.

**PAHs contribute nothing** (`f_align = 0`, an HD23 assumption), and **DL07 and Zubko have
no polarized optics** at all; `build_dl07` frees any polarized arrays left by a previous
astrodust initialization, so `lamI_pol` returns zeros rather than stale values.

---

## 5. Where to start next

Work units, in order:

1. Get a magnetic field direction into every cell. This is the largest single gap and
   everything else waits on it.
2. Derive the two angles for each cell and viewing direction: `gamma` between the field
   and the line of sight (not a sky-plane angle), and `phi` for the projected field.
   Check them against a uniform field before trusting anything downstream.
3. Build the Stokes emission vector cell by cell from `lamI_total` and `lamI_pol`,
   applying `sin^2(gamma)`, the turbulent factor, and the position angle. With no
   extinction this already produces a polarization map worth inspecting on its own.
4. Apply the extinction matrix along the line of sight using the Peest closed form, with
   `C_ext` and `C_polext` from the table and `C_circ = 0`.
5. Validate. Run the Peest analytic benchmarks first, since they isolate the transfer
   solver from the dust model, then compare against the HAWC+ 154 um maps of NGC 891.
   Kim et al. (2023) needed turbulent depolarization in the disk, so a model that gets
   the field geometry right but omits it will over-polarize there.

---

## 6. Which document to read

| Question | Document |
|---|---|
| Why do this, what is the alignment physics, what is the transfer equation, what are the acceptance targets | `docs/aligned_grain_polarization.pdf` |
| How an RT code calls the library: signatures, units, what is already integrated | `docs/SEDust_user_manual.pdf` |
| Implementation detail of the polarization update | `docs/sedust_polarization_implementation.pdf` |
| What was decided this session and how far the numbers go | this file |

`aligned_grain_polarization.tex` was written before any of this existed, so where it says
SEDust would need to add `Cpol`, read that SEDust has it. Its §6.1, §7 and §8.1 have been
updated to the current state; §§1-5 and §§9-11 are physics and remain valid as written.

---

## 7. References

- Draine, B. T., & Hensley, B. S. 2021b, ApJ, 919, 65 (DH21b). `2021ApJ...919...65D`
- Hensley, B. S., & Draine, B. T. 2023, ApJ, 948, 55 (HD23). `2023ApJ...948...55H`
- Jones, T. J., et al. 2020, AJ, 160, 167. `2020AJ....160..167J`
- Kim, J., et al. 2023, AJ, 165, 223. `2023AJ....165..223K`
- Peest, C., et al. 2023, A&A, 673, A112. `2023A&A...673A.112P`
- Planck Collaboration 2020, A&A, 641, A12. `2020A&A...641A..12P`
- Reissl, S., Wolf, S., & Brauer, R. 2016, A&A, 593, A87 (POLARIS I). `2016A&A...593A..87R`
- Seon, K.-I. 2018, ApJ, 862, 87. `2018ApJ...862...87S`
