# SEDust

A self-contained library for the optical properties and thermal emission of
interstellar dust: dielectric functions -> T-matrix / Mie cross sections ->
grain temperature distributions -> emergent infrared SED.

SEDust is **model-agnostic**. The HD23 astrodust+PAH model, the Draine & Li
(2007) carbonaceous+silicate model, the Mathis, Rumpl & Nordsieck (1977)
graphite+silicate model, and the Zubko et al. (2004) BARE-GR-S model are handled
as peers through one derived type (`dust_model_t`) and one emission call
(`dust_emission`). The solver core does not know which model it is running.

The whole package can also be linked into a Fortran 3D radiative-transfer code
as a static library, `libsedust.a`, with a two-step API: initialize once, then
solve one cell at a time. Alongside the emission call, `dust_extinction` returns
the model's size-integrated extinction, scattering, and asymmetry cross sections
per H on the same wavelength grid, so a transfer code takes its opacity and its
emission from one model object. Those cross sections come from the precomputed
`data/<model>/kext_*.dat` table the builder loaded, which is the same product
`calc_kext.x` writes; the size integral behind it is `size_integrated_extinction`.
Every model builder and solver call takes an
optional `status` argument, so a missing input file or an invalid model is
reported back to the host instead of stopping the process.

## Which version

SEDust is split into two versions because most radiative-transfer codes either
do not carry polarization at all, or carry it only for spherical grains.

This is **version 1.00**, the scalar branch: unpolarized cross sections, thermal
emission, and the `dust_emission` / `dust_extinction` API for a radiative-transfer
host. It covers that common case.

Use **version 1.20** instead for polarized transfer that accounts for the
non-spherical (**spheroidal**) grain shape -- dichroic extinction, birefringence,
and scattering by aligned spheroidal grains. It adds the direction-dependent
extinction matrix and the aligned-grain scattering matrix on top of everything
here.

The two versions are branches of one repository, and `v1.00` is the default
branch, so a plain clone lands here:

```
git clone git@github.com:seoncafe/SEDust.git          # this version, 1.00
git clone -b v1.20 git@github.com:seoncafe/SEDust.git # the polarized version
```

An existing clone switches with `git checkout v1.20`.

## Layout

```
SEDust/
  sed/          the SED solver: cross sections, enthalpy, P(T), emission
    src/        library modules + drivers
    rt_example/ minimal example of linking libsedust.a into an RT code
  mc/           Draine & Anderson (1985) Monte Carlo solver (independent check)
  tmatrix/      Mishchenko T-matrix engine + driver; writes the Q table
  data/         one directory per dust model, plus what the models share
    astrodust/  everything that model owns: `sedust_astrodust.h5` (the
    dl07/       primary form -- its wavelength axis, cross-section tables
    mrn/        and extinction curve in one file), the same as text
    zubko/
                (`q_*.dat`, `kext_*.dat`), and, where the model IS a set of
                files, its definition (the ZDA config, optics, calorimetry)
    dielectric/ shared material data: the D03 / DH21 / D16 optical
                constants and the PAH cross sections.  A dielectric function
                is not one model's -- DL07, MRN and Zubko read the same D03
                astrosilicate -- so it does not live inside a model directory
    release/    published reference tables (HD23, Draine) and the HD23 size
                distribution
  docs/         technical reports and the library user manual
  pyutil/       small Python helpers (radiation fields, SED from Cabs, and
                sedust_h5 -- the reader for data/<model>/sedust_<model>.h5)
```

Everything the code reads at run time ships with the package. There are no
paths outside this directory.

## What it computes

- **Four grain models as peers** — HD23 astrodust+PAH, Draine & Li (2007)
  carbonaceous+silicate, Mathis, Rumpl & Nordsieck (1977) graphite+silicate,
  Zubko et al. (2004) BARE-GR-S — behind one `dust_model_t` and one
  `dust_emission` call, plus a fifth builder that takes a model defined by a
  descriptor file.
- **Optics from first principles** — dielectric functions through Mie and the
  Mishchenko T-matrix (spheroids) to stored `Q` tables, and the size-integrated
  `C_ext` / `C_abs` / `C_sca` / `<cos>` per H that a transfer code takes as its
  opacity.
- **The ionizing band on the same grid** — one wavelength axis reaching to
  1.0e-4 um, with `include_euv` choosing the view rather than the file, so a
  host that transports ionizing radiation reads the band off the table.
  See [README_HOWTO.md](README_HOWTO.md#the-ionizing-band).
- **Four stochastic-heating solvers** — two Guhathakurta & Draine matrix
  solvers, an energy-space transition matrix with three cooling kernels
  (thermal-discrete, thermal-continuous, exact-statistical), and a pure
  equilibrium mode; plus an independent Monte Carlo solver in `mc/` that shares
  no code with them.
  See [README_HOWTO.md](README_HOWTO.md#stochastic-heating).
- **One command-line vocabulary** for every driver, where a word names the same
  physics whichever program reads it and a program without a referent for it
  refuses by name rather than ignoring it.
  See [README_HOWTO.md](README_HOWTO.md#build-and-run).
- **HDF5 products** — one `data/<model>/sedust_<model>.h5` per model carrying
  its wavelength axis, cross-section tables and extinction curve, with the same
  content as text beside it.

## Documentation

- **[README_HOWTO.md](README_HOWTO.md)** — build, run, every run setting the
  drivers take, the ionizing band, and the stochastic-heating solvers.
- [`docs/astrodust_sed_report.pdf`](docs/astrodust_sed_report.pdf) — the astrodust+PAH pipeline, its validation
  against the HD23 release, and the resolution of the far-infrared offset.
- [`docs/SEDust_user_manual.pdf`](docs/SEDust_user_manual.pdf) — the `libsedust.a` API: model builders,
  channels, solver options, and how to link it into an RT code.
- [`docs/mc_pT_report.pdf`](docs/mc_pT_report.pdf) — the Monte Carlo algorithm, its adaptive-grid
  engines, and its validation against the matrix solvers.

## References

- [Mathis, Rumpl & Nordsieck 1977, ApJ, 217, 425](https://ui.adsabs.harvard.edu/abs/1977ApJ...217..425M/abstract) — the a^-3.5 power-law graphite+silicate model.
- [Draine & Lee 1984, ApJ, 285, 89](https://ui.adsabs.harvard.edu/abs/1984ApJ...285...89D/abstract) — the graphite and astrosilicate dielectric functions, and the MRN abundances adopted here.
- [Draine & Anderson 1985, ApJ, 292, 494](https://ui.adsabs.harvard.edu/abs/1985ApJ...292..494D/abstract) — the Monte Carlo temperature-history method the `mc/` solver follows.
- [Guhathakurta & Draine 1989, ApJ, 345, 230](https://ui.adsabs.harvard.edu/abs/1989ApJ...345..230G/abstract) — the transition-matrix solution for `P(T)` of stochastically heated grains.
- [Draine & Li 2001, ApJ, 551, 807](https://ui.adsabs.harvard.edu/abs/2001ApJ...551..807D/abstract) — the grain enthalpy used here, and the LD01 PAH cross-section vintage.
- [Weingartner & Draine 2001, ApJ, 548, 296](https://ui.adsabs.harvard.edu/abs/2001ApJ...548..296W/abstract) — the Milky Way size distributions the DL07 model is built on.
- [Zubko, Dwek & Arendt 2004, ApJS, 152, 211](https://ui.adsabs.harvard.edu/abs/2004ApJS..152..211Z/abstract) — the ZDA BARE-GR-S composition, size distributions, and calorimetry.
- [Draine & Li 2007, ApJ, 657, 810](https://ui.adsabs.harvard.edu/abs/2007ApJ...657..810D/abstract) — the DL07 PAH cross sections and emission model.
- [Mishchenko & Travis 1998, JQSRT, 60, 309](https://ui.adsabs.harvard.edu/abs/1998JQSRT..60..309M/abstract) — the T-matrix method used for the spheroidal astrodust optics.
- [Draine & Hensley 2021, ApJ, 909, 94](https://ui.adsabs.harvard.edu/abs/2021ApJ...909...94D/abstract) — the astrodust dielectric function.
- [Hensley & Draine 2023, ApJ, 948, 55](https://ui.adsabs.harvard.edu/abs/2023ApJ...948...55H/abstract) — the astrodust+PAH model this package reproduces.

## Author

Kwang-il Seon (KASI/UST)

---

Last updated: 2026-08-21 20:09 KST
