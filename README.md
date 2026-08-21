# SEDust

Version 1.20.

A self-contained library for the optical properties and thermal emission of
interstellar dust: dielectric functions -> T-matrix / Mie cross sections ->
grain temperature distributions -> emergent infrared SED.

SEDust owns the dust physics. A radiative-transfer code that links it is
responsible for photon transport and geometry only: everything about the
grains, from cross sections through temperature distributions to emission and
polarization, is computed here and handed over as wavelength-resolved
quantities per H atom.

SEDust is **model-agnostic**. The HD23 astrodust+PAH model, the Draine & Li
(2007) carbonaceous+silicate model, and the Zubko et al. (2004) BARE-GR-S model
are handled as peers through one derived type (`dust_model_t`) and one emission
call (`dust_emission`). The solver core does not know which model it is running.

The whole package can also be linked into a Fortran 3D radiative-transfer code
as a static library, `libsedust.a`, with a two-step API: initialize once, then
solve one cell at a time. Every model builder and solver call takes an optional
`status` argument, so a missing input file or an invalid model is reported back
to the host instead of stopping the process.

## Which version

SEDust is split into two versions because most radiative-transfer codes either
do not carry polarization at all, or carry it only for spherical grains.

This is **version 1.20**, the polarized branch. Use it for polarized transfer
that accounts for the non-spherical (**spheroidal**) grain shape -- dichroic
extinction, birefringence, and scattering by aligned spheroidal grains -- for
which it provides the direction-dependent extinction matrix `K(theta_i)` and the
aligned-grain scattering matrix, in addition to the scalar optics and thermal
emission.

For transfer that does not carry polarization, use **version 1.00**, the scalar
branch: it computes the same unpolarized cross sections and emission through the
same `dust_emission` / `dust_extinction` API without the polarized optics or
their data tables. (Version 1.20 can also run scalar-only, by building a model
with `load_polarized_optics = .false.`, when a host wants both from one build.)

**1.20 is a superset of 1.00.** Every scalar capability of the scalar branch is
present here and computes the same numbers: the same builders, the same
`dust_emission` / `dust_extinction`, the same `lam_min` extension into the
ionizing band, the same `calc_kext.x` products. The differences are the added
polarized quantities and one consequence of adding them — the `status` codes of
`sed_init` / `build_astrodust` are numbered differently, because codes 3-5 are
taken by the polarized inputs, so what is `3` and `4` in v1.00 is `6` and `7`
here. Code that only checks `status /= 0` is unaffected.

The two versions are branches of one repository. `v1.00` is the default branch,
so this version has to be asked for by name:

```
git clone -b v1.20 git@github.com:seoncafe/SEDust.git # this version, 1.20
git clone git@github.com:seoncafe/SEDust.git          # the scalar version
```

An existing clone switches with `git checkout v1.20`.

## Layout

```
SEDust/
  sed/          the SED solver: cross sections, enthalpy, P(T), emission
    src/        library modules + drivers
    rt_example/ examples of linking libsedust.a into an RT code, one
                minimal and two showing the polarized quantities
  mc/           Draine & Anderson (1985) Monte Carlo solver (independent check)
  tmatrix/      Mishchenko T-matrix engine + driver; writes the Q table
  data/         one directory per dust model, plus what the models share
    astrodust/  everything that model owns: `sedust_astrodust.h5` (the
    dl07/       primary form -- its wavelength axis, cross-section tables
    zubko/      and extinction curve in one file), the same as text
                (`q_*.dat`, `kext_*.dat`), and, where the model IS a set of
                files, its definition (the ZDA config, optics, calorimetry)
    dielectric/ shared material data: the D03 / DH21 / D16 optical
                constants and the PAH cross sections.  A dielectric function
                is not one model's -- DL07 and Zubko read the same D03
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

- **Three grain models as peers** — HD23 astrodust+PAH, Draine & Li (2007)
  carbonaceous+silicate, Zubko et al. (2004) BARE-GR-S — behind one
  `dust_model_t` and one `dust_emission` call, plus a fourth builder that takes
  a model defined by a descriptor file.
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
- **Polarized optics of spheroidal grains** — dichroic and birefringent
  extinction, polarized emission, and the scattering matrix of aligned grains,
  with a runtime-settable alignment efficiency.
  See [README_HOWTO.md](README_HOWTO.md#polarization).

## Documentation

- **[README_HOWTO.md](README_HOWTO.md)** — build, run, every run setting the
  drivers take, the ionizing band, the stochastic-heating solvers, and the
  polarized optics.
- [`docs/astrodust_sed_report.pdf`](docs/astrodust_sed_report.pdf) — the astrodust+PAH pipeline, its validation
  against the HD23 release, and the resolution of the far-infrared offset.
- [`docs/SEDust_user_manual.pdf`](docs/SEDust_user_manual.pdf) — the `libsedust.a` API: model builders,
  channels, solver options, and how to link it into an RT code.
- [`docs/mc_pT_report.pdf`](docs/mc_pT_report.pdf) — the Monte Carlo algorithm, its adaptive-grid
  engines, and its validation against the matrix solvers.
- [`docs/sedust_polarization_implementation.pdf`](docs/sedust_polarization_implementation.pdf) — how the polarized optics are
  built: the orientation-resolved table, the derived cross sections, the
  aligned-grain scattering matrix and its two-layer API, the implementation
  decisions and their reasons, and the verification.
- [`docs/aligned_grain_polarization.pdf`](docs/aligned_grain_polarization.pdf) — background on grain alignment and
  polarized radiative transfer, and what a radiative-transfer code would need
  in order to use the polarized optics.

## References

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

Last updated: 2026-08-21 16:58 KST
