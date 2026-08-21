"""Reader for the SEDust optics products, data/<model>/sedust_<model>.h5.

One file per dust model -- astrodust, dl07, mrn, zubko -- holding the wavelength
axis, the (lambda, a_eff) cross-section tables of each grain population, the
size-integrated extinction curve, and (in the polarized branch) the
orientation-resolved and scattering-matrix products.  `sed/calc_qtable.x` and
`sed/calc_kext.x` write them; `sed/src/sedust_product.f90` is the Fortran
counterpart of this module, and the two are meant to return the same numbers
from the same file.

THE IONIZING BAND.  Each file carries ONE wavelength axis, the widest the model
has, and /grid records `i_lyman`, the 1-based index of the LAST wavelength at
or below the Lyman limit (0.0912 um) -- the covering rule of the Fortran
writer's `lyman_index` -- so the non-ionizing view starts at or just below the
limit and a host whose transport floor is the limit interpolates inside the
table.  Every reader here takes

    include_euv=False   lambda[i_lyman:] and the same slice of every
                        wavelength-indexed array -- the non-ionizing product
    include_euv=True    the whole axis

and takes the cut from the file's own `i_lyman`, never from a wavelength
comparison of its own, so that Python and Fortran slice at the same node.  The
paired *_euv text products are views of one array here, and cannot disagree.

ARRAY SHAPES.  Fortran declares the cross sections (n_lam, n_a) and HDF5 stores
them in that order, so h5py reads them back as **(n_a, n_lam)**: each grain
radius is one contiguous spectrum, and the wavelength is the LAST axis.  The
readers below keep that orientation.

    from pyutil.sedust_h5 import read_kext, read_qtable, describe

    k = read_kext('data/zubko/sedust_zubko.h5')      # non-ionizing, 865 points
    q = read_qtable('data/dl07/sedust_dl07.h5', 'sil', include_euv=True)
    describe('data/astrodust/sedust_astrodust.h5')

Nothing here computes anything: a value that is not in the file comes back as
None, never as a zero that could be mistaken for a computed one.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from typing import Dict, List, Optional

import numpy as np
import h5py

__all__ = [
    'model_dir',
    'model_file',
    'read_grid',
    'read_kext',
    'read_qtable',
    'read_polarized',
    'describe',
    'Grid',
    'Kext',
    'QTable',
    'Polarized',
]

# The Lyman limit [um].  Only used to check what a file says about itself; the
# cut always comes from the file's own i_lyman.
LAM_LYMAN_UM = 0.0912

MODELS = ('astrodust', 'dl07', 'mrn', 'zubko')


def model_dir(data_dir: str, model: str) -> str:
    """The directory a model owns.

    Everything a model owns lives in one place, ``data_dir/<model>/``: its HDF5
    product, its cross-section tables, its extinction curves and, where the
    model IS a set of files (ZDA), its definition. What is shared between
    models -- the dielectric functions, the published reference tables -- sits
    beside those directories rather than inside one of them. Mirrors the
    Fortran ``sedust_dir``.
    """
    return os.path.join(data_dir, model)


def model_file(data_dir: str, model: str) -> str:
    """Path of a model's product.  Mirrors Fortran `sedust_h5_file`."""
    return os.path.join(model_dir(data_dir, model), f'sedust_{model}.h5')


def _attrs(obj) -> Dict[str, object]:
    """Attributes with HDF5's fixed-length strings decoded to str."""
    out = {}
    for k, v in obj.attrs.items():
        if isinstance(v, bytes):
            v = v.decode('utf-8', 'replace')
        elif isinstance(v, np.ndarray) and v.dtype.kind == 'S':
            v = [s.decode('utf-8', 'replace') for s in v]
        elif isinstance(v, np.generic):
            v = v.item()
        out[k] = v
    return out


def _lam_start(f: h5py.File, include_euv: bool) -> int:
    """0-based index the wavelength axis is cut at.

    A file that states no i_lyman is read whole, which is the reading under
    which a narrow request returns everything rather than silently dropping
    the short end of an axis whose Lyman node is unknown.
    """
    if include_euv:
        return 0
    i_lyman = int(f['grid'].attrs.get('i_lyman', 1))
    return max(0, i_lyman - 1)


# ---------------------------------------------------------------------------
# /grid
# ---------------------------------------------------------------------------

@dataclass
class Grid:
    lam: np.ndarray                 # [um], strictly increasing
    i_lyman: int                    # 1-based, as the file states it
    n_full: int                     # length of the uncut axis
    include_euv: bool
    attrs: Dict[str, object] = field(default_factory=dict)

    @property
    def n(self) -> int:
        return len(self.lam)


def read_grid(path: str, include_euv: bool = False) -> Grid:
    """The model's wavelength axis, cut as `include_euv` asks."""
    with h5py.File(path, 'r') as f:
        i0 = _lam_start(f, include_euv)
        d = f['grid/lambda']
        return Grid(lam=d[i0:], i_lyman=int(f['grid'].attrs.get('i_lyman', 1)),
                    n_full=d.shape[0], include_euv=include_euv,
                    attrs=_attrs(f['grid']))


# ---------------------------------------------------------------------------
# /kext
# ---------------------------------------------------------------------------

@dataclass
class Kext:
    lam: np.ndarray                 # [um]
    albedo: np.ndarray              # C_sca/C_ext
    g: np.ndarray                   # <cos>
    C_ext: np.ndarray               # [cm^2/H]
    C_abs: np.ndarray
    C_sca: np.ndarray
    K_abs: Optional[np.ndarray]     # [cm^2/g]; None when the model states no mass
    C_polext: Optional[np.ndarray]  # polarized branch only
    M_dust_per_H: Optional[float]   # [g/H]
    attrs: Dict[str, object] = field(default_factory=dict)


def read_kext(path: str, include_euv: bool = False) -> Kext:
    """The size-integrated extinction curve.

    This is the same curve `calc_kext.x` writes as text beside the file, and
    the one an RT host gets back from `dust_extinction`.

    Raises KeyError when the file carries no /kext -- which is what a product
    written by `calc_qtable.x` and not yet extended by `calc_kext.x` looks
    like.  Absence is an error rather than an empty array on purpose: a curve
    of zeros is a physical statement, and this is not one.
    """
    with h5py.File(path, 'r') as f:
        if 'kext' not in f:
            raise KeyError(f'{path} carries no /kext; run  ./calc_kext.x '
                           f'{f.attrs.get("model", "<model>")} euv  after calc_qtable.x')
        i0 = _lam_start(f, include_euv)
        g = f['kext']
        opt = lambda n: g[n][i0:] if n in g else None
        return Kext(lam=f['grid/lambda'][i0:],
                    albedo=g['albedo'][i0:], g=g['g'][i0:],
                    C_ext=g['C_ext'][i0:], C_abs=g['C_abs'][i0:],
                    C_sca=g['C_sca'][i0:],
                    K_abs=opt('K_abs'), C_polext=opt('C_polext'),
                    M_dust_per_H=g.attrs.get('M_dust_per_H'),
                    attrs=_attrs(g))


# ---------------------------------------------------------------------------
# /qtable
# ---------------------------------------------------------------------------

@dataclass
class QTable:
    component: str
    lam: np.ndarray                 # (n_lam) [um]
    a_eff: np.ndarray               # (n_a) [um]
    Q_abs: np.ndarray               # (n_a, n_lam)
    Q_ext: Optional[np.ndarray]     # None for a population that does not scatter
    Q_sca: Optional[np.ndarray]
    g: Optional[np.ndarray]
    regime: Optional[np.ndarray]    # astrodust: 0 T-matrix, 10 Rayleigh, 20 geometric
    rho_bulk: Optional[float]       # [g/cm^3]
    scatters: bool
    attrs: Dict[str, object] = field(default_factory=dict)


def components(path: str) -> List[str]:
    """Grain populations the file carries, in file order."""
    with h5py.File(path, 'r') as f:
        return sorted(f['qtable']) if 'qtable' in f else []


def read_qtable(path: str, component: str, include_euv: bool = False) -> QTable:
    """One grain population's cross sections, (n_a, n_lam) in this orientation.

    A population whose cross sections are an absorption prescription rather
    than a Mie solution -- the PAHs -- carries Q_abs alone, and Q_ext, Q_sca
    and g come back None with `scatters` False.  The Fortran size integral
    substitutes Q_ext = Q_abs and Q_sca = g = 0 for such a population; that
    substitution is the solver's, not the file's, so it is not made here.
    """
    with h5py.File(path, 'r') as f:
        key = f'qtable/{component}'
        if key not in f:
            raise KeyError(f'{path} has no {key}; it carries '
                           f'{sorted(f["qtable"]) if "qtable" in f else "no /qtable"}')
        i0 = _lam_start(f, include_euv)
        g = f[key]
        opt = lambda n: g[n][:, i0:] if n in g else None
        reg = g['regime'][:, i0:].astype(int) if 'regime' in g else None
        return QTable(component=component,
                      lam=f['grid/lambda'][i0:], a_eff=g['a_eff'][:],
                      Q_abs=g['Q_abs'][:, i0:],
                      Q_ext=opt('Q_ext'), Q_sca=opt('Q_sca'), g=opt('g'),
                      regime=reg,
                      rho_bulk=g.attrs.get('rho_bulk_g_cm3'),
                      scatters='Q_sca' in g,
                      attrs=_attrs(g))


# ---------------------------------------------------------------------------
# /polarized  (v1.20)
# ---------------------------------------------------------------------------

@dataclass
class Polarized:
    lam: np.ndarray                       # its own axis, shorter than /grid/lambda
    covers_euv: bool
    pol_valid_from: Optional[float]       # [um]; shortward of it there is no table
    qjori: Dict[str, np.ndarray] = field(default_factory=dict)
    scatmat_random: Dict[str, np.ndarray] = field(default_factory=dict)
    scatmat_aligned: Dict[str, np.ndarray] = field(default_factory=dict)
    attrs: Dict[str, object] = field(default_factory=dict)


def read_polarized(path: str) -> Polarized:
    """The orientation-resolved and scattering-matrix products of the
    polarized branch.

    These have their OWN wavelength axis, /polarized/lambda, which stops at the
    Lyman limit -- the orientation-resolved table is not computed shortward of
    it.  That is why they are not sliced by /grid's i_lyman and why the group
    records `covers_euv` and `pol_valid_from`: a reader has to be able to tell
    a band that was not computed from one that came out zero.

    Raises KeyError on a scalar-branch file, which carries no /polarized.
    """
    with h5py.File(path, 'r') as f:
        if 'polarized' not in f:
            raise KeyError(f'{path} carries no /polarized (scalar branch)')
        p = f['polarized']
        a = _attrs(p)
        out = Polarized(lam=p['lambda'][:] if 'lambda' in p else np.empty(0),
                        covers_euv=bool(a.get('covers_euv', 0)),
                        pol_valid_from=a.get('pol_valid_from'),
                        attrs=a)
        for sub, into in (('qjori', out.qjori),
                          ('scatmat_random', out.scatmat_random),
                          ('scatmat_aligned', out.scatmat_aligned)):
            if sub in p:
                for name in p[sub]:
                    into[name] = p[f'{sub}/{name}'][...]
                into['_attrs'] = _attrs(p[sub])
        return out


# ---------------------------------------------------------------------------
# describe
# ---------------------------------------------------------------------------

def describe(path: str, show: bool = True) -> str:
    """A readable summary of one product: axes, groups, shapes, provenance.

    The provenance attributes are printed in full rather than summarized.  The
    format carries them because a shipped optics table in this tree once named
    a dielectric function it had demonstrably not been computed from, and a
    stale label of that kind costs real time to unpick.
    """
    L: List[str] = []
    with h5py.File(path, 'r') as f:
        L.append(f'{path}  ({os.path.getsize(path) / 1e6:.1f} MB)')
        for k, v in _attrs(f).items():
            L.append(f'  {k:<12s} {v}')

        if 'grid' in f:
            lam = f['grid/lambda']
            ga = _attrs(f['grid'])
            ily = int(ga.get('i_lyman', 1))
            L.append(f'  /grid/lambda  {lam.shape[0]} points, '
                     f'{lam[0]:.5e} .. {lam[-1]:.5e} um')
            L.append(f'    i_lyman = {ily}  ->  lambda = {lam[ily-1]:.6e} um; '
                     f'include_euv=False gives {lam.shape[0] - ily + 1} points')

        if 'qtable' in f:
            L.append('  /qtable')
            for c in sorted(f['qtable']):
                g = f[f'qtable/{c}']
                ca = _attrs(g)
                sets = ', '.join(f'{n}{f[f"qtable/{c}/{n}"].shape}' for n in sorted(g))
                L.append(f'    {c}')
                L.append(f'      {sets}')
                L.append(f'      rho = {ca.get("rho_bulk_g_cm3")} g/cm^3, '
                         f'scatters: {ca.get("scatters")}')
                L.append(f'      method: {ca.get("method")}')
                L.append(f'      source: {ca.get("source")}')

        if 'kext' in f:
            g = f['kext']
            ka = _attrs(g)
            L.append('  /kext   ' + ', '.join(sorted(g)))
            L.append(f'      M_dust/N_H = {ka.get("M_dust_per_H")} g/H')
            L.append(f'      size_dist_file: {ka.get("size_dist_file")}')
            L.append(f'      text_product:   {ka.get("text_product")}')
        else:
            L.append('  /kext   ABSENT -- run ./calc_kext.x <model> euv')

        if 'polarized' in f:
            p = f['polarized']
            pa = _attrs(p)
            L.append('  /polarized')
            if 'lambda' in p:
                pl = p['lambda']
                L.append(f'      lambda {pl.shape[0]} points, '
                         f'{pl[0]:.5e} .. {pl[-1]:.5e} um')
            L.append(f'      covers_euv = {pa.get("covers_euv")}, '
                     f'pol_valid_from = {pa.get("pol_valid_from")}')
            for sub in ('qjori', 'scatmat_random', 'scatmat_aligned'):
                if sub in p:
                    sets = ', '.join(f'{n}{p[f"{sub}/{n}"].shape}' for n in sorted(p[sub]))
                    L.append(f'      {sub}: {sets}')

    s = '\n'.join(L)
    if show:
        print(s)
    return s


if __name__ == '__main__':
    import sys
    if len(sys.argv) < 2:
        print('usage: python -m pyutil.sedust_h5 <sedust_<model>.h5> [...]')
        raise SystemExit(1)
    for p in sys.argv[1:]:
        describe(p)
        print()
