"""Does pyutil.sedust_h5 read the products the way the Fortran does?

Run from the repository root:

    python3 pyutil/tests/test_sedust_h5.py [data_dir]

Three things are checked for every shipped model, against the text products
`calc_qtable.x` and `calc_kext.x` write beside the HDF5 file:

  * the CUT.  include_euv=False must return exactly the wavelengths the narrow
    text product carries -- same count, same values.  That is the claim the
    format rests on: the paired *_euv files are views of one array, so the two
    cannot disagree.  A reader that cut on its own wavelength comparison
    instead of on the file's i_lyman would pass a count check and fail this.
  * the VALUES.  Every cross-section table and the extinction curve must equal
    the text product to the precision that product is written with -- seven
    significant digits for the Q tables, thirteen for the kext curve.
  * the SHAPE.  Cross sections come back (n_a, n_lam): each grain radius one
    contiguous spectrum, wavelength last.
"""

from __future__ import annotations

import os
import sys

import numpy as np
import h5py

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__)))))

from pyutil.sedust_h5 import (read_grid, read_kext, read_qtable, read_polarized,
                              components, model_file, describe)

# model -> (component -> text basename stem), and the kext product stems.
# zubko: the q_zubko_* text products are the D03 Mie recomputation, which the
# file stores as the *_mie_d03 groups; the sil/gra/pah groups hold the ZDA
# tables as distributed, whose text form (suvSil_121_1201.dat, ...) is a
# different format this test does not parse.
QTEXT = {
    'astrodust': {'pah_neu': 'q_astrodust_pah_neu',
                  'pah_ion': 'q_astrodust_pah_ion'},
    'dl07':      {'sil': 'q_dl07_sil', 'gra': 'q_dl07_gra',
                  'pah_neu': 'q_dl07_pah_neu', 'pah_ion': 'q_dl07_pah_ion'},
    'zubko':     {'sil_mie_d03': 'q_zubko_sil', 'gra_mie_d03': 'q_zubko_gra',
                  'pah_mie_d03': 'q_zubko_pah'},
}
# The full component set each file carries (astrodust's Q table itself lives
# in /qtable/astrodust with no text counterpart on this grid).
COMPONENTS = {
    'astrodust': ['astrodust', 'pah_ion', 'pah_neu'],
    'dl07':      ['gra', 'pah_ion', 'pah_neu', 'sil'],
    'zubko':     ['gra', 'gra_mie_d03', 'pah', 'pah_mie_d03',
                  'sil', 'sil_mie_d03'],
}
KTEXT = {
    'astrodust': ('kext_astrodust_MW_euv.dat', 'kext_astrodust_MW.dat'),
    'dl07':      ('kext_dl07_MW_euv.dat',      'kext_dl07_MW.dat'),
    'zubko':     ('kext_zubko_BARE_GR_S_euv.dat', 'kext_zubko_BARE_GR_S.dat'),
}
# Every product a model owns is in data/<model>/, so both the Q tables and the
# extinction curves are found there, beside the HDF5 file they were written
# with.  There is no shared qtable/ directory any more.
# The astrodust Q table is written on one wavelength set only; its narrow
# counterpart is the HDF5 cut and has no text file of its own.
QTEXT_HAS_NARROW = {'astrodust': False, 'dl07': True, 'zubko': True}

TOL_Q = 5.0e-7      # seven written digits
TOL_K = 5.0e-12     # thirteen written digits
TOL_LAM = 1.0e-12   # the axis is the same numbers, not a re-derivation

nbad = 0


def check(cond: bool, what: str, detail: str = '') -> None:
    global nbad
    mark = 'ok  ' if cond else '*** '
    print(f'    {mark}{what}{("  " + detail) if detail else ""}')
    if not cond:
        nbad += 1


def read_q_text(path: str):
    """lambda-major rows -> (lam, a_eff, columns as (n_a, n_lam))."""
    d = np.loadtxt(path, comments='#')
    na = len(np.unique(d[:, 1]))
    nlam = d.shape[0] // na
    lam = d[::na, 0]
    aeff = d[:na, 1]
    cols = [d[:, k].reshape(nlam, na).T for k in range(2, d.shape[1])]
    return lam, aeff, cols


def maxrel(a, b):
    den = np.where(np.abs(b) > 0, np.abs(b), 1.0)
    return float((np.abs(a - b) / den).max())


def main(data_dir: str) -> int:
    for model in ('astrodust', 'dl07', 'zubko'):
        path = model_file(data_dir, model)
        print(f'\n=== {path}')
        if not os.path.exists(path):
            check(False, 'file exists')
            continue

        wide = read_grid(path, include_euv=True)
        narrow = read_grid(path, include_euv=False)
        print(f'  grid: {wide.n} wide, {narrow.n} narrow, i_lyman = {wide.i_lyman}')
        check(narrow.n == wide.n - wide.i_lyman + 1, 'narrow count = n - i_lyman + 1')
        check(np.array_equal(narrow.lam, wide.lam[wide.i_lyman - 1:]),
              'narrow axis is the wide axis sliced at i_lyman')
        # The writer's covering rule (lyman_index in sedust_product.f90):
        # i_lyman is the LAST node at or below the limit, so the narrow view
        # starts at or just below 0.0912 um and the next node is above it.
        lim = 0.0912 * (1.0 + 1.0e-9)
        check(narrow.lam[0] <= lim < narrow.lam[1],
              'i_lyman is the last node at or below the Lyman limit',
              f'{narrow.lam[0]:.6e} <= 0.0912 < {narrow.lam[1]:.6e}')

        # ---- Q tables vs the text products --------------------------------
        for comp, stem in QTEXT[model].items():
            for euv, suffix in ((True, '_euv'), (False, '')):
                if not euv and not QTEXT_HAS_NARROW[model]:
                    continue
                tp = os.path.join(data_dir, model, stem + suffix + '.dat')
                if not os.path.exists(tp):
                    # astrodust's PAH tables carry no suffix at all.
                    tp = os.path.join(data_dir, model, stem + '.dat')
                    if not os.path.exists(tp) or not euv:
                        continue
                q = read_qtable(path, comp, include_euv=euv)
                tl, ta, tc = read_q_text(tp)
                tag = f'{comp}{suffix or " (narrow)"}'
                if len(tl) != len(q.lam) or len(ta) != len(q.a_eff):
                    check(False, f'{tag}: grid vs {os.path.basename(tp)}',
                          f'{len(q.lam)}x{len(q.a_eff)} vs {len(tl)}x{len(ta)}')
                    continue
                check(q.Q_abs.shape == (len(q.a_eff), len(q.lam)),
                      f'{tag}: shape is (n_a, n_lam)', str(q.Q_abs.shape))
                worst = max(maxrel(q.lam, tl), maxrel(q.a_eff, ta))
                # 3-column text is Q_abs alone; 7-column is
                # Q_ext, Q_abs, Q_sca, albedo, g.
                pairs = [(q.Q_abs, tc[0])] if not q.scatters else \
                        [(q.Q_ext, tc[0]), (q.Q_abs, tc[1]),
                         (q.Q_sca, tc[2]), (q.g, tc[4])]
                for h, t in pairs:
                    worst = max(worst, maxrel(h, t))
                check(worst < TOL_Q, f'{tag}: vs {os.path.basename(tp)}',
                      f'max rel {worst:.1e}')

        # ---- kext vs the text products ------------------------------------
        wide_txt, narrow_txt = KTEXT[model]
        for euv, name in ((True, wide_txt), (False, narrow_txt)):
            tp = os.path.join(data_dir, model, name)
            if not os.path.exists(tp):
                continue
            k = read_kext(path, include_euv=euv)
            t = np.loadtxt(tp, comments='#')
            if len(t) != len(k.lam):
                check(False, f'kext{"" if euv else " (narrow)"}: row count',
                      f'{len(k.lam)} vs {len(t)}')
                continue
            # kext_astrodust_MW.dat is the frozen legacy format: its
            # cross-section columns carry 8 significant digits but albedo and
            # <cos> are f10.6, so where those are of order 1e-7 the file holds
            # 0.000000 and a RELATIVE comparison against them is meaningless.
            # Those two columns are therefore checked absolutely, to half a
            # unit in the last decimal the format writes.
            legacy = (not euv) and model == 'astrodust'
            tol = 5.0e-8 if legacy else TOL_K
            tag = f'kext{"" if euv else " (narrow)"}: vs {name}'
            worst = maxrel(k.lam, t[:, 0])
            for h, col in ((k.C_ext, 3), (k.C_abs, 4), (k.C_sca, 5)):
                worst = max(worst, maxrel(h, t[:, col]))
            if k.K_abs is not None and t.shape[1] > 6:
                worst = max(worst, maxrel(k.K_abs, t[:, 6]))
            check(worst < tol, tag + '  [lambda, C_*]', f'max rel {worst:.1e}')
            if legacy:
                dabs = max(float(np.abs(k.albedo - t[:, 1]).max()),
                           float(np.abs(k.g - t[:, 2]).max()))
                check(dabs <= 5.0e-7, tag + '  [albedo, <cos>]',
                      f'max abs {dabs:.1e} (the format writes 6 decimals)')
            else:
                w2 = max(maxrel(k.albedo, t[:, 1]), maxrel(k.g, t[:, 2]))
                check(w2 < tol, tag + '  [albedo, <cos>]', f'max rel {w2:.1e}')
            if euv and k.M_dust_per_H:
                check(abs(k.C_abs / k.M_dust_per_H - k.K_abs).max() <
                      1e-12 * np.abs(k.K_abs).max(),
                      'K_abs = C_abs / M_dust_per_H')

        # ---- what the file says it holds ----------------------------------
        check(sorted(components(path)) == COMPONENTS[model],
              'components as expected', str(sorted(components(path))))

    check_polarized(data_dir)

    print()
    if nbad == 0:
        print('ALL CHECKS PASSED')
        return 0
    print(f'{nbad} CHECK(S) FAILED')
    return 1


def check_polarized(data_dir: str) -> None:
    """read_polarized on a scalar product, and on a synthetic polarized one.

    No v1.20 file exists yet -- those groups are written in step 5 -- so the
    second half writes the layout of section 2 of the migration note with h5py
    and reads it back.  That exercises the reader, and says nothing about
    whether the writer will agree with it; when the writer lands, this test
    should be pointed at a real file instead.
    """
    import tempfile
    print('\n=== /polarized')
    # The astrodust product carries /polarized in the polarized branch and not
    # in the scalar one, and the same file is copied between the two trees, so
    # this test has to accept either and check what it finds.
    p = model_file(data_dir, 'astrodust')
    if os.path.exists(p):
        try:
            pol = read_polarized(p)
        except KeyError:
            pol = None
            check(True, 'scalar branch: no /polarized, read_polarized raises KeyError')
        if pol is not None:
            check(len(pol.lam) > 0, 'polarized branch: /polarized/lambda present',
                  f'{len(pol.lam)} points, {pol.lam[0]:.4e} .. {pol.lam[-1]:.4e} um')
            check(pol.pol_valid_from is not None and
                  abs(pol.pol_valid_from - pol.lam[0]) < 1e-12,
                  'pol_valid_from is where the axis starts',
                  f'{pol.pol_valid_from}')
            check(pol.covers_euv == (pol.lam[0] < 0.0912),
                  'covers_euv agrees with the axis', str(pol.covers_euv))
            q = pol.qjori
            check('Q_ext' in q and q['Q_ext'].shape[0] == 3,
                  'qjori carries all three orientations',
                  str(q.get('Q_ext', np.empty(0)).shape))
            check(bool(np.all(np.isfinite(q['Q_ext']))), 'qjori is finite')
            for sub, key in (('scatmat_random', 'F_tot'),
                             ('scatmat_aligned', 'Z')):
                d = getattr(pol, sub)
                check(key in d, f'{sub} carries {key}',
                      str(d.get(key, np.empty(0)).shape))

    n_lam, n_a = 1129, 169
    with tempfile.TemporaryDirectory() as td:
        fp = os.path.join(td, 'sedust_astrodust.h5')
        with h5py.File(fp, 'w') as f:
            g = f.create_group('polarized')
            g.attrs['covers_euv'] = 0
            g.attrs['pol_valid_from'] = 0.0912
            g.create_dataset('lambda', data=np.logspace(-1, 4, n_lam))
            q = g.create_group('qjori')
            q.attrs['jori_convention'] = np.bytes_(b'1 = k parallel to a')
            q.create_dataset('a_eff', data=np.logspace(-3, 0, n_a))
            for name in ('Q_ext', 'Q_abs', 'Q_sca', 'Q_re'):
                q.create_dataset(name, data=np.zeros((3, n_a, n_lam)))
            r = g.create_group('scatmat_random')
            r.create_dataset('theta', data=np.linspace(0, 180, 37))
            a = g.create_group('scatmat_aligned')
            a.attrs['profile_name'] = np.bytes_(b'test')
            for name, shape in (('theta_i', 7), ('theta_s', 37), ('phi', 13)):
                a.create_dataset(name, data=np.zeros(shape))
        pol = read_polarized(fp)
        check(len(pol.lam) == n_lam, 'polarized: its own lambda axis',
              f'{len(pol.lam)} points')
        check(pol.covers_euv is False and pol.pol_valid_from == 0.0912,
              'polarized: covers_euv and pol_valid_from',
              f'{pol.covers_euv}, {pol.pol_valid_from}')
        check(pol.qjori['Q_ext'].shape == (3, n_a, n_lam),
              'polarized: qjori is (3, n_a, n_lam)',
              str(pol.qjori['Q_ext'].shape))
        check(pol.qjori['_attrs'].get('jori_convention') == '1 = k parallel to a',
              'polarized: group attributes come back decoded')
        check(len(pol.scatmat_aligned['theta_i']) == 7 and
              len(pol.scatmat_random['theta']) == 37,
              'polarized: both scattering-matrix groups read')


if __name__ == '__main__':
    ddir = sys.argv[1] if len(sys.argv) > 1 else 'data'
    raise SystemExit(main(ddir))
