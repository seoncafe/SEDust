#!/usr/bin/env python3
"""THEMIS: this work against the DustEM products the model is distributed with.

Extinction  -> data/themis/reference/EXT_J13.RES   (the file HD23 cite)
Emission    -> data/themis/reference/SED_J13.RES

Both references are DustEM's own output for GRAIN_J13.DAT, so every input is
shared except the solver: the optics, the calorimetry, the size distribution and
the radiation field are the same files or the same formula on both sides.
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'text.usetex': True, 'font.family': 'serif', 'font.size': 9})

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
REF  = os.path.join(ROOT, 'data', 'themis', 'reference')
FOURPI = 4.0 * np.pi


def read_res(path):
    """DustEM .RES: '# ...' header, then 'ntype nwave', then nwave rows."""
    txt = [l for l in open(path) if not l.lstrip().startswith('#')]
    ntype, nwave = (int(x) for x in txt[0].split()[:2])
    d = np.array([[float(x) for x in l.split()] for l in txt[1:1 + nwave]])
    return ntype, d


def main():
    fig, ax = plt.subplots(2, 2, figsize=(7.2, 5.0),
                           gridspec_kw={'height_ratios': [2.4, 1]})

    # ---- extinction ---------------------------------------------------
    nt, d = read_res(os.path.join(REF, 'EXT_J13.RES'))
    lref, ext_ref = d[:, 0], d[:, -1]
    ours = np.loadtxt(os.path.join(ROOT, 'data', 'themis',
                                   'kext_themis_euv.dat'))
    lam, cext = ours[:, 0], ours[:, 3] * 1e21      # cm^2/H -> the RES scaling
    o = np.interp(np.log(lref), np.log(lam), cext)

    a = ax[0, 0]
    a.loglog(lref, ext_ref, 'k-', lw=1.4, label=r'DustEM \texttt{EXT\_J13.RES}')
    a.loglog(lam, cext, 'r--', lw=0.9, label=r'this work')
    for i in range(nt):
        a.loglog(lref, d[:, 1 + i] + d[:, 1 + nt + i], color='0.6', lw=0.4)
    a.set_ylabel(r'$\tau_{\rm ext}/N_{\rm H}\ \times10^{21}$ [cm$^2$/H]')
    a.set_ylim(1e-8, 5)
    a.set_title(r'Extinction (grey: the four populations)')
    a.legend(fontsize=7, frameon=False, loc='lower left')

    a = ax[1, 0]
    a.semilogx(lref, (o / ext_ref - 1) * 1e6, 'r-', lw=0.8)
    a.axhline(0, color='0.5', lw=0.5)
    a.set_ylabel(r'residual [ppm]')
    a.set_xlabel(r'$\lambda$ [$\mu$m]')
    a.set_ylim(-1, 1)

    # ---- emission -----------------------------------------------------
    nt, d = read_res(os.path.join(REF, 'SED_J13.RES'))
    lref, sed_ref = d[:, 0], d[:, -1]
    ours = np.loadtxt(os.path.join(ROOT, 'sed', 'output',
                                   'sed_themis_morig.dat'), comments='#')
    lam, sed = ours[:, 0], FOURPI * ours[:, 1]
    o = np.interp(np.log(lref), np.log(lam), sed)
    m = (sed_ref > 0) & (lref >= 1) & (lref <= 3000)

    a = ax[0, 1]
    a.loglog(lref[m], sed_ref[m], 'k-', lw=1.4, label=r'DustEM \texttt{SED\_J13.RES}')
    a.loglog(lref[m], o[m], 'r--', lw=0.9, label=r'this work')
    a.set_ylabel(r'$4\pi\nu I_\nu/N_{\rm H}$ [erg s$^{-1}$ H$^{-1}$]')
    a.set_title(r'Emission, original Mathis field, $U=1$')
    a.legend(fontsize=7, frameon=False, loc='lower center')

    a = ax[1, 1]
    a.semilogx(lref[m], (o[m] / sed_ref[m] - 1) * 100, 'r-', lw=0.8)
    a.axhline(0, color='0.5', lw=0.5)
    a.set_ylabel(r'residual [\%]')
    a.set_xlabel(r'$\lambda$ [$\mu$m]')
    a.set_ylim(-10, 5)

    for row in ax:
        for a in row:
            a.grid(True, which='both', ls=':', alpha=0.3)
    fig.tight_layout()
    out = os.path.join(HERE, 'figs', 'dustem_themis_comparison.pdf')
    fig.savefig(out)
    print('wrote', out)



if __name__ == '__main__':
    main()
