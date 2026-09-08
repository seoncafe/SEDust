#!/usr/bin/env python3
"""Figures for the G18 Model D scattering asymmetry SEDust computes for itself.

The DustEM distribution ships no G_<gtype>.DAT for the two spheroid
populations of Guillet et al. (2018) Model D, so SEDust computes them with
tmatrix/driver/spheroid_asymmetry_table.f90.  Four checks stand behind those
tables, and this script draws them:

  1  g18d_mie_amCBEx.pdf         our Mie with the BE a-C index against
                                 DustEM's own amCBE SPHERE tables
  2  g18d_mie_suvSil.pdf         our Mie with the WD01 astrosilicate index
                                 against Draine's suvSil_81 sphere table
  3  g18d_spheroid_vs_dustem.pdf our random-orientation T-matrix Q against the
                                 published Q of the two spheroid populations
  4  g18d_sphere_vs_spheroid.pdf how far the volume-equivalent sphere sits
                                 from the spheroid in g, versus x
  5  g18d_gbar.pdf               the size-integrated <cos theta> of g18d that
                                 comes out of the pipeline, beside THEMIS

Inputs are the files those programs write:
  tmatrix/output/mie_vs_amCBEx.dat        ./compare_mie_spheres.x
  tmatrix/output/mie_vs_suvSil_81.dat     ./compare_mie_spheres.x
  tmatrix/output/g18d_<gtype>_regimes.dat ./spheroid_asymmetry_table.x
  data/g18d/kext_g18d.dat                 sed/calc_kext.x
  data/themis/kext_themis.dat             sed/calc_kext.x
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams.update({'text.usetex': True, 'font.family': 'serif', 'font.size': 9})

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
TM   = os.path.join(ROOT, 'tmatrix', 'output')
FIGS = os.path.join(HERE, 'figs')

AMC = 'amCBE_0.3333x'
SIL = 'aSil2001BE6pctG_0.4x'


def nearest_sizes(a_all, wanted):
    """The tabulated radii closest to the ones asked for, without repeats."""
    uniq = np.unique(a_all)
    out = []
    for w in wanted:
        k = uniq[np.argmin(np.abs(np.log(uniq) - np.log(w)))]
        if k not in out:
            out.append(k)
    return out


# ----------------------------------------------------------------------
# 1 and 2: the material tables and our Mie, on spheres
# ----------------------------------------------------------------------
def sphere_check(datafile, outfile, refname, title, sizes, band=None):
    d = np.loadtxt(os.path.join(TM, datafile))
    lam, a, x = d[:, 0], d[:, 1], d[:, 2]
    qa, qa_r, qs, qs_r, g, g_r = (d[:, i] for i in range(3, 9))
    use = d[:, 9] > 0 if d.shape[1] > 9 else np.ones(len(lam), bool)

    fig, ax = plt.subplots(2, 3, figsize=(9.6, 4.8),
                           gridspec_kw={'height_ratios': [2.2, 1]})
    colors = plt.cm.viridis(np.linspace(0, 0.9, len(sizes)))

    for (col, (ours, ref, name, logy)) in enumerate(
            [(qa, qa_r, r'$Q_{\rm abs}$', True),
             (qs, qs_r, r'$Q_{\rm sca}$', True),
             (g,  g_r,  r'$\langle\cos\theta\rangle$', False)]):
        top, bot = ax[0, col], ax[1, col]
        for c, a0 in zip(colors, sizes):
            s = (a == a0) & use
            o = np.argsort(lam[s])
            lab = r'$a = ' + ('%.3g' % a0) + r'\ \mu$m'
            if logy:
                top.loglog(lam[s][o], ref[s][o], '-', color=c, lw=1.3, label=lab)
                top.loglog(lam[s][o], ours[s][o], '--', color='k', lw=0.6)
            else:
                top.semilogx(lam[s][o], ref[s][o], '-', color=c, lw=1.3, label=lab)
                top.semilogx(lam[s][o], ours[s][o], '--', color='k', lw=0.6)
            if logy:
                r = np.where(ref[s][o] > 0, ours[s][o] / np.where(ref[s][o] > 0,
                             ref[s][o], 1.0) - 1.0, np.nan) * 100.0
                bot.semilogx(lam[s][o], r, '-', color=c, lw=0.8)
            else:
                bot.semilogx(lam[s][o], ours[s][o] - ref[s][o], '-', color=c, lw=0.8)
        top.set_ylabel(name)
        top.set_xlabel('')
        bot.axhline(0, color='0.5', lw=0.5)
        bot.set_xlabel(r'$\lambda$ [$\mu$m]')
        bot.set_ylabel(r'ours $-$ ref [\%]' if logy else r'ours $-$ ref')
        if col == 0:
            top.legend(fontsize=6, frameon=False, loc='lower left')
    ax[0, 0].set_title(title + '\n' + r'solid: ' + refname + r',  dashed: this work',
                       fontsize=8, loc='left')
    if band is not None:
        for k in range(3):
            for row in range(2):
                ax[row, k].axvspan(band[0], band[1], color='0.92', zorder=-5)
    fig.tight_layout()
    out = os.path.join(FIGS, outfile)
    fig.savefig(out)
    plt.close(fig)
    print('wrote', out)


# ----------------------------------------------------------------------
# 3: the spheroid populations against the published Q tables
# ----------------------------------------------------------------------
def spheroid_check(sizes_amc, sizes_sil):
    fig, ax = plt.subplots(2, 2, figsize=(8.0, 5.4))
    for row, (gtype, label, sizes) in enumerate(
            [(AMC, r'BE a-C, $b/a = 1/3$', sizes_amc),
             (SIL, r'astrosilicate $+$ 6\% a-C, $b/a = 0.4$', sizes_sil)]):
        d = np.loadtxt(os.path.join(TM, 'g18d_%s_regimes.dat' % gtype))
        lam, a, x, avl = d[:, 0], d[:, 1], d[:, 2], d[:, 3]
        qa, qs = d[:, 4], d[:, 5]
        qa_r, qs_r = d[:, 8], d[:, 9]
        reg = d[:, 10].astype(int)
        ok = (reg == 2) & (qa_r > 0) & (qs_r > 0)
        colors = plt.cm.viridis(np.linspace(0, 0.9, len(sizes)))
        for col, (ours, ref, name) in enumerate(
                [(qa, qa_r, r'$Q_{\rm abs}$'), (qs, qs_r, r'$Q_{\rm sca}$')]):
            b = ax[row, col]
            for c, a0 in zip(colors, sizes):
                s = ok & (a == a0)
                if s.sum() == 0:
                    continue
                o = np.argsort(x[s])
                b.semilogx(x[s][o], (ours[s][o] / ref[s][o] - 1.0) * 100.0, '-',
                           color=c, lw=0.9,
                           label=r'$a_V = ' + ('%.3g' % a0) + r'\ \mu$m')
            b.axhline(0, color='0.5', lw=0.5)
            b.axhspan(-1, 1, color='0.92', zorder=-5)
            b.set_xlabel(r'$x = 2\pi a_V/\lambda$')
            b.set_ylabel(name + r':  ours $-$ DustEM [\%]')
            b.set_ylim(-5, 5)
            if col == 0:
                b.set_title(label, fontsize=8, loc='left')
                b.legend(fontsize=6, frameon=False, ncol=2)
    fig.tight_layout()
    out = os.path.join(FIGS, 'g18d_spheroid_vs_dustem.pdf')
    fig.savefig(out)
    plt.close(fig)
    print('wrote', out)


# ----------------------------------------------------------------------
# 4: how good the volume-equivalent sphere is as a stand-in for g
# ----------------------------------------------------------------------
def sphere_vs_spheroid():
    fig, ax = plt.subplots(1, 2, figsize=(8.0, 3.2))
    for k, (gtype, label) in enumerate(
            [(AMC, r'BE a-C, $b/a = 1/3$'),
             (SIL, r'astrosilicate $+$ 6\% a-C, $b/a = 0.4$')]):
        d = np.loadtxt(os.path.join(TM, 'g18d_%s_regimes.dat' % gtype))
        x, g, gs, reg = d[:, 2], d[:, 6], d[:, 7], d[:, 10].astype(int)
        s = reg == 2
        b = ax[k]
        b.semilogx(x[s], gs[s] - g[s], '.', ms=0.8, color='0.4', alpha=0.4)
        # running maximum of |dg| in log-spaced bins
        edges = np.logspace(np.log10(0.1), np.log10(max(x[s].max(), 0.2)), 40)
        mid = np.sqrt(edges[:-1] * edges[1:])
        hi = []
        for lo, up in zip(edges[:-1], edges[1:]):
            m = s & (x >= lo) & (x < up)
            hi.append(np.abs(gs[m] - g[m]).max() if m.sum() else np.nan)
        b.semilogx(mid, hi, 'r-', lw=1.2, label=r'$\max |\Delta g|$ in bin')
        b.axhline(0, color='0.5', lw=0.5)
        b.set_xlabel(r'$x = 2\pi a_V/\lambda$')
        b.set_ylabel(r'$g_{\rm sphere} - g_{\rm spheroid}$')
        b.set_title(label, fontsize=8, loc='left')
        b.legend(fontsize=7, frameon=False)
    fig.tight_layout()
    out = os.path.join(FIGS, 'g18d_sphere_vs_spheroid.pdf')
    fig.savefig(out)
    plt.close(fig)
    print('wrote', out)


# ----------------------------------------------------------------------
# 5: what the pipeline delivers
# ----------------------------------------------------------------------
def gbar_end_to_end():
    fig, ax = plt.subplots(1, 2, figsize=(8.0, 3.2))
    for path, lab, style in [
            (os.path.join(ROOT, 'data', 'g18d', 'kext_g18d.dat'),
             r'\texttt{g18d} (Guillet et al. 2018 model D)', 'r-'),
            (os.path.join(ROOT, 'data', 'themis', 'kext_themis.dat'),
             r'\texttt{themis} (Jones et al. 2017)', 'b--')]:
        d = np.loadtxt(path)
        lam, alb, gbar = d[:, 0], d[:, 1], d[:, 2]
        ax[0].semilogx(lam, gbar, style, lw=1.1, label=lab)
        ax[1].semilogx(lam, alb, style, lw=1.1, label=lab)
    for b, name in zip(ax, [r'$\langle\cos\theta\rangle$', 'albedo']):
        b.set_xlabel(r'$\lambda$ [$\mu$m]')
        b.set_ylabel(name)
        b.set_xlim(0.09, 30)
        b.set_ylim(0, 1)
        b.legend(fontsize=7, frameon=False, loc='upper right')
    fig.tight_layout()
    out = os.path.join(FIGS, 'g18d_gbar.pdf')
    fig.savefig(out)
    plt.close(fig)
    print('wrote', out)


def main():
    os.makedirs(FIGS, exist_ok=True)

    d = np.loadtxt(os.path.join(TM, 'mie_vs_amCBEx.dat'))
    sizes = nearest_sizes(d[:, 1], [0.003, 0.01, 0.05, 0.2, 1.0])
    sphere_check('mie_vs_amCBEx.dat', 'g18d_mie_amCBEx.pdf',
                 r'DustEM \texttt{Q\_amCBEx} / \texttt{G\_amCBEx}',
                 r'BE amorphous carbon spheres (Zubko et al. 1996)', sizes,
                 band=(800.0, 1.0e5))

    d = np.loadtxt(os.path.join(TM, 'mie_vs_suvSil_81.dat'))
    sizes = nearest_sizes(d[:, 1], [0.003, 0.01, 0.05, 0.2, 1.0])
    sphere_check('mie_vs_suvSil_81.dat', 'g18d_mie_suvSil.pdf',
                 r"Draine \texttt{suvSil\_81}",
                 r'WD01 astrosilicate spheres', sizes)

    d = np.loadtxt(os.path.join(TM, 'g18d_%s_regimes.dat' % AMC))
    s_amc = nearest_sizes(d[:, 1], [0.01, 0.03, 0.1, 0.3, 1.0])
    d = np.loadtxt(os.path.join(TM, 'g18d_%s_regimes.dat' % SIL))
    s_sil = nearest_sizes(d[:, 1], [0.01, 0.03, 0.1, 0.2, 0.29])
    spheroid_check(s_amc, s_sil)
    sphere_vs_spheroid()
    gbar_end_to_end()


if __name__ == '__main__':
    main()
