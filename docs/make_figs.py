#!/usr/bin/env python3
"""Generate the figures used by astrodust_sed_report.tex.

   figs/q_vs_lambda.pdf             Q_ext and Q_abs at four grain sizes
   figs/tau_residual.pdf            our tau_Ad/N_H vs HD23 extinction.dat
   figs/d03_comparison.pdf          HD23 Ad+PAH vs WD01/D03 (2003, 2009)
   figs/sed_stages.pdf              three enthalpy stages vs astrodust_irem
   figs/cabs_per_grain_vs_hd23.pdf  single-grain C_abs vs the release Q table
   figs/sed_total_with_pah.pdf      Ad + PAH vs HD23 model_irem
   figs/sed_total_d16_threeway.pdf  graphite source of the xi blend (see below)
   figs/polarization.pdf            p_max/N_H vs HD23 polarized_extinction

Every path is resolved relative to this file, so the script reads only the
tree it sits in and writes only into that tree's figs/.

Each figure prints the numbers the report quotes for it, so a run of this
script is also the measurement the text is checked against.

sed_total_d16_threeway needs the PAH SED solved with two non-production
graphite sources for the xi blend, which the driver writes as

    ./calc_sed.x astrodust gra_d03_sphere      -> ..._d03gra_PAH.dat
    ./calc_sed.x astrodust gra_d16_spheroid    -> ..._sphdgra_PAH.dat

alongside the production run. That figure is skipped, with a message,
until both files exist.

Matplotlib renders all text through LaTeX on this machine, so labels carry
no non-ASCII character.
"""
from __future__ import annotations
import os
import gzip
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams.update({'text.usetex': True, 'font.family': 'serif'})

# --------- paths
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, '..'))
DATA = os.path.join(ROOT, 'data')
AD   = os.path.join(DATA, 'astrodust')
REL  = os.path.join(DATA, 'release')
SED  = os.path.join(ROOT, 'sed', 'output')
FIGS = os.path.join(HERE, 'figs')

# --------- the DH21 orientation-resolved table shipped by HD23
Q_HD23  = os.path.join(AD, 'q_DH21Ad_P0.20_Fe0.00_1.400.dat.gz')
Q_OURS  = os.path.join(AD, 'q_astrodust_P0.20_Fe0.00_1.400.dat')
NA, NW, HEADER = 169, 1129, 12

UM2CM = 1.0e-4

# The three enthalpy stages of the report: the label used in the text, and
# the driver output that carries it. C2 is the density-corrected Stage-1
# prefactor (`./calc_sed.x astrodust c2`).
STAGES = [('S1_C1', 'sed_astrodust_S1.dat'),
          ('S1_C2', 'sed_astrodust_c2_S1.dat'),
          ('S2',    'sed_astrodust_S2.dat')]


# --------------------------------------------------------------------------
# loaders
# --------------------------------------------------------------------------
def load_wav():
    return np.loadtxt(os.path.join(REL, 'wav.dat'), skiprows=1)


def load_logU():
    return np.loadtxt(os.path.join(REL, 'logU.dat'), skiprows=1)


def load_release_col(name, jcol):
    arr = np.loadtxt(os.path.join(REL, name), skiprows=4)
    return arr[:, jcol]


def load_ours_sed(name):
    out = np.loadtxt(os.path.join(SED, name), comments='#')
    return out[:, 0], out[:, 1]


def jcol_for(logU_target=0.2):
    logU = load_logU()
    return int(np.argmin(np.abs(logU - logU_target)))


def load_size_dist():
    """a [um], dn_Ad/n_H, dn_PAH/n_H, f_ion, f_align on the HD23 size grid."""
    sd = np.loadtxt(os.path.join(REL, 'size_distribution.dat'), skiprows=4)
    return sd[:, 0], sd[:, 1], sd[:, 2], sd[:, 3], sd[:, 4]


def load_q_ours():
    """Our random-orientation Q table, reshaped to (n_lam, n_aeff).

    Columns: lambda[um] a_eff[um] Q_ext Q_abs Q_sca albedo g flag, with
    a_eff varying fastest.
    """
    q = np.loadtxt(Q_OURS)
    lam = np.unique(q[:, 0])
    aef = np.unique(q[:, 1])
    nL, na = len(lam), len(aef)
    if nL * na != q.shape[0]:
        raise ValueError(f'{Q_OURS}: {q.shape[0]} rows is not {nL} x {na}')
    out = {k: q[:, j].reshape(nL, na)
           for k, j in (('ext', 2), ('abs', 3), ('sca', 4), ('g', 6))}
    return lam, aef, out


def _read_free_format(path, n, skip=2):
    with open(path) as f:
        for _ in range(skip):
            f.readline()
        vals = []
        for line in f:
            vals.extend(float(s) for s in line.split())
            if len(vals) >= n:
                break
    return np.asarray(vals[:n])


def load_q_hd23():
    """Orientation-resolved HD23 Q table.

    Returns (aeff, wave, blocks) with blocks[(kq, jori)] of shape (NW, NA);
    kq = 0, 1, 2 for ext, abs, sca and jori = 0, 1, 2 for the file's
    jori = 1 (k || a), 2 (k perp a, E || a), 3 (k perp a, E perp a).
    """
    aeff = _read_free_format(os.path.join(AD, 'DH21_aeff'), NA)
    wave = _read_free_format(os.path.join(AD, 'DH21_wave'), NW)
    blocks = {}
    with gzip.open(Q_HD23, 'rt') as f:
        head = [f.readline() for _ in range(HEADER)]
        if 'Q_ext' not in ''.join(head):
            raise ValueError(f'{Q_HD23}: unexpected header, check HEADER={HEADER}')
        for kq in range(3):
            for jori in range(3):
                rows = [[float(s) for s in f.readline().split()]
                        for _ in range(NW)]
                arr = np.asarray(rows)
                if arr.shape != (NW, NA):
                    raise ValueError(f'{Q_HD23}: block ({kq},{jori}) is {arr.shape}')
                blocks[(kq, jori)] = arr
    return aeff, wave, blocks


def interp_in_a(q_lam, a_in, a_out):
    """Log-linear in a_eff, linear in Q, clamped at the grid edges.

    This is what sed/src/q_table.f90 `interp_q_in_a` does, so the size
    integrals here reproduce the Fortran pipeline's.
    """
    return np.interp(np.log(a_out), np.log(a_in), q_lam)


def loglog_interp(x_in, y_in, x_out):
    return np.exp(np.interp(np.log(x_out), np.log(x_in), np.log(y_in)))


def band_stats(lam, rel, lo, hi, name):
    m = (lam >= lo) & (lam <= hi)
    r = rel[m]
    i = int(np.argmax(np.abs(r)))
    print(f'    {name:16s} median {np.median(r)*100:+7.2f}%   '
          f'mean|rel| {np.mean(np.abs(r))*100:6.2f}%   '
          f'max {np.abs(r[i])*100:6.2f}% at {lam[m][i]:9.2f} um')


# --------------------------------------------------------------------------
# Figure: Q vs lambda at four sizes
# --------------------------------------------------------------------------
def fig_q_vs_lambda():
    lam_all, a_all, Q = load_q_ours()

    sizes_target = [0.001, 0.01, 0.1, 1.0]      # um
    sizes_use = [int(np.argmin(np.abs(a_all - t))) for t in sizes_target]

    fig, axes = plt.subplots(1, 2, figsize=(9, 3.5), sharex=True)
    colors = ['#1f77b4', '#2ca02c', '#d62728', '#7f7f7f']
    for col, key, label in [(0, 'ext', r'$Q_{\rm ext}$'),
                            (1, 'abs', r'$Q_{\rm abs}$')]:
        ax = axes[col]
        for ia, c in zip(sizes_use, colors):
            ax.loglog(lam_all, Q[key][:, ia], color=c, lw=1.0,
                      label=r'$a_{\rm eff}=' + f'{a_all[ia]:.3g}' + r'\,\mu{\rm m}$')
        ax.set_xlabel(r'$\lambda$ [$\mu$m]')
        ax.set_ylabel(label)
        ax.grid(True, which='both', ls=':', alpha=0.5)
        ax.set_xlim(0.09, 4e4)
    axes[0].legend(loc='lower left', fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGS, 'q_vs_lambda.pdf'))
    plt.close(fig)
    print(f'    sizes used [um]: {[float(a_all[i]) for i in sizes_use]}')


# --------------------------------------------------------------------------
# Figure: tau_Ad/N_H residual against the HD23 release
# --------------------------------------------------------------------------
def fig_tau_residual():
    """tau_Ad/N_H = sum_a (dn_Ad/n_H)(a) pi a^2 Q_ext(lambda, a), and the
    same for the scattering channel, against extinction.dat column 2 and
    scattering.dat column 2 -- both astrodust alone, not the Ad+PAH total.

    The size sum is done here rather than read from a product file. The
    shipped kext_astrodust_MW.dat carries the Ad+PAH total, a different
    quantity from the astrodust-only tau the report compares, and the model
    has no astrodust-only extinction product. Column 4 of that file is
    checked against extinction.dat column 4 below, which is the unit check
    on it: both are cm^2 per H, so the ratio must sit at 1.
    """
    lam_q, aef_q, Q = load_q_ours()
    a_sd, dn_ad = load_size_dist()[:2]
    pia2 = np.pi * (a_sd * UM2CM)**2

    tau_e = np.array([np.sum(dn_ad * pia2 * interp_in_a(Q['ext'][i], aef_q, a_sd))
                      for i in range(len(lam_q))])
    tau_s = np.array([np.sum(dn_ad * pia2 * interp_in_a(Q['sca'][i], aef_q, a_sd))
                      for i in range(len(lam_q))])

    ext = np.loadtxt(os.path.join(REL, 'extinction.dat'), skiprows=2)
    sca = np.loadtxt(os.path.join(REL, 'scattering.dat'), skiprows=2)
    lam = ext[:, 0]
    ref_e, ref_s = ext[:, 1], sca[:, 1]
    ours_e = loglog_interp(lam_q, tau_e, lam)
    ours_s = loglog_interp(lam_q, tau_s, lam)
    rel_e = (ours_e - ref_e) / ref_e
    rel_s = (ours_s - ref_s) / ref_s

    ie, iss = int(np.argmax(np.abs(rel_e))), int(np.argmax(np.abs(rel_s)))
    print(f'    tau_Ad  max |rel| = {abs(rel_e[ie])*100:.2f}% at {lam[ie]:.2f} um')
    print(f'    tau_sca max |rel| = {abs(rel_s[iss])*100:.2f}% at {lam[iss]:.2f} um')

    kext = np.loadtxt(os.path.join(AD, 'kext_astrodust_MW.dat'), comments='#')
    r = loglog_interp(kext[:, 0], kext[:, 3], lam) / ext[:, 3]
    print(f'    unit check, kext_astrodust_MW C_ext/H over extinction.dat '
          f'total: min {r.min():.4f} max {r.max():.4f} median {np.median(r):.4f}')

    fig, ax = plt.subplots(2, 1, figsize=(7, 4), sharex=True,
                           gridspec_kw={'height_ratios': [2, 1]})
    ax[0].loglog(lam, ref_e, 'k-', lw=1.2, label=r'HD23 $\tau_{\rm Ad}/N_H$')
    ax[0].loglog(lam, ours_e, color='#d62728', lw=0.8, ls='--', label=r'ours')
    ax[0].loglog(lam, ref_s, color='#1f77b4', lw=0.6, alpha=0.6,
                 label=r'HD23 $\tau^{\rm sca}_{\rm Ad}/N_H$')
    ax[0].set_ylabel(r'$\tau/N_H$ [cm$^2$/H]')
    ax[0].grid(True, which='both', ls=':', alpha=0.4)
    ax[0].legend(loc='lower left', fontsize=8)

    ax[1].semilogx(lam, rel_e * 100.0, 'k-', lw=0.7)
    ax[1].axhline(0, color='gray', lw=0.5)
    ax[1].axhline(+1, color='r', lw=0.5, ls=':')
    ax[1].axhline(-1, color='r', lw=0.5, ls=':')
    ax[1].set_ylabel(r'$(\tau_{\rm ours}-\tau_{\rm HD23})/\tau_{\rm HD23}$ [\%]')
    ax[1].set_xlabel(r'$\lambda$ [$\mu$m]')
    ax[1].grid(True, which='both', ls=':', alpha=0.4)
    ax[1].set_ylim(-1.5, 1.5)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGS, 'tau_residual.pdf'))
    plt.close(fig)


# --------------------------------------------------------------------------
# Figure: SED, three enthalpy stages against HD23 astrodust_irem
# --------------------------------------------------------------------------
def fig_sed_stages():
    wav = load_wav()
    ref = load_release_col('astrodust_irem.dat', jcol_for(0.2))
    ipk = int(np.argmax(ref))

    fig, ax = plt.subplots(2, 1, figsize=(7, 4.5), sharex=True,
                           gridspec_kw={'height_ratios': [2, 1]})
    ax[0].loglog(wav, ref, 'k-', lw=1.2, label=r'HD23 Ad')
    colors = {'S1_C1': '#d62728', 'S1_C2': '#2ca02c', 'S2': '#1f77b4'}
    for tag, fname in STAGES:
        lam, si = load_ours_sed(fname)
        ax[0].loglog(lam, si, color=colors[tag], lw=0.7, ls='--',
                     label=tag.replace('_', '\\_'))
    ax[0].set_ylabel(r'$\lambda I_\lambda^{\rm Ad}/N_H$ [erg\,s$^{-1}$\,sr$^{-1}$\,H$^{-1}$]')
    ax[0].set_xlim(0.5, 4e4)
    ax[0].set_ylim(1e-30, 1e-23)
    ax[0].legend(loc='lower center', fontsize=8, ncol=4)
    ax[0].grid(True, which='both', ls=':', alpha=0.4)

    for tag, fname in STAGES:
        lam, si = load_ours_sed(fname)
        si_on_ref = np.interp(np.log10(wav), np.log10(lam), si)
        rel = (si_on_ref - ref) / np.where(ref > 0, ref, 1.0)
        ax[1].semilogx(wav, rel * 100.0, color=colors[tag], lw=0.7)
        print(f'    {tag}:  peak ({wav[ipk]:.1f} um) {rel[ipk]*100:+.2f}%')
        band_stats(wav, rel, 30.0, 300.0, '30-300 um')
        band_stats(wav, rel, 300.0, 3000.0, '300-3000 um')
        band_stats(wav, rel, 5.0, 30.0, '5-30 um')
    ax[1].axhline(0, color='gray', lw=0.5)
    ax[1].set_xlabel(r'$\lambda$ [$\mu$m]')
    ax[1].set_ylabel(r'$(\rm{ours}-\rm{HD23})/\rm{HD23}$ [\%]')
    ax[1].set_ylim(-50, 80)
    ax[1].set_xlim(0.5, 4e4)
    ax[1].grid(True, which='both', ls=':', alpha=0.4)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGS, 'sed_stages.pdf'))
    plt.close(fig)


# --------------------------------------------------------------------------
# Figure: single-grain C_abs against the HD23 orientation-resolved Q table
# --------------------------------------------------------------------------
def fig_cabs_single_grain():
    """C_abs(lambda) of one astrodust grain, ours against the release table.

    Our Q table is the random-orientation average; the release table's
    jori = 1 block is k || a, a single orientation in which E is always
    perpendicular to the symmetry axis, so it is Q_perp alone. The ratio of
    the two is the flat plateau the report quotes. The random-orientation
    combination (Q^{jori=2} + 2 Q^{jori=3})/3 built from the same release
    file is printed as the control.

    Both tables are on the DH21 grid, so no wavelength interpolation is
    needed; the grids are checked against each other before use.
    """
    aeff, wave, blocks = load_q_hd23()
    lam_q, aef_q, Q = load_q_ours()
    if not (np.allclose(lam_q, wave, rtol=1e-6) and
            np.allclose(aef_q, aeff, rtol=1e-6)):
        raise ValueError('our Q table and the release table are on different grids')

    ja = int(np.argmin(np.abs(aeff - 0.1259)))
    a_um = aeff[ja]
    area = np.pi * (a_um * UM2CM)**2

    Cabs_hd = area * blocks[(1, 0)][:, ja]                     # jori = 1
    Cabs_rand = area * (blocks[(1, 1)][:, ja]
                        + 2.0 * blocks[(1, 2)][:, ja]) / 3.0   # (1/3, 2/3)
    Cabs_ours = area * Q['abs'][:, ja]

    ratio = Cabs_ours / Cabs_hd
    m = (wave >= 30.0) & (wave <= 30000.0)
    print(f'    a_eff = {a_um:.5f} um (DH21 index {ja})')
    print(f'    30-30000 um, ours / jori=1 : mean {ratio[m].mean():.5f}  '
          f'sigma {ratio[m].std():.2e}  -> {(ratio[m].mean()-1)*100:+.2f}%  '
          f'(1/mean = {1.0/ratio[m].mean():.4f})')
    rc = (Cabs_ours / Cabs_rand)[m]
    print(f'    30-30000 um, ours / (1/3 jori2 + 2/3 jori3) : '
          f'mean {rc.mean():.6f}  sigma {rc.std():.2e}')

    fig, ax = plt.subplots(2, 1, figsize=(7, 5), sharex=True,
                           gridspec_kw={'height_ratios': [2, 1]})
    ax[0].loglog(wave, Cabs_hd, 'k-', lw=1.1,
                 label=r'HD23 release $C_{\rm abs}$ (jori=1)')
    ax[0].loglog(wave, Cabs_ours, color='#d62728', lw=0.8, ls='--',
                 label=r'ours $C_{\rm abs}$ (random orientation)')
    ax[0].set_ylabel(r'$C_{\rm abs}$ [cm$^2$]')
    ax[0].set_xlim(0.08, 4e4)
    ax[0].set_ylim(1e-19, 1e-9)
    ax[0].grid(True, which='both', ls=':', alpha=0.4)
    ax[0].set_title(r'Single-grain $C_{\rm abs}$ at $a = ' + f'{a_um:.3f}'
                    + r'\,\mu$m (astrodust)')
    ax[0].legend(loc='lower left', fontsize=8)

    ax[1].semilogx(wave, (ratio - 1.0) * 100.0, 'k-', lw=0.7)
    ax[1].axhline(0, color='gray', lw=0.5)
    ax[1].axhline((ratio[m].mean() - 1.0) * 100.0, color='#d62728',
                  lw=0.5, ls=':')
    ax[1].set_xlabel(r'$\lambda$ [$\mu$m]')
    ax[1].set_ylabel(r'(ours $-$ HD23)/HD23 [\%]')
    ax[1].set_xlim(0.08, 4e4)
    ax[1].set_ylim(-20, 20)
    ax[1].grid(True, which='both', ls=':', alpha=0.4)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGS, 'cabs_per_grain_vs_hd23.pdf'))
    plt.close(fig)


# --------------------------------------------------------------------------
# Figure: total SED (Ad + PAH) against HD23 model_irem
# --------------------------------------------------------------------------
def fig_sed_total_with_pah():
    wav = load_wav()
    j = jcol_for(0.2)
    mod_ref = load_release_col('model_irem.dat', j)
    pah_ref = load_release_col('PAH_irem.dat',   j)

    lam_ad,  si_ad  = load_ours_sed(STAGES[0][1])
    lam_pah, si_pah = load_ours_sed('sed_astrodust_PAH.dat')
    si_tot = si_ad + si_pah

    fig, ax = plt.subplots(2, 1, figsize=(7, 4.5), sharex=True,
                           gridspec_kw={'height_ratios': [2, 1]})
    ax[0].loglog(wav,    mod_ref, 'k-',  lw=1.2, label='HD23 total')
    ax[0].loglog(wav,    pah_ref, color='#7f7f7f', lw=0.6, ls=':', label='HD23 PAH')
    ax[0].loglog(lam_ad, si_tot,  color='#d62728', lw=0.8, ls='--',
                 label=r'ours Ad(S1\_C1) + PAH')
    ax[0].loglog(lam_pah, si_pah, color='#1f77b4', lw=0.6, ls='--',
                 label='ours PAH only')
    ax[0].set_ylabel(r'$\lambda I_\lambda/N_H$ [erg\,s$^{-1}$\,sr$^{-1}$\,H$^{-1}$]')
    ax[0].set_xlim(0.5, 4e4)
    ax[0].set_ylim(1e-30, 1e-23)
    ax[0].legend(loc='lower center', fontsize=8, ncol=2)
    ax[0].grid(True, which='both', ls=':', alpha=0.4)

    si_tot_on_ref = np.interp(np.log10(wav), np.log10(lam_ad), si_tot)
    rel = (si_tot_on_ref - mod_ref) / np.where(mod_ref > 0, mod_ref, 1.0)
    ax[1].semilogx(wav, rel * 100.0, color='#d62728', lw=0.7)
    ax[1].axhline(0, color='gray', lw=0.5)
    ax[1].set_xlabel(r'$\lambda$ [$\mu$m]')
    ax[1].set_ylabel(r'$(\rm ours-HD23)/\rm HD23$ [\%]')
    ax[1].set_xlim(0.5, 4e4)
    ax[1].set_ylim(-100, 100)
    ax[1].grid(True, which='both', ls=':', alpha=0.4)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGS, 'sed_total_with_pah.pdf'))
    plt.close(fig)

    for name, lo, hi in [('1-5 um', 1.0, 5.0), ('5-30 um', 5.0, 30.0),
                         ('30-300 um', 30.0, 300.0),
                         ('300-3000 um', 300.0, 3000.0)]:
        band_stats(wav, rel, lo, hi, name)
    # Band-integrated PAH ratio, the quantity the report's tables quote.
    si_pah_on_ref = np.interp(np.log10(wav), np.log10(lam_pah), si_pah)
    for name, lo, hi in [('NIR 1-5', 1.0, 5.0), ('MIR 5-30', 5.0, 30.0),
                         ('FIR 30-300', 30.0, 300.0),
                         ('sub-mm 300-3000', 300.0, 3000.0)]:
        m = (wav >= lo) & (wav <= hi)
        num = np.trapz(si_pah_on_ref[m] / wav[m], wav[m])
        den = np.trapz(pah_ref[m] / wav[m], wav[m])
        print(f'    PAH band-integrated ours/HD23, {name:16s} {num/den:6.3f}')
    pah_feature_ratios(si_pah_on_ref, pah_ref, wav)


# DL07 Table 1 mode centres and widths, as qpah.f90 carries them. A Drude
# profile written in (lambda/lambda_j - lambda_j/lambda) has half-maximum at
# that variable = +-gamma_j, i.e. FWHM = gamma_j*lambda_j in wavelength.
PAH_LAM_J = [0.0722, 0.2175, 1.050, 1.260, 1.905, 3.300, 5.270, 5.700,
             6.220, 6.690, 7.417, 7.598, 7.850, 8.330, 8.610, 10.68,
             11.23, 11.33, 11.99, 12.62, 12.69, 13.48, 14.19, 15.90,
             16.45, 17.04, 17.375, 17.87, 18.92, 15.0]
PAH_GAMMA_J = [0.195, 0.217, 0.055, 0.11, 0.09, 0.012, 0.034, 0.035,
               0.030, 0.070, 0.126, 0.044, 0.053, 0.052, 0.039, 0.020,
               0.012, 0.032, 0.045, 0.042, 0.013, 0.040, 0.025, 0.020,
               0.014, 0.065, 0.012, 0.016, 0.10, 0.8]
# Each named feature and the DL07 modes that make it up.
PAH_FEATURES = [('3.3 C-H stretch',    [6]),
                ('6.2 C-C stretch',    [9]),
                ('7.7 C-C complex',    [11, 12, 13]),
                ('8.6 C-H in-plane',   [15]),
                ('11.3 C-H oop',       [17, 18]),
                ('12.0-12.7 duo/trio', [19, 20, 21]),
                ('14.19 quartet',      [23]),
                ('17.0 complex',       [26, 27, 28])]


def pah_feature_ratios(ours_on_ref, ref, wav):
    """Band-integrated ours/HD23 over each PAH feature.

    The band of a feature is lambda_j (1 +- gamma_j) -- one Drude FWHM to
    either side of the mode centre -- united over the modes that make the
    feature up, and widened until it holds at least three release-grid
    points where the FWHM is narrower than the grid (only the 3.3 um C-H
    stretch, whose FWHM is 0.040 um against a 0.042 um grid spacing).
    The bands include the underlying continuum; these are band ratios, not
    continuum-subtracted feature strengths.
    """
    lnw = np.log(wav)
    print('    PAH feature bands, ours/HD23 (band = lambda_j +- FWHM_j)')
    for name, modes in PAH_FEATURES:
        lo = min(PAH_LAM_J[m-1] * (1.0 - PAH_GAMMA_J[m-1]) for m in modes)
        hi = max(PAH_LAM_J[m-1] * (1.0 + PAH_GAMMA_J[m-1]) for m in modes)
        c, h = 0.5*(np.log(lo)+np.log(hi)), 0.5*(np.log(hi)-np.log(lo))
        while ((wav >= np.exp(c-h)) & (wav <= np.exp(c+h))).sum() < 3:
            h *= 1.05
        m = (wav >= np.exp(c-h)) & (wav <= np.exp(c+h))
        r = np.trapz(ours_on_ref[m], lnw[m]) / np.trapz(ref[m], lnw[m])
        print(f'      {name:<22}{np.exp(c-h):7.3f}-{np.exp(c+h):<8.3f}'
              f'{m.sum():>3d} pts{r:>9.3f}{(r-1)*100:>+8.1f}%')


# --------------------------------------------------------------------------
# Figure: HD23 astrodust+PAH against WD01/D03, UV through NIR
# --------------------------------------------------------------------------
def fig_d03_comparison():
    """Extinction per H, albedo and asymmetry of the HD23 astrodust+PAH
    model against the older WD01/D03 model, on the D03 wavelength grid.

    Cross sections per H come from the release files extinction.dat and
    scattering.dat (column 4 = Ad+PAH total). g is built from our own
    T-matrix Q table integrated over the HD23 size distribution; PAH g is
    taken as 0 (Rayleigh limit over UV-NIR), so PAHs enter only through
    their small scattering cross section in the denominator.
    """
    def _load_d03(fname):
        # Six numerical columns from row 81 on; some rows carry a trailing
        # comment (e.g. "out-of-plane C-H bend"), so read only those six.
        arr = np.loadtxt(os.path.join(REL, fname), skiprows=80,
                         usecols=(0, 1, 2, 3, 4, 5))
        msk = (arr[:, 0] >= 0.09) & (arr[:, 0] <= 5.0)
        arr = arr[msk]
        return arr[np.argsort(arr[:, 0])]

    d03_03 = _load_d03('kext_albedo_WD_MW_3.1_60_D03.all_2003')
    d03_09 = _load_d03('kext_albedo_WD_MW_3.1_60_D03.all_2009')
    lam_D,  albedo_D,  g_D,  Cext_D  = (d03_03[:, 0], d03_03[:, 1],
                                        d03_03[:, 2], d03_03[:, 3])
    lam_D9, albedo_D9, g_D9, Cext_D9 = (d03_09[:, 0], d03_09[:, 1],
                                        d03_09[:, 2], d03_09[:, 3])

    ext = np.loadtxt(os.path.join(REL, 'extinction.dat'), skiprows=2)
    sca = np.loadtxt(os.path.join(REL, 'scattering.dat'), skiprows=2)
    lam_H  = ext[:, 0]
    Cext_H = ext[:, 3]
    Csca_H = sca[:, 3]
    albedo_H = Csca_H / Cext_H

    lam_q, aef_q, Q = load_q_ours()
    a_sd, dn_ad = load_size_dist()[:2]
    pia2 = np.pi * (a_sd * UM2CM)**2
    Csca_g_ad_q = np.empty(len(lam_q))
    for i in range(len(lam_q)):
        Qs = interp_in_a(Q['sca'][i], aef_q, a_sd)
        gi = interp_in_a(Q['g'][i],   aef_q, a_sd)
        Csca_g_ad_q[i] = np.sum(dn_ad * pia2 * Qs * gi)
    Csca_g_ad_H = np.interp(np.log(lam_H), np.log(lam_q), Csca_g_ad_q)
    g_H = np.where(Csca_H > 0, Csca_g_ad_H / Csca_H, 0.0)

    msk = (lam_H >= 0.09) & (lam_H <= 5.0)
    fig, ax = plt.subplots(1, 3, figsize=(11, 3.3))
    for a in ax:
        a.set_xscale('log')

    ax[0].loglog(lam_D,  Cext_D,  'k-',  lw=1.3, label='WD01 / D03 (2003)')
    ax[0].loglog(lam_D9, Cext_D9, color='0.4', linestyle='-.', lw=1.3,
                 label='WD01 / D03 (2009)')
    ax[0].loglog(lam_H[msk], Cext_H[msk], 'r--', lw=1.3, label='HD23 (Ad+PAH)')
    ax[0].set_xlabel(r'$\lambda~[\mu \mathrm{m}]$')
    ax[0].set_ylabel(r'$C_\mathrm{ext}/N_H~[\mathrm{cm}^2~\mathrm{H}^{-1}]$')
    ax[0].set_title('Extinction')
    ax[0].legend(loc='best', fontsize=7, frameon=False)

    ax[1].semilogx(lam_D,  albedo_D,  'k-',  lw=1.3, label='WD01 / D03 (2003)')
    ax[1].semilogx(lam_D9, albedo_D9, color='0.4', linestyle='-.', lw=1.3,
                   label='WD01 / D03 (2009)')
    ax[1].semilogx(lam_H[msk], albedo_H[msk], 'r--', lw=1.3, label='HD23')
    ax[1].set_xlabel(r'$\lambda~[\mu \mathrm{m}]$')
    ax[1].set_ylabel(r'albedo $= C_\mathrm{sca}/C_\mathrm{ext}$')
    ax[1].set_title('Albedo')
    ax[1].set_ylim(0.0, 0.85)
    ax[1].legend(loc='best', fontsize=7, frameon=False)

    ax[2].semilogx(lam_D,  g_D,  'k-',  lw=1.3, label='WD01 / D03 (2003)')
    ax[2].semilogx(lam_D9, g_D9, color='0.4', linestyle='-.', lw=1.3,
                   label='WD01 / D03 (2009)')
    ax[2].semilogx(lam_H[msk], g_H[msk], 'r--', lw=1.3, label='HD23')
    ax[2].set_xlabel(r'$\lambda~[\mu \mathrm{m}]$')
    ax[2].set_ylabel(r'asymmetry $g = \langle \cos\theta \rangle$')
    ax[2].set_title('Asymmetry parameter')
    ax[2].set_ylim(0.0, 0.85)
    ax[2].legend(loc='best', fontsize=7, frameon=False)

    fig.tight_layout()
    fig.savefig(os.path.join(FIGS, 'd03_comparison.pdf'))
    plt.close(fig)

    iv = int(np.argmin(np.abs(lam_H - 0.55)))
    print(f'    at 0.55 um: HD23 Cext/H {Cext_H[iv]:.4e}, albedo '
          f'{albedo_H[iv]:.4f}, g {g_H[iv]:.4f}')


# --------------------------------------------------------------------------
# Figure: graphite source of the xi blend, three ways
# --------------------------------------------------------------------------
def fig_sed_total_d16_threeway():
    variants = [
        ('D03 sphere',   '#d62728', '--', 'sed_astrodust_d03gra_PAH.dat'),
        ('D16 sphere',   '#1f77b4', '-',  'sed_astrodust_PAH.dat'),
        ('D16 spheroid', '#2ca02c', '-.', 'sed_astrodust_sphdgra_PAH.dat'),
    ]
    missing = [f for _, _, _, f in variants
               if not os.path.exists(os.path.join(SED, f))]
    if missing:
        print('    skipped: not yet run -> ' + ', '.join(missing))
        print('    (./calc_sed.x astrodust gra_d03_sphere and'
              ' ./calc_sed.x astrodust gra_d16_spheroid write them)')
        return

    wav = load_wav()
    j   = jcol_for(0.2)
    mod_ref = load_release_col('model_irem.dat', j)
    pah_ref = load_release_col('PAH_irem.dat',   j)
    lam, si_ad = load_ours_sed(STAGES[0][1])

    fig, ax = plt.subplots(2, 1, figsize=(7, 4.8), sharex=True,
                           gridspec_kw={'height_ratios': [2, 1]})
    ax[0].loglog(wav, mod_ref, 'k-', lw=1.3, label='HD23 total')
    ax[0].loglog(wav, pah_ref, color='#7f7f7f', lw=0.6, ls=':',
                 label='HD23 PAH only')
    for label, color, ls, fname in variants:
        _, pah = load_ours_sed(fname)
        ax[0].loglog(lam, si_ad + pah, color=color, lw=0.8, ls=ls,
                     label=r'ours Ad(S1\_C1)+PAH, ' + label)
    ax[0].set_ylabel(r'$\lambda I_\lambda/N_H$ [erg\,s$^{-1}$\,sr$^{-1}$\,H$^{-1}$]')
    ax[0].set_xlim(0.5, 4e4)
    ax[0].set_ylim(1e-30, 1e-23)
    ax[0].legend(loc='lower center', fontsize=7, ncol=2)
    ax[0].grid(True, which='both', ls=':', alpha=0.4)

    for label, color, ls, fname in variants:
        _, pah = load_ours_sed(fname)
        si_on = np.interp(np.log10(wav), np.log10(lam), si_ad + pah)
        rel   = (si_on - mod_ref) / np.where(mod_ref > 0, mod_ref, 1.0)
        ax[1].semilogx(wav, rel * 100.0, color=color, lw=0.8, ls=ls, label=label)
    ax[1].axhline(0, color='gray', lw=0.5)
    ax[1].set_xlabel(r'$\lambda$ [$\mu$m]')
    ax[1].set_ylabel(r'$(\rm ours-HD23)/\rm HD23$ [\%]')
    ax[1].set_xlim(0.5, 4e4)
    ax[1].set_ylim(-80, 80)
    ax[1].grid(True, which='both', ls=':', alpha=0.4)
    ax[1].legend(loc='lower right', fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGS, 'sed_total_d16_threeway.pdf'))
    plt.close(fig)

    # The PAH-only band-integrated ratios the report tabulates, and the two
    # residuals of the total its caption quotes. Integration is
    # int (lambda I_lambda) dln(lambda), trapezoid on the release grid.
    lnw = np.log(wav)

    def band(y, lo, hi):
        m = (wav >= lo) & (wav <= hi)
        return np.trapz(y[m], lnw[m])

    pah_on_ref = {}
    for label, _, _, fname in variants:
        _, pah = load_ours_sed(fname)
        pah_on_ref[label] = np.interp(np.log10(wav), np.log10(lam), pah)
    labels = [v[0] for v in variants]
    print('    PAH-only band ratio ours/HD23   ' +
          ''.join(f'{l:>15}' for l in labels))
    for name, lo, hi in [('<1 um', 0.0, 1.0), ('1-5 um', 1.0, 5.0),
                         ('5-15 um', 5.0, 15.0), ('15-60 um', 15.0, 60.0),
                         ('60-300 um', 60.0, 300.0),
                         ('300-3000 um', 300.0, 3000.0),
                         ('0.1-3000 um', 0.1, 3000.0)]:
        ref = band(pah_ref, lo, hi)
        print(f'      {name:<28}' +
              ''.join(f'{band(pah_on_ref[l], lo, hi)/ref:>15.4f}'
                      for l in labels))
    si_on = {l: np.interp(np.log10(wav), np.log10(lam), si_ad) + pah_on_ref[l]
             for l in labels}
    m_fir = (wav >= 30.0) & (wav <= 300.0)
    i_fir = int(np.argmax(np.where(m_fir, mod_ref, -np.inf)))
    print(f"    {'total resid. at the %.1f um FIR peak' % wav[i_fir]:<38}" +
          ''.join(f'{(si_on[l][i_fir]/mod_ref[i_fir]-1)*100:>+14.1f}%'
                  for l in labels))
    print(f"    {'total resid., median 3.3-30 um':<38}" +
          ''.join(f'{np.median((si_on[l]/mod_ref-1)[(wav>=3.3)&(wav<=30.0)])*100:>+14.1f}%'
                  for l in labels))


# --------------------------------------------------------------------------
# Figure: polarized extinction
# --------------------------------------------------------------------------
def fig_polarization():
    """p_max/N_H = sum_a dn_Ad(a) f_align(a) C_pol^ext(lambda, a), with
    C_pol^ext = 0.5 pi a^2 (Q_ext^{E perp a} - Q_ext^{E || a}) at k perp a,
    against the release polarized_extinction.dat.

    The sign is the release file's, the same one the polarized branch's own
    C_polext product carries: p_max comes out positive over the optical.
    Only |p_max| is plotted, so the sign matters solely for the residual.
    """
    aeff_q, wave_q, blocks = load_q_hd23()
    Qext_par  = blocks[(0, 1)]      # jori = 2: k perp a, E || a
    Qext_perp = blocks[(0, 2)]      # jori = 3: k perp a, E perp a

    a_dist, dn_ad, _, _, f_align = load_size_dist()

    log_a_in, log_a_out = np.log(aeff_q), np.log(a_dist)
    Qpar_d  = np.empty((NW, len(a_dist)))
    Qperp_d = np.empty_like(Qpar_d)
    for jw in range(NW):
        Qpar_d[jw]  = np.interp(log_a_out, log_a_in, Qext_par[jw])
        Qperp_d[jw] = np.interp(log_a_out, log_a_in, Qext_perp[jw])

    Cpol = 0.5 * np.pi * (a_dist * UM2CM)**2 * (Qperp_d - Qpar_d)
    pmax_ours = (Cpol * (dn_ad * f_align)[None, :]).sum(axis=1)

    ref = np.loadtxt(os.path.join(REL, 'polarized_extinction.dat'), comments='#')
    lam_ref, pmax_ref = ref[:, 0], ref[:, 1]
    pmax_on_ref = np.interp(np.log10(lam_ref), np.log10(wave_q), pmax_ours)

    rel = (pmax_on_ref - pmax_ref) / np.where(np.abs(pmax_ref) > 0,
                                              np.abs(pmax_ref), 1.0)
    ipk_o = int(np.argmax(np.abs(pmax_on_ref)))
    ipk_r = int(np.argmax(np.abs(pmax_ref)))
    print(f'    peak |p_max|/N_H: ours {abs(pmax_on_ref[ipk_o]):.4e} at '
          f'{lam_ref[ipk_o]:.2f} um, HD23 {abs(pmax_ref[ipk_r]):.4e} at '
          f'{lam_ref[ipk_r]:.2f} um')
    # p_max changes sign near 0.107 um, where a relative residual carries no
    # information; the report calls that out as the one outlier. Quote the
    # bands away from the crossing.
    for name, lo, hi in [('0.3-1 um', 0.3, 1.0), ('1-5 um', 1.0, 5.0),
                         ('5-30 um', 5.0, 30.0), ('30-300 um', 30.0, 300.0),
                         ('300-3000 um', 300.0, 3000.0),
                         ('3000-30000 um', 3000.0, 30000.0)]:
        band_stats(lam_ref, rel, lo, hi, name)

    fig, ax = plt.subplots(2, 1, figsize=(7, 4.5), sharex=True,
                           gridspec_kw={'height_ratios': [2, 1]})
    ax[0].loglog(lam_ref, np.abs(pmax_ref), 'k-', lw=1.2, label='HD23')
    ax[0].loglog(lam_ref, np.abs(pmax_on_ref), color='#d62728', lw=0.7,
                 ls='--', label='ours')
    ax[0].set_xlim(0.09, 4e4)
    ax[0].set_ylabel(r'$|p^{\rm max}_\lambda|/N_H$ [cm$^2$/H]')
    ax[0].legend(loc='upper right', fontsize=8)
    ax[0].grid(True, which='both', ls=':', alpha=0.4)

    ax[1].semilogx(lam_ref, rel * 100.0, 'k-', lw=0.7)
    ax[1].axhline(0, color='gray', lw=0.5)
    ax[1].set_ylim(-1, 1)
    ax[1].set_xlim(0.09, 4e4)
    ax[1].set_xlabel(r'$\lambda$ [$\mu$m]')
    ax[1].set_ylabel(r'$(\rm ours-HD23)/|\rm HD23|$ [\%]')
    ax[1].grid(True, which='both', ls=':', alpha=0.4)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGS, 'polarization.pdf'))
    plt.close(fig)


if __name__ == '__main__':
    os.makedirs(FIGS, exist_ok=True)
    print(f'tree: {ROOT}')
    for msg, fn in [('Q vs lambda',             fig_q_vs_lambda),
                    ('tau residual',            fig_tau_residual),
                    ('D03 comparison',          fig_d03_comparison),
                    ('SED stages',              fig_sed_stages),
                    ('single-grain C_abs',      fig_cabs_single_grain),
                    ('SED Ad + PAH',            fig_sed_total_with_pah),
                    ('xi-blend graphite 3-way', fig_sed_total_d16_threeway),
                    ('polarization',            fig_polarization)]:
        print(f'--- {msg}')
        fn()
    print('Done.')
