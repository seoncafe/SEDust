#!/usr/bin/env python3
"""
Compare the two MC-pT solver modes against the production astrodust pipeline:
  - default (use_mc_for_large_grain=F): calc_Teq for large + MC for small
  - mc_all (use_mc_for_large_grain=T):  MC for all grains
  - reference: astrodust pipeline (calc_Teq + calc_P)
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
SED_REF_DIR = os.path.join(ROOT, "..", "sed", "output")

MC_DEFAULT = os.path.join(ROOT, "output", "sed_default_S2_U1.585_irem_mc.dat")
MC_ALL     = os.path.join(ROOT, "output", "sed_mc_all_S2_U1.585_irem_mc.dat")
REF_S2     = os.path.join(SED_REF_DIR, "astrodust_irem_ours_S2.dat")
REF_PAH    = os.path.join(SED_REF_DIR, "astrodust_irem_ours_PAH.dat")


def load_mc(path):
    arr = np.loadtxt(path, comments="#")
    return arr[:, 0], arr[:, 1], arr[:, 2], arr[:, 3]


def load_ref():
    s2  = np.loadtxt(REF_S2,  comments="#")
    pah = np.loadtxt(REF_PAH, comments="#")
    lam = s2[:, 0]
    return lam, s2[:, 1], np.interp(lam, pah[:, 0], pah[:, 1])


def main():
    lam_ref, li_ref_ad, li_ref_pah = load_ref()
    li_ref_tot = li_ref_ad + li_ref_pah

    lam_d, _, _, _              = load_mc(MC_DEFAULT)
    _,     li_def_tot, li_def_ad, li_def_pah = load_mc(MC_DEFAULT)
    _,     li_mca_tot, li_mca_ad, li_mca_pah = load_mc(MC_ALL)

    # Reference on the MC lambda grid
    li_ref_tot_g = np.interp(lam_d, lam_ref, li_ref_tot)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7.5, 7.0),
                                   gridspec_kw={"height_ratios": [3, 1]},
                                   sharex=True)
    # Top: SEDs
    ax1.plot(lam_ref, li_ref_tot,  color="k",  lw=1.2, label="Reference (calc_Teq + calc_P)")
    ax1.plot(lam_d,   li_def_tot,  color="C0", lw=0.9, ls="--",
             label="Default: Teq(large) + MC(small)")
    ax1.plot(lam_d,   li_mca_tot,  color="C3", lw=0.9, ls=":",
             label=r"$\texttt{use\_mc\_for\_large\_grain=T}$ (MC all)")
    ax1.set_xscale("log"); ax1.set_yscale("log")
    ax1.set_xlim(0.5, 5.0e3); ax1.set_ylim(1e-30, 5e-23)
    ax1.set_ylabel(r"$\lambda I_\lambda / N_H$ (erg s$^{-1}$ cm$^{-2}$ sr$^{-1}$)",
                   fontsize=9)
    ax1.legend(loc="lower center", fontsize=8, frameon=True)
    ax1.set_title("Two MC-pT solver modes vs production astrodust pipeline, $\\log U=0.20$, Stage S2",
                  fontsize=9)
    ax1.tick_params(labelsize=8)

    # Bottom: MC/ref ratio for both modes
    pos = li_ref_tot_g > 0
    r_def = np.full_like(li_def_tot, np.nan)
    r_mca = np.full_like(li_mca_tot, np.nan)
    r_def[pos] = li_def_tot[pos] / li_ref_tot_g[pos]
    r_mca[pos] = li_mca_tot[pos] / li_ref_tot_g[pos]
    ax2.plot(lam_d[pos], r_def[pos], color="C0", lw=0.9, ls="--",
             label="Default / ref")
    ax2.plot(lam_d[pos], r_mca[pos], color="C3", lw=0.9, ls=":",
             label="MC-all / ref")
    ax2.axhline(1.0, color="0.6", lw=0.5)
    for h in [0.97, 1.03]:
        ax2.axhline(h, color="0.7", lw=0.4, ls=":")
    ax2.set_xscale("log"); ax2.set_xlim(0.5, 5.0e3); ax2.set_ylim(0.8, 1.2)
    ax2.set_ylabel("MC / ref", fontsize=9)
    ax2.set_xlabel(r"$\lambda$ ($\mu$m)", fontsize=9)
    ax2.legend(loc="lower right", fontsize=8, frameon=True)
    ax2.tick_params(labelsize=8)

    fig.tight_layout()
    out = os.path.join(ROOT, "output", "sed_mc_vs_ref.pdf")
    fig.savefig(out, dpi=140)
    print("wrote", out)


if __name__ == "__main__":
    main()
