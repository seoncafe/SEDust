#!/usr/bin/env python3
"""
Compare the two MC SED modes with the production astrodust solver:
  - default (use_mc_for_large_grain=F): calc_Teq for grains that pass the
    150 eV equilibrium gate, MC for the others
  - all-MC  (use_mc_for_large_grain=T): MC for every grain
  - reference: calc_sed.x astrodust (default solver 'heuristic': the
    Guhathakurta-Draine matrix solution with look-ahead temperature-window
    narrowing, and calc_Teq for equilibrium grains)
All three use the Mathis ISRF (Draine-corrected) at U = 1.585, Stage 2.

Usage: plot_sed_compare.py [--mc DIR] [--ref DIR] [--out FILE.pdf]
"""
import argparse
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["text.usetex"] = True
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    p.add_argument("--mc", default=os.path.join(ROOT, "output"))
    p.add_argument("--ref", default=os.path.join(ROOT, "..", "sed", "output"))
    p.add_argument("--out", default=os.path.join(ROOT, "output", "sed_mc_vs_ref.pdf"))
    args = p.parse_args()

    s2 = np.loadtxt(os.path.join(args.ref, "sed_astrodust_S2.dat"), comments="#")
    pah = np.loadtxt(os.path.join(args.ref, "sed_astrodust_PAH.dat"), comments="#")
    lam_ref = s2[:, 0]
    li_ref_tot = s2[:, 1] + np.interp(lam_ref, pah[:, 0], pah[:, 1])

    dflt = np.loadtxt(os.path.join(args.mc, "sed_default_S2_U1.585_irem_mc.dat"), comments="#")
    mall = np.loadtxt(os.path.join(args.mc, "sed_mc_all_S2_U1.585_irem_mc.dat"), comments="#")
    lam = dflt[:, 0]
    if not np.allclose(mall[:, 0], lam):
        raise SystemExit("default and all-MC spectra use different wavelength grids")
    li_ref_g = np.interp(lam, lam_ref, li_ref_tot)
    x_limits = (0.5, 5.0e3)
    ref_visible = (lam_ref >= x_limits[0]) & (lam_ref <= x_limits[1])
    mc_visible = (lam >= x_limits[0]) & (lam <= x_limits[1])
    flux = np.concatenate((li_ref_tot[ref_visible], dflt[mc_visible, 1],
                           mall[mc_visible, 1]))
    flux = flux[np.isfinite(flux) & (flux > 0)]
    if not len(flux):
        raise SystemExit("no positive finite SED values in the plotted wavelength range")
    flux_limits = (0.75 * np.min(flux), 1.35 * np.max(flux))

    ratio_visible = mc_visible & np.isfinite(li_ref_g) & (li_ref_g > 0)
    ratio = np.concatenate((dflt[ratio_visible, 1] / li_ref_g[ratio_visible],
                            mall[ratio_visible, 1] / li_ref_g[ratio_visible]))
    ratio = ratio[np.isfinite(ratio)]
    if not len(ratio):
        raise SystemExit("no finite MC/reference ratios in the plotted wavelength range")
    ratio_padding = max(0.05, 0.1 * (np.max(ratio) - np.min(ratio)))
    ratio_limits = (max(0.0, np.min(ratio) - ratio_padding),
                    np.max(ratio) + ratio_padding)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7.5, 7.0),
                                   gridspec_kw={"height_ratios": [3, 1]},
                                   sharex=True)
    ax1.plot(lam_ref, li_ref_tot, color="k", lw=1.2,
             label=r"Reference: \texttt{calc\_sed.x astrodust} (GD matrix)")
    ax1.plot(lam, dflt[:, 1], color="C0", lw=0.9, ls="--",
             label=r"Default: $T_{\rm eq}$ (equilibrium gate) + MC")
    ax1.plot(lam, mall[:, 1], color="C3", lw=0.9, ls=":",
             label=r"All-MC: \texttt{use\_mc\_for\_large\_grain=T}")
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.set_xlim(*x_limits)
    ax1.set_ylim(*flux_limits)
    ax1.set_ylabel(r"$\lambda I_\lambda / N_{\rm H}$ (erg s$^{-1}$ sr$^{-1}$ H$^{-1}$)",
                   fontsize=9)
    ax1.legend(loc="lower center", fontsize=8, frameon=True)
    ax1.set_title(r"MC SED modes vs.\ production astrodust solver, Mathis ISRF "
                  r"$\log U=0.20$, Stage 2", fontsize=9)
    ax1.tick_params(labelsize=8)

    pos = li_ref_g > 0
    ax2.plot(lam[pos], dflt[pos, 1] / li_ref_g[pos], color="C0", lw=0.9, ls="--",
             label="Default / ref")
    ax2.plot(lam[pos], mall[pos, 1] / li_ref_g[pos], color="C3", lw=0.9, ls=":",
             label="All-MC / ref")
    ax2.axhline(1.0, color="0.6", lw=0.5)
    for h in (0.97, 1.03):
        ax2.axhline(h, color="0.7", lw=0.4, ls=":")
    ax2.set_xscale("log")
    ax2.set_xlim(*x_limits)
    ax2.set_ylim(*ratio_limits)
    ax2.set_ylabel("MC / ref", fontsize=9)
    ax2.set_xlabel(r"$\lambda$ ($\mu$m)", fontsize=9)
    ax2.legend(loc="lower right", fontsize=8, frameon=True)
    ax2.tick_params(labelsize=8)

    fig.tight_layout()
    fig.savefig(args.out)
    print("wrote", args.out)


if __name__ == "__main__":
    main()
