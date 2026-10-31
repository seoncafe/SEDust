#!/usr/bin/env python3
"""
Plot the radiation-equivalent MC temperature T_rad = <T^4>^(1/4) of each
equilibrium-gate grain against its calc_Teq temperature, versus grain radius.

Reads <data>/diag_Teq_vs_Trad.dat with columns
  a_eff[um]  Teq[K]  Trad_MC[K]  Trad/Teq  hist_weight_sum  population
written by mc_sed_loop when use_mc_for_large_grain=.true. The three
populations are plotted separately in increasing grain radius.

Usage: plot_Teq_vs_Trad.py [--data DIR] [--out FILE.pdf]
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
POPULATIONS = ("astrodust", "PAH_neutral", "PAH_ionized")
LABELS = ("Astrodust", "Neutral PAH", "Ionized PAH")


def load_diagnostic(path):
    rows = []
    with open(path, encoding="utf-8") as stream:
        for line_number, line in enumerate(stream, 1):
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            fields = line.split()
            if len(fields) != 6 or fields[5] not in POPULATIONS:
                raise ValueError(f"{path}:{line_number}: expected five values and a population")
            values = np.asarray([float(value) for value in fields[:5]])
            if not np.all(np.isfinite(values)) or values[0] <= 0 or values[1] <= 0:
                raise ValueError(f"{path}:{line_number}: invalid radius or temperature")
            rows.append((fields[5], values))
    if not rows:
        raise ValueError(f"{path}: no diagnostic rows")
    return rows


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    p.add_argument("--data", default=os.path.join(ROOT, "output"))
    p.add_argument("--out", default=os.path.join(ROOT, "output", "diag_Teq_vs_Trad.pdf"))
    args = p.parse_args()
    rows = load_diagnostic(os.path.join(args.data, "diag_Teq_vs_Trad.dat"))

    fig, axes = plt.subplots(2, 1, figsize=(7, 6), sharex=True,
                             gridspec_kw={"height_ratios": [2, 1]})
    ax = axes[0]
    for index, (population, label) in enumerate(zip(POPULATIONS, LABELS)):
        subset = np.asarray([values for name, values in rows if name == population])
        if not len(subset):
            continue
        subset = subset[np.argsort(subset[:, 0])]
        a_um, Teq, Trad, ratio, _ = subset.T
        color = f"C{index}"
        ax.semilogx(a_um, Teq, "o", ms=3, mfc="none", color=color,
                    label=f"{label}: " + r"$T_{\rm eq}$")
        ax.semilogx(a_um, Trad, "x", ms=3, color=color,
                    label=f"{label}: " + r"$\langle T^4\rangle^{1/4}$")
        axes[1].semilogx(a_um, ratio, ".", ms=3, color=color, label=label)
    ax.set_ylabel(r"$T$ (K)")
    ax.legend(loc="upper right", fontsize=7, ncol=2)
    ax.grid(True, ls=":", alpha=0.5)
    ax.set_title(r"Equilibrium-gate grains in all-MC mode (two-pass engine)", fontsize=10)

    ax = axes[1]
    ax.axhline(1.0, color="0.5", lw=0.5)
    ax.set_xlabel(r"$a_{\rm eff}$ ($\mu$m)")
    ax.set_ylabel(r"$T_{\rm rad}^{\rm MC} / T_{\rm eq}$")
    ax.grid(True, ls=":", alpha=0.5)
    ax.legend(loc="best", fontsize=7)

    fig.tight_layout()
    fig.savefig(args.out)
    print("wrote", args.out)


if __name__ == "__main__":
    main()
