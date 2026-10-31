#!/usr/bin/env python3
"""
Plot Monte Carlo dP/dlnT (Fig 24.6-style) for graphite grains of various
radii heated by the Mathis ISRF (U=1).

Reads <data>/fig24_6/<aLabel>_U1_PT.dat (T_lo, T_hi, T_mid, dP/dT, dP/dlnT).

Usage: plot_fig24_6.py [--data DIR] [--out FILE.pdf]
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

SIZES = [
    ("200A", "0.02"),
    ("100A", "0.01"),
    ("50A",  "0.005"),
    ("20A",  "0.002"),
    ("10A",  "0.001"),
]
COLORS = ["C0", "C1", "C2", "C3", "C4"]


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    p.add_argument("--data", default=os.path.join(ROOT, "output"))
    p.add_argument("--out", default=os.path.join(ROOT, "output", "fig24_6.pdf"))
    args = p.parse_args()

    fig, ax = plt.subplots(figsize=(6.0, 5.0))
    highest_temperature = []
    smallest_density = []
    largest_density = []
    for (lbl, atxt), col in zip(SIZES, COLORS):
        path = os.path.join(args.data, "fig24_6", f"{lbl}_U1_PT.dat")
        values = np.loadtxt(path, comments="#", ndmin=2)
        if values.shape[1] != 5 or not np.isfinite(values).all():
            raise ValueError(f"Invalid temperature probability table: {path}")
        Tlo, Thi, Tmid, dPdT, dPdlnT = values.T
        if np.any(Tlo <= 0) or np.any(Thi <= Tlo) or np.any(dPdT < 0) or np.any(dPdlnT < 0):
            raise ValueError(f"Invalid temperature bins or probability: {path}")
        sel = dPdlnT > 0
        if not np.any(sel):
            raise ValueError(f"No positive temperature probability in {path}")
        highest_temperature.append(float(np.max(Thi[sel])))
        smallest_density.append(float(np.min(dPdlnT[sel])))
        largest_density.append(float(np.max(dPdlnT[sel])))
        # Empty bins remain gaps in the plotted distribution.
        ax.plot(Tmid, np.where(sel, dPdlnT, np.nan), color=col, lw=1.2,
                label=rf"$a = {atxt}\,\mu$m")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(0.5, 1.05 * max(highest_temperature))
    ax.set_ylim(0.8 * min(smallest_density), 10.0 * max(largest_density))
    ax.set_xlabel(r"$T$ (K)", fontsize=10)
    ax.set_ylabel(r"$dP/d\ln T$", fontsize=10)
    ax.set_title(r"MC $P(T)$ for graphite, Mathis ISRF ($U=1$)", fontsize=10)
    ax.legend(loc="upper right", fontsize=8, frameon=True)
    ax.grid(True, which="both", ls=":", lw=0.4, alpha=0.5)
    fig.tight_layout()
    fig.savefig(args.out)
    print("wrote", args.out)


if __name__ == "__main__":
    main()
