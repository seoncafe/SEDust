#!/usr/bin/env python3
"""
Plot Monte Carlo dP/dlnT (Fig 24.6-style) for graphite grains of various
radii heated by Mathis 1983 ISRF (U=1).

Reads output/fig24_6/<aLabel>_U1_PT.dat (T_lo, T_hi, T_mid, dP/dT, dP/dlnT).
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT  = os.path.join(ROOT, "output", "fig24_6")

SIZES = [
    ("200A", 0.020, "0.02"),
    ("100A", 0.010, "0.01"),
    ("50A",  0.005, "0.005"),
    ("20A",  0.002, "0.002"),
    ("10A",  0.001, "0.001"),
]
COLORS = ["C0", "C1", "C2", "C3", "C4"]


def load_PT(path):
    arr = np.loadtxt(path, comments="#")
    Tlo, Thi, Tmid, dPdT, dPdlnT = arr.T
    return Tmid, dPdT, dPdlnT


def main():
    fig, ax = plt.subplots(figsize=(6.0, 5.0))
    for (lbl, aum, atxt), col in zip(SIZES, COLORS):
        path = os.path.join(OUT, f"{lbl}_U1_PT.dat")
        Tmid, dPdT, dPdlnT = load_PT(path)
        sel = dPdlnT > 0
        ax.plot(Tmid[sel], dPdlnT[sel], color=col, lw=1.2,
                label=fr"$a = {atxt}\,\mu$m")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(3.0, 1.5e3)
    ax.set_ylim(1.0e-8, 5.0)
    ax.set_xlabel(r"$T$ (K)", fontsize=10)
    ax.set_ylabel(r"$dP/d\ln T$", fontsize=10)
    ax.set_title("MC P(T) for graphite, Mathis ISRF (U=1)", fontsize=10)
    ax.legend(loc="upper right", fontsize=8, frameon=True)
    ax.grid(True, which="both", ls=":", lw=0.4, alpha=0.5)
    fig.tight_layout()
    out = os.path.join(ROOT, "output", "fig24_6.pdf")
    fig.savefig(out, dpi=140)
    print("wrote", out)


if __name__ == "__main__":
    main()
