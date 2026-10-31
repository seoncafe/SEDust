#!/usr/bin/env python3
"""
Single-grain P(T) comparison: adaptive-grid Monte Carlo
(mc_run_engine_2pass) vs the Guhathakurta-Draine matrix solver (calc_P).
Reads pT_compare_<pop>_NN_aXXXXA.dat files written by main_pT_compare.x
and produces one four-panel figure for each population present.

The data file has four blocks separated by '# ===' lines:
  Block 1 (GD)          : T [K]   P(T_i)   dP/dlnT
  Block 2 (MC fixed)    : T_mid   dP/dT    dP/dlnT   (log grid 0.5-5000 K, 600 bins)
  Block 3 (MC 2pass)    : T_mid   dP/dT    dP/dlnT   <-- plotted here
  Block 4 (MC buffered) : T_mid   dP/dT    dP/dlnT

The two-pass block uses a grid fitted to the recorded trajectory, with a
0.5 K lower limit (the engine's numerical temperature floor).  GD uses the
fixed T_first grid of main_pT_compare (NT_init bins from T_lo to T_hi).
For pop='pah' main_pT_compare uses the f_ion-blended carbonaceous cross
section (Cabs_pah), not the separate neutral and cation grains of the SED.

Usage: plot_pT_compare.py [--data DIR] [--outdir DIR]
"""
import argparse
import glob
import os
import re
import numpy as np
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["text.usetex"] = True
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)


def parse_rows(block):
    rows = [[float(x) for x in line.split()] for line in block.splitlines()
            if line.strip() and not line.strip().startswith("#")]
    return np.array(rows)


def load_one(path):
    """Return (T_gd, P_gd_per_bin, T_mc, dPdlnT_mc, label); MC from block 3."""
    text = open(path).read()
    blocks = re.split(r"^# ===.*$", text, flags=re.MULTILINE)
    if len(blocks) < 5:
        raise ValueError(f"expected 5 segments (preamble + 4 blocks) in {path}, "
                         f"got {len(blocks)}")
    gd = parse_rows(blocks[1])
    mc_2p = parse_rows(blocks[3])
    m = re.search(r"_a(\d+)A\.dat$", path)
    label = rf"$a = {int(m.group(1))}$\,\AA" if m else os.path.basename(path)
    return gd[:, 0], gd[:, 1], mc_2p[:, 0], mc_2p[:, 2], label


def gd_dpdlnt(T_gd, P_gd_bin):
    """GD probability of each T_first bin converted to dP/dlnT."""
    lnT = np.log(T_gd)
    dlnT = np.empty_like(T_gd)
    dlnT[1:-1] = 0.5 * (lnT[2:] - lnT[:-2])
    dlnT[0] = lnT[1] - lnT[0]
    dlnT[-1] = lnT[-1] - lnT[-2]
    return P_gd_bin / dlnT


def plot_pop(files, title, outpath, ymin=1e-6, ymax=30.0):
    fig, axes = plt.subplots(2, 2, figsize=(8.5, 6.5))
    axes = axes.flatten()
    for ax, fpath in zip(axes, files):
        T_gd, P_gd_bin, T_mc, dPdlnT_mc, lbl = load_one(fpath)
        dPdlnT_gd = gd_dpdlnt(T_gd, P_gd_bin)
        gd_pos = dPdlnT_gd > 0
        mc_pos = dPdlnT_mc > 0
        ax.plot(T_gd[gd_pos], dPdlnT_gd[gd_pos], color="k", lw=1.4,
                label=r"GD (calc\_P)")
        ax.plot(T_mc[mc_pos], dPdlnT_mc[mc_pos], color="C3", lw=0.9, ls="--",
                label="MC (two-pass, adaptive)")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(0.4, 4.0e3)
        ax.set_ylim(ymin, ymax)
        ax.set_xlabel(r"$T$ (K)", fontsize=9)
        ax.set_ylabel(r"$dP/d\ln T$", fontsize=9)
        ax.set_title(lbl, fontsize=10)
        ax.legend(loc="upper right", fontsize=8, frameon=True)
        ax.grid(True, which="major", ls=":", lw=0.4, alpha=0.5)
        ax.tick_params(labelsize=8)
    for ax in axes[len(files):]:
        ax.set_visible(False)
    fig.suptitle(title, fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(outpath)
    print("wrote", outpath)


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    p.add_argument("--data", default=os.path.join(ROOT, "output"))
    p.add_argument("--outdir", default=os.path.join(ROOT, "output"))
    args = p.parse_args()
    sets = [
        ("ad", r"Single-grain $P(T)$: astrodust S2, Mathis ISRF ($\log U=0.20$), "
               r"$N_{\rm evt}=5\times10^5$", "pT_compare_astrodust.pdf"),
        ("pah", r"Single-grain $P(T)$: PAH ($f_{\rm ion}$-blended), Mathis ISRF "
                r"($\log U=0.20$), $N_{\rm evt}=5\times10^5$", "pT_compare_pah.pdf"),
    ]
    for tag, title, name in sets:
        files = sorted(glob.glob(os.path.join(args.data, f"pT_compare_{tag}_[0-9][0-9]_a*A.dat")))
        if len(files) > 4:
            raise SystemExit(f"{len(files)} pT_compare_{tag} files in {args.data}; expected at most 4")
        if files:
            plot_pop(files, title, os.path.join(args.outdir, name))


if __name__ == "__main__":
    main()
