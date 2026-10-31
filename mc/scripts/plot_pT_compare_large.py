#!/usr/bin/env python3
"""
Single-grain P(T) for large (equilibrium-regime) astrodust grains: the three
MC engines (fixed log grid 0.5-5000 K with 600 bins, adaptive two-pass,
adaptive buffered) side by side, plus the Guhathakurta-Draine calc_P
solution on the fixed T_first grid of main_pT_compare (NT_init log bins from
T_lo to T_hi), whose bins are much wider than these distributions.

Reads pT_compare_large_NN_aXXXXXA.dat written by main_pT_compare.x with the
4-block layout:
  # === GD block            (T_first grid)
  # === MC fixed-grid block (log 0.5-5000 K, 600 bins)
  # === MC 2pass block      (adaptive)
  # === MC buffered block   (adaptive)

Usage: plot_pT_compare_large.py [--data DIR] [--out FILE.pdf]
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
    return np.array(rows) if rows else np.zeros((0, 3))


def load_one(path):
    text = open(path).read()
    blocks = re.split(r"^# ===.*$", text, flags=re.MULTILINE)
    if len(blocks) < 5:
        raise ValueError(f"expected 5 segments (preamble + 4 blocks) in {path}, "
                         f"got {len(blocks)}")
    gd, fixed, twop, buf = (parse_rows(b) for b in blocks[1:5])
    m = re.search(r"_a(\d+)A\.dat$", path)
    aA = int(m.group(1)) if m else 0
    m_teq = re.search(r"Teq\[K\]=\s*([0-9.eE+-]+)", blocks[0])
    Teq = float(m_teq.group(1)) if m_teq else None
    return gd, fixed, twop, buf, aA, Teq


def gd_dpdlnt(T_gd, P_per_bin):
    """GD probability of each T_first bin converted to dP/dlnT."""
    lnT = np.log(T_gd)
    dlnT = np.empty_like(T_gd)
    dlnT[1:-1] = 0.5 * (lnT[2:] - lnT[:-2])
    dlnT[0] = lnT[1] - lnT[0]
    dlnT[-1] = lnT[-1] - lnT[-2]
    return P_per_bin / dlnT


def fmt_size(aA):
    """Radius label from the file's radius in Angstrom."""
    if aA >= 2000:
        return rf"$a = {aA/10000:.2f}\,\mu$m"
    return rf"$a = {aA}$\,\AA"


def plot_set(files, title, outpath, nrows=4, ncols=2):
    fig, axes = plt.subplots(nrows, ncols, figsize=(8.6, 2.6*nrows))
    axes = axes.flatten()
    if len(files) > len(axes):
        raise SystemExit(f"{len(files)} input files for {len(axes)} panels")
    for ax, fpath in zip(axes, files):
        gd, fixed, twop, buf, aA, Teq = load_one(fpath)
        span = [Teq] if Teq else []
        if gd.size > 0:
            y_gd = gd_dpdlnt(gd[:, 0], gd[:, 1])
            pos = y_gd > 0
            ax.plot(gd[pos, 0], y_gd[pos], color="k", lw=1.1, label=r"GD calc\_P")
            span += list(gd[y_gd > 1e-2, 0])
        for arr, col, ls, lbl in [
                (fixed, "C0", "-",  "MC fixed grid"),
                (twop,  "C3", "--", "MC two-pass (adaptive)"),
                (buf,   "C2", ":",  "MC buffered (adaptive)"),
        ]:
            if arr.size == 0:
                continue
            pos = (arr[:, 2] > 0) & np.isfinite(arr[:, 2])
            ax.plot(arr[pos, 0], arr[pos, 2], color=col, ls=ls, lw=1.0, label=lbl,
                    marker="o" if pos.sum() < 3 else None, ms=3)
            span += list(arr[arr[:, 2] > 1e-2, 0])
        if Teq is not None and Teq > 0:
            ax.axvline(Teq, color="0.6", lw=0.6, ls="-.", zorder=0,
                       label=rf"$T_{{\rm eq}} = {Teq:.2f}$ K")
        lo, hi = min(span), max(span)
        pad = 0.05 * (hi - lo)
        ax.set_xlim(lo - pad, hi + pad)
        ax.set_yscale("log")
        ax.set_ylim(1e-2, 1e4)
        ax.set_xlabel(r"$T$ (K)", fontsize=8)
        ax.set_ylabel(r"$dP/d\ln T$", fontsize=8)
        ax.set_title(fmt_size(aA), fontsize=9)
        ax.legend(loc="upper left", fontsize=6, frameon=True)
        ax.grid(True, which="major", ls=":", lw=0.4, alpha=0.5)
        ax.tick_params(labelsize=7)
    for ax in axes[len(files):]:
        ax.set_visible(False)
    fig.suptitle(title, fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.98))
    fig.savefig(outpath)
    print("wrote", outpath)


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    p.add_argument("--data", default=os.path.join(ROOT, "output"))
    p.add_argument("--out", default=os.path.join(ROOT, "output", "pT_compare_large.pdf"))
    args = p.parse_args()
    files = sorted(glob.glob(os.path.join(args.data, "pT_compare_large_[0-9][0-9]_a*A.dat")))
    if not files:
        raise SystemExit("No pT_compare_large_*.dat files found; run "
                         "./main_pT_compare.x run_pT_compare_large.nml first.")
    plot_set(files,
             r"Single-grain $P(T)$ for large astrodust grains: three MC engines and GD "
             r"($N_{\rm evt}=2\times10^5$, $\log U=0.20$)",
             args.out)


if __name__ == "__main__":
    main()
