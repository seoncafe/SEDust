#!/usr/bin/env python3
"""
Single-grain P(T) for *large* (equilibrium-regime) grains: shows the three
MC engines (fixed-grid log[2,5000K], adaptive 2pass, adaptive buffered)
side by side, plus the Guhathakurta-Draine calc_P attempt on the wide
T_first grid (which collapses to ~ delta at Teq, or P==0 if the grid is
too wide to resolve).

Reads pT_compare_*_NN_aXXXXA.dat written by main_pT_compare.x with the
4-block layout:
  # === GD block  (T_first grid)
  # === MC fixed-grid block       (log 2-5000 K, NHIST=600)
  # === MC 2pass block            (adaptive)
  # === MC buffered block         (adaptive)
"""
import os, glob, re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)


def parse_rows(block):
    rows = []
    for line in block.splitlines():
        line = line.strip()
        if (not line) or line.startswith("#"):
            continue
        rows.append([float(x) for x in line.split()])
    return np.array(rows) if rows else np.zeros((0, 3))


def load_one(path):
    text = open(path).read()
    blocks = re.split(r"^# ===.*$", text, flags=re.MULTILINE)
    if len(blocks) < 5:
        raise ValueError(f"expected 5 segments (preamble + 4 blocks) in {path}, "
                         f"got {len(blocks)}")
    gd    = parse_rows(blocks[1])
    fixed = parse_rows(blocks[2])
    twop  = parse_rows(blocks[3])
    buf   = parse_rows(blocks[4])
    m = re.search(r"_a(\d+)A\.dat$", path)
    aA = int(m.group(1)) if m else 0
    # Parse Teq from the preamble (first header block)
    Teq = None
    m_teq = re.search(r"Teq\[K\]=\s*([0-9.eE+-]+)", blocks[0])
    if m_teq:
        Teq = float(m_teq.group(1))
    return gd, fixed, twop, buf, aA, Teq


def gd_dpdlnt(T_gd, P_per_bin):
    """Convert GD per-bin probability into dP/dlnT on T_first."""
    lnT = np.log(np.maximum(T_gd, 1e-30))
    dlnT = np.empty_like(T_gd)
    dlnT[1:-1] = 0.5 * (lnT[2:] - lnT[:-2])
    dlnT[0]    = lnT[1] - lnT[0]
    dlnT[-1]   = lnT[-1] - lnT[-2]
    return P_per_bin / np.maximum(dlnT, 1e-30)


def fmt_size(aA):
    """Format aA (in Angstroms) as a human-friendly label."""
    if aA >= 10000:
        return fr"$a = {aA/10000:.1f}\,\mu{{\rm m}}$"
    if aA >= 1000:
        return fr"$a = {aA/10000:.2f}\,\mu{{\rm m}}\ ({aA}\,\mathrm{{\AA}})$"
    return fr"$a = {aA}\,\mathrm{{\AA}}$"


def plot_set(files, title, outpath, nrows=4, ncols=2):
    fig, axes = plt.subplots(nrows, ncols, figsize=(8.6, 2.6*nrows))
    axes = axes.flatten()
    for i, fpath in enumerate(files):
        if i >= len(axes):
            break
        gd, fixed, twop, buf, aA, Teq = load_one(fpath)
        ax = axes[i]

        # GD
        if gd.size > 0:
            T_gd = gd[:, 0]
            P_gd = gd[:, 1]
            y_gd = gd_dpdlnt(T_gd, P_gd)
            pos = y_gd > 0
            ax.plot(T_gd[pos], y_gd[pos], color="k", lw=1.1,
                    label="GD calc\\_P")

        for arr, col, ls, lbl in [
                (fixed, "C0", "-",  "MC fixed-grid"),
                (twop,  "C3", "--", "MC 2pass (adaptive)"),
                (buf,   "C2", ":",  "MC buffered (adaptive)"),
        ]:
            if arr.size == 0:
                continue
            Tmid = arr[:, 0]
            y    = arr[:, 2]
            pos  = (y > 0) & np.isfinite(y)
            ax.plot(Tmid[pos], y[pos], color=col, ls=ls, lw=1.0, label=lbl)

        # Teq vertical line
        if Teq is not None and Teq > 0:
            ax.axvline(Teq, color="purple", lw=0.8, ls="-.",
                       label=fr"$T_{{\rm eq}} = {Teq:.2f}$ K")

        ax.set_xscale("log")
        ax.set_yscale("log")
        # Auto-fit x range around adaptive engines' active region
        all_T = []
        for arr in (twop, buf):
            if arr.size > 0:
                y = arr[:, 2]
                m = (y > 0) & np.isfinite(y)
                if m.any():
                    all_T.append(arr[m, 0])
        if all_T:
            Tcat = np.concatenate(all_T)
            T_lo = max(Tcat.min() * 0.5, 1.0)
            T_hi = Tcat.max() * 2.0
            # Ensure Teq stays inside the visible range, even when the MC
            # trajectory collapses far from it.
            if Teq is not None and Teq > 0:
                T_lo = min(T_lo, Teq * 0.5)
                T_hi = max(T_hi, Teq * 2.0)
            ax.set_xlim(T_lo, max(T_hi, 5e1))
        else:
            ax.set_xlim(2, 2e3)
        ax.set_ylim(1e-3, 1e4)
        ax.set_xlabel(r"$T$ (K)", fontsize=8)
        ax.set_ylabel(r"$dP/d\ln T$", fontsize=8)
        ax.set_title(fmt_size(aA), fontsize=9)
        ax.legend(loc="upper left", fontsize=6, frameon=True)
        ax.grid(True, which="both", ls=":", lw=0.4, alpha=0.5)
        ax.tick_params(labelsize=7)
    # Hide any unused axes
    for j in range(len(files), len(axes)):
        axes[j].set_visible(False)
    fig.suptitle(title, fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.98))
    fig.savefig(outpath, dpi=140)
    print("wrote", outpath)


def main():
    files = sorted(glob.glob(os.path.join(
        ROOT, "output", "pT_compare_large_*_a*A.dat")))
    if not files:
        raise SystemExit("No pT_compare_large_*.dat files found; run "
                         "./main_pT_compare.x run_pT_compare_large.nml first.")
    plot_set(files,
             r"Single-grain $P(T)$ for large astrodust grains: "
             r"3 MC engines + GD ($N_{\mathrm{evt}}=2\times10^5$, $\log U=0.20$)",
             os.path.join(ROOT, "output", "pT_compare_large.pdf"))


if __name__ == "__main__":
    main()
