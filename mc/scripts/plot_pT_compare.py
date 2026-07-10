#!/usr/bin/env python3
"""
Single-grain P(T) comparison: adaptive-grid Monte Carlo
(mc_run_engine_2pass) vs Guhathakurta-Draine matrix solver (calc_P).
Reads pT_compare_<pop>_NN_aXXXXA.dat files written by main_pT_compare.x
and produces a multi-panel figure.

The data file has four blocks separated by '# ===' lines (see
mc_pT/CLAUDE.md session log 2026-05-17):
  Block 1 (GD)          : T [K]   P(T_i)   dP/dlnT_approx
  Block 2 (MC fixed)    : T_mid   dP/dT    dP/dlnT
  Block 3 (MC 2pass)    : T_mid   dP/dT    dP/dlnT     <-- used here
  Block 4 (MC buffered) : T_mid   dP/dT    dP/dlnT

The adaptive 2pass block is used because the fixed-grid log[2,5000]K
NHIST=600 histogram collapses to a single populated bin once Teq narrows
the trajectory (a >~ 80 A in the astrodust population).  2pass adaptively
centers the grid on the observed [Tmin, Tmax] of the trajectory replay
and gives ~50x sharper resolution.  The integrated SED is invariant to
engine choice to better than 0.5%; only the P(T) shape differs.
"""
import os, glob, re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)


def load_one(path):
    """Return (T_gd, P_gd_per_bin, T_mc, dPdlnT_mc) and the radius label.
    MC values come from block 3 (mc_run_engine_2pass adaptive grid)."""
    text = open(path).read()
    blocks = re.split(r"^# ===.*$", text, flags=re.MULTILINE)
    if len(blocks) < 5:
        raise ValueError(f"expected 5 segments (preamble + 4 blocks) in {path}, "
                         f"got {len(blocks)}")
    def parse_rows(block):
        rows = []
        for line in block.splitlines():
            line = line.strip()
            if (not line) or line.startswith("#"):
                continue
            rows.append([float(x) for x in line.split()])
        return np.array(rows)
    gd    = parse_rows(blocks[1])
    mc_2p = parse_rows(blocks[3])         # 2pass adaptive
    m = re.search(r"_a(\d+)A\.dat$", path)
    if m:
        aA = int(m.group(1))
        label = fr"$a = {aA}\,\mathrm{{\AA}}$"
    else:
        label = os.path.basename(path)
    return gd[:, 0], gd[:, 1], mc_2p[:, 0], mc_2p[:, 2], label


def plot_pop(files, title, outpath, ymin=1e-6):
    n = len(files)
    fig, axes = plt.subplots(2, 2, figsize=(8.5, 6.5), sharex=False)
    axes = axes.flatten()
    for i, fpath in enumerate(files):
        T_gd, P_gd_bin, T_mc, dPdlnT_mc, lbl = load_one(fpath)
        # Convert the GD probability in each bin to dP/dlnT.  T_first is log-spaced
        # so dlnT_i = (1/2) ln(T_{i+1}/T_{i-1}).
        lnT = np.log(np.maximum(T_gd, 1e-30))
        dlnT = np.zeros_like(T_gd)
        dlnT[1:-1] = 0.5 * (lnT[2:] - lnT[:-2])
        dlnT[0]    = lnT[1] - lnT[0]
        dlnT[-1]   = lnT[-1] - lnT[-2]
        dPdlnT_gd = P_gd_bin / np.maximum(dlnT, 1e-30)
        ax = axes[i]
        gd_pos = dPdlnT_gd > 0
        mc_pos = dPdlnT_mc > 0
        ax.plot(T_gd[gd_pos], dPdlnT_gd[gd_pos], color="k", lw=1.4,
                label="GD (calc\\_P)")
        ax.plot(T_mc[mc_pos], dPdlnT_mc[mc_pos], color="C3", lw=0.9, ls="--",
                label="MC (2pass adaptive)")
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlim(2, 2e3); ax.set_ylim(ymin, 5)
        ax.set_xlabel(r"$T$ (K)", fontsize=9)
        ax.set_ylabel(r"$dP/d\ln T$", fontsize=9)
        ax.set_title(lbl, fontsize=10)
        ax.legend(loc="lower left", fontsize=8, frameon=True)
        ax.grid(True, which="both", ls=":", lw=0.4, alpha=0.5)
        ax.tick_params(labelsize=8)
    fig.suptitle(title, fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(outpath, dpi=140)
    print("wrote", outpath)


def newest_per_index(pattern):
    """Glob `pattern`, then keep only the newest file for each panel index
    NN (the `_NN_` field).  Guards against stale duplicates left over from
    earlier runs that used a different radius-label naming convention."""
    by_idx = {}
    for path in glob.glob(pattern):
        m = re.search(r"_(\d\d)_a", os.path.basename(path))
        if not m:
            continue
        idx = m.group(1)
        if idx not in by_idx or os.path.getmtime(path) > os.path.getmtime(by_idx[idx]):
            by_idx[idx] = path
    return [by_idx[k] for k in sorted(by_idx)]


def main():
    ad_files = newest_per_index(os.path.join(ROOT, "output", "pT_compare_ad_*_a*A.dat"))
    if ad_files:
        plot_pop(ad_files,
                 r"Single-grain $P(T)$: Astrodust S2, Mathis ISRF ($\log U=0.20$), $N_{\mathrm{evt}}=5\times10^5$",
                 os.path.join(ROOT, "output", "pT_compare_astrodust.pdf"))
    pah_files = newest_per_index(os.path.join(ROOT, "output", "pT_compare_pah_*_a*A.dat"))
    if pah_files:
        plot_pop(pah_files,
                 r"Single-grain $P(T)$: PAH, Mathis ISRF ($\log U=0.20$), $N_{\mathrm{evt}}=5\times10^5$",
                 os.path.join(ROOT, "output", "pT_compare_pah.pdf"))


if __name__ == "__main__":
    main()
