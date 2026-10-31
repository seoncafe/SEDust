#!/usr/bin/env python3
"""
Plot Monte Carlo T(t) histories for PIIM Figure 24.5.

For each (a, U) panel, reads <data>/fig24_5/<aLabel>_<ULabel>_Tt.dat, the
sub-step trajectory written by main_mc.x (stride-sampled to limit file size),
and <data>/fig24_5/<...>_evt.dat, the event-boundary records.  The first
_evt.dat record is the start of the recorded interval and the last is its end
at t_max; the rows in between are photon absorptions (T_post > T_pre).

The plotted curve is the raw sub-step trajectory connected with line
segments -- no analytic reconstruction.  Each absorption is drawn as a
vertical segment from T_pre to T_post at its event time.

Usage: plot_fig24_5.py [--data DIR] [--out FILE.pdf]
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

SIZES    = ["200A", "100A", "50A", "20A", "10A"]   # PIIM Fig 24.5 top-to-bottom
U_LABELS = ["U1", "U100"]
U_TITLES = [r"$U = 1$", r"$U = 10^{2}$"]


def load_table(path, ncol):
    if not os.path.exists(path):
        raise SystemExit(f"missing input: {path}")
    arr = np.loadtxt(path, comments="#")
    return arr.reshape(-1, ncol)


def plot_panel(ax, traj, evt, title):
    ax.plot(traj[:, 0], traj[:, 1], color="k", lw=0.5)
    jumps = evt[1:-1]
    if jumps.size:
        nan = np.full(len(jumps), np.nan)
        x = np.column_stack([jumps[:, 0], jumps[:, 0], nan]).ravel()
        y = np.column_stack([jumps[:, 1], jumps[:, 2], nan]).ravel()
        ax.plot(x, y, color="k", lw=0.5)
    ax.set_yscale("log")
    ax.set_ylim(2.0, 2.0e3)
    ax.set_xlim(0, 1.0e5)
    ax.text(0.97, 0.92, title, transform=ax.transAxes,
            ha="right", va="top", fontsize=8)
    ax.tick_params(labelsize=8)


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    p.add_argument("--data", default=os.path.join(ROOT, "output"))
    p.add_argument("--out", default=os.path.join(ROOT, "output", "fig24_5.pdf"))
    args = p.parse_args()
    out_dir = os.path.join(args.data, "fig24_5")

    fig, axes = plt.subplots(5, 2, figsize=(8.0, 7.5), sharex="col")
    for j, (Ulbl, Utitle) in enumerate(zip(U_LABELS, U_TITLES)):
        for i, sz in enumerate(SIZES):
            traj = load_table(os.path.join(out_dir, f"{sz}_{Ulbl}_Tt.dat"), 2)
            evt = load_table(os.path.join(out_dir, f"{sz}_{Ulbl}_evt.dat"), 4)
            n_evt = int(np.sum(evt[1:-1, 2] > evt[1:-1, 1]))
            title = (rf"$a = {sz[:-1]}$\,\AA, " + Utitle + ", "
                     + rf"$N_{{\rm evt}} = {n_evt}$")
            plot_panel(axes[i, j], traj, evt, title)
    for ax in axes[-1, :]:
        ax.set_xlabel(r"$t$ (s)", fontsize=9)
    for ax in axes[:, 0]:
        ax.set_ylabel(r"$T$ (K)", fontsize=9)
    fig.suptitle(r"MC reproduction of PIIM Fig.~24.5 -- graphite, DL01 enthalpy",
                 fontsize=10, y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.985))
    fig.savefig(args.out)
    print("wrote", args.out)


if __name__ == "__main__":
    main()
