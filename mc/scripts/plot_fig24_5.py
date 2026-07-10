#!/usr/bin/env python3
"""
Plot Monte Carlo T(t) histories for PIIM Figure 24.5.

For each (a, U) panel, reads output/fig24_5/<aLabel>_<ULabel>_Tt.dat which
contains the actual MC sub-step trajectory (one row per cool_segment sub-step,
stride-sampled to limit file size), plus output/fig24_5/<...>_evt.dat which
holds the event-boundary (T_pre, T_post) records.

The plotted curve is the raw sub-step trajectory connected with line
segments -- no analytic reconstruction. Vertical jumps at photon-absorption
events are shown by adding a (T_pre -> T_post) segment at each event time
from the _evt.dat record.
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT  = os.path.join(ROOT, "output", "fig24_5")

SIZES    = ["200A", "100A", "50A", "20A", "10A"]   # PIIM Fig 24.5 top-to-bottom
U_LABELS = ["U1", "U100"]
U_TITLES = ["U = 1", r"U = $10^{2}$"]


def load_traj(path):
    """Sub-step trajectory: (t[s], T[K]) per row."""
    if not os.path.exists(path):
        return np.empty(0), np.empty(0)
    arr = np.loadtxt(path, comments="#")
    if arr.size == 0:
        return np.empty(0), np.empty(0)
    if arr.ndim == 1:
        arr = arr[None, :]
    return arr[:, 0], arr[:, 1]


def load_evt(path):
    """Event-boundary record: (t, T_pre, T_post, dt_event) per row."""
    if not os.path.exists(path):
        return np.empty(0), np.empty(0), np.empty(0)
    arr = np.loadtxt(path, comments="#")
    if arr.size == 0:
        return np.empty(0), np.empty(0), np.empty(0)
    if arr.ndim == 1:
        arr = arr[None, :]
    return arr[:, 0], arr[:, 1], arr[:, 2]


def plot_panel(ax, t_traj, T_traj, t_evt, T_pre, T_post, title):
    # Plot the actual sub-step trajectory
    if t_traj.size:
        ax.plot(t_traj, T_traj, color="k", lw=0.5)
    # Draw vertical jumps at each event from T_pre to T_post (skip the
    # initial record where T_pre == T_post == T_init).
    for k in range(1, len(t_evt)):
        if T_post[k] > T_pre[k]:
            ax.plot([t_evt[k], t_evt[k]], [T_pre[k], T_post[k]],
                    color="k", lw=0.5)
    ax.set_yscale("log")
    ax.set_ylim(2.0, 2.0e3)            # extends below T_CMB = 2.725 K
    ax.set_xlim(0, 1.0e5)
    ax.text(0.97, 0.92, title, transform=ax.transAxes,
            ha="right", va="top", fontsize=8,
            bbox=dict(boxstyle="round,pad=0.2", fc="w", ec="0.7", alpha=0.85))
    ax.tick_params(labelsize=8)


def main():
    fig, axes = plt.subplots(5, 2, figsize=(8.0, 7.5), sharex="col")
    for j, (Ulbl, Utitle) in enumerate(zip(U_LABELS, U_TITLES)):
        for i, sz in enumerate(SIZES):
            t_traj, T_traj = load_traj(os.path.join(OUT, f"{sz}_{Ulbl}_Tt.dat"))
            t_evt, T_pre, T_post = load_evt(os.path.join(OUT, f"{sz}_{Ulbl}_evt.dat"))
            title = f"a = {sz},  {Utitle},  N_evt = {max(len(t_evt) - 1, 0)}"
            plot_panel(axes[i, j], t_traj, T_traj, t_evt, T_pre, T_post, title)
    for ax in axes[-1, :]:
        ax.set_xlabel(r"$t$ (s)", fontsize=9)
    for ax in axes[:, 0]:
        ax.set_ylabel(r"$T$ (K)", fontsize=9)
    fig.suptitle(r"MC reproduction of PIIM Fig.\ 24.5 -- graphite, DL01 enthalpy",
                 fontsize=10, y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.985))
    outpdf = os.path.join(ROOT, "output", "fig24_5.pdf")
    fig.savefig(outpdf, dpi=140)
    print("wrote", outpdf)


if __name__ == "__main__":
    main()
