#!/usr/bin/env python3
"""
Plot the Trad_MC / Teq ratio of each grain versus grain radius, exposing the
integrator-induced bias in MC-all mode.

Reads output/diag_Teq_vs_Trad.dat with columns:
  a_eff[um]  Teq[K]  Trad_MC[K]  Trad/Teq  hist_weight_sum
written by mc_sed_loop when use_mc_for_large_grain=.true.
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DAT  = os.path.join(ROOT, "output", "diag_Teq_vs_Trad.dat")

arr = np.loadtxt(DAT)
arr = arr[np.argsort(arr[:, 0])]
a_um, Teq, Trad, ratio, wsum = arr.T

fig, axes = plt.subplots(2, 1, figsize=(7, 6), sharex=True,
                         gridspec_kw={"height_ratios": [2, 1]})
ax = axes[0]
ax.semilogx(a_um, Teq,  color="C0", lw=1.5, label=r"calc\_Teq")
ax.semilogx(a_um, Trad, color="C3", lw=1.0, ls="--", label=r"MC $\langle T^4\rangle^{1/4}$")
ax.set_ylabel(r"$T$ (K)")
ax.legend(loc="upper right")
ax.grid(True, ls=":", alpha=0.5)
ax.set_title("Single-grain Teq vs MC $T_{\\rm rad}$ "
             "(\\texttt{use\\_mc\\_for\\_large\\_grain}=T, 2pass engine)",
             fontsize=10)

ax = axes[1]
ax.semilogx(a_um, ratio, color="k", lw=1.0)
ax.axhline(1.0, color="0.5", lw=0.5)
for h in (0.9, 0.7):
    ax.axhline(h, color="0.7", lw=0.4, ls=":")
ax.set_xlabel(r"$a_{\rm eff}$ ($\mu$m)")
ax.set_ylabel(r"$T_{\rm rad}^{\rm MC} / T_{\rm eq}$")
ax.set_ylim(0.0, 1.1)
ax.grid(True, ls=":", alpha=0.5)

out = os.path.join(ROOT, "output", "diag_Teq_vs_Trad.pdf")
fig.tight_layout()
fig.savefig(out, dpi=140)
print("wrote", out)
