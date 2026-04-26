#!/usr/bin/env python3
"""
Plot C++ vs Fortran surfdisp96 comparison results.
Run from the build directory where test_compare wrote the CSV files.
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

CASES = [
    ("love_fund_phase.csv",  "Love fund. phase vel."),
    ("love_fund_group.csv",  "Love fund. group vel."),
    ("love_1st_phase.csv",   "Love 1st-higher phase vel."),
    ("ray_fund_phase.csv",   "Rayleigh fund. phase vel."),
    ("ray_fund_group.csv",   "Rayleigh fund. group vel."),
    ("ray_1st_phase.csv",    "Rayleigh 1st-higher phase vel."),
]

def load(fname):
    if not os.path.exists(fname):
        print(f"  [skip] {fname} not found")
        return None, None, None
    data = np.genfromtxt(fname, delimiter=",", skip_header=1)
    if data.ndim == 1:
        data = data[np.newaxis, :]
    mask = (data[:, 1] > 0) & (data[:, 2] > 0)
    return data[mask, 0], data[mask, 1], data[mask, 2]   # t, cpp, fortran

fig = plt.figure(figsize=(16, 10))
fig.suptitle(
    "surfdisp96: C++ vs Fortran comparison\n"
    "(10-layer crustal/mantle model, flat Earth)",
    fontsize=13, fontweight="bold")

gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.55, wspace=0.35)

for idx, (fname, title) in enumerate(CASES):
    row, col = divmod(idx, 3)
    t, cpp, fort = load(fname)

    inner = gridspec.GridSpecFromSubplotSpec(
        2, 1, subplot_spec=gs[row, col], height_ratios=[3, 1], hspace=0.08)
    ax_main = fig.add_subplot(inner[0])
    ax_res  = fig.add_subplot(inner[1], sharex=ax_main)

    ax_main.set_title(title, fontsize=9)

    if t is not None and len(t) > 0:
        ax_main.semilogx(t, fort, "b-",  lw=2,   label="Fortran", alpha=0.7)
        ax_main.semilogx(t, cpp,  "r--", lw=1.5, label="C++",     alpha=0.9)
        ax_main.set_ylabel("Velocity (km/s)", fontsize=8)
        ax_main.tick_params(labelbottom=False, labelsize=7)
        ax_main.legend(fontsize=7, loc="best")
        ax_main.grid(True, which="both", ls=":", alpha=0.4)

        rel = np.abs(cpp - fort) / (np.abs(fort) + 1e-30)
        ax_res.semilogx(t, rel * 100, "k-", lw=1)
        ax_res.set_ylabel("|Δ| %", fontsize=7)
        ax_res.set_xlabel("Period (s)", fontsize=8)
        ax_res.tick_params(labelsize=7)
        ax_res.grid(True, which="both", ls=":", alpha=0.4)
        ax_res.axhline(0, color="k", lw=0.5)

        print(f"{title:40s}  max={rel.max()*100:.4f}%  mean={rel.mean()*100:.4f}%")
    else:
        ax_main.text(0.5, 0.5, "no data", ha="center", va="center",
                     transform=ax_main.transAxes, fontsize=9, color="gray")
        ax_res.set_visible(False)

out = "comparison.png"
plt.savefig(out, dpi=150, bbox_inches="tight")
print(f"\nSaved {out}")
plt.show()
