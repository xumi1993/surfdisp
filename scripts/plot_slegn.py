#!/usr/bin/env python3
"""
plot_slegn.py
Plot C++ vs Fortran comparison for slegn96 and slegnpu outputs.
Run from the build directory where the CSV files were written.
"""

import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def read_csv(fname):
    """Read depth_km, cpp_value, fortran_value columns."""
    if not os.path.exists(fname):
        print(f"WARNING: {fname} not found, skipping.")
        return None, None, None
    data = np.loadtxt(fname, delimiter=',', skiprows=1)
    depth   = data[:, 0]
    cpp_val = data[:, 1]
    frt_val = data[:, 2]
    return depth, cpp_val, frt_val


def pct_resid(cpp, frt):
    """Percentage residual |cpp - frt| / max(|frt|, 1e-30) * 100."""
    denom = np.where(np.abs(frt) > 1e-30, np.abs(frt), 1e-30)
    return np.abs(cpp - frt) / denom * 100.0


def plot_panel(ax, ax_res, depth, cpp, frt, label, color_cpp='C0', color_frt='C1'):
    """Plot main comparison + residual sub-panel."""
    ax.plot(frt, depth, color=color_frt, lw=2,   label='Fortran', linestyle='--')
    ax.plot(cpp, depth, color=color_cpp, lw=1.5, label='C++',     linestyle='-')
    ax.set_ylabel('Depth (km)')
    ax.set_title(label, fontsize=10)
    ax.invert_yaxis()
    ax.legend(fontsize=7)
    ax.grid(True, alpha=0.3)

    resid = pct_resid(cpp, frt)
    ax_res.plot(resid, depth, color='C2', lw=1.2)
    ax_res.set_xlabel('|Δ| %')
    ax_res.invert_yaxis()
    ax_res.set_xscale('log')
    ax_res.grid(True, alpha=0.3)
    ax_res.set_xlim(left=1e-12)


# =============================================================
# Figure 1: slegn96
# =============================================================
files96 = {
    'disp':   'slegn96_disp.csv',
    'stress': 'slegn96_stress.csv',
    'dc2db':  'slegn96_dc2db.csv',
    'dc2dh':  'slegn96_dc2dh.csv',
    'dc2dr':  'slegn96_dc2dr.csv',
}

labels96 = {
    'disp':   'Displacement',
    'stress': 'Stress',
    'dc2db':  'dc/dVs',
    'dc2dh':  'dc/dh',
    'dc2dr':  'dc/drho',
}

n_panels = len(files96)
fig1, axes = plt.subplots(2, n_panels,
                           figsize=(4 * n_panels, 10),
                           gridspec_kw={'height_ratios': [3, 1]})
fig1.suptitle('slegn96: C++ vs Fortran (T=10 s, flat earth)', fontsize=12)

for col, key in enumerate(files96):
    depth, cpp, frt = read_csv(files96[key])
    ax_main = axes[0, col]
    ax_res  = axes[1, col]
    if depth is None:
        ax_main.set_visible(False)
        ax_res.set_visible(False)
        continue
    plot_panel(ax_main, ax_res, depth, cpp, frt, labels96[key])

plt.tight_layout()
fig1.savefig('slegn96_comparison.png', dpi=150, bbox_inches='tight')
print("Saved slegn96_comparison.png")
plt.close(fig1)


# =============================================================
# Figure 2: slegnpu group kernels
# =============================================================
filespu = {
    'du2db': 'slegnpu_du2db.csv',
    'du2dh': 'slegnpu_du2dh.csv',
    'du2dr': 'slegnpu_du2dr.csv',
}

labelspu = {
    'du2db': 'dU/dVs',
    'du2dh': 'dU/dh',
    'du2dr': 'dU/drho',
}

n_panels2 = len(filespu)
fig2, axes2 = plt.subplots(2, n_panels2,
                            figsize=(4 * n_panels2, 10),
                            gridspec_kw={'height_ratios': [3, 1]})
fig2.suptitle('slegnpu: C++ vs Fortran group kernels (T=20 s, flat earth)', fontsize=12)

for col, key in enumerate(filespu):
    depth, cpp, frt = read_csv(filespu[key])
    ax_main = axes2[0, col]
    ax_res  = axes2[1, col]
    if depth is None:
        ax_main.set_visible(False)
        ax_res.set_visible(False)
        continue
    plot_panel(ax_main, ax_res, depth, cpp, frt, labelspu[key])

plt.tight_layout()
fig2.savefig('slegnpu_comparison.png', dpi=150, bbox_inches='tight')
print("Saved slegnpu_comparison.png")
plt.close(fig2)

print("Done.")
