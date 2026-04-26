#!/usr/bin/env python3
"""
Plot sregn comparison: C++ vs Fortran (no pandas)
"""

import matplotlib.pyplot as plt
import numpy as np
import os

def read_csv(csv_file):
    """Read CSV file, return layer, cpp_value, fortran_value."""
    data = []
    with open(csv_file, 'r') as f:
        header = f.readline().strip()
        for line in f:
            parts = line.strip().split(',')
            if len(parts) >= 3:
                try:
                    data.append((float(parts[0]), float(parts[1]), float(parts[2])))
                except ValueError:
                    pass
    if not data:
        return None, None, None
    data = np.array(data)
    return data[:, 0], data[:, 1], data[:, 2]

def plot_comparison(csv_file, title, ylabel):
    """Plot C++ vs Fortran comparison from CSV."""
    if not os.path.exists(csv_file):
        return None

    layer, cpp_val, fort_val = read_csv(csv_file)
    if layer is None:
        return None

    # Relative error
    rel_err = np.abs(cpp_val - fort_val) / (np.abs(fort_val) + 1e-30)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # Values comparison
    ax1.plot(layer, cpp_val, 'o-', label='C++', alpha=0.7)
    ax1.plot(layer, fort_val, 's--', label='Fortran', alpha=0.7)
    ax1.set_xlabel('Layer')
    ax1.set_ylabel(ylabel)
    ax1.set_title(f'{title} - Values')
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # Relative error
    ax2.semilogy(layer, rel_err, 'o-', color='red', alpha=0.7)
    ax2.set_xlabel('Layer')
    ax2.set_ylabel('Relative Error')
    ax2.set_title(f'{title} - Relative Error')
    ax2.grid(True, alpha=0.3)
    ax2.axhline(1e-4, color='orange', linestyle='--', label='Threshold 1e-4')
    ax2.legend()

    plt.tight_layout()
    return fig

# Create output directory
os.makedirs('plots_sregn', exist_ok=True)

print("Generating sregn comparison plots...")

# sregn96 test
plots = [
    ('sregn96_dispu.csv', 'sregn96 Horizontal Displacement', 'ur (km/s)'),
    ('sregn96_dispw.csv', 'sregn96 Vertical Displacement', 'uz (km/s)'),
    ('sregn96_dc2da.csv', 'sregn96 dc/dVp Kernel', 'dc/dVp'),
    ('sregn96_dc2db.csv', 'sregn96 dc/dVs Kernel', 'dc/dVs'),
    ('sregn96_dc2dh.csv', 'sregn96 dc/dh Kernel', 'dc/dh'),
    ('sregn96_dc2dr.csv', 'sregn96 dc/drho Kernel', 'dc/drho'),
]

for csv_file, title, ylabel in plots:
    print(f"  {csv_file}...", end='', flush=True)
    fig = plot_comparison(csv_file, title, ylabel)
    if fig:
        fig.savefig(f'plots_sregn/{csv_file[:-4]}.png', dpi=150, bbox_inches='tight')
        plt.close(fig)
        print(" OK")
    else:
        print(" SKIP")

# sregn96_hti test
hti_plots = [
    ('sregn96_hti_dc2dgc.csv', 'sregn96_hti dc/dgc Kernel (HTI)', 'dc/dgc'),
    ('sregn96_hti_dc2dgs.csv', 'sregn96_hti dc/dgs Kernel (HTI)', 'dc/dgs'),
]

for csv_file, title, ylabel in hti_plots:
    print(f"  {csv_file}...", end='', flush=True)
    fig = plot_comparison(csv_file, title, ylabel)
    if fig:
        fig.savefig(f'plots_sregn/{csv_file[:-4]}.png', dpi=150, bbox_inches='tight')
        plt.close(fig)
        print(" OK")
    else:
        print(" SKIP")

# sregnpu test (group velocity kernels)
pu_plots = [
    ('sregnpu_du2da.csv', 'sregnpu dU/dVp Kernel', 'dU/dVp'),
    ('sregnpu_du2db.csv', 'sregnpu dU/dVs Kernel', 'dU/dVs'),
    ('sregnpu_du2dh.csv', 'sregnpu dU/dh Kernel', 'dU/dh'),
    ('sregnpu_du2dr.csv', 'sregnpu dU/drho Kernel', 'dU/drho'),
]

for csv_file, title, ylabel in pu_plots:
    print(f"  {csv_file}...", end='', flush=True)
    fig = plot_comparison(csv_file, title, ylabel)
    if fig:
        fig.savefig(f'plots_sregn/{csv_file[:-4]}.png', dpi=150, bbox_inches='tight')
        plt.close(fig)
        print(" OK")
    else:
        print(" SKIP")

print("\nPlots saved to plots_sregn/")
