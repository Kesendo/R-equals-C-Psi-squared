"""Hardware version of bellplus_trajectory_on_mandelbrot.png.

Same Mandelbrot background and trajectory style as the April original
(critical_slowing_mandelbrot_overlay.py) but the points along the
trajectory are the actual hardware measurements from
`run_cusp_precision.py` on Kingston pair (14, 15).

Reads:
    cusp_precision_ibm_kingston_*.json  (pass path as 1st arg)

Writes:
    visualizations/bellplus_trajectory_on_mandelbrot_hardware.png
    visualizations/bellplus_trajectory_on_mandelbrot_hardware_zoom.png
"""
import json
import math
import os
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

if len(sys.argv) < 2:
    print("Usage: python _plot_cusp_mandelbrot_hardware.py <cusp_precision_*.json>")
    sys.exit(1)

json_path = Path(sys.argv[1])
with open(json_path, encoding='utf-8') as f:
    data = json.load(f)

t_meas = np.array([d['t_us'] for d in data['cpsi_data']])
cpsi_meas = np.array([d['cpsi'] for d in data['cpsi_data']])
gamma_fit = data['gamma_fit']
backend = data.get('backend', 'unknown')
pair = data['pair']['qubits']

# ====================================================================
#  Mandelbrot set (vectorised escape time)
# ====================================================================
print("Computing Mandelbrot set...")


def compute_mandelbrot(x_range, y_range, nx, ny, max_iter=200):
    x = np.linspace(*x_range, nx)
    y = np.linspace(*y_range, ny)
    X, Y = np.meshgrid(x, y)
    C = X + 1j * Y
    escape = np.zeros_like(X, dtype=int)
    Z = np.zeros_like(C)
    mask = np.ones_like(X, dtype=bool)
    for i in range(max_iter):
        Z[mask] = Z[mask] ** 2 + C[mask]
        escaped = mask & (np.abs(Z) > 2)
        escape[escaped] = i
        mask = mask & ~escaped
    escape[mask] = max_iter
    return X, Y, escape


# ====================================================================
#  Theoretical trajectory at γ_fit
# ====================================================================

t_theory = np.linspace(0, 8.0, 2000)
f_theory = np.exp(-4 * gamma_fit * t_theory)
cpsi_theory = f_theory * (1 + f_theory ** 2) / 6

VIS_DIR = Path(__file__).parent.parent / "visualizations"
VIS_DIR.mkdir(exist_ok=True)


def make_plot(x_range, y_range, nx, ny, suffix, figsize=(14, 10)):
    X, Y, escape = compute_mandelbrot(x_range, y_range, nx, ny)
    log_escape = np.log(escape.astype(float) + 1)

    fig, ax = plt.subplots(figsize=figsize)

    # Mandelbrot background — same style as critical_slowing_mandelbrot_overlay.py
    ax.imshow(log_escape,
              extent=[x_range[0], x_range[1], y_range[0], y_range[1]],
              origin='lower', cmap='bone',
              aspect='equal', interpolation='bilinear')

    # Theoretical trajectory (faint line on the real axis)
    ax.plot(cpsi_theory, np.zeros_like(cpsi_theory), color='lightblue',
            linewidth=1.0, alpha=0.7, label=fr'F25 theory at $\gamma^{{fit}} = {gamma_fit*1e3:.2f}$/ms')

    # Hardware points colored by time
    sc = ax.scatter(cpsi_meas, np.zeros_like(cpsi_meas),
                    c=t_meas, cmap='viridis', s=80, zorder=5,
                    edgecolor='black', linewidth=0.6,
                    label=f'Hardware ({backend}, qubits {pair})')

    # Cusp marker
    ax.plot(0.25, 0, 'X', color='lime', markersize=18, mew=2.5,
            label=r'$C\Psi = 1/4$  (cardioid cusp)', zorder=6)

    # Start and end markers
    ax.plot(1/3, 0, 'o', color='red', markersize=12,
            label=r'$t = 0$:  $C\Psi = 1/3$', zorder=6)
    ax.plot(0, 0, 's', color='cyan', markersize=10,
            label=r'$t \to \infty$:  $C\Psi \to 0$', zorder=6)

    # Annotations
    if suffix == "":
        ax.annotate("Cardioid cusp\n$c = 1/4$",
                    xy=(0.25, 0), xytext=(0.45, 0.45),
                    fontsize=10, color='green', ha='center',
                    arrowprops=dict(arrowstyle='->', color='green', lw=1.2))
        ax.annotate("Mandelbrot set\n(bound orbits)",
                    xy=(-0.7, -0.6), fontsize=11, ha='center', color='dimgray')
        ax.annotate("Divergence\nzone", xy=(0.42, 0.15), fontsize=9, color='orange')

    ax.set_xlabel(r"$\mathrm{Re}(c)$", fontsize=12)
    ax.set_ylabel(r"$\mathrm{Im}(c)$", fontsize=12)
    ax.set_title(
        fr"Bell$^+$ decoherence trajectory on the Mandelbrot set"
        fr"  —  $\mathtt{{{backend}}}$  (qubits {pair})",
        fontsize=12,
    )
    ax.legend(loc='upper left', fontsize=9)
    ax.grid(alpha=0.2)

    cbar = plt.colorbar(sc, ax=ax, label=r"$t$  $[\mu s]$", shrink=0.7)

    out = VIS_DIR / f"bellplus_trajectory_on_mandelbrot_hardware{suffix}.png"
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"Saved: {out}")


# Full view
make_plot((-2.2, 0.8), (-1.2, 1.2), 1200, 900, "", figsize=(14, 10))

# Zoomed near cardioid cusp — taller aspect to keep equal axes readable
make_plot((-0.10, 0.50), (-0.20, 0.20), 1200, 800, "_zoom", figsize=(12, 6))
