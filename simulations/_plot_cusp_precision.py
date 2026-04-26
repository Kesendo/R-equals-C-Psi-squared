#!/usr/bin/env python3
"""Plot CΨ(t) hardware trajectory + F25 prediction overlay.

Reads the JSON output from `run_cusp_precision.py` (in the IBM tomography
directory) and produces two figures:

  visualizations/cusp_precision_trajectory.png  — CΨ vs t with F25 fit line
  visualizations/cusp_precision_residuals.png   — pointwise (data − fit) with σ shot noise

The F25 closed form is CΨ(t) = f·(1+f²)/6 with f = exp(−4γt). The hardware
γ is extracted from the trajectory itself (γ_fit), and a dashed line shows
what the calibration's T2-echo would have predicted (γ_calib) for contrast.
"""
import json
import math
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

if len(sys.argv) < 2:
    print("Usage: python _plot_cusp_precision.py <path_to_cusp_precision_*.json> [gamma_t2echo_per_us]")
    sys.exit(1)

json_path = Path(sys.argv[1])
with open(json_path, encoding='utf-8') as f:
    data = json.load(f)

t = np.array([d['t_us'] for d in data['cpsi_data']])
cpsi = np.array([d['cpsi'] for d in data['cpsi_data']])
# If the run used --gamma-override, JSON's gamma_calib equals the override.
# Pass the actual T2-echo γ as 2nd CLI arg to compare hardware to that prediction.
gamma_calib = float(sys.argv[2]) if len(sys.argv) >= 3 else data['gamma_calib']
gamma_fit = data['gamma_fit']
F_CROSS = data['F_CROSS']
K_CROSS = data['K_CROSS']
backend = data.get('backend', 'unknown')
pair = data['pair']['qubits']
shots = data['shots']

t_cross_calib = K_CROSS / gamma_calib
t_cross_fit = K_CROSS / gamma_fit


def cpsi_bellplus(t_us, gamma):
    f = np.exp(-4.0 * gamma * t_us)
    return f * (1.0 + f * f) / 6.0


# Smooth curve for plotting the prediction
t_fine = np.linspace(0, t.max() * 1.05, 600)
cpsi_calib = cpsi_bellplus(t_fine, gamma_calib)
cpsi_fitline = cpsi_bellplus(t_fine, gamma_fit)

# Shot noise for error bars: each Pauli expectation has variance ~1/N_shots,
# CΨ is a function of multiple expectations; we approximate σ_CΨ ≈ 0.5/√shots.
sigma = 0.5 / math.sqrt(shots)

# ════════════════════════════════════════════════════════════════════
#  FIGURE 1: trajectory
# ════════════════════════════════════════════════════════════════════

fig, ax = plt.subplots(figsize=(8, 5.5))

ax.axhline(0.25, color='gray', linestyle=':', linewidth=1.0,
           label=r'$C\Psi = 1/4$  (cusp)')
ax.axvline(t_cross_calib, color='C1', linestyle=':', alpha=0.5,
           label=fr'$t_{{cross}}^{{calib}} = {t_cross_calib:.1f}\,\mu s$  ($\gamma^{{calib}}/\,$ms = {gamma_calib*1e3:.2f})')
ax.axvline(t_cross_fit, color='C2', linestyle=':', alpha=0.7,
           label=fr'$t_{{cross}}^{{fit}} = {t_cross_fit:.2f}\,\mu s$  ($\gamma^{{fit}}/\,$ms = {gamma_fit*1e3:.2f})')

ax.plot(t_fine, cpsi_calib, color='C1', linestyle='--', linewidth=1.2, alpha=0.6,
        label=r'F25 with $\gamma^{calib}$  (T2-echo)')
ax.plot(t_fine, cpsi_fitline, color='C2', linewidth=1.6,
        label=r'F25 with $\gamma^{fit}$  (in-situ)')

ax.errorbar(t, cpsi, yerr=sigma, fmt='o', color='black', markersize=5,
            capsize=2, elinewidth=0.8, label=f'Hardware ({backend}, qubits {pair})')

ax.set_xlabel(r'$t$  $[\mu s]$', fontsize=11)
ax.set_ylabel(r'$C\Psi$', fontsize=11)
ax.set_title(rf'Bell$^+$ trajectory through $C\Psi = 1/4$  on $\mathtt{{{backend}}}$'
             f'\nF25 fit: $\\gamma^{{fit}}/\\gamma^{{calib}}$ = {gamma_fit/gamma_calib:.2f}', fontsize=11)
ax.legend(loc='upper right', fontsize=8.5)
ax.grid(alpha=0.3)
ax.set_xlim(left=0)
ax.set_ylim(bottom=0)

out1 = Path(__file__).parent.parent / 'visualizations' / 'cusp_precision_trajectory.png'
fig.tight_layout()
fig.savefig(out1, dpi=150)
print(f"Saved: {out1}")

# ════════════════════════════════════════════════════════════════════
#  FIGURE 2: residuals
# ════════════════════════════════════════════════════════════════════

fig2, ax2 = plt.subplots(figsize=(8, 4.5))

residual_calib = cpsi - cpsi_bellplus(t, gamma_calib)
residual_fit = cpsi - cpsi_bellplus(t, gamma_fit)

ax2.axhline(0, color='gray', linewidth=0.8, alpha=0.6)
ax2.axhspan(-sigma, sigma, color='gray', alpha=0.15,
            label=fr'$\pm\sigma_{{shot}} \approx 0.5/\sqrt{{{shots}}} = {sigma:.4f}$')

ax2.plot(t, residual_calib, 'o--', color='C1', markersize=5,
         label='residual vs $\\gamma^{calib}$ (T2-echo prediction)')
ax2.plot(t, residual_fit, 'o-', color='C2', markersize=5,
         label='residual vs $\\gamma^{fit}$ (in-situ fit)')

ax2.axvline(t_cross_fit, color='C2', linestyle=':', alpha=0.5,
            label=fr'$t_{{cross}}^{{fit}}$')

ax2.set_xlabel(r'$t$  $[\mu s]$', fontsize=11)
ax2.set_ylabel(r'$C\Psi_{measured} - C\Psi_{F25}$', fontsize=11)
ax2.set_title(f'Pointwise residuals on {backend}, qubits {pair}', fontsize=11)
ax2.legend(loc='best', fontsize=9)
ax2.grid(alpha=0.3)

out2 = Path(__file__).parent.parent / 'visualizations' / 'cusp_precision_residuals.png'
fig2.tight_layout()
fig2.savefig(out2, dpi=150)
print(f"Saved: {out2}")

print(f"\n  γ_calib / ms = {gamma_calib*1e3:.4f}")
print(f"  γ_fit  / ms = {gamma_fit*1e3:.4f}")
print(f"  Ratio       = {gamma_fit/gamma_calib:.3f}")
print(f"  RMS fit res = {np.sqrt(np.mean(residual_fit**2)):.5f}")
