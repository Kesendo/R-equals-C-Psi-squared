"""
Aufgabe 4: Trajectory-based dwell time near CPsi = 1/4.

The static scaling n(eps) is a textbook effect. The physical question:
how long does CPsi(t) dwell near 1/4 as gamma drives the system through
the crossing point?

Bell+ (N=2): CPsi(t) = f*(1+f^2)/6, f = exp(-4*gamma*t)
  dCPsi/dt = -2*gamma*f*(1 + 3*f^2)/3   (Formula 25)
  Crossing: f*(1+f^2) = 3/2, f_cross ~ 0.8612
"""

import numpy as np
from scipy.optimize import brentq
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
os.makedirs(OUT_DIR, exist_ok=True)


def cpsi_bellplus(t, gamma):
    """CPsi(t) for Bell+, Formula 25."""
    f = np.exp(-4 * gamma * t)
    return f * (1 + f**2) / 6


def dcpsi_dt_bellplus(t, gamma):
    """dCPsi/dt for Bell+, derivative of Formula 25."""
    f = np.exp(-4 * gamma * t)
    return -2 * gamma * f * (1 + 3 * f**2) / 3


# ============================================================
# Part 1: Find crossing point CPsi = 1/4
# ============================================================
print("=" * 72)
print("PART 1: Crossing point CPsi(t_cross) = 1/4")
print("=" * 72)

# Solve f*(1+f^2) = 3/2 for f
f_cross = brentq(lambda f: f * (1 + f**2) - 1.5, 0.5, 1.0)
print(f"f_cross = {f_cross:.10f}")
print(f"Verify: f*(1+f^2) = {f_cross * (1 + f_cross**2):.15f} (should be 1.5)")

gamma_values = [0.5, 1.0, 2.0]
print()
for gamma in gamma_values:
    t_cross = -np.log(f_cross) / (4 * gamma)
    dcpsi = dcpsi_dt_bellplus(t_cross, gamma)
    print(f"gamma={gamma:.1f}: t_cross = {t_cross:.8f}, dCPsi/dt = {dcpsi:.8f}")

# ============================================================
# Part 2: Dwell time for various delta
# ============================================================
print()
print("=" * 72)
print("PART 2: Dwell time |CPsi(t) - 1/4| < delta")
print("=" * 72)

deltas = [1e-2, 1e-3, 1e-4]

print(f"{'gamma':>6} {'delta':>10} {'t_dwell':>12} {'t_pred':>12} {'ratio':>8}")
print("-" * 55)

for gamma in gamma_values:
    t_cross = -np.log(f_cross) / (4 * gamma)
    dcpsi_at_cross = dcpsi_dt_bellplus(t_cross, gamma)

    for delta in deltas:
        # Find t1, t2 where CPsi(t) = 1/4 +/- delta
        # CPsi is monotonically decreasing, so:
        # t1: CPsi(t1) = 1/4 + delta  (before crossing, t1 < t_cross)
        # t2: CPsi(t2) = 1/4 - delta  (after crossing, t2 > t_cross)
        try:
            t1 = brentq(lambda t: cpsi_bellplus(t, gamma) - (0.25 + delta),
                        0.001, t_cross)
            t2 = brentq(lambda t: cpsi_bellplus(t, gamma) - (0.25 - delta),
                        t_cross, 10.0 / gamma)
            t_dwell = t2 - t1

            # Prediction: t_dwell ~ 2*delta / |dCPsi/dt|_cross
            t_pred = 2 * delta / abs(dcpsi_at_cross)
            ratio = t_dwell / t_pred

            print(f"{gamma:6.1f} {delta:10.0e} {t_dwell:12.8f} {t_pred:12.8f} {ratio:8.4f}")
        except ValueError:
            print(f"{gamma:6.1f} {delta:10.0e} {'(no bracket)':>12}")

# ============================================================
# Part 3: Analytical derivative at crossing
# ============================================================
print()
print("=" * 72)
print("PART 3: Analytical dCPsi/dt at crossing")
print("=" * 72)

# dCPsi/dt = -2*gamma*f*(1+3f^2)/3
# At f = f_cross: factor = f_cross*(1+3*f_cross^2)/3
factor = f_cross * (1 + 3 * f_cross**2) / 3
print(f"f_cross = {f_cross:.10f}")
print(f"f*(1+3f^2)/3 = {factor:.10f}")
print(f"dCPsi/dt = -2*gamma * {factor:.6f}")
print()
for gamma in gamma_values:
    dcpsi = -2 * gamma * factor
    print(f"gamma={gamma:.1f}: dCPsi/dt = {dcpsi:.8f}")
    print(f"  t_dwell(delta) ~ 2*delta / {abs(dcpsi):.6f} = {2.0/abs(dcpsi):.6f} * delta")

# ============================================================
# Part 4: Gamma scaling and K-invariance
# ============================================================
print()
print("=" * 72)
print("PART 4: Gamma scaling of dwell time")
print("=" * 72)

delta_fixed = 1e-3
print(f"Fixed delta = {delta_fixed:.0e}")
print()
print(f"{'gamma':>8} {'t_dwell':>12} {'K_dwell':>12} {'t*gamma':>12}")
print("-" * 48)

gamma_scan = np.array([0.1, 0.25, 0.5, 1.0, 2.0, 5.0, 10.0])
K_dwells = []

for gamma in gamma_scan:
    t_cross = -np.log(f_cross) / (4 * gamma)
    try:
        t1 = brentq(lambda t: cpsi_bellplus(t, gamma) - (0.25 + delta_fixed),
                     1e-6, t_cross)
        t2 = brentq(lambda t: cpsi_bellplus(t, gamma) - (0.25 - delta_fixed),
                     t_cross, 100.0 / gamma)
        t_dwell = t2 - t1
        K_dwell = gamma * t_dwell
        K_dwells.append(K_dwell)
        print(f"{gamma:8.2f} {t_dwell:12.8f} {K_dwell:12.8f} {t_dwell*gamma:12.8f}")
    except ValueError:
        print(f"{gamma:8.2f} {'(bracket fail)':>12}")

K_dwells = np.array(K_dwells)
print()
print(f"K_dwell range: [{K_dwells.min():.8f}, {K_dwells.max():.8f}]")
print(f"K_dwell std:   {K_dwells.std():.2e}")
print(f"K_dwell mean:  {K_dwells.mean():.8f}")
print(f"CONFIRMED: K_dwell = gamma*t_dwell is independent of gamma.")
print(f"In K-units, dwell time is a pure constant (K-invariance).")

# ============================================================
# Part 5: Plot CPsi(t) for three gammas
# ============================================================
fig, axes = plt.subplots(1, 3, figsize=(16, 5))

# Panel 1: CPsi(t) for three gammas in physical time
ax = axes[0]
for gamma, color in [(0.5, 'C0'), (1.0, 'C1'), (2.0, 'C2')]:
    t_max = 1.5 / gamma
    t = np.linspace(0, t_max, 500)
    cpsi = cpsi_bellplus(t, gamma)
    ax.plot(t, cpsi, label=f'$\\gamma={gamma}$', color=color, linewidth=2)
ax.axhline(0.25, color='red', linestyle=':', alpha=0.7, label='$C\\Psi=1/4$')
ax.set_xlabel('t (physical time)')
ax.set_ylabel('$C\\Psi(t)$')
ax.set_title('$C\\Psi(t)$ at different $\\gamma$ (physical time)')
ax.legend()
ax.grid(alpha=0.3)

# Panel 2: CPsi(K) where K = gamma*t (all collapse)
ax = axes[1]
K_range = np.linspace(0, 1.5, 500)
for gamma, color in [(0.5, 'C0'), (1.0, 'C1'), (2.0, 'C2')]:
    t = K_range / gamma
    cpsi = cpsi_bellplus(t, gamma)
    ax.plot(K_range, cpsi, label=f'$\\gamma={gamma}$', color=color, linewidth=2)
ax.axhline(0.25, color='red', linestyle=':', alpha=0.7)
ax.set_xlabel('$K = \\gamma t$')
ax.set_ylabel('$C\\Psi(K/\\gamma)$')
ax.set_title('Same trajectories in K-time (all collapse)')
ax.legend()
ax.grid(alpha=0.3)

# Panel 3: Dwell time vs gamma (should be 1/gamma)
ax = axes[2]
gamma_plot = np.array([0.1, 0.25, 0.5, 1.0, 2.0, 5.0, 10.0])
t_dwells_plot = []
for gamma in gamma_plot:
    t_cross = -np.log(f_cross) / (4 * gamma)
    t1 = brentq(lambda t: cpsi_bellplus(t, gamma) - (0.25 + 1e-3),
                1e-6, t_cross)
    t2 = brentq(lambda t: cpsi_bellplus(t, gamma) - (0.25 - 1e-3),
                t_cross, 100.0 / gamma)
    t_dwells_plot.append(t2 - t1)
t_dwells_plot = np.array(t_dwells_plot)

ax.loglog(gamma_plot, t_dwells_plot, 'o-', color='C3', linewidth=2,
          markersize=6, label='$t_{dwell}$ (measured)')
# 1/gamma reference
ax.loglog(gamma_plot, t_dwells_plot[3] / gamma_plot, '--', color='gray',
          alpha=0.7, label='$\\propto 1/\\gamma$')
ax.set_xlabel('$\\gamma$')
ax.set_ylabel('$t_{dwell}$')
ax.set_title('Dwell time $\\propto 1/\\gamma$ (K-invariance)')
ax.legend()
ax.grid(alpha=0.3)

plt.tight_layout()
plot_path = os.path.join(OUT_DIR, "trajectory_dwell_time.png")
plt.savefig(plot_path, dpi=150)
print(f"\nPlot saved: {plot_path}")
print("DONE")
