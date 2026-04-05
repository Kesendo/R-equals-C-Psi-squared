"""
Aufgabe 5: Mandelbrot trajectory overlay.

Plot the Bell+ decoherence trajectory as a curve on the Mandelbrot set.
The trajectory is a 1D path along the real axis in c-space:
  c(t) = CPsi(t) = f*(1+f^2)/6, f = exp(-4*gamma*t)
  from c = 1/3 (outside cardioid, divergence) through c = 1/4 (cusp)
  to c = 0 (center, deep classical regime).
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.collections import LineCollection
import os

VIS_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       "visualizations")
os.makedirs(VIS_DIR, exist_ok=True)


def mandelbrot_escape(c, max_iter=200):
    """Return escape iteration count for Mandelbrot set."""
    z = 0 + 0j
    for i in range(max_iter):
        z = z * z + c
        if abs(z) > 2:
            return i
    return max_iter


# ============================================================
# Compute Mandelbrot set
# ============================================================
print("Computing Mandelbrot set...")
nx, ny = 1200, 900
x_range = (-2.2, 0.8)
y_range = (-1.2, 1.2)

x = np.linspace(x_range[0], x_range[1], nx)
y = np.linspace(y_range[0], y_range[1], ny)
X, Y = np.meshgrid(x, y)
C = X + 1j * Y

# Vectorized escape time
escape = np.zeros_like(X, dtype=int)
Z = np.zeros_like(C)
mask = np.ones_like(X, dtype=bool)

max_iter = 200
for i in range(max_iter):
    Z[mask] = Z[mask] ** 2 + C[mask]
    escaped = mask & (np.abs(Z) > 2)
    escape[escaped] = i
    mask = mask & ~escaped

escape[mask] = max_iter  # points that never escaped

# ============================================================
# Bell+ trajectory along real axis
# ============================================================
print("Computing Bell+ trajectory...")
gamma = 1.0
t_traj = np.linspace(0, 5.0, 2000)
f_traj = np.exp(-4 * gamma * t_traj)
cpsi_traj = f_traj * (1 + f_traj**2) / 6

# Trajectory is along real axis: Im(c) = 0
# c goes from ~1/3 (t=0) to ~0 (t=inf)

# ============================================================
# Plot
# ============================================================
fig, ax = plt.subplots(figsize=(14, 10))

# Mandelbrot set: black interior, smooth exterior
# Use log coloring for exterior
log_escape = np.log(escape.astype(float) + 1)
ax.imshow(log_escape, extent=[x_range[0], x_range[1], y_range[0], y_range[1]],
          cmap='bone', origin='lower', aspect='equal', interpolation='bilinear')

# Overlay: Bell+ trajectory as colored line (color = time)
points = np.array([cpsi_traj, np.zeros_like(cpsi_traj)]).T.reshape(-1, 1, 2)
segments = np.concatenate([points[:-1], points[1:]], axis=1)
norm = plt.Normalize(t_traj.min(), t_traj.max())
lc = LineCollection(segments, cmap='plasma', norm=norm, linewidths=3.5, zorder=5)
lc.set_array(t_traj[:-1])
line = ax.add_collection(lc)
cbar = fig.colorbar(line, ax=ax, label='Time t ($\\gamma=1$)', shrink=0.6,
                    pad=0.02)

# Mark crossing point at (1/4, 0)
ax.plot(0.25, 0, 'x', color='lime', markersize=14, markeredgewidth=3,
        zorder=10, label='$C\\Psi = 1/4$ (cusp)')

# Mark start and end
ax.plot(cpsi_traj[0], 0, 'o', color='red', markersize=8, zorder=10,
        label=f'$t=0$: $C\\Psi = {cpsi_traj[0]:.3f}$')
ax.plot(cpsi_traj[-1], 0, 's', color='cyan', markersize=8, zorder=10,
        label=f'$t\\to\\infty$: $C\\Psi \\to 0$')

# Annotation for the cardioid cusp
ax.annotate('Cardioid cusp\n$c = 1/4$',
            xy=(0.25, 0), xytext=(0.5, 0.5),
            fontsize=11, color='lime',
            arrowprops=dict(arrowstyle='->', color='lime', lw=1.5),
            zorder=11)

# Zone labels
ax.text(-0.5, -0.9, 'Mandelbrot Set\n(bound orbits)', fontsize=12,
        color='white', ha='center', style='italic')
ax.text(0.15, 0.15, 'Classical\nzone', fontsize=10, color='yellow',
        ha='center', alpha=0.8)
ax.text(0.32, 0.15, 'Divergent\nzone', fontsize=10, color='orange',
        ha='center', alpha=0.8)

ax.set_xlabel('Re(c)', fontsize=12)
ax.set_ylabel('Im(c)', fontsize=12)
ax.set_title('Bell+ decoherence trajectory on the Mandelbrot set',
             fontsize=14, fontweight='bold')
ax.legend(loc='upper left', fontsize=10)

plt.tight_layout()
out_path = os.path.join(VIS_DIR, "bellplus_trajectory_on_mandelbrot.png")
plt.savefig(out_path, dpi=150)
print(f"Saved: {out_path}")

# ============================================================
# Zoomed view near the cusp
# ============================================================
fig2, ax2 = plt.subplots(figsize=(10, 8))

# Recompute at higher resolution near cusp
nx2, ny2 = 800, 600
x2_range = (-0.1, 0.5)
y2_range = (-0.3, 0.3)

x2 = np.linspace(x2_range[0], x2_range[1], nx2)
y2 = np.linspace(y2_range[0], y2_range[1], ny2)
X2, Y2 = np.meshgrid(x2, y2)
C2 = X2 + 1j * Y2

escape2 = np.zeros_like(X2, dtype=int)
Z2 = np.zeros_like(C2)
mask2 = np.ones_like(X2, dtype=bool)

for i in range(max_iter):
    Z2[mask2] = Z2[mask2] ** 2 + C2[mask2]
    escaped2 = mask2 & (np.abs(Z2) > 2)
    escape2[escaped2] = i
    mask2 = mask2 & ~escaped2
escape2[mask2] = max_iter

log_escape2 = np.log(escape2.astype(float) + 1)
ax2.imshow(log_escape2,
           extent=[x2_range[0], x2_range[1], y2_range[0], y2_range[1]],
           cmap='bone', origin='lower', aspect='equal', interpolation='bilinear')

# Trajectory (zoomed)
mask_zoom = (cpsi_traj >= x2_range[0]) & (cpsi_traj <= x2_range[1])
t_zoom = t_traj[mask_zoom]
c_zoom = cpsi_traj[mask_zoom]

points2 = np.array([c_zoom, np.zeros_like(c_zoom)]).T.reshape(-1, 1, 2)
segments2 = np.concatenate([points2[:-1], points2[1:]], axis=1)
lc2 = LineCollection(segments2, cmap='plasma', norm=norm, linewidths=4, zorder=5)
lc2.set_array(t_zoom[:-1])
ax2.add_collection(lc2)
fig2.colorbar(lc2, ax=ax2, label='Time t', shrink=0.6)

ax2.plot(0.25, 0, 'x', color='lime', markersize=18, markeredgewidth=4, zorder=10)
ax2.plot(1/3, 0, 'o', color='red', markersize=10, zorder=10,
         label='$t=0$: $C\\Psi=1/3$')

# Cardioid boundary (parametric)
theta = np.linspace(0, 2 * np.pi, 1000)
c_cardioid = 0.5 * np.exp(1j * theta) - 0.25 * np.exp(2j * theta)
ax2.plot(c_cardioid.real, c_cardioid.imag, '--', color='lime', alpha=0.4,
         linewidth=1, label='Cardioid boundary')

ax2.set_xlabel('Re(c)', fontsize=12)
ax2.set_ylabel('Im(c)', fontsize=12)
ax2.set_title('Zoomed: Bell+ trajectory through the cusp at $c = 1/4$',
              fontsize=13, fontweight='bold')
ax2.legend(loc='upper left')

plt.tight_layout()
out_path2 = os.path.join(VIS_DIR, "bellplus_trajectory_on_mandelbrot_zoom.png")
plt.savefig(out_path2, dpi=150)
print(f"Saved: {out_path2}")

print("DONE")
