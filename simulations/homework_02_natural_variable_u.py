"""
Homework #2: Die natuerliche Variable u(t) = C(Psi + R)
Beide Varianten testen.

Variante (a) statisch: Fuer feste (C, Psi), iteriere u_{n+1} = u_n^2 + c
  mit c = C*Psi, u_0 = c. Beobachte Konvergenz/Divergenz.

Variante (b) dynamisch: Entlang Bell+ unter Z-Dephasing,
  CPsi(t) = f*(1+f^2)/6, f = exp(-4*gamma*t).
  Zu jedem t: Iteriere bis Konvergenz, extrahiere R_inf(t).
  Plotte R(t) gegen t UND gegen xi = ln(Psi).
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
os.makedirs(OUT_DIR, exist_ok=True)


def iterate_mandelbrot(c, max_iter=2000, bailout=1e6, tol=1e-14):
    """Iteriere u_{n+1} = u_n^2 + c mit u_0 = c."""
    u = complex(c)
    for n in range(max_iter):
        u_new = u * u + c
        if abs(u_new) > bailout:
            return False, u_new, n
        if abs(u_new - u) < tol:
            return True, u_new, n
        u = u_new
    return False, u, max_iter


print("=" * 72)
print("VARIANTE (a): Mandelbrot-Dynamik bei festen (C, Psi)")
print("=" * 72)

test_points = [
    (1.0, 0.10, "weit unter Rand"),
    (1.0, 0.20, "nahe Rand"),
    (1.0, 0.24, "knapp unter"),
    (1.0, 0.249, "ultra-kritisch"),
    (1.0, 0.25, "genau am Rand"),
    (1.0, 0.26, "knapp ueber"),
    (1.0, 0.30, "ueber Rand"),
    (0.875, 0.289, "Bell+ Kreuzung"),
    (0.625, 0.167, "Bell+ nach Kreuzung"),
]

print(f"{'C':>6} {'Psi':>6} {'CPsi':>8} {'conv':>5} {'n':>5} {'R_iter':>14} {'R_analytic':>14} {'diff':>10}  label")
print("-" * 100)
for C, Psi, label in test_points:
    c = C * Psi
    converged, u_inf, n = iterate_mandelbrot(c)
    if converged:
        R_iter = u_inf.real / C - Psi
        if c <= 0.25:
            disc = 1 - 4 * c
            R_an = (1 - 2 * c - np.sqrt(disc)) / (2 * C)
            diff = abs(R_iter - R_an)
            print(f"{C:6.3f} {Psi:6.3f} {c:8.5f} {'yes':>5} {n:5d} {R_iter:14.8f} {R_an:14.8f} {diff:10.2e}  {label}")
        else:
            print(f"{C:6.3f} {Psi:6.3f} {c:8.5f} {'yes':>5} {n:5d} {R_iter:14.8f} {'complex':>14} {'-':>10}  {label}")
    else:
        print(f"{C:6.3f} {Psi:6.3f} {c:8.5f} {'NO':>5} {n:5d} {'DIVERGED':>14} {'-':>14} {'-':>10}  {label}")

print()
print("=" * 72)
print("VARIANTE (b): R(t) entlang Bell+ Dekohaerenzpfad")
print("=" * 72)

# Bell+ unter Z-Dephasing, N=2:
#   f(t) = exp(-4*gamma*t)
#   rho(t) hat Diagonal-Werte 1/2 und Off-Diagonal f/2
#   C(t) = Tr(rho^2) = (1 + f^2)/2
#   Psi(t) = L1/(d-1) = f/3 (d=4, L1 = |f|)
#   CPsi(t) = f*(1+f^2)/6
#
# Das ist Formel 25 aus ANALYTICAL_FORMULAS.md

gamma = 1.0
t_array = np.linspace(0.001, 3.0, 300)  # vermeide t=0 wegen log(Psi)
f_array = np.exp(-4 * gamma * t_array)

C_t = (1 + f_array**2) / 2
Psi_t = f_array / 3
CPsi_t = C_t * Psi_t  # = f*(1+f^2)/6
xi_t = np.log(Psi_t)

# Fuer jeden Zeitpunkt: Mandelbrot iterieren, R_inf extrahieren
R_t = np.zeros_like(t_array)
converged_flags = np.zeros_like(t_array, dtype=bool)
iter_counts = np.zeros_like(t_array, dtype=int)

for i, (C, Psi, cpsi) in enumerate(zip(C_t, Psi_t, CPsi_t)):
    c = cpsi  # c = CPsi fuer Mandelbrot-Map
    converged, u_inf, n = iterate_mandelbrot(c)
    converged_flags[i] = converged
    iter_counts[i] = n
    if converged:
        R_t[i] = u_inf.real / C - Psi
    else:
        R_t[i] = np.nan

# Analytische R_- Kurve zum Vergleich (wo c <= 0.25)
R_analytic = np.full_like(t_array, np.nan)
mask_conv = CPsi_t <= 0.25
disc = 1 - 4 * CPsi_t[mask_conv]
R_analytic[mask_conv] = (1 - 2 * CPsi_t[mask_conv] - np.sqrt(disc)) / (2 * C_t[mask_conv])

# u(t) selbst, die natuerliche Variable
u_t = C_t * (Psi_t + R_t)

print(f"Zeitraster: {len(t_array)} Punkte, t in [{t_array[0]:.3f}, {t_array[-1]:.3f}]")
print(f"CPsi(t=0)  = {CPsi_t[0]:.6f}   (Startwert, nahe Maximum 1/6)")
print(f"CPsi(t=end)= {CPsi_t[-1]:.6e}   (sehr klein)")
print(f"Konvergenz: {converged_flags.sum()}/{len(t_array)} Iterationen")
print(f"max Iter-Count: {iter_counts.max()}, min: {iter_counts.min()}, median: {int(np.median(iter_counts))}")

# Finde Kreuzungspunkt CPsi=1/4
idx_cross = None
for i in range(len(CPsi_t) - 1):
    if CPsi_t[i] >= 0.25 >= CPsi_t[i+1] or CPsi_t[i] <= 0.25 <= CPsi_t[i+1]:
        idx_cross = i
        break

if idx_cross is not None:
    print(f"CPsi = 1/4 Kreuzung bei t ~ {t_array[idx_cross]:.4f}")
else:
    print(f"CPsi(0) = {CPsi_t[0]:.6f} < 1/4  (keine Kreuzung, startet bereits unter Rand)")
    print("  Bell+ N=2 startet bei CPsi_max = 1/6 = 0.1667")
    print("  Das ist UNTER 1/4, also sind wir immer im konvergenten Regime.")

# Diff iterated vs analytic
if mask_conv.any():
    diff = np.abs(R_t[mask_conv] - R_analytic[mask_conv])
    print(f"max |R_iter - R_analytic| (conv. Regime): {np.nanmax(diff):.2e}")
    print(f"mean |R_iter - R_analytic|:              {np.nanmean(diff):.2e}")

# ============================================================
# Plots
# ============================================================
fig, axes = plt.subplots(2, 2, figsize=(12, 9))

ax = axes[0, 0]
ax.plot(t_array, C_t, label='C(t) = purity', color='C0')
ax.plot(t_array, Psi_t, label='Psi(t) = L1/(d-1)', color='C1')
ax.plot(t_array, CPsi_t, label='CPsi(t)', color='C2', linewidth=2)
ax.axhline(0.25, color='red', linestyle='--', alpha=0.5, label='CPsi = 1/4')
ax.set_xlabel('t')
ax.set_ylabel('value')
ax.set_title('Bell+ decoherence path (Z-dephasing, gamma=1)')
ax.legend()
ax.grid(alpha=0.3)

ax = axes[0, 1]
ax.plot(t_array, R_t, label='R(t) iterated', color='C3', linewidth=2)
ax.plot(t_array, R_analytic, label='R_-(t) analytic', color='black', linestyle=':', linewidth=2)
ax.set_xlabel('t')
ax.set_ylabel('R')
ax.set_title('R(t) = Mandelbrot fixed point')
ax.legend()
ax.grid(alpha=0.3)

ax = axes[1, 0]
ax.plot(xi_t, R_t, label='R(xi)', color='C3', linewidth=2)
ax.set_xlabel('xi = ln(Psi)')
ax.set_ylabel('R')
ax.set_title('R vs decoherence clock xi = ln(Psi)')
ax.legend()
ax.grid(alpha=0.3)

ax = axes[1, 1]
ax.plot(t_array, u_t, label='u(t) = C(Psi+R)', color='C4', linewidth=2)
ax.plot(t_array, CPsi_t, label='CPsi(t) = c', color='C2', linestyle='--')
ax.set_xlabel('t')
ax.set_ylabel('u, c')
ax.set_title('Natural variable u(t) vs c(t)=CPsi(t)')
ax.legend()
ax.grid(alpha=0.3)

plt.tight_layout()
out_path = os.path.join(OUT_DIR, "homework_02_natural_variable_u.png")
plt.savefig(out_path, dpi=120)
print(f"\nPlot saved: {out_path}")

# ============================================================
# Die wichtigste Auswertung: verhaelt sich u(t) linear in xi?
# ============================================================
print()
print("=" * 72)
print("Struktur-Check: u(t) vs xi(t)")
print("=" * 72)

# Linearer Fit u vs xi
valid = ~np.isnan(u_t)
slope, intercept = np.polyfit(xi_t[valid], u_t[valid], 1)
u_fit = slope * xi_t + intercept
residuals = u_t - u_fit
rms = np.sqrt(np.nanmean(residuals**2))
print(f"Linear fit u = a*xi + b:  a = {slope:.6f}, b = {intercept:.6f}")
print(f"RMS residual: {rms:.4e}  (relative to |u| range {np.nanmax(u_t)-np.nanmin(u_t):.4f})")

# Exponentieller Fit: u ~ A*exp(B*xi)?
log_u_valid = u_t > 0
slope_log, intercept_log = np.polyfit(xi_t[log_u_valid], np.log(u_t[log_u_valid]), 1)
print(f"Log fit ln(u) = a*xi + b: a = {slope_log:.6f}, b = {intercept_log:.6f}")
print(f"  -> u ~ {np.exp(intercept_log):.4f} * Psi^{slope_log:.4f}")

# Direkter Vergleich: u/Psi und u/CPsi als Funktion von t
ratio_u_psi = u_t / Psi_t
ratio_u_cpsi = u_t / CPsi_t
print()
print(f"u/Psi:   min={np.nanmin(ratio_u_psi):.4f}, max={np.nanmax(ratio_u_psi):.4f}")
print(f"u/CPsi:  min={np.nanmin(ratio_u_cpsi):.4f}, max={np.nanmax(ratio_u_cpsi):.4f}")
print()
print("DONE")
