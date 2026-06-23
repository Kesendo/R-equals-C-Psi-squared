#!/usr/bin/env python3
"""EQ-030 last try: Bures velocity v_B(t) at the cusp boundary.

Tom's instinct: the rest lies in the in-between. Probably true. But let's
check whether Bures velocity at t* (the boundary crossing) tells us anything
independent of the θ-descent shape.

Bures distance between density matrices:
    D_B(ρ_1, ρ_2)² = 2 · (1 − F(ρ_1, ρ_2))
    F = Tr(√(√ρ_1 · ρ_2 · √ρ_1))²       (Uhlmann fidelity)

Bures velocity:
    v_B(t) = lim_{dt→0} D_B(ρ(t), ρ(t+dt)) / dt

For each soft case under L_Z and L_Z + T1:
  1. Evolve ρ(t) at fine resolution.
  2. Compute v_B(t) by finite difference.
  3. Report v_B at the boundary crossing t* and in a window around it.

If v_B is just |dθ/dt| in different units, this adds nothing. If v_B
captures a separate motion-direction in state space (orthogonal to the
CΨ direction), the discrimination might sharpen further.
"""
import math
import sys
from pathlib import Path

import numpy as np
from scipy.linalg import expm, sqrtm

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass

SCRIPT_DIR = Path(__file__).parent
sys.path.insert(0, str(SCRIPT_DIR))
import framework as fw


def purity(rho):
    return float(np.real(np.trace(rho @ rho)))


def psi_norm(rho):
    d = rho.shape[0]
    l1 = float(np.sum(np.abs(rho)) - np.sum(np.abs(np.diag(rho))))
    return l1 / (d - 1)


def cpsi(rho):
    return purity(rho) * psi_norm(rho)


def theta_from_cpsi(c):
    if c <= 0.25:
        return 0.0
    return math.degrees(math.atan(math.sqrt(4 * c - 1)))


def bures_distance(rho_a, rho_b):
    """D_B(ρ_a, ρ_b)² = 2(1 - F), F = Tr(√(√ρ_a ρ_b √ρ_a))²."""
    sqrt_a = sqrtm(rho_a)
    inner = sqrt_a @ rho_b @ sqrt_a
    sqrt_inner = sqrtm(inner)
    fidelity = float(np.real(np.trace(sqrt_inner)) ** 2)
    fidelity = max(0.0, min(1.0, fidelity))
    return math.sqrt(max(0.0, 2.0 * (1.0 - fidelity)))


def evolve(L, rho_0, times):
    d = rho_0.shape[0]
    rho_vec = rho_0.T.reshape(-1).copy()
    out = []
    for t in times:
        rho_t_vec = expm(L * t) @ rho_vec
        out.append(rho_t_vec.reshape(d, d).T)
    return out


def bures_velocity_trajectory(traj, times):
    """v_B(t_i) ≈ D_B(ρ_i, ρ_{i+1}) / Δt, evaluated at midpoints."""
    v = []
    for i in range(len(traj) - 1):
        D = bures_distance(traj[i], traj[i + 1])
        dt = times[i + 1] - times[i]
        v.append(D / dt if dt > 0 else 0.0)
    return np.array(v)


def find_first_cross_idx(cpsi_t, threshold=0.25):
    below = cpsi_t <= threshold
    if not below.any():
        return None
    return int(np.argmax(below))


def find_last_above_idx(cpsi_t, threshold=0.25):
    above = cpsi_t > threshold
    if not above.any():
        return None
    return len(cpsi_t) - 1 - int(np.argmax(above[::-1]))


def main():
    N = 3
    GAMMA_DEPH = 0.1
    GAMMA_T1 = 0.01
    J = 1.0
    bonds = [(i, i + 1) for i in range(N - 1)]

    plus = np.array([1, 1], dtype=complex) / math.sqrt(2)
    minus = np.array([1, -1], dtype=complex) / math.sqrt(2)
    psi = np.kron(plus, np.kron(minus, plus))
    rho_0 = np.outer(psi, psi.conj())

    cases = [
        ('truly XX+YY', [('X', 'X', J), ('Y', 'Y', J)], 31),
        ('IY+YI',       [('I', 'Y', J), ('Y', 'I', J)], 0),
        ('XY+YX',       [('X', 'Y', J), ('Y', 'X', J)], 1),
        ('YZ+ZY',       [('Y', 'Z', J), ('Z', 'Y', J)], 28),
        ('XZ+ZX',       [('X', 'Z', J), ('Z', 'X', J)], 29),
        ('XZ+XZ',       [('X', 'Z', J), ('X', 'Z', J)], 40),
    ]

    t_max = 6.0  # near the early cusp window only — saves time
    dt = 0.01
    n_steps = int(t_max / dt)
    times = np.linspace(0, t_max, n_steps + 1)
    mid_times = (times[:-1] + times[1:]) / 2

    print(f"Bures-velocity-at-cusp analysis")
    print(f"  N={N}, |+−+⟩, γ_deph={GAMMA_DEPH}, γ_T1={GAMMA_T1}, "
          f"dt={dt}, t_max={t_max}")
    print()

    print(f"  {'case':<14s}  {'drop':>4s}  {'channel':<8s}  "
          f"{'t* (last)':>10s}  {'v_B@t*':>10s}  "
          f"{'v_B mean[t*-0.3, t*]':>22s}  {'|dθ/dt|@t*':>12s}")
    print('-' * 100)
    rows = []

    for label, terms, drop in cases:
        H = fw._build_bilinear(N, bonds, terms)
        for ch_label, gamma_t1_val in [('pure-Z', 0.0), ('+T1', GAMMA_T1)]:
            if gamma_t1_val > 0:
                L = fw.lindbladian_z_plus_t1(
                    H, [GAMMA_DEPH] * N, [gamma_t1_val] * N
                )
            else:
                L = fw.lindbladian_z_dephasing(H, [GAMMA_DEPH] * N)
            traj = evolve(L, rho_0, times)
            cpsi_t = np.array([cpsi(r) for r in traj])
            theta_t = np.array([theta_from_cpsi(c) for c in cpsi_t])
            v_B = bures_velocity_trajectory(traj, times)

            # Last index above the boundary
            la_idx = find_last_above_idx(cpsi_t)
            if la_idx is None or la_idx == len(times) - 1:
                rows.append({
                    'label': label, 'drop': drop, 'ch': ch_label,
                    't_cross': None, 'vB_at_cross': None,
                    'vB_mean_window': None, 'dtheta_dt_at_cross': None,
                })
                print(f"  {label:<14s}  {drop:>4d}  {ch_label:<8s}  "
                      f"{'never':>10s}  {'—':>10s}  {'—':>22s}  {'—':>12s}")
                continue

            t_cross = times[la_idx]
            # v_B at the crossing (the v_B element straddling the boundary)
            vB_at_cross = float(v_B[la_idx]) if la_idx < len(v_B) else None

            # v_B mean in [t_cross - 0.3, t_cross]
            window_mask = (mid_times >= t_cross - 0.3) & (mid_times <= t_cross)
            vB_mean_window = float(np.mean(v_B[window_mask])) if window_mask.any() else None

            # |dθ/dt| at the crossing (degrees/unit-time)
            if la_idx > 0:
                dtheta = (theta_t[la_idx + 1] - theta_t[la_idx - 1]) / (2 * dt)
                dtheta_at_cross = float(abs(dtheta))
            else:
                dtheta_at_cross = None

            rows.append({
                'label': label, 'drop': drop, 'ch': ch_label,
                't_cross': t_cross, 'vB_at_cross': vB_at_cross,
                'vB_mean_window': vB_mean_window,
                'dtheta_dt_at_cross': dtheta_at_cross,
            })
            tcs = f"{t_cross:.3f}"
            vBs = f"{vB_at_cross:.4f}" if vB_at_cross is not None else "—"
            vBws = f"{vB_mean_window:.4f}" if vB_mean_window is not None else "—"
            dts = f"{dtheta_at_cross:.2f}" if dtheta_at_cross is not None else "—"
            print(f"  {label:<14s}  {drop:>4d}  {ch_label:<8s}  "
                  f"{tcs:>10s}  {vBs:>10s}  {vBws:>22s}  {dts:>12s}")
        print()

    # Pearson(drop, vB_at_cross) under +T1
    t1_rows = [r for r in rows if r['ch'] == '+T1' and r['vB_at_cross'] is not None]
    if len(t1_rows) >= 3:
        drops = np.array([r['drop'] for r in t1_rows])
        vBs = np.array([r['vB_at_cross'] for r in t1_rows])
        if drops.std() > 0 and vBs.std() > 0:
            rho = float(np.corrcoef(drops, vBs)[0, 1])
            print()
            print(f"  Pearson(drop, v_B at cusp [+T1]) = {rho:+.4f}")

        vBws = np.array([r['vB_mean_window'] for r in t1_rows
                         if r['vB_mean_window'] is not None])
        drops_w = np.array([r['drop'] for r in t1_rows
                            if r['vB_mean_window'] is not None])
        if len(drops_w) >= 3 and drops_w.std() > 0 and vBws.std() > 0:
            rho_w = float(np.corrcoef(drops_w, vBws)[0, 1])
            print(f"  Pearson(drop, mean v_B in [t*-0.3, t*] [+T1]) = {rho_w:+.4f}")

    print()
    print("Reading guide:")
    print("  v_B at t* = state-space speed at the moment of crossing CΨ=1/4.")
    print("  - Small v_B = state lingers near the boundary (long fragile tail).")
    print("  - Large v_B = state crosses briskly (no tail).")
    print("  If v_B and |dθ/dt|/2 differ significantly, v_B captures motion")
    print("  outside the CΨ-direction (orthogonal in state space).")


if __name__ == "__main__":
    main()
