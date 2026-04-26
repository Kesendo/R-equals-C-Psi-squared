#!/usr/bin/env python3
"""EQ-030 angular Winkel: θ = arctan(√(4·CΨ − 1)) trajectories under +T1.

Tom's framing:
  CΨ = Purity × Ψ-norm  (Purity = Tr(ρ²), Ψ-norm = L₁/(d−1))
  θ = arctan(√(4·CΨ − 1))  defined for CΨ > 1/4
  At |+−+⟩ initial: CΨ = 1, θ = arctan(√3) = 60°
  At CΨ = 1/4 boundary: θ = 0°

The protected-set drop under +T1 is an *algebraic* (static) measure. θ(t) is
a *state-trajectory* measure. Both are framework-grounded (Π-protection on
the operator side; CΨ = ¼ as the cusp from WHAT_WE_FOUND).

This script:
  1. Evolves ρ(t) for the 4 bond-flipped soft Hamiltonians under (a) pure-Z,
     (b) +T1.
  2. Adds truly J(XX+YY) and a fragile soft (XZ+XZ, drop=40) for reference.
  3. Computes CΨ(t) and θ(t) at each step.
  4. Reports: when does θ → 0 (CΨ crosses ¼)?

Hypothesis: T1-robust softs (XY+YX, IY+YI) keep θ above 0 for far longer
under +T1 than T1-fragile softs (YZ+ZY, XZ+ZX). This is the state-level
signature of the algebraic robustness.
"""
import math
import sys
from pathlib import Path

import numpy as np

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


def evolve_trajectory(L, rho_0, times):
    """Compute ρ(t) at each time via matrix exponential of L on vec(ρ)."""
    d = rho_0.shape[0]
    rho_vec = rho_0.reshape(-1)  # column-major (vec) flatten
    # vec(ρ) convention: vec[i + d*j] = ρ_ij in row-major Python indexing
    # but our framework's L uses Kronecker-product convention. Match it.
    # Standard: dρ/dt = L · vec(ρ) where vec(ρ)_{i+d*j} = ρ_{ij} (column-major)
    # numpy's reshape is row-major; convert via .T then .ravel().
    rho_vec = rho_0.T.reshape(-1).copy()
    trajectory = []
    for t in times:
        Ut = np.linalg.matrix_power if False else None  # no, use expm
        from scipy.linalg import expm
        rho_t_vec = expm(L * t) @ rho_vec
        rho_t = rho_t_vec.reshape(d, d).T  # invert the .T from before
        trajectory.append(rho_t)
    return trajectory


def main():
    N = 3
    GAMMA_DEPH = 0.1
    GAMMA_T1 = 0.01  # 0.1 * γ_deph
    J = 1.0
    bonds = [(i, i + 1) for i in range(N - 1)]

    plus = np.array([1, 1], dtype=complex) / math.sqrt(2)
    minus = np.array([1, -1], dtype=complex) / math.sqrt(2)
    psi = np.kron(plus, np.kron(minus, plus))
    rho_0 = np.outer(psi, psi.conj())

    cases = [
        ('truly XX+YY', [('X', 'X', J), ('Y', 'Y', J)]),
        ('IY+YI',       [('I', 'Y', J), ('Y', 'I', J)]),
        ('XY+YX',       [('X', 'Y', J), ('Y', 'X', J)]),
        ('YZ+ZY',       [('Y', 'Z', J), ('Z', 'Y', J)]),
        ('XZ+ZX',       [('X', 'Z', J), ('Z', 'X', J)]),
        ('XZ+XZ',       [('X', 'Z', J), ('X', 'Z', J)]),
    ]

    t_max = 30.0
    dt = 0.02
    n_steps = int(t_max / dt)
    times = np.linspace(0, t_max, n_steps + 1)

    # Initial sanity
    p0 = purity(rho_0)
    psi0 = psi_norm(rho_0)
    c0 = cpsi(rho_0)
    theta0 = theta_from_cpsi(c0)
    print(f"Initial state |+−+⟩, N={N}, d={2**N}")
    print(f"  Purity = {p0:.4f}, Ψ-norm = {psi0:.4f}, "
          f"CΨ = {c0:.4f}, θ = {theta0:.2f}°")
    print(f"  (CΨ = 1, θ = 60°)  — coherent product state")
    print()

    print(f"Evolution: γ_deph = {GAMMA_DEPH}, γ_T1 = {GAMMA_T1}, "
          f"t_max = {t_max}, |+−+⟩")
    print()
    print(f"  For each H, report θ(t) trajectory at sampled times under (a)")
    print(f"  pure-Z and (b) +T1. Mark: t* where CΨ crosses 1/4 (θ → 0).")
    print()

    print(f"{'case':<14s}  {'channel':<8s}  "
          f"{'θ_max':>8s}  {'t*_first':>9s}  {'t*_last':>9s}  "
          f"{'∫θ dt':>10s}  {'#crossings':>11s}")
    print('-' * 90)

    summary = []
    for label, terms in cases:
        H = fw._build_bilinear(N, bonds, terms)
        L_pureZ = fw.lindbladian_z_dephasing(H, [GAMMA_DEPH] * N)
        L_T1ext = fw.lindbladian_z_plus_t1(
            H, [GAMMA_DEPH] * N, [GAMMA_T1] * N
        )

        for ch_label, L in [('pure-Z', L_pureZ), ('+T1', L_T1ext)]:
            traj = evolve_trajectory(L, rho_0, times)
            theta_t = np.array([theta_from_cpsi(cpsi(r)) for r in traj])
            cpsi_t = np.array([cpsi(r) for r in traj])

            theta_max = float(theta_t.max())

            # First crossing: when CΨ first drops to ≤ 1/4
            below = cpsi_t <= 0.25
            t_first = float(times[np.argmax(below)]) if below.any() else None

            # Last instance of being above 1/4 (after first crossing)
            above = cpsi_t > 0.25
            if above.any():
                last_above_idx = len(times) - 1 - int(np.argmax(above[::-1]))
                t_last = float(times[last_above_idx])
            else:
                t_last = None

            # Number of zero-crossings of (CΨ - 1/4)
            sign_seq = np.sign(cpsi_t - 0.25)
            crossings = int(np.sum(np.abs(np.diff(sign_seq)) > 0.5))

            # Integrated dwell ∫θ dt (degrees · time units)
            integrated = float(np.trapezoid(np.maximum(theta_t, 0.0), times))

            t_first_str = f"{t_first:.2f}" if t_first is not None else "—"
            t_last_str = f"{t_last:.2f}" if t_last is not None else "—"
            print(f"{label:<14s}  {ch_label:<8s}  "
                  f"{theta_max:>7.2f}°  {t_first_str:>9s}  {t_last_str:>9s}  "
                  f"{integrated:>10.2f}  {crossings:>11d}")

            summary.append({
                'label': label,
                'channel': ch_label,
                'theta_max': theta_max,
                't_first': t_first,
                't_last': t_last,
                'integrated': integrated,
                'crossings': crossings,
                'theta_t': theta_t,
                'cpsi_t': cpsi_t,
            })
        print()

    print()
    print("Summary: pure-Z vs +T1 ∫θ dt comparison (cusp dwell loss under T1)")
    print()
    print(f"  {'case':<14s}  {'∫θ pure-Z':>11s}  {'∫θ +T1':>11s}  "
          f"{'Δ∫θ':>9s}  {'rel %':>7s}  {'algebraic drop':>14s}")
    drops_for_label = {
        'truly XX+YY': 31,  # truly drops 32→1 = 31
        'IY+YI': 0,
        'XY+YX': 1,
        'YZ+ZY': 28,
        'XZ+ZX': 29,
        'XZ+XZ': 40,
    }
    rows_for_corr = []
    for label, _ in cases:
        s_pure = next((s for s in summary if s['label'] == label and s['channel'] == 'pure-Z'), None)
        s_t1 = next((s for s in summary if s['label'] == label and s['channel'] == '+T1'), None)
        if s_pure is None or s_t1 is None:
            continue
        ip = s_pure['integrated']
        it = s_t1['integrated']
        delta = it - ip
        rel = (delta / ip * 100) if ip > 0 else 0.0
        drop = drops_for_label.get(label, '—')
        rows_for_corr.append((drop, delta) if isinstance(drop, int) else None)
        drop_str = str(drop) if isinstance(drop, int) else drop
        print(f"  {label:<14s}  {ip:>11.2f}  {it:>11.2f}  "
              f"{delta:>+9.2f}  {rel:>+6.1f}%  {drop_str:>14s}")

    valid = [r for r in rows_for_corr if r is not None]
    if len(valid) >= 3:
        drops_arr = np.array([r[0] for r in valid])
        deltas_arr = np.array([r[1] for r in valid])
        if drops_arr.std() > 0 and deltas_arr.std() > 0:
            rho = float(np.corrcoef(drops_arr, deltas_arr)[0, 1])
            print()
            print(f"  Pearson(algebraic drop, Δ∫θ) = {rho:+.4f}")
            print(f"  (negative → T1 closes the cusp more for cases with bigger drops)")


if __name__ == "__main__":
    main()
