#!/usr/bin/env python3
"""Multi-lens tour extended: PTF α_i + closure, plus Universal Carrier readings.

PTF baseline: L_A = truly XX+YY chain (canonical XY chain). For each other case,
L_B = canonical L_A with the case's H replacing the bilinear. Per-site purity
trajectory under L_A vs L_B → fit α_i such that P_B(i, t) ≈ P_A(i, α_i·t).
Closure check: Σ_i ln(α_i) ≈ 0 (PTF prediction, Tier 2 in perturbative window).

Universal Carrier constants (Q = J/γ₀, t_peak = 1/(4γ₀), AbsorptionQuantum = 2γ₀)
are chain-level invariants — printed once.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import minimize_scalar

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass

SCRIPT_DIR = Path(__file__).parent
sys.path.insert(0, str(SCRIPT_DIR))
import framework as fw
from framework.lindblad import lindbladian_z_dephasing, lindbladian_z_plus_t1


def per_site_purity_trajectory(L, rho_0, N, t_grid):
    """Propagate ρ_0 under L, sample per-site purity at t_grid.

    P_i(t) = ½(1 + <X_i>² + <Y_i>² + <Z_i>²).
    """
    d = 2 ** N
    evals, R = np.linalg.eig(L)
    R_inv = np.linalg.inv(R)
    rho_vec = rho_0.flatten()
    c0 = R_inv @ rho_vec
    # site Paulis (operators on full d×d space)
    paulis = {}
    sx = np.array([[0, 1], [1, 0]], dtype=complex)
    sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
    sz = np.array([[1, 0], [0, -1]], dtype=complex)
    eye = np.eye(2, dtype=complex)
    for i in range(N):
        for label, op in [('X', sx), ('Y', sy), ('Z', sz)]:
            full = np.array([[1]], dtype=complex)
            for j in range(N):
                full = np.kron(full, op if j == i else eye)
            paulis[(i, label)] = full
    P = np.zeros((N, len(t_grid)))
    for k, t in enumerate(t_grid):
        ct = c0 * np.exp(evals * t)
        rho_t = (R @ ct).reshape(d, d)
        for i in range(N):
            sx_i = np.real(np.trace(paulis[(i, 'X')] @ rho_t))
            sy_i = np.real(np.trace(paulis[(i, 'Y')] @ rho_t))
            sz_i = np.real(np.trace(paulis[(i, 'Z')] @ rho_t))
            P[i, k] = 0.5 * (1.0 + sx_i**2 + sy_i**2 + sz_i**2)
    return P


def fit_alpha(t_grid, P_A_site, P_B_site, bounds=(0.1, 10.0)):
    """Fit α s.t. P_A(α·t) ≈ P_B(t). Returns (α, rmse, on_boundary)."""
    interp = interp1d(t_grid, P_A_site, bounds_error=False,
                      fill_value=(float(P_A_site[0]), float(P_A_site[-1])),
                      kind='cubic')

    def mse(alpha):
        d = interp(alpha * t_grid) - P_B_site
        return float(np.mean(d * d))

    res = minimize_scalar(mse, bounds=bounds, method='bounded', options={'xatol': 1e-6})
    alpha = float(res.x)
    rmse = float(np.sqrt(res.fun))
    on_boundary = (abs(alpha - bounds[0]) < 1e-3 or abs(alpha - bounds[1]) < 1e-3)
    return alpha, rmse, on_boundary


def build_L(N, bonds, terms, gamma_l, gamma_t1_l, J=1.0):
    """Build full Lindbladian for given bilinear terms + Z-dephasing + T1."""
    terms_j = [(t[0], t[1], J) for t in terms]
    H = fw._build_bilinear(N, bonds, terms_j)
    if all(g == 0.0 for g in gamma_t1_l):
        return fw.lindbladian_z_dephasing(H, gamma_l)
    return fw.lindbladian_z_plus_t1(H, gamma_l, gamma_t1_l)


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
        ('truly XX+YY', [('X', 'X'), ('Y', 'Y')]),
        ('XY+YX',       [('X', 'Y'), ('Y', 'X')]),
        ('IY+YI',       [('I', 'Y'), ('Y', 'I')]),
        ('YZ+ZY',       [('Y', 'Z'), ('Z', 'Y')]),
        ('XZ+ZX',       [('X', 'Z'), ('Z', 'X')]),
        ('XZ+XZ',       [('X', 'Z'), ('X', 'Z')]),
    ]

    # Universal Carrier constants (chain-level invariants)
    Q_chain = J / GAMMA_DEPH
    t_peak = 1.0 / (4.0 * GAMMA_DEPH)
    absorption_quantum = 2.0 * GAMMA_DEPH

    print(f"PTF + Universal Carrier tour")
    print(f"  N={N}, |+−+⟩, γ_deph={GAMMA_DEPH}, γ_T1={GAMMA_T1}, J={J}, topology=chain")
    print()
    print(f"Universal Carrier (chain-level invariants, depend on γ₀ only):")
    print(f"  Q = J/γ₀          = {Q_chain:>8.4f}")
    print(f"  t_peak = 1/(4γ₀)  = {t_peak:>8.4f}")
    print(f"  Abs quantum = 2γ₀ = {absorption_quantum:>8.4f}")
    print()

    # Baseline L_A = truly XX+YY case
    gamma_l = [GAMMA_DEPH] * N
    gamma_t1_l = [GAMMA_T1] * N
    L_A = build_L(N, bonds, cases[0][1], gamma_l, gamma_t1_l, J)

    t_max = 8.0
    n_t = 200
    t_grid = np.linspace(0, t_max, n_t)
    P_A = per_site_purity_trajectory(L_A, rho_0, N, t_grid)

    # PTF readings per case (case-H as L_B, truly as L_A)
    print(f"PTF readings (L_A = truly XX+YY baseline; α_i = per-site time rescaling):")
    print()
    cols = [
        ('case',          14),
        ('α[0]',           9),
        ('α[1]',           9),
        ('α[2]',           9),
        ('Σ ln α',         10),
        ('max RMSE',       10),
        ('all in (0.1,10)?', 18),
    ]
    print(''.join(f"{name:<{w}s}" for name, w in cols))
    print('-' * sum(w for _, w in cols))

    for label, terms in cases:
        L_B = build_L(N, bonds, terms, gamma_l, gamma_t1_l, J)
        P_B = per_site_purity_trajectory(L_B, rho_0, N, t_grid)
        alphas = np.zeros(N)
        rmses = np.zeros(N)
        on_b = np.zeros(N, dtype=bool)
        for i in range(N):
            a, r, b = fit_alpha(t_grid, P_A[i], P_B[i])
            alphas[i] = a
            rmses[i] = r
            on_b[i] = b
        sigma_log = float(np.sum(np.log(np.clip(alphas, 1e-30, None))))

        row = [
            f"{label:<14s}",
            f"{alphas[0]:>8.4f} ",
            f"{alphas[1]:>8.4f} ",
            f"{alphas[2]:>8.4f} ",
            f"{sigma_log:>+9.4f} ",
            f"{rmses.max():>9.4f} ",
            f"{'yes' if not any(on_b) else 'NO':<18s}",
        ]
        print(''.join(row))

    print()
    print("Reading guide:")
    print("  α[i]            = time-rescaling factor s.t. P_A(α[i]·t) ≈ P_B(i, t)")
    print("                    α=1 → site evolves like baseline; α<1 → slower; α>1 → faster.")
    print("  Σ ln α          = PTF closure law (Tier 2). PTF predicts ≈ 0 in perturbative window.")
    print("  max RMSE        = worst per-site fit residual (mismatch with L_B not actually a")
    print("                    rescaling of L_A → fit is forced; large RMSE means PTF model")
    print("                    breaks for this case).")
    print()
    print("Observation: PTF was designed for SINGLE-BOND defects (small δJ on XY chain).")
    print("Forcing it onto wildly-different H (XY+YX, YZ+ZY, etc.) tests how far the")
    print("perspectival-time-field reading extends. Large RMSE = the case is outside the")
    print("perturbative regime PTF was tuned for; the α values are then 'best L_A-rescaling")
    print("approximation' rather than a structurally honest factorisation.")


if __name__ == "__main__":
    main()
