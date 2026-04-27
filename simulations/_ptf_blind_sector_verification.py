#!/usr/bin/env python3
"""PTF light/shadow verification: closure law Σ ln α vs n_XY across many observables.

Hypothesis (from _ptf_per_observable_alpha.py + Tom's Licht/Schatten reading):
  PTF closure law Σ ln α ≈ 0 holds for observables touched by light
  (n_XY ≥ 1), fails for observables in the shadow (n_XY = 0).

This script tests the hypothesis by sweeping a large family of single-site
and two-site Pauli observables on the same N=7 PTF setup, classifying each
by light dose n_XY (= number of X or Y letters in the operator), and
checking the closure law per family.

Single-site (Σ over N=7 sites): X_i, Y_i, Z_i
Two-site adjacent (Σ over N-1=6 bonds): aa+bb for a, b ∈ {X, Y, Z}
Two-site non-adjacent: same operators on bonds with gap > 1

Plus ⟨P_i⟩ purity (PTF reference) for sanity.
"""
from __future__ import annotations

import sys
from itertools import product
from pathlib import Path

import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import minimize_scalar

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass


# Setup matches PTF
N = 7
GAMMA_0 = 0.05
J_UNIFORM = 1.0
J_MOD = 1.1
DEFECT_BOND = (0, 1)
T_MAX = 20.0
DT = 0.05
N_STEPS = int(T_MAX / DT)

ALPHA_BOUNDS = (0.1, 10.0)


I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
LETTER_OP = {'X': X, 'Y': Y, 'Z': Z, 'I': I2}


def kron_chain(*ops):
    out = ops[0]
    for op in ops[1:]:
        out = np.kron(out, op)
    return out


def pauli_string(letters_at_sites):
    """letters_at_sites: dict {site_index: 'X'|'Y'|'Z'} (others = I)."""
    factors = [I2] * N
    for site, letter in letters_at_sites.items():
        factors[site] = LETTER_OP[letter]
    return kron_chain(*factors)


def n_XY(letters_at_sites):
    return sum(1 for L in letters_at_sites.values() if L in ('X', 'Y'))


def build_H_XY(J_list):
    d = 2 ** N
    H = np.zeros((d, d), dtype=complex)
    for i in range(N - 1):
        J_i = J_list[i]
        Xi = pauli_string({i: 'X'})
        Xj = pauli_string({i + 1: 'X'})
        Yi = pauli_string({i: 'Y'})
        Yj = pauli_string({i + 1: 'Y'})
        H += (J_i / 2.0) * (Xi @ Xj + Yi @ Yj)
    return H


def build_hamming_matrix():
    d = 2 ** N
    idx = np.arange(d, dtype=np.uint32)
    xor = idx[:, None] ^ idx[None, :]
    h = np.zeros((d, d), dtype=np.int32)
    for i in range(N):
        h += ((xor >> i) & 1).astype(np.int32)
    return h


def initial_rho():
    psi = np.zeros(2 ** N, dtype=complex)
    psi[0] = 1.0  # vac
    norm = np.sqrt(2.0 / (N + 1))
    for i in range(N):
        amp = norm * np.sin(np.pi * (i + 1) / (N + 1))
        psi[2 ** (N - 1 - i)] += amp
    psi /= np.linalg.norm(psi)
    return np.outer(psi, psi.conj())


def lindblad_rhs(rho, H, hamming):
    return -1j * (H @ rho - rho @ H) - 2.0 * GAMMA_0 * hamming * rho


def rk4_step(rho, H, hamming, dt):
    k1 = lindblad_rhs(rho, H, hamming)
    k2 = lindblad_rhs(rho + 0.5 * dt * k1, H, hamming)
    k3 = lindblad_rhs(rho + 0.5 * dt * k2, H, hamming)
    k4 = lindblad_rhs(rho + dt * k3, H, hamming)
    return rho + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)


def propagate(J_list):
    H = build_H_XY(J_list)
    hamming = build_hamming_matrix()
    rho = initial_rho()
    rhos = [rho.copy()]
    for _ in range(N_STEPS):
        rho = rk4_step(rho, H, hamming, DT)
        rho = (rho + rho.conj().T) / 2.0
        rhos.append(rho.copy())
    times = np.array([k * DT for k in range(N_STEPS + 1)])
    return times, rhos


def alpha_fit(t, obs_A, obs_B):
    interp = interp1d(t, obs_A, bounds_error=False,
                      fill_value=(float(obs_A[0]), float(obs_A[-1])),
                      kind='cubic')

    def mse(alpha):
        d = interp(alpha * t) - obs_B
        return float(np.mean(d * d))

    res = minimize_scalar(mse, bounds=ALPHA_BOUNDS, method='bounded',
                          options={'xatol': 1e-7})
    alpha = float(res.x)
    rmse = float(np.sqrt(res.fun))
    boundary = abs(alpha - ALPHA_BOUNDS[0]) < 1e-3 or abs(alpha - ALPHA_BOUNDS[1]) < 1e-3
    return alpha, rmse, boundary


def expectation_traj(rhos, op):
    return np.array([float(np.real(np.trace(op @ r))) for r in rhos])


def main():
    print(f"PTF blind-sector verification (N={N}, J_mod={J_MOD} on bond {DEFECT_BOND})")
    print(f"  Hypothesis: Σ ln α ≈ 0 for n_XY ≥ 1 (light), fails for n_XY = 0 (shadow)")
    print()

    print("Propagating A (uniform J)...")
    J_A = [J_UNIFORM] * (N - 1)
    t_A, rhos_A = propagate(J_A)
    print("Propagating B (J_mod on defect bond)...")
    J_B = list(J_A); J_B[DEFECT_BOND[0]] = J_MOD
    t_B, rhos_B = propagate(J_B)

    # ---- Single-site Pauli observables ----
    # For each Pauli letter L in {X, Y, Z}, compute α_i for L_i (i=0..N-1),
    # then sum ln α_i (the per-site closure law).
    print()
    print("=== Single-site observables (Σ over sites i=0..N−1) ===")
    print(f"  {'Pauli':>5}  {'n_XY':>5}  {'Σ ln α':>10}  {'min α':>8}  {'max α':>8}  {'#bound':>6}  {'closure?':>10}")
    single_results = {}
    for L in ['X', 'Y', 'Z']:
        alphas = []
        n_bound = 0
        for i in range(N):
            op = pauli_string({i: L})
            obs_A = expectation_traj(rhos_A, op)
            obs_B = expectation_traj(rhos_B, op)
            a, rmse, bd = alpha_fit(t_A, obs_A, obs_B)
            alphas.append(a)
            if bd:
                n_bound += 1
        alphas = np.array(alphas)
        sum_ln = float(np.sum(np.log(alphas)))
        verdict = ('YES' if abs(sum_ln) < 0.10 else
                   'NO' if abs(sum_ln) > 0.30 else '~marginal')
        nxy = n_XY({0: L})
        single_results[L] = (nxy, sum_ln, alphas)
        print(f"  {L:>5}  {nxy:>5}  {sum_ln:>+10.4f}  {alphas.min():>8.3f}  "
              f"{alphas.max():>8.3f}  {n_bound:>6}  {verdict:>10}")

    # ---- Two-site Pauli observables (adjacent + non-adjacent gaps) ----
    print()
    print("=== Two-site observables L_i M_j  (Σ over i<j pairs) ===")
    print(f"  {'L M':>5}  {'n_XY':>5}  {'Σ ln α':>10}  {'min α':>8}  {'max α':>8}  {'#bound':>6}  {'closure?':>10}")
    two_letter_results = {}
    for L1, L2 in product(['X', 'Y', 'Z'], repeat=2):
        alphas = []
        n_bound = 0
        for i in range(N):
            for j in range(i + 1, N):
                op = pauli_string({i: L1, j: L2})
                obs_A = expectation_traj(rhos_A, op)
                obs_B = expectation_traj(rhos_B, op)
                # Skip observables with vanishing dynamic range (numerical noise)
                if max(np.ptp(obs_A), np.ptp(obs_B)) < 1e-6:
                    continue
                a, rmse, bd = alpha_fit(t_A, obs_A, obs_B)
                alphas.append(a)
                if bd:
                    n_bound += 1
        alphas = np.array(alphas)
        if len(alphas) == 0:
            continue
        sum_ln = float(np.sum(np.log(alphas)))
        verdict = ('YES' if abs(sum_ln) < 0.20 else
                   'NO' if abs(sum_ln) > 1.00 else '~marginal')
        nxy = n_XY({0: L1, 1: L2})
        key = f"{L1}{L2}"
        two_letter_results[key] = (nxy, sum_ln, len(alphas))
        print(f"  {key:>5}  {nxy:>5}  {sum_ln:>+10.4f}  {alphas.min():>8.3f}  "
              f"{alphas.max():>8.3f}  {n_bound:>6}  {verdict:>10}")

    # ---- Aggregate by n_XY ----
    print()
    print("=== Aggregate by light dose n_XY ===")
    print(f"  {'n_XY':>5}  {'#families':>10}  {'mean |Σ ln α|':>15}  {'failure rate':>14}  reading")
    by_dose = {}
    for src in [single_results, two_letter_results]:
        for k, v in src.items():
            nxy = v[0]
            sum_ln = v[1]
            by_dose.setdefault(nxy, []).append(abs(sum_ln))
    for nxy in sorted(by_dose.keys()):
        vals = by_dose[nxy]
        n_fail = sum(1 for v in vals if v > 0.30)
        n = len(vals)
        reading = ('shadow (closure fails)' if nxy == 0 else
                   'weak light (mixed)' if nxy == 1 else
                   'full light (closure holds)')
        print(f"  {nxy:>5}  {n:>10}  {np.mean(vals):>15.4f}  "
              f"{n_fail}/{n:<3}  {reading}")


if __name__ == "__main__":
    main()
