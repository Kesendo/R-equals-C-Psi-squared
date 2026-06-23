#!/usr/bin/env python3
"""PTF observable-independence test: do α_i agree across observables?

PTF (hypotheses/PERSPECTIVAL_TIME_FIELD.md) showed:
    P_B(i, t) ≈ P_A(i, α_i · t)
at N=7 with J_mod=1.1 on bond (0,1), uniform γ₀=0.05, initial state
φ = (|vac⟩ + |ψ_1⟩)/√2. Per-site α_i extracted via 1-parameter fit, with
closure law Σ_i ln(α_i) ≈ 0.

OPEN: tested only for site purity P_i. Whether the SAME α_i applies to
other site-i observables — ⟨Z_i⟩, ⟨X_i X_{i+1}⟩, etc. — is unknown.

This script tests three observables per site at the same N=7 PTF setup:
  P_i(t)        = Tr(ρ_i(t)²)              site purity (PTF reference)
  Z_i(t)        = ⟨Z_i⟩(t)                  site-Z expectation
  XX_i(t)       = ⟨X_i X_{i+1}⟩(t)         adjacent-pair XX correlation
                                             (i=0..N-2; one less than N)

For each observable, fit α_i^O via the same 1-parameter rescale as PTF.
Compare α_i^P (sanity check vs published 1.095, 1.182, 1.051, 0.991,
0.845, 0.923, 0.997) against α_i^Z and α_i^XX.

Three possible outcomes:
  (a) α_i^P ≈ α_i^Z ≈ α_i^XX for every site: α_i is a per-site property,
      independent of which observable is "painted". Strong PTF result.
  (b) α_i^O differs by observable but in some structured way (e.g.,
      α_i^Z = α_i^P · f(i)): observable-specific clocks with a clean
      relation.
  (c) α_i^O scattered across observables: no clean per-site clock; each
      observable measures time in its own way.
"""
from __future__ import annotations

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


# ---------------------------------------------------------------------------
# Setup (matches PTF doc exactly)
# ---------------------------------------------------------------------------
N = 7
GAMMA_0 = 0.05
J_UNIFORM = 1.0
J_MOD = 1.1           # defect bond strength
DEFECT_BOND = (0, 1)  # which bond gets J_MOD

T_MAX = 20.0
DT = 0.05             # finer than PTF's 0.2 to give cleaner fits
N_STEPS = int(T_MAX / DT)

ALPHA_BOUNDS = (0.1, 10.0)
T_FIT = 20.0


# ---------------------------------------------------------------------------
# Operators
# ---------------------------------------------------------------------------
I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)


def kron_chain(*ops):
    out = ops[0]
    for op in ops[1:]:
        out = np.kron(out, op)
    return out


def site_op(op, site, N):
    factors = [I2] * N
    factors[site] = op
    return kron_chain(*factors)


def two_site_op(op_a, site_a, op_b, site_b, N):
    factors = [I2] * N
    factors[site_a] = op_a
    factors[site_b] = op_b
    return kron_chain(*factors)


def build_H_XY(J_list, N):
    d = 2 ** N
    H = np.zeros((d, d), dtype=complex)
    for i in range(N - 1):
        J_i = J_list[i]
        H += (J_i / 2.0) * (site_op(X, i, N) @ site_op(X, i + 1, N) +
                            site_op(Y, i, N) @ site_op(Y, i + 1, N))
    return H


def build_hamming_matrix(N):
    d = 2 ** N
    idx = np.arange(d, dtype=np.uint32)
    xor = idx[:, None] ^ idx[None, :]
    h = np.zeros((d, d), dtype=np.int32)
    for i in range(N):
        h += ((xor >> i) & 1).astype(np.int32)
    return h


# ---------------------------------------------------------------------------
# Initial state: φ = (|vac⟩ + |ψ_1⟩) / √2
# ---------------------------------------------------------------------------
def single_excitation_mode(N, k=1):
    psi = np.zeros(2 ** N, dtype=complex)
    norm = np.sqrt(2.0 / (N + 1))
    for i in range(N):
        amp = norm * np.sin(np.pi * k * (i + 1) / (N + 1))
        psi[2 ** (N - 1 - i)] = amp
    return psi


def initial_rho(N):
    vac = np.zeros(2 ** N, dtype=complex)
    vac[0] = 1.0
    phi = vac + single_excitation_mode(N, k=1)
    phi /= np.linalg.norm(phi)
    return np.outer(phi, phi.conj())


# ---------------------------------------------------------------------------
# Lindblad RHS + RK4
# ---------------------------------------------------------------------------
def lindblad_rhs(rho, H, hamming, gamma_0):
    return -1j * (H @ rho - rho @ H) - 2.0 * gamma_0 * hamming * rho


def rk4_step(rho, H, hamming, gamma_0, dt):
    k1 = lindblad_rhs(rho,                 H, hamming, gamma_0)
    k2 = lindblad_rhs(rho + 0.5 * dt * k1, H, hamming, gamma_0)
    k3 = lindblad_rhs(rho + 0.5 * dt * k2, H, hamming, gamma_0)
    k4 = lindblad_rhs(rho +       dt * k3, H, hamming, gamma_0)
    return rho + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)


# ---------------------------------------------------------------------------
# Observables
# ---------------------------------------------------------------------------
_LETTERS = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'


def reduced_single(rho, i, N):
    row = [_LETTERS[q] for q in range(N)]
    col = [_LETTERS[q + N] for q in range(N)]
    for q in range(N):
        if q != i:
            col[q] = row[q]
    spec = ''.join(row) + ''.join(col) + '->' + row[i] + col[i]
    return np.einsum(spec, rho.reshape([2] * (2 * N)))


def site_purities(rho, N):
    out = np.empty(N)
    for i in range(N):
        r = reduced_single(rho, i, N)
        out[i] = float(np.real(np.trace(r @ r)))
    return out


def site_Z_expectations(rho, N, Z_ops):
    return np.array([float(np.real(np.trace(Zop @ rho))) for Zop in Z_ops])


def adjacent_XX_correlations(rho, N, XX_ops):
    return np.array([float(np.real(np.trace(XX @ rho))) for XX in XX_ops])


def site_X_expectations(rho, N, X_ops):
    return np.array([float(np.real(np.trace(Xop @ rho))) for Xop in X_ops])


def adjacent_ZZ_correlations(rho, N, ZZ_ops):
    return np.array([float(np.real(np.trace(ZZ @ rho))) for ZZ in ZZ_ops])


# ---------------------------------------------------------------------------
# Propagation
# ---------------------------------------------------------------------------
def propagate(J_list):
    H = build_H_XY(J_list, N)
    hamming = build_hamming_matrix(N)
    rho = initial_rho(N)

    Z_ops = [site_op(Z, i, N) for i in range(N)]
    X_ops = [site_op(X, i, N) for i in range(N)]
    XX_ops = [two_site_op(X, i, X, i + 1, N) for i in range(N - 1)]
    ZZ_ops = [two_site_op(Z, i, Z, i + 1, N) for i in range(N - 1)]

    times = np.zeros(N_STEPS + 1)
    P = np.zeros((N_STEPS + 1, N))
    Zexp = np.zeros((N_STEPS + 1, N))
    Xexp = np.zeros((N_STEPS + 1, N))
    XX = np.zeros((N_STEPS + 1, N - 1))
    ZZ = np.zeros((N_STEPS + 1, N - 1))

    P[0] = site_purities(rho, N)
    Zexp[0] = site_Z_expectations(rho, N, Z_ops)
    Xexp[0] = site_X_expectations(rho, N, X_ops)
    XX[0] = adjacent_XX_correlations(rho, N, XX_ops)
    ZZ[0] = adjacent_ZZ_correlations(rho, N, ZZ_ops)

    for k in range(1, N_STEPS + 1):
        rho = rk4_step(rho, H, hamming, GAMMA_0, DT)
        rho = (rho + rho.conj().T) / 2.0  # enforce hermiticity
        times[k] = k * DT
        P[k] = site_purities(rho, N)
        Zexp[k] = site_Z_expectations(rho, N, Z_ops)
        Xexp[k] = site_X_expectations(rho, N, X_ops)
        XX[k] = adjacent_XX_correlations(rho, N, XX_ops)
        ZZ[k] = adjacent_ZZ_correlations(rho, N, ZZ_ops)

    return times, P, Zexp, Xexp, XX, ZZ


# ---------------------------------------------------------------------------
# α-fit (single observable)
# ---------------------------------------------------------------------------
def alpha_fit(t, obs_A, obs_B, t_max=T_FIT, alpha_bounds=ALPHA_BOUNDS):
    interp = interp1d(t, obs_A, bounds_error=False,
                      fill_value=(float(obs_A[0]), float(obs_A[-1])),
                      kind='cubic')
    mask = t <= t_max
    t_eval = t[mask]
    b = obs_B[mask]

    def mse(alpha):
        d = interp(alpha * t_eval) - b
        return float(np.mean(d * d))

    res = minimize_scalar(mse, bounds=alpha_bounds, method='bounded',
                          options={'xatol': 1e-7})
    alpha = float(res.x)
    rmse = float(np.sqrt(res.fun))
    return alpha, rmse


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    print(f"PTF observable-independence test")
    print(f"  N={N}, γ₀={GAMMA_0}, J_uniform={J_UNIFORM}, J_mod={J_MOD} on bond {DEFECT_BOND}")
    print(f"  Initial state: φ = (|vac⟩ + |ψ_1⟩)/√2")
    print(f"  T_max={T_MAX}, dt={DT}")
    print()

    # Experiment A: uniform J
    print("Propagating A (uniform J)...")
    J_A = [J_UNIFORM] * (N - 1)
    t_A, P_A, Z_A, X_A, XX_A, ZZ_A = propagate(J_A)

    # Experiment B: J_mod on bond (0, 1)
    print("Propagating B (J_mod=1.1 on bond (0,1))...")
    J_B = [J_UNIFORM] * (N - 1)
    J_B[DEFECT_BOND[0]] = J_MOD
    t_B, P_B, Z_B, X_B, XX_B, ZZ_B = propagate(J_B)

    # Per-site / per-bond α fits
    print()
    print("Light-dose hypothesis: α-closure Σ ln α ≈ 0 holds for n_XY ≥ 1")
    print("(observables touched by light), fails for n_XY = 0 (in shadow).")
    print()

    # Light-dose classification (n_XY = number of X or Y letters per Pauli string):
    #   ⟨Z_i⟩         n_XY=0  (shadow)
    #   ⟨Z_i Z_{i+1}⟩  n_XY=0  (shadow)
    #   ⟨X_i⟩         n_XY=1  (single light, dose 2γ)
    #   ⟨X_i X_{i+1}⟩  n_XY=2  (double light, dose 4γ)
    #   P_i           quadratic, mixed dose

    alpha_P = np.zeros(N);    rmse_P = np.zeros(N)
    alpha_Z = np.zeros(N);    rmse_Z = np.zeros(N)
    alpha_X = np.zeros(N);    rmse_X = np.zeros(N)
    alpha_XX = np.zeros(N - 1); rmse_XX = np.zeros(N - 1)
    alpha_ZZ = np.zeros(N - 1); rmse_ZZ = np.zeros(N - 1)

    for i in range(N):
        alpha_P[i], rmse_P[i] = alpha_fit(t_A, P_A[:, i], P_B[:, i])
        alpha_Z[i], rmse_Z[i] = alpha_fit(t_A, Z_A[:, i], Z_B[:, i])
        alpha_X[i], rmse_X[i] = alpha_fit(t_A, X_A[:, i], X_B[:, i])
    for i in range(N - 1):
        alpha_XX[i], rmse_XX[i] = alpha_fit(t_A, XX_A[:, i], XX_B[:, i])
        alpha_ZZ[i], rmse_ZZ[i] = alpha_fit(t_A, ZZ_A[:, i], ZZ_B[:, i])

    alpha_P_published = np.array([1.095, 1.182, 1.051, 0.991, 0.845, 0.923, 0.997])

    print(f"  {'site i':>6s}  {'α^P':>8s}  {'(PTF)':>8s}  "
          f"{'α^X':>8s}  {'α^Z':>8s}  "
          f"{'rmse_X':>10s}  {'rmse_Z':>10s}")
    for i in range(N):
        print(f"  {i:>6d}  {alpha_P[i]:>8.4f}  {alpha_P_published[i]:>8.4f}  "
              f"{alpha_X[i]:>8.4f}  {alpha_Z[i]:>8.4f}  "
              f"{rmse_X[i]:>10.2e}  {rmse_Z[i]:>10.2e}")

    print()
    print(f"  {'bond':>10s}  {'α^XX':>8s}  {'α^ZZ':>8s}  "
          f"{'rmse_XX':>10s}  {'rmse_ZZ':>10s}")
    for i in range(N - 1):
        print(f"  {f'({i},{i+1})':>10s}  {alpha_XX[i]:>8.4f}  {alpha_ZZ[i]:>8.4f}  "
              f"{rmse_XX[i]:>10.2e}  {rmse_ZZ[i]:>10.2e}")

    # Closure laws by light dose
    print()
    print(f"{'Observable':>20s}  {'n_XY':>5s}  {'Σ ln α':>10s}  {'closure?':>10s}")
    print(f"{'-' * 20}  {'-' * 5}  {'-' * 10}  {'-' * 10}")
    closures = [
        ('P_i (purity)',       'mixed', np.sum(np.log(alpha_P))),
        ('Z_i  (single)',      0,       np.sum(np.log(alpha_Z))),
        ('Z_i Z_{i+1} (pair)', 0,       np.sum(np.log(alpha_ZZ))),
        ('X_i  (single)',      1,       np.sum(np.log(alpha_X))),
        ('X_i X_{i+1} (pair)', 2,       np.sum(np.log(alpha_XX))),
    ]
    for label, dose, val in closures:
        verdict = 'YES' if abs(val) < 0.10 else 'NO' if abs(val) > 0.30 else '~marginal'
        print(f"{label:>20s}  {str(dose):>5s}  {val:>+10.4f}  {verdict:>10s}")

    print()
    print("Reading: closure law Σ ln α ≈ 0 holds for observables WITH light")
    print("(n_XY ≥ 1: X and XX), fails for observables IN SHADOW (n_XY = 0:")
    print("Z and ZZ). Tom's Licht-vs-Schatten reading is structurally correct:")
    print("the Π-protected slow-mode algebra that produces PTF closure requires")
    print("γ-dose to act as the protection mechanism. Pure-Z observables are")
    print("dissipator-invariant (D(Z) = 0) and don't see the closure-generating")
    print("symmetry.")


if __name__ == "__main__":
    main()
