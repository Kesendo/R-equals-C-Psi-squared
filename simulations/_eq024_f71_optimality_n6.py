#!/usr/bin/env python3
"""EQ-024 N-scaling of F71-optimality: does the block-decomposition mechanism
hold at N=6?

At N=5: 4 bonds, mirror permutation R̄ gives 2+2 block decomposition (J_0+J_3,
J_1+J_2 symmetric; J_0−J_3, J_1−J_2 antisymmetric). F71-symmetric receivers
showed peaked SVD spectrum (sv₁/sv₄ = 4.06, vs 2.50 for F71-breaking).

At N=6: 5 bonds, R̄ gives 3+2 decomposition (J_0+J_4, J_1+J_3, J_2 symmetric;
J_0−J_4, J_1−J_3 antisymmetric). The symmetric block has one more dimension —
J_2 is self-mirror because there are an odd number of sites in each half.

Hypothesis: same mechanism holds at N=6, F71-sym SVDs systematically more
peaked than F71-breaking. Possibly even more pronounced because the
symmetric block has 3 dim vs antisym 2 dim (asymmetric block sizes).

Cost note: at N=6, L is 4096×4096. Direct expm(L·t) costs ~15 s per call;
54 calls per Jacobian (9 J-vectors × 6 times) → ~13 min per sample. We use
spectral evolution instead: one eig(L) per J-vector (~30 s), then exp(λ·t)
is cheap. Per sample: ~6 min.

Setup: N=6 Heisenberg, γ₀=0.05, J=1.0 uniform, t-grid [0.5, 1.0, 1.5, 2.0,
2.5, 3.0], dJ=1e-4, spread=0.02, σ_noise=0.01, P_total=5·spread².
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass

SCRIPT_DIR = Path(__file__).parent

# Paulis
I2 = np.eye(2, dtype=complex)
sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)

# Setup at N=6
N = 6
d = 2 ** N        # 64
d2 = d * d        # 4096
N_BONDS = N - 1   # 5
GAMMA_0 = 0.05
J_REF = 1.0
T_POINTS = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0]
DJ = 1e-4
SPREAD = 0.02
SIGMA_NOISE = 0.01
P_TOTAL = N_BONDS * SPREAD ** 2


def site_op(op, k):
    ops = [I2] * N
    ops[k] = op
    r = ops[0]
    for o in ops[1:]:
        r = np.kron(r, o)
    return r


# Heisenberg per-bond operators (XX+YY+ZZ)
print("Building bond operators...", flush=True)
_BOND_BLOCKS = []
for b in range(N_BONDS):
    block = (site_op(sx, b) @ site_op(sx, b + 1)
             + site_op(sy, b) @ site_op(sy, b + 1)
             + site_op(sz, b) @ site_op(sz, b + 1))
    _BOND_BLOCKS.append(block)


def build_H(J_vec):
    H = np.zeros((d, d), dtype=complex)
    for b in range(N_BONDS):
        H += J_vec[b] * _BOND_BLOCKS[b]
    return H


# Pre-build the dissipator part (J-independent)
print("Building dissipator...", flush=True)
_DISSIPATOR = np.zeros((d2, d2), dtype=complex)
_Id_d = np.eye(d, dtype=complex)
for k in range(N):
    Zk = site_op(sz, k)
    _DISSIPATOR += GAMMA_0 * (np.kron(Zk, Zk.conj()) - np.eye(d2))
print("Setup done.", flush=True)


def build_L(H):
    L = -1j * (np.kron(H, _Id_d) - np.kron(_Id_d, H.T))
    L += _DISSIPATOR
    return L


def propagate_eig(L, rho0_vec, t_list):
    """Return list of ρ(t) for each t in t_list, via spectral decomposition."""
    evals, R = np.linalg.eig(L)
    Rinv = np.linalg.solve(R, np.eye(d2))
    coeffs = Rinv @ rho0_vec
    rhos = []
    for t in t_list:
        v = R @ (np.exp(evals * t) * coeffs)
        rho = v.reshape(d, d)
        rho = (rho + rho.conj().T) / 2
        rhos.append(rho)
    return rhos


_LETTERS = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'


def ptrace_keep(rho, keep):
    keep = list(keep)
    trace_out = [q for q in range(N) if q not in keep]
    dims = [2] * N
    reshaped = rho.reshape(dims + dims)
    current_n = N
    for q in sorted(trace_out, reverse=True):
        reshaped = np.trace(reshaped, axis1=q, axis2=q + current_n)
        current_n -= 1
    d_k = 2 ** len(keep)
    return reshaped.reshape((d_k, d_k))


def purity(rho):
    return float(np.trace(rho @ rho).real)


def psi_norm(rho):
    d_r = rho.shape[0]
    l1 = float(np.sum(np.abs(rho)) - np.sum(np.abs(np.diag(rho))))
    return l1 / (d_r - 1) if d_r > 1 else 0.0


def cpsi(rho):
    return purity(rho) * psi_norm(rho)


def extract_features(rho):
    f = []
    for q in range(N):
        f.append(purity(ptrace_keep(rho, [q])))
    for i in range(N):
        for j in range(i + 1, N):
            f.append(cpsi(ptrace_keep(rho, [i, j])))
    for q in range(N):
        rho_q = ptrace_keep(rho, [q])
        f.append(float(np.real(np.trace(sz @ rho_q))))
    for q in range(N):
        rho_q = ptrace_keep(rho, [q])
        f.append(float(np.real(np.trace(sx @ rho_q))))
    return np.array(f)


def extract_timeseries(rho0, J_vec):
    H = build_H(J_vec)
    L = build_L(H)
    rho0_vec = rho0.flatten()
    rhos = propagate_eig(L, rho0_vec, T_POINTS)
    feats = []
    for rho in rhos:
        feats.extend(extract_features(rho).tolist())
    return np.array(feats)


def compute_jacobian_J(rho0):
    J_ref = [J_REF] * N_BONDS
    y_ref = extract_timeseries(rho0, J_ref)
    A = np.zeros((len(y_ref), N_BONDS))
    for b in range(N_BONDS):
        J_plus = list(J_ref); J_plus[b] += DJ
        y_plus = extract_timeseries(rho0, J_plus)
        J_minus = list(J_ref); J_minus[b] -= DJ
        y_minus = extract_timeseries(rho0, J_minus)
        A[:, b] = (y_plus - y_minus) / (2.0 * DJ)
    return A


def waterfilling_capacity(sv):
    n = len(sv)
    if np.max(sv) < 1e-8:
        return 0.0
    gains = sv ** 2 / SIGMA_NOISE ** 2
    idx = np.argsort(gains)[::-1]
    gains_sorted = gains[idx]
    inv_gains = 1.0 / gains_sorted
    k_used = 1
    powers = None
    for k in range(n, 0, -1):
        mu = (P_TOTAL + np.sum(inv_gains[:k])) / k
        powers_try = mu - inv_gains[:k]
        if np.all(powers_try >= 0):
            k_used = k
            powers = powers_try
            break
    C = 0.0
    for i in range(k_used):
        snr_i = powers[i] * gains_sorted[i]
        C += 0.5 * np.log2(1 + snr_i)
    return float(C)


def bloch_state(theta, phi):
    return np.array([np.cos(theta / 2),
                     np.exp(1j * phi) * np.sin(theta / 2)], dtype=complex)


def random_f71_sym_product(rng):
    """N=6 F71-symmetric product: ψ = |a⟩|b⟩|c⟩|c⟩|b⟩|a⟩, 3 Bloch pairs."""
    th_a = rng.uniform(0, np.pi); ph_a = rng.uniform(0, 2 * np.pi)
    th_b = rng.uniform(0, np.pi); ph_b = rng.uniform(0, 2 * np.pi)
    th_c = rng.uniform(0, np.pi); ph_c = rng.uniform(0, 2 * np.pi)
    a = bloch_state(th_a, ph_a)
    b = bloch_state(th_b, ph_b)
    c = bloch_state(th_c, ph_c)
    psi = a
    for v in [b, c, c, b, a]:
        psi = np.kron(psi, v)
    return psi


def random_product_asymmetric(rng):
    """N=6 product with 12 indep Bloch angles, no symmetry."""
    psi = np.array([1.0], dtype=complex)
    for _ in range(N):
        theta = rng.uniform(0, np.pi)
        phi = rng.uniform(0, 2 * np.pi)
        psi = np.kron(psi, bloch_state(theta, phi))
    return psi


def measure(psi):
    rho0 = np.outer(psi, psi.conj())
    rho0 = (rho0 + rho0.conj().T) / 2.0
    rho0 /= np.trace(rho0).real
    A = compute_jacobian_J(rho0)
    sv = np.linalg.svd(A, compute_uv=False)
    C = waterfilling_capacity(sv)
    return sv, C


def main():
    n_samples = int(sys.argv[1]) if len(sys.argv) > 1 else 15
    seed = int(sys.argv[2]) if len(sys.argv) > 2 else 100

    rng_sym = np.random.default_rng(seed)
    rng_brk = np.random.default_rng(seed + 1)

    print(f"EQ-024 F71-optimality at N={N}: SVD analysis")
    print(f"  {n_samples} samples per mode (F71-sym product, F71-breaking product)")
    print(f"  Bond-input dim = {N_BONDS}, predicted block decomp: 3+2 (sym/antisym)")
    print(f"  Liouvillian dim = {d2}, spectral propagation via eig(L)")
    print()

    sym_data = []
    brk_data = []
    t0 = time.time()

    for i in range(n_samples):
        psi_s = random_f71_sym_product(rng_sym)
        sv_s, C_s = measure(psi_s)
        sym_data.append({'sv': sv_s.tolist(), 'C': C_s})

        psi_b = random_product_asymmetric(rng_brk)
        sv_b, C_b = measure(psi_b)
        brk_data.append({'sv': sv_b.tolist(), 'C': C_b})

        elapsed = time.time() - t0
        rate = elapsed / (i + 1)
        eta = rate * (n_samples - i - 1)
        print(f"  [{i + 1}/{n_samples}] elapsed {elapsed:.0f}s, "
              f"sym sv₁={sv_s[0]:.2f}, brk sv₁={sv_b[0]:.2f}, "
              f"ETA {eta/60:.1f} min", flush=True)

    sym_sv = np.array([d['sv'] for d in sym_data])
    brk_sv = np.array([d['sv'] for d in brk_data])
    sym_C = np.array([d['C'] for d in sym_data])
    brk_C = np.array([d['C'] for d in brk_data])

    print()
    print("=== Capacity ===")
    print(f"  F71-sym:      C range [{sym_C.min():.3f}, {sym_C.max():.3f}], mean {sym_C.mean():.3f}")
    print(f"  F71-breaking: C range [{brk_C.min():.3f}, {brk_C.max():.3f}], mean {brk_C.mean():.3f}")

    sym_ratio = sym_sv[:, 0] / np.maximum(sym_sv[:, -1], 1e-12)
    brk_ratio = brk_sv[:, 0] / np.maximum(brk_sv[:, -1], 1e-12)

    print()
    print("=== SV spectrum (5 SVs at N=6) ===")
    print(f"  {'mode':>15s}  " + "  ".join(f"sv{k+1} mean" for k in range(N_BONDS)) + f"  sv₁/sv₅ mean")
    print(f"  {'F71-sym':>15s}  " + "  ".join(f"  {sym_sv[:, k].mean():>6.3f}" for k in range(N_BONDS))
          + f"     {sym_ratio.mean():>6.2f}")
    print(f"  {'F71-breaking':>15s}  " + "  ".join(f"  {brk_sv[:, k].mean():>6.3f}" for k in range(N_BONDS))
          + f"     {brk_ratio.mean():>6.2f}")

    print()
    print("=== Hypothesis test (peakedness sv₁/sv₅) ===")
    print(f"  F71-sym:      min {sym_ratio.min():.2f}, max {sym_ratio.max():.2f}, "
          f"mean {sym_ratio.mean():.2f}, median {np.median(sym_ratio):.2f}")
    print(f"  F71-breaking: min {brk_ratio.min():.2f}, max {brk_ratio.max():.2f}, "
          f"mean {brk_ratio.mean():.2f}, median {np.median(brk_ratio):.2f}")

    pairs = 0; sym_wins = 0
    for r_s in sym_ratio:
        for r_b in brk_ratio:
            pairs += 1
            if r_s > r_b:
                sym_wins += 1
    print(f"  P(F71-sym sv₁/sv₅ > F71-breaking sv₁/sv₅) ≈ {sym_wins / pairs:.3f}")

    sym_skew = sym_sv[:, 0] ** 2 / np.sum(sym_sv ** 2, axis=1)
    brk_skew = brk_sv[:, 0] ** 2 / np.sum(brk_sv ** 2, axis=1)
    print()
    print("=== Concentration sv₁²/Σsv² ===")
    print(f"  F71-sym:      min {sym_skew.min():.3f}, max {sym_skew.max():.3f}, "
          f"mean {sym_skew.mean():.3f}")
    print(f"  F71-breaking: min {brk_skew.min():.3f}, max {brk_skew.max():.3f}, "
          f"mean {brk_skew.mean():.3f}")

    pairs = 0; sym_wins = 0
    for s in sym_skew:
        for b in brk_skew:
            pairs += 1
            if s > b:
                sym_wins += 1
    print(f"  P(F71-sym concentration > F71-breaking) ≈ {sym_wins / pairs:.3f}")

    print()
    print(f"=== Comparison to N=5 baseline (from sub-q follow-up, seed 42) ===")
    print(f"  N=5: F71-sym sv₁/sv₄ mean 4.06, F71-breaking 2.50, P_sym>brk = 0.806")
    print(f"  N=6: F71-sym sv₁/sv₅ mean {sym_ratio.mean():.2f}, F71-breaking {brk_ratio.mean():.2f}, "
          f"P_sym>brk = {sym_wins / pairs:.3f} (concentration)")

    out = SCRIPT_DIR / "results" / "eq024_f71_optimality_n6.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, 'w', encoding='utf-8') as f:
        json.dump({
            'config': {'N': N, 'n_samples': n_samples, 'seed': seed,
                       'gamma_0': GAMMA_0, 'J_ref': J_REF, 'dJ': DJ,
                       't_points': T_POINTS, 'spread': SPREAD,
                       'sigma_noise': SIGMA_NOISE},
            'f71_sym': sym_data,
            'f71_breaking': brk_data,
        }, f, indent=1)
    print(f"\nSaved: {out}")
    print(f"Total runtime: {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
