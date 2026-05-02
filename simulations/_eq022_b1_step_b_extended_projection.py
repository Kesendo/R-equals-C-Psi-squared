#!/usr/bin/env python3
"""_eq022_b1_step_b_extended_projection.py — EQ-022 (b1) Step (b).

Extended channel projection: 2c-dim effective model.

Motivation (from step a)
------------------------
The Dicke probe |S_n⟩⟨S_{n+1}|/2 lives entirely in the c-dim channel-
uniform subspace span(|c_1⟩, ..., |c_c⟩) where |c_k⟩ is the orthonormal
uniform projector on HD=2k-1 states. M_H_eff in this basis is diagonal,
so the c×c effective L lacks Q-dependent rate mixing. But probe DYNAMICS
under L is Q-dependent because L_H|c_k⟩ has components OUTSIDE the
channel-uniform subspace.

Strategy
--------
For each k, define |c_k^⊥⟩ as the unit vector along the orthogonal
complement of L_H|c_k⟩ relative to span(|c_1⟩, ..., |c_c⟩).
Together {|c_k⟩, |c_k^⊥⟩} span a 2c-dim subspace. Build L_eff in this
basis. If L_eff is closed (L_H · vector ∈ subspace, no further leakage),
the Dicke probe dynamics is exactly captured by the 2c×2c model and
Q_peak follows from its eigenstructure.

If 2c is not closed, iterate: compute next-generation orthogonals
|c_k^{⊥⊥}⟩, etc., until closure.

Verification path
-----------------
- Build at multiple (N, n) with same c. Test L_eff N-invariance.
- Compute Q_peak from the 2c×2c (or larger) effective L.
- Compare to F86 values: 1.6 (c=3), 1.8 (c=4, c=5).
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

sys.path.insert(0, str(Path(__file__).resolve().parent))
import framework as fw  # noqa: E402


def _block_L_split_total(N, n, gamma_0):
    """Aggregate the per-bond M_H from `block_L_split_xy` to (D, M_H_total)."""
    D, M_H_per_bond, _, _ = fw.block_L_split_xy(N, n, gamma_0)
    return D, sum(M_H_per_bond)


def gram_schmidt_extend(P_existing, candidates, atol=1e-12):
    """Append columns from `candidates` to `P_existing`, orthonormalising via
    modified Gram-Schmidt. Returns the augmented basis matrix.

    P_existing: (M, k) matrix with orthonormal columns.
    candidates: (M, m) matrix of vectors to orthogonalise against existing
                AND against each other.

    Skips candidates that are linearly dependent on the current basis (norm
    after projection-out < atol).
    """
    Q = P_existing.copy()
    for j in range(candidates.shape[1]):
        v = candidates[:, j].copy()
        # Project out current basis
        for i in range(Q.shape[1]):
            v = v - (Q[:, i].conj() @ v) * Q[:, i]
        nv = np.linalg.norm(v)
        if nv > atol:
            v /= nv
            Q = np.column_stack([Q, v])
    return Q


def build_extended_subspace(N, n, gamma_0, max_generations=5, atol=1e-12):
    """Iteratively build the smallest L_H-invariant subspace containing the
    channel-uniform projectors.

    Generation 0: |c_1⟩, ..., |c_c⟩ (channel uniforms).
    Generation g+1: append L_H · (gen g vectors), orthogonalised.

    Stops when no new linearly-independent vectors emerge (subspace closed
    under L_H).
    """
    D, M_H = _block_L_split_total(N, n, gamma_0)
    P_chan, HDs = fw.hd_channel_basis(N, n)
    c = P_chan.shape[1]

    P = P_chan.copy()
    gen_sizes = [c]
    for g in range(max_generations):
        last_block = P[:, sum(gen_sizes[:-1]):]  # most recent generation
        # Apply M_H (J-coefficient) to last generation
        cands = M_H @ last_block
        P_new = gram_schmidt_extend(P, cands, atol=atol)
        added = P_new.shape[1] - P.shape[1]
        if added == 0:
            return P, HDs, gen_sizes
        gen_sizes.append(added)
        P = P_new
    return P, HDs, gen_sizes


def effective_L(N, n, gamma_0, P):
    """Project block L = D + J·M_H onto subspace spanned by columns of P.
    Returns (D_eff, M_H_eff)."""
    D, M_H = _block_L_split_total(N, n, gamma_0)
    D_eff = P.conj().T @ D @ P
    M_H_eff = P.conj().T @ M_H @ P
    return D_eff, M_H_eff


def find_Q_peak_from_split(D_eff, M_H_eff, gamma_0, J_grid, mode_observable=None):
    """Scan J on c×c (or larger) effective; find Q_peak via |dRe(λ_slow)/dQ| max.

    `mode_observable`, if given, is a c-vector encoding the Dicke probe in
    the effective basis. We use it to weight which eigenvalue's Re(λ)
    matters; otherwise we use the slowest mode.
    """
    Qs = J_grid / gamma_0
    proxy = np.zeros_like(J_grid)
    for i, J in enumerate(J_grid):
        L = D_eff + J * M_H_eff
        if mode_observable is None:
            evals = np.linalg.eigvals(L)
            proxy[i] = float(np.max(evals.real))
        else:
            evals, V = np.linalg.eig(L)
            # Weight Re(λ) by |⟨probe | v_j⟩|² (right-eigenvector overlap with probe)
            W = np.linalg.inv(V)
            cs = W @ mode_observable
            weights = np.abs(cs) ** 2
            proxy[i] = float(np.sum(weights * evals.real) / np.sum(weights))
    dQ = Qs[1] - Qs[0]
    deriv = np.gradient(proxy, dQ)
    abs_deriv = np.abs(deriv)
    peak_idx = int(np.argmax(abs_deriv))
    return Qs[peak_idx], abs_deriv[peak_idx], proxy


def main():
    gamma_0 = 0.05

    print(f"# EQ-022 (b1) Step (b): extended projection")
    print(f"# gamma_0 = {gamma_0}, J grid: Q in [0.05, 5.0]")
    print()

    test_cases = [
        # c = 2
        (4, 1), (5, 1), (6, 1), (7, 1),
        # c = 3
        (5, 2), (6, 2), (7, 2),
        # c = 4
        (7, 3), (8, 3),
    ]

    by_c = {}
    for (N, n) in test_cases:
        c = fw.chromaticity(N, n)
        t0 = time.time()
        P, HDs, gen_sizes = build_extended_subspace(N, n, gamma_0)
        D_eff, M_H_eff = effective_L(N, n, gamma_0, P)
        elapsed = time.time() - t0
        by_c.setdefault(c, []).append({
            "N": N, "n": n,
            "P": P, "HDs": HDs, "gen_sizes": gen_sizes,
            "D_eff": D_eff, "M_H_eff": M_H_eff,
            "elapsed": elapsed,
        })

    for c, runs in sorted(by_c.items()):
        print(f"## c = {c}")
        for r in runs:
            N, n = r["N"], r["n"]
            print(f"  (N={N}, n={n}): generations = {r['gen_sizes']}, "
                  f"total dim = {sum(r['gen_sizes'])}, built in {r['elapsed']:.2f} s")
        print()

    # Inspect c×c → 2c×2c → ... structure for c=4 case
    print("## L_eff structure at c=4, N=8, n=3")
    r = by_c[4][1]  # N=8, n=3
    print(f"  basis dim = {sum(r['gen_sizes'])}, generations = {r['gen_sizes']}")
    with np.printoptions(precision=4, suppress=True, linewidth=200):
        print("  D_eff / gamma_0 =")
        print(r["D_eff"] / gamma_0)
        print()
        print("  M_H_eff (J-coefficient) =")
        print(r["M_H_eff"])
    print()

    # Q_peak scan from the extended L_eff
    print("## Q_peak from extended L_eff (J-grid scan, slowest-eigenvalue proxy)")
    J_grid = np.linspace(0.005, 0.25, 246)  # Q in [0.1, 5.0], dQ ~ 0.02
    print(f"  J grid: Q in [{J_grid[0]/gamma_0:.2f}, {J_grid[-1]/gamma_0:.2f}], "
          f"dQ ≈ {(J_grid[1]-J_grid[0])/gamma_0:.3f}")
    print()
    for c, runs in sorted(by_c.items()):
        for r in runs:
            N, n = r["N"], r["n"]
            Q_peak, peak_val, _ = find_Q_peak_from_split(
                r["D_eff"], r["M_H_eff"], gamma_0, J_grid)
            print(f"  c={c}, N={N}, n={n}, dim={sum(r['gen_sizes'])}: "
                  f"Q_peak = {Q_peak:.3f}, peak |dRe(λ_slow)/dQ| = {peak_val:.4f}")
    print()

    # Q_peak with Dicke probe weighted observable
    print("## Q_peak weighted by Dicke probe (probe-weighted Re(λ))")
    for c, runs in sorted(by_c.items()):
        for r in runs:
            N, n = r["N"], r["n"]
            P = r["P"]
            rho_probe = fw.dicke_block_probe(N, n)
            probe_in_basis = P.conj().T @ rho_probe
            Q_peak, peak_val, _ = find_Q_peak_from_split(
                r["D_eff"], r["M_H_eff"], gamma_0, J_grid,
                mode_observable=probe_in_basis)
            probe_norm_in_basis = float(np.linalg.norm(probe_in_basis))
            print(f"  c={c}, N={N}, n={n}, dim={sum(r['gen_sizes'])}: "
                  f"Q_peak = {Q_peak:.3f}, |peak deriv| = {peak_val:.4f}, "
                  f"probe norm in basis = {probe_norm_in_basis:.4f} "
                  f"(full norm = {np.linalg.norm(rho_probe):.4f})")
    print()

    # external reference values (F86 in docs/ANALYTICAL_FORMULAS.md)
    print("## F86 reference: Q_peak(c=3) = 1.6, Q_peak(c=4) = Q_peak(c=5) = 1.8")
    print("## c=2 wobbles 1.4-1.6 (finite-size sensitive)")


if __name__ == "__main__":
    main()
