#!/usr/bin/env python3
"""q_scale_n_scaling.py

Extends Q_SCALE_THREE_BANDS findings to N in {7, 8} via block-restricted
Liouvillian construction.

Background
----------
The Liouvillian L preserves popcount-block structure (H is U(1)-invariant,
Z-dephasing is popcount-diagonal). For the coherence probe

    rho_cc_init = (|S_n><S_{n+1}| + h.c.) / 2

only the (popcount=n, popcount=n+1) CC block is populated, and L maps
this block into itself. We can eigendecompose L restricted to this block
WITHOUT constructing the full d^2 x d^2 Liouvillian.

Block dimension: C(N, n) * C(N, n+1). At N=8, max block is 56 * 70 = 3920
(at n=3, n=4). This fits trivially in memory, while the full L at N=8
would be 65536 x 65536 = 64 GB.

Task context: TASK_Q_SCALE_N_SCALING.md (Sub-A: band tightening, Sub-B:
inner-richness-quench at c=4). This script runs the Q-scan at N=7 and
N=8, tabulates Q_onset / Q_peak / W_plateau, and outputs JSON.

Construction (algebra)
----------------------
For CC block |p><q| with popcount(p)=n, popcount(q)=n+1:

  H|p><q| - |p><q|H decomposes to sum over nearest-neighbor bonds where
  either p or q has (0,1) or (1,0) bit pattern at that bond:
    H|p><q|: bond b acts on p giving |p'><q| with coefficient J_b, where
             p' = p XOR (mask of bond b) iff p has (0,1) or (1,0) at bond.
    |p><q|H: symmetric for q.

  Dephasing: Z_i |p><q| Z_i = z_i(p)*z_i(q) |p><q| where z_i(x) = 1-2*bit_i(x).
    Sum over i: N - 2*D where D = Hamming(p, q).
    Net dephasing contribution: -2*gamma_0*D on diagonal.

So L_block[(p, q), (p, q)] = -2*gamma_0*D(p, q),
   L_block[(p', q), (p, q)] = -1j*J_b if p' = p^mask_b and p has (0,1)/(1,0) at b,
   L_block[(p, q'), (p, q)] = +1j*J_b if q' = q^mask_b and q has (0,1)/(1,0) at b.

Rules: em-dashes forbidden. UTF-8 stdout.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from math import comb
from pathlib import Path

import numpy as np
from scipy.linalg import eig, expm

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")


RESULTS_DIR = Path(__file__).parent / "results" / "q_scale_n_scaling"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# Block-restricted L construction
# ---------------------------------------------------------------------------
def popcount_states(N, n):
    """List of computational basis states with popcount n, big-endian."""
    return [x for x in range(2 ** N) if bin(x).count("1") == n]


def build_L_cc_block(N, n, J_list, gamma_0):
    """Construct L restricted to the (popcount=n, popcount=n+1) CC block.

    Basis: flat index k = p_idx * |P_{n+1}| + q_idx, where p in popcount-n
    states and q in popcount-(n+1) states. Returns (L, P_n, P_np1).
    """
    assert 0 <= n <= N - 1
    P_n = popcount_states(N, n)
    P_np1 = popcount_states(N, n + 1)
    Mn = len(P_n)
    Mnp1 = len(P_np1)
    M = Mn * Mnp1
    p_to_idx = {p: i for i, p in enumerate(P_n)}
    q_to_idx = {q: i for i, q in enumerate(P_np1)}

    def basis_idx(p, q):
        return p_to_idx[p] * Mnp1 + q_to_idx[q]

    L = np.zeros((M, M), dtype=complex)

    for p in P_n:
        for q in P_np1:
            i = basis_idx(p, q)

            # Dephasing (diagonal): -2 * gamma_0 * Hamming(p, q)
            D = bin(p ^ q).count("1")
            L[i, i] += -2.0 * gamma_0 * D

            # H acting from left: -i * H|p><q|
            # H_{bond b} sends |p> -> J_b |p'> where p' = p XOR mask(b) iff bits at bond differ.
            for bond in range(N - 1):
                # bit_b and bit_{b+1} of p, with site 0 = MSB (big-endian)
                bit_b = (p >> (N - 1 - bond)) & 1
                bit_bp1 = (p >> (N - 2 - bond)) & 1
                if bit_b != bit_bp1:
                    mask = (1 << (N - 1 - bond)) | (1 << (N - 2 - bond))
                    p_new = p ^ mask
                    # p_new has same popcount as p (bit swap)
                    L[basis_idx(p_new, q), i] += -1j * J_list[bond]

            # H acting from right: +i * |p><q|H
            # Symmetric: bond b acts on q.
            for bond in range(N - 1):
                bit_b = (q >> (N - 1 - bond)) & 1
                bit_bp1 = (q >> (N - 2 - bond)) & 1
                if bit_b != bit_bp1:
                    mask = (1 << (N - 1 - bond)) | (1 << (N - 2 - bond))
                    q_new = q ^ mask
                    L[basis_idx(p, q_new), i] += 1j * J_list[bond]

    return L, P_n, P_np1


# ---------------------------------------------------------------------------
# Dicke coherence probe in block basis
# ---------------------------------------------------------------------------
def dicke_cc_probe_flat(N, n):
    """Return rho_cc = (|S_n><S_{n+1}| + h.c.) / 2 restricted to the (n, n+1)
    CC block, as a flat vector in the block basis.

    Note: the Hermitian operator has support in both (n, n+1) and (n+1, n)
    blocks. For the spatial-sum S observable we work with the (n, n+1) block
    alone and treat the (n+1, n) part separately if needed. For W extraction
    we only need the (n, n+1) block part.

    |S_n> = (1/sqrt(C(N,n))) * Sum_{p: popcount=n} |p>
    |S_n><S_{n+1}| = (1 / sqrt(C(N,n) * C(N,n+1))) * Sum_{p, q} |p><q|

    So rho_cc restricted to (n, n+1) block has every basis element coefficient
    = 1/(2 * sqrt(C(N,n) * C(N,n+1))).
    """
    Mn = comb(N, n)
    Mnp1 = comb(N, n + 1)
    coeff = 1.0 / (2.0 * np.sqrt(Mn * Mnp1))
    M = Mn * Mnp1
    return np.full(M, coeff, dtype=complex)


def site_local_cc_probe_flat(N, n, site, spectator_bits):
    """Site-local spectator probe in (n, n+1) CC block.
    spectator_bits: list of n bits on the other N-1 sites (in ascending site index).
    """
    assert len(spectator_bits) == N - 1
    assert sum(spectator_bits) == n
    # Build the ket_0 (site bit = 0) and ket_1 (site bit = 1) computational basis states
    ket_0_idx = 0
    ket_1_idx = 0
    spec_iter = iter(spectator_bits)
    for s in range(N):
        if s == site:
            bit_0 = 0
            bit_1 = 1
        else:
            bit = next(spec_iter)
            bit_0 = bit
            bit_1 = bit
        ket_0_idx |= bit_0 << (N - 1 - s)
        ket_1_idx |= bit_1 << (N - 1 - s)

    # In (n, n+1) block: p = ket_0_idx (popcount n), q = ket_1_idx (popcount n+1)
    P_n = popcount_states(N, n)
    P_np1 = popcount_states(N, n + 1)
    p_to_idx = {p: i for i, p in enumerate(P_n)}
    q_to_idx = {q: i for i, q in enumerate(P_np1)}
    Mnp1 = len(P_np1)

    vec = np.zeros(len(P_n) * Mnp1, dtype=complex)
    idx = p_to_idx[ket_0_idx] * Mnp1 + q_to_idx[ket_1_idx]
    vec[idx] = 0.5  # coefficient of |p><q| in the Hermitian operator, one element only
    return vec


# ---------------------------------------------------------------------------
# Observable: S(t) spatial-sum coherence purity, restricted to block
# ---------------------------------------------------------------------------
def spatial_sum_coh_purity_from_block(rho_flat, N, n, P_n, P_np1):
    """S = Sum_i 2 * |(rho_i)_{0,1}|^2 for the given (n, n+1) block content.

    rho_flat is the (n, n+1) CC block coefficients in basis |p><q| order.
    After Tr_{~i}, the coherence (rho_i)_{0,1} picks up contributions only
    from |p><q| pairs where p and q differ exactly at site i (1-site-
    differing at site i, with p_i = 0 and q_i = 1).
    """
    Mnp1 = len(P_np1)
    p_to_idx = {p: i for i, p in enumerate(P_n)}
    q_to_idx = {q: i for i, q in enumerate(P_np1)}
    S = 0.0
    for i in range(N):
        # Coherence at site i: sum over (p, q) 1-site-differing at i
        coh = 0j
        for p in P_n:
            # q = p with bit i flipped (p_i = 0 -> q_i = 1)
            bit_i_p = (p >> (N - 1 - i)) & 1
            if bit_i_p != 0:
                continue  # p_i must be 0 for q_i = 1
            q = p | (1 << (N - 1 - i))
            # q must have popcount n+1 which it does (flipped 0 to 1)
            idx = p_to_idx[p] * Mnp1 + q_to_idx[q]
            coh += rho_flat[idx]
        # The h.c. part (|q><p| in the (n+1, n) block) contributes complex conjugate
        # to (rho_i)_{1,0}, which is the conjugate of (rho_i)_{0,1}.
        # For rho Hermitian, (rho_i)_{0,1} * + (rho_i)_{1,0} already accounts for both.
        # But we only store (n, n+1) block; the h.c. is in (n+1, n) block which gives
        # (rho_i)_{1,0} = coh^* in this convention.
        # For |rho_i)_{0,1}|^2 we take |coh|^2. Multiply by 2 for Hermitian factor.
        S += 2.0 * abs(coh) ** 2
    return S


# ---------------------------------------------------------------------------
# Propagation in block basis
# ---------------------------------------------------------------------------
def propagate_block(L_block, rho_init, times):
    """Propagate rho_init through L_block via step-wise expm(L*dt)."""
    if len(times) < 2:
        raise ValueError("Need at least 2 times.")
    dt = float(times[1] - times[0])
    assert np.allclose(np.diff(times), dt), "times must be uniform"
    U_dt = expm(L_block * dt)
    out = np.empty((len(times), len(rho_init)), dtype=complex)
    rho = rho_init.astype(complex).copy()
    for k in range(len(times)):
        if k > 0:
            rho = U_dt @ rho
        out[k] = rho
    return out


# ---------------------------------------------------------------------------
# W computation (dressed-mode weight) in block basis
# ---------------------------------------------------------------------------
def compute_W_block(L_block, rho_init, gamma_0, N_qubits, chromaticity,
                    pure_rate_tol_rel=0.015):
    """Eigendecompose L_block, project rho_init onto left-eigenvector basis,
    identify modes at pure dephasing rates 2*gamma_0*(2k+1) for k=0..c-1,
    return W = fraction of weight on non-pure (dressed) modes.
    """
    eigvals, V_R = eig(L_block)
    # Solve V_R @ c = rho_init for c (avoids explicit inverse)
    c_k = np.linalg.solve(V_R, rho_init)
    w = np.abs(c_k) ** 2
    total_w = w.sum()
    if total_w < 1e-20:
        return 0.0, eigvals, V_R, c_k

    rates = -eigvals.real
    pure_rates = [2 * gamma_0 * (2 * k + 1) for k in range(chromaticity)]
    is_pure = np.zeros_like(rates, dtype=bool)
    for pr in pure_rates:
        is_pure |= np.abs(rates - pr) < pure_rate_tol_rel * pr
    w_dressed = w[~is_pure & (w > 1e-12)].sum()
    return float(w_dressed / total_w), eigvals, V_R, c_k


# ---------------------------------------------------------------------------
# Main Q-scan entry points
# ---------------------------------------------------------------------------
DEFAULT_Q_GRID = [
    0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0,
    1.2, 1.4, 1.5, 1.6, 1.8, 2.0, 3.0, 5.0, 10.0, 20.0, 50.0,
]


def q_scan_block(N, n, gamma_0, J_base_q_one, Q_grid, verbose=True):
    """For each Q in Q_grid, build L_block at J = Q*gamma_0, compute W and
    abs(K_CC_pr)_max (via finite difference in J on S(t) trajectory).

    Returns dict with per-Q: W, abs_K_max, t_peak.
    """
    c = min(n, N - 1 - n) + 1
    chrom = c
    out = {}
    # Pre-compute probe
    rho_init = dicke_cc_probe_flat(N, n)
    # Time grid for K via propagation
    T_max = 100.0
    n_times = 201
    times = np.linspace(0.0, T_max, n_times)
    # dJ relative to J (0.05*J as in Q_SCALE_THREE_BANDS)
    for Q in Q_grid:
        J = Q * gamma_0
        J_list = [J] * (N - 1)
        t0 = time.time()

        # Build L_block for W
        L_block, P_n, P_np1 = build_L_cc_block(N, n, J_list, gamma_0)

        # W via eigendecomposition
        W, eigvals, V_R, c_k = compute_W_block(L_block, rho_init, gamma_0,
                                                 N, chrom)

        # abs(K_CC_pr)_max: finite difference of S(t) at bond = middle (interior)
        bond = (N - 1) // 2  # middle bond (N even: left-of-middle; N odd: exact)
        dJ = 0.05 * J if J > 0 else 1e-5
        J_plus = list(J_list); J_plus[bond] = J + dJ
        J_minus = list(J_list); J_minus[bond] = J - dJ
        L_plus, _, _ = build_L_cc_block(N, n, J_plus, gamma_0)
        L_minus, _, _ = build_L_cc_block(N, n, J_minus, gamma_0)

        rho_traj_plus = propagate_block(L_plus, rho_init, times)
        rho_traj_minus = propagate_block(L_minus, rho_init, times)
        S_plus = np.array([spatial_sum_coh_purity_from_block(r, N, n, P_n, P_np1)
                           for r in rho_traj_plus])
        S_minus = np.array([spatial_sum_coh_purity_from_block(r, N, n, P_n, P_np1)
                            for r in rho_traj_minus])
        K_traj = (S_plus - S_minus) / (2 * dJ)
        # Skip boundary indices
        abs_K_int = np.abs(K_traj[5:-5])
        idx_max = int(np.argmax(abs_K_int)) + 5
        abs_K_max = float(np.abs(K_traj[idx_max]))
        t_peak = float(times[idx_max])

        elapsed = time.time() - t0
        out[Q] = {
            "Q": Q, "J": J,
            "W": W,
            "abs_K_max": abs_K_max,
            "t_peak": t_peak,
            "elapsed_s": elapsed,
        }
        if verbose:
            print(f"  Q={Q:>6.3f}  J={J:>8.5f}  W={W:>.4f}  |K|_max={abs_K_max:>.4e}  "
                  f"t_peak={t_peak:>5.2f}  ({elapsed:.1f}s)")

    return out


def extract_characteristic_points(q_data):
    """From q_scan_block output, extract Q_onset, Q_peak (primary), W_peak,
    W_plateau.
    """
    Qs = sorted(q_data.keys())
    Ws = [q_data[q]["W"] for q in Qs]
    Ks = [q_data[q]["abs_K_max"] for q in Qs]

    # Q_onset: first Q with W > 0.05
    Q_onset = next((q for q, w in zip(Qs, Ws) if w > 0.05), None)

    # Q_peak: Q at first local maximum of W
    # Simplest: argmax over the grid
    idx_peak_w = int(np.argmax(Ws))
    Q_peak_w = Qs[idx_peak_w]
    W_peak = Ws[idx_peak_w]

    # Q where abs_K_max is maximal
    idx_peak_k = int(np.argmax(Ks))
    Q_peak_k = Qs[idx_peak_k]
    K_peak = Ks[idx_peak_k]

    # Plateau: W at largest Q
    W_plateau = Ws[-1]

    return {
        "Q_onset": Q_onset,
        "Q_peak_w": Q_peak_w,
        "W_peak": W_peak,
        "Q_peak_k": Q_peak_k,
        "abs_K_peak": K_peak,
        "W_plateau_at_Q_max": W_plateau,
        "Q_max": Qs[-1],
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--N", type=int, default=7)
    parser.add_argument("--n", type=str, default="all",
                        help="Comma-separated n values or 'all'")
    parser.add_argument("--gamma-0", type=float, default=0.05)
    parser.add_argument("--Q-grid", type=str, default=None,
                        help="Comma-separated Q values, overrides default grid")
    parser.add_argument("--tag", type=str, default=None)
    parser.add_argument("--sanity", action="store_true",
                        help="Run sanity check against N=5 n=1 existing data")
    args = parser.parse_args()

    if args.sanity:
        print("=" * 78)
        print("SANITY CHECK: N=5 n=1 via block-L vs existing data")
        print("=" * 78)
        N = 5
        n = 1
        gamma_0 = 0.05
        # Compare to existing values (from Q_SCALE_THREE_BANDS at gamma_0=0.05):
        # Q=1.0: W=0.666; Q=1.6: W=0.864; Q=20: W=0.490
        spot_checks = [(1.0, 0.666), (1.6, 0.864), (5.0, 0.672), (20.0, 0.490)]
        for Q, W_expected in spot_checks:
            J = Q * gamma_0
            J_list = [J] * (N - 1)
            L_block, P_n, P_np1 = build_L_cc_block(N, n, J_list, gamma_0)
            rho_init = dicke_cc_probe_flat(N, n)
            W, _, _, _ = compute_W_block(L_block, rho_init, gamma_0, N, chromaticity=2)
            err = abs(W - W_expected)
            status = "OK" if err < 0.005 else "MISMATCH"
            print(f"  Q={Q:>5.1f}  W_block={W:.4f}  W_expected={W_expected:.4f}  "
                  f"err={err:.4f}  [{status}]")

        # Also check the pure-rate structure at J=0
        L0, _, _ = build_L_cc_block(N, n, [0.0] * (N - 1), gamma_0)
        ev0, _ = eig(L0)
        rates0 = -ev0.real
        unique_rates = sorted(set(np.round(rates0, 4)))
        print(f"\n  J=0 active rates at N=5 n=1: {unique_rates}")
        print(f"  Expected (c=2): 2*gamma_0 = {2*gamma_0:.2f}, 6*gamma_0 = {6*gamma_0:.2f}")
        return 0

    N = args.N
    if args.n == "all":
        n_values = list(range(N))  # 0 to N-1
    else:
        n_values = [int(x) for x in args.n.split(",")]

    if args.Q_grid is None:
        Q_grid = DEFAULT_Q_GRID
    else:
        Q_grid = [float(x) for x in args.Q_grid.split(",")]

    tag = args.tag if args.tag else f"N{N}_gamma0_{args.gamma_0:.2f}".replace(".", "p")

    print("=" * 78)
    print(f"Q-scale scan at N={N}, gamma_0={args.gamma_0}, n in {n_values}")
    print(f"Q grid: {Q_grid}")
    print(f"Output tag: {tag}")
    print("=" * 78)

    all_results = {}
    t_start = time.time()
    for n in n_values:
        print(f"\n--- Block (n={n}, n+1={n+1}), c={min(n, N-1-n)+1} ---")
        q_data = q_scan_block(N, n, args.gamma_0, args.gamma_0, Q_grid, verbose=True)
        chars = extract_characteristic_points(q_data)
        print(f"  Summary: Q_onset={chars['Q_onset']}, Q_peak(W)={chars['Q_peak_w']}, "
              f"W_peak={chars['W_peak']:.4f}, Q_peak(|K|)={chars['Q_peak_k']}, "
              f"|K|_peak={chars['abs_K_peak']:.4f}, W_plateau(Q={chars['Q_max']})={chars['W_plateau_at_Q_max']:.4f}")
        all_results[f"n_{n}"] = {
            "q_data": q_data,
            "characteristic": chars,
            "chromaticity": min(n, N - 1 - n) + 1,
        }

    out = {
        "config": {
            "N": N,
            "gamma_0": args.gamma_0,
            "Q_grid": Q_grid,
            "n_values": n_values,
            "tag": tag,
        },
        "blocks": all_results,
        "elapsed_s": time.time() - t_start,
    }
    out_path = RESULTS_DIR / f"q_scan_{tag}.json"
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2, default=str)
    print(f"\nSaved: {out_path}")
    print(f"Total walltime: {out['elapsed_s']:.1f} s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
