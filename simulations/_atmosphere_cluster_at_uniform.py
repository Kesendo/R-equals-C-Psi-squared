"""Test 1: does the 16+16 (N=6) / 12+12 (N=5) cluster exist at uniform gamma?

If YES (a 16+16 F1-mirror-pair degenerate L-eigenspace exists already at eps=0),
the cluster is an intrinsic Lebensader-protected family; non-uniform palindromic
gamma only SWEEPS it. The atmosphere knob doesn't create the structure, it
displaces it.

If NO (no such cluster at eps=0), the cluster is itself a non-uniform-gamma
construction — the F1 algebra alone doesn't pre-organize it.

Method: at eps=0 build L, eigendecompose, find F1-mirror pairs (eigenvalues
(lambda, -2sigma-lambda)), then look for clusters of 16+16 (N=6) or 12+12 (N=5)
F1-pairs at common |Im|.

For comparison, also report the eps=eps* case.

Investigation only.
"""
import sys
from collections import Counter

import numpy as np

sys.stdout.reconfigure(encoding="utf-8")

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
GAMMA0 = 0.05
J = 0.075


def site_op(op, k, N):
    m = np.array([[1.0 + 0j]])
    for i in range(N):
        m = np.kron(m, op if i == k else I2)
    return m


def chain_H(N):
    d = 2 ** N
    H = np.zeros((d, d), dtype=complex)
    for b in range(N - 1):
        for P in (X, Y):
            H += site_op(P, b, N) @ site_op(P, b + 1, N)
    return H


def build_L(H, gamma, N):
    d = 2 ** N
    Id = np.eye(d, dtype=complex)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for k in range(N):
        Zk = site_op(Z, k, N)
        L += gamma[k] * (np.kron(Zk, Zk.conj()) - np.eye(d * d, dtype=complex))
    return L


def sym_shape(N):
    i = np.arange(N, dtype=float)
    u = (i - (N - 1) / 2.0) ** 2
    u = u - u.mean()
    return u / np.max(np.abs(u))


def analyze(N, gamma, label):
    """Build L, eigendecompose, characterize cluster structure."""
    print(f"\n--- {label} ---")
    print(f"gamma = {gamma}")
    sigma = float(np.sum(gamma))
    print(f"sigma = sum(gamma) = {sigma:.5f}, -2*sigma = {-2*sigma:.5f}")

    H = J * chain_H(N)
    L = build_L(H, gamma, N)
    ev = np.linalg.eigvals(L)

    # Find F1-mirror pairs: for each lambda, look for -2sigma - lambda
    target_sum = -2.0 * sigma
    n_modes = len(ev)
    paired = np.zeros(n_modes, dtype=bool)
    pairs = []
    for i in range(n_modes):
        if paired[i]:
            continue
        partner_target = target_sum - ev[i]
        # search for closest unpaired partner
        diffs = np.abs(ev - partner_target)
        diffs[paired] = np.inf
        j = int(np.argmin(diffs))
        if j != i and diffs[j] < 1e-6:
            pairs.append((i, j, ev[i], ev[j]))
            paired[i] = True
            paired[j] = True

    # Self-paired modes: lambda = -sigma (own F1 partner)
    self_paired = [i for i in range(n_modes) if not paired[i]]
    print(f"F1-pair structure: {len(pairs)} mirror pairs, {len(self_paired)} self-paired (lambda = -sigma)")

    # Look for clusters: pairs whose |Im(lambda)| coincides
    # Group pairs by quantized |Im|, count
    pair_aims = np.array([abs(ev[i].imag) for i, j, _, _ in pairs])
    # quantize
    aim_round = np.round(pair_aims, 6)
    aim_count = Counter(aim_round.tolist())
    # Find big clusters
    print(f"\nlargest |Im|-clusters (pair count, |Im|):")
    sorted_clusters = sorted(aim_count.items(), key=lambda kv: -kv[1])
    for aim_val, count in sorted_clusters[:15]:
        if count >= 4:
            # Get the Re values of these pairs
            relevant_pairs = [p for p in pairs if abs(abs(p[2].imag) - aim_val) < 5e-7]
            re_vals = sorted(set(round(p[2].real, 4) for p in relevant_pairs))
            re_str = ', '.join(f'{r:+.4f}' for r in re_vals[:6])
            print(f"  {count:3d} pairs at |Im|={aim_val:.6e}, Re values: {re_str}")
    return ev, sigma, pairs


def find_F1_re_pairs(ev, sigma, n_decimals=4):
    """Group eigenvalues by (round(|Im|, n_dec), round(Re, n_dec)) and find
    F1 Re-pairs."""
    aim_round = np.round(np.abs(ev.imag), n_decimals)
    re_round = np.round(ev.real, n_decimals)
    # All (Re, |Im|) eigenvalue classes with multiplicity
    classes = Counter(zip(re_round.tolist(), aim_round.tolist()))
    # Find Re pairs summing to -2*sigma
    target_sum = round(-2.0 * sigma, n_decimals - 1)
    re_pairs = []
    for (re_a, aim), n_a in classes.items():
        re_b = round(target_sum - re_a, n_decimals)
        if (re_b, aim) in classes and re_a < re_b:
            n_b = classes[(re_b, aim)]
            re_pairs.append(((re_a, re_b), aim, n_a, n_b))
    return re_pairs


def main():
    for N, eps_star in ((5, -0.3357), (6, -0.83)):
        print(f"\n========== N={N} ==========")
        d = 2 ** N

        # Case 1: uniform gamma (eps=0)
        gamma_uniform = np.full(N, GAMMA0)
        ev0, sigma0, pairs0 = analyze(N, gamma_uniform, f"eps=0 (uniform gamma)")

        # Look for symmetric Re-pair clusters around -sigma
        print(f"\n  symmetric Re-pair clusters around -sigma = {-sigma0:.5f}:")
        re_pairs0 = find_F1_re_pairs(ev0, sigma0, n_decimals=4)
        # Sort by count
        re_pairs0.sort(key=lambda x: -(x[2] + x[3]))
        for (re_a, re_b), aim, n_a, n_b in re_pairs0[:10]:
            if n_a + n_b >= 8:
                print(f"    Re=[{re_a:+.4f}, {re_b:+.4f}] (delta={(re_b-re_a)/2:+.4f}), "
                      f"|Im|={aim:.5f}, counts: {n_a}+{n_b}")

        # Case 2: eps=eps* (dip)
        u = sym_shape(N)
        gamma_dip = GAMMA0 * (1.0 + eps_star * u)
        ev1, sigma1, pairs1 = analyze(N, gamma_dip, f"eps={eps_star} (dip)")

        print(f"\n  symmetric Re-pair clusters around -sigma = {-sigma1:.5f}:")
        re_pairs1 = find_F1_re_pairs(ev1, sigma1, n_decimals=4)
        re_pairs1.sort(key=lambda x: -(x[2] + x[3]))
        for (re_a, re_b), aim, n_a, n_b in re_pairs1[:10]:
            if n_a + n_b >= 8:
                print(f"    Re=[{re_a:+.4f}, {re_b:+.4f}] (delta={(re_b-re_a)/2:+.4f}), "
                      f"|Im|={aim:.5f}, counts: {n_a}+{n_b}")


if __name__ == "__main__":
    main()
