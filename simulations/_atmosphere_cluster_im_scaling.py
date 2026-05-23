"""Phase A: scaling sweep of cluster |Im| at uniform gamma.

Hypothesis (dimensional analysis): for uniform gamma, L = J*(-i*H_dim + (1/Q)*D_dim)
so Im(L_eigenvalue) = J * Im(L_dim(Q)) where Q = J/gamma. Equivalently,
|Im|_cluster / gamma is a function of Q and N only.

Test:
  1. Vary gamma at fixed J: see how |Im| varies with Q.
  2. Vary J at fixed gamma: same Q sweep, should give consistent function f(Q, N).
  3. Vary N at fixed Q: see N-dependence of |Im|/gamma at given Q.

Identify the cluster as the 32-mode (N=6) / 24-mode (N=5) / smaller family at the
smallest non-real |Im| with F1-paired Re sub-clusters around -sigma.

Investigation only.
"""
import sys
import math
from collections import Counter

import numpy as np

sys.stdout.reconfigure(encoding="utf-8")

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)


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


def find_cluster_im(N, J, gamma0, expected_cluster_size=None):
    """Identify the cluster |Im| at uniform gamma.

    The cluster is the family of F1-mirror-paired modes at the smallest non-zero
    |Im| with Re sub-clusters at -sigma +/- delta. We find it by:
      1. Build L at uniform gamma.
      2. Group eigenvalues by |Im| (rounded), find clusters of 8+ pairs at given |Im|.
      3. Among clusters with two distinct Re values summing to -2*sigma, find the
         one with smallest |Im|.
    """
    H = J * chain_H(N)
    gamma = np.full(N, gamma0)
    sigma = float(np.sum(gamma))
    L = build_L(H, gamma, N)
    ev = np.linalg.eigvals(L)
    aim = np.abs(ev.imag)
    re = ev.real

    # Filter: non-real (|Im| > 1e-10)
    nonreal_mask = aim > 1e-10
    aim_nz = aim[nonreal_mask]
    re_nz = re[nonreal_mask]

    # Sort by |Im|, group
    aim_round = np.round(aim_nz, 6)
    # For each unique |Im|, find Re values and count
    unique_im = sorted(set(aim_round.tolist()))
    for im_val in unique_im:
        same_im_mask = aim_round == im_val
        re_at_im = sorted(set(np.round(re_nz[same_im_mask], 4).tolist()))
        if len(re_at_im) >= 2:
            # Check F1-mirror pair
            # find two Re values symmetric around -sigma
            for r_lo in re_at_im:
                r_hi = round(-2 * sigma - r_lo, 4)
                if r_hi in re_at_im and r_lo < r_hi:
                    count_lo = int(np.sum((aim_round == im_val) & (np.round(re_nz, 4) == r_lo)))
                    count_hi = int(np.sum((aim_round == im_val) & (np.round(re_nz, 4) == r_hi)))
                    total = count_lo + count_hi
                    # Skip the trivial |Im|=0 cluster (we want non-real cluster)
                    # We accept clusters of count >= 8 (smaller than our 32-mode family
                    # for smaller N or different cluster families)
                    if total >= 8:
                        return im_val, sigma, r_lo, r_hi, total
    return None, sigma, None, None, 0


def main():
    print("Phase A: cluster |Im| scaling sweep at uniform gamma\n")

    # Test 1: vary gamma at fixed J, fixed N=6
    print("=== Test 1: vary gamma0 (J=0.075 fixed), N=6 ===")
    J_fixed = 0.075
    gamma_list = [0.025, 0.05, 0.075, 0.1, 0.15, 0.2]
    print(f"{'gamma0':>8} {'Q=J/g0':>8} {'|Im|':>14} {'|Im|/gamma':>14} {'|Im|/J':>14} {'sigma':>10} {'(Re_lo, Re_hi)':>22} {'count':>6}")
    print("-" * 110)
    results_t1 = []
    for g0 in gamma_list:
        im_val, sigma, r_lo, r_hi, count = find_cluster_im(6, J_fixed, g0)
        if im_val is None:
            print(f"{g0:>8.4f}  no cluster found")
            continue
        Q = J_fixed / g0
        print(f"{g0:>8.4f} {Q:>8.4f} {im_val:>14.6e} {im_val/g0:>14.6f} {im_val/J_fixed:>14.6f} "
              f"{sigma:>10.4f} {f'({r_lo:+.4f}, {r_hi:+.4f})':>22} {count:>6d}")
        results_t1.append((g0, Q, im_val, sigma, r_lo, r_hi, count))
        sys.stdout.flush()

    # Test 2: vary N at moderate Q (Q=1.5) with fixed gamma=0.05
    print(f"\n=== Test 2: vary N (J chosen so Q=1.5, gamma0=0.05 fixed) ===")
    gamma_fixed = 0.05
    Q_fixed = 1.5
    print(f"{'N':>3} {'J':>8} {'|Im|':>14} {'|Im|/gamma':>14} {'|Im|/J':>14} {'sigma':>10} {'(Re_lo, Re_hi)':>22} {'count':>6}")
    print("-" * 100)
    for N in [4, 5, 6]:
        J = Q_fixed * gamma_fixed
        im_val, sigma, r_lo, r_hi, count = find_cluster_im(N, J, gamma_fixed)
        if im_val is None:
            print(f"N={N:>2d}: no cluster found")
            continue
        print(f"{N:>3d} {J:>8.4f} {im_val:>14.6e} {im_val/gamma_fixed:>14.6f} {im_val/J:>14.6f} "
              f"{sigma:>10.4f} {f'({r_lo:+.4f}, {r_hi:+.4f})':>22} {count:>6d}")
        sys.stdout.flush()

    # Test 3: vary J at fixed gamma=0.05, N=6 (same Q-sweep but along J-axis)
    print(f"\n=== Test 3: vary J (gamma0=0.05 fixed), N=6 ===")
    J_list = [0.025, 0.05, 0.075, 0.1, 0.15, 0.2]
    gamma_fixed = 0.05
    print(f"{'J':>8} {'Q=J/g':>8} {'|Im|':>14} {'|Im|/gamma':>14} {'|Im|/J':>14} {'sigma':>10} {'(Re_lo, Re_hi)':>22} {'count':>6}")
    print("-" * 110)
    for J in J_list:
        im_val, sigma, r_lo, r_hi, count = find_cluster_im(6, J, gamma_fixed)
        if im_val is None:
            print(f"J={J:>8.4f}: no cluster found")
            continue
        Q = J / gamma_fixed
        print(f"{J:>8.4f} {Q:>8.4f} {im_val:>14.6e} {im_val/gamma_fixed:>14.6f} {im_val/J:>14.6f} "
              f"{sigma:>10.4f} {f'({r_lo:+.4f}, {r_hi:+.4f})':>22} {count:>6d}")
        sys.stdout.flush()

    # Sanity check: does (|Im|/gamma) computed via Test 1 (var gamma, fixed J) match
    # the (|Im|/gamma) at same Q via Test 3 (var J, fixed gamma)?
    print(f"\n=== Sanity: |Im|/gamma is Q-function only? ===")
    print("If dimensional analysis holds, |Im|/gamma at same Q must match across Tests 1 & 3.")
    print("Compare Test 1 (g0=0.05, Q=1.5) vs Test 3 (J=0.075, Q=1.5): should be identical (already same point).")
    # Already covered: g0=0.05/J=0.075 is the same point in both tests; Q=1.5 in both.

    # Try a few Q values that appear in both tests:
    # Test 1: gamma in {0.025, 0.05, 0.075, 0.1, 0.15, 0.2} with J=0.075 -> Q in {3, 1.5, 1, 0.75, 0.5, 0.375}
    # Test 3: J in {0.025, 0.05, 0.075, 0.1, 0.15, 0.2} with gamma=0.05 -> Q in {0.5, 1, 1.5, 2, 3, 4}
    # Common Q values: 0.5, 1, 1.5, 3 -- check at these Q values that |Im|/gamma matches.


if __name__ == "__main__":
    main()
