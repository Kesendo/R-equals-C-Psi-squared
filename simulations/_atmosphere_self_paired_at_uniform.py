"""(B) Sezieren the 24 self-F1-paired Re=-sigma modes at uniform gamma N=6.

At uniform gamma (eps=0), N=6: 24 self-F1-paired pairs at Re=-sigma=-0.3,
|Im|=0.15 (from Test 1). These are the F1-fixed-point set: lambda = -sigma
means the F1 partner -2sigma-lambda = -sigma = lambda itself.

For F1 (popcount-flip) (m_L, m_R) <-> (N-m_L, N-m_R), self-fixed points have
m_L = N-m_L and m_R = N-m_R, i.e. (m_L, m_R) = (N/2, N/2) = (3, 3) for N=6.
The 24 self-paired modes should live in the (3, 3) popcount block.

Connect to the Admixture-Lebensader reading: that doc says the slow mode at
(ceil(N/2), ceil(N/2)) is the physical Lebensader. Are these 24 modes the
slow-mode neighborhood? Their Re=-sigma is NOT the slow-mode position (slow
modes have small |Re|, not large |Re|). But they ARE in the same popcount block.

Method:
  1. Build L at eps=0, eigendecompose.
  2. Find 24 self-paired modes at Re=-sigma, |Im|=0.15.
  3. For each, project eigenvector onto:
     - Pauli basis (per-mode top Pauli strings)
     - joint-popcount basis ((m_L, m_R) breakdown)
  4. Same for the slow modes (smallest non-zero |Re|) for comparison.

Investigation only.
"""
import sys
import itertools
from collections import Counter

import numpy as np

sys.stdout.reconfigure(encoding="utf-8")

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
PAULIS = {'I': I2, 'X': X, 'Y': Y, 'Z': Z}
LABELS = ('I', 'X', 'Y', 'Z')
GAMMA0 = 0.05
J = 0.075
N = 6


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


def popcount(i):
    return bin(i).count('1')


def joint_popcount_distribution(v, N):
    """For vec(rho)-like vector v of length 4^N=d^2 with d=2^N, return dict
    {(m_L, m_R): |weight|^2} where m_L = popcount(bra index), m_R = popcount(ket index)."""
    d = 2 ** N
    v = np.asarray(v).reshape(d, d)
    dist = Counter()
    for i in range(d):
        for j in range(d):
            m_L = popcount(i)
            m_R = popcount(j)
            dist[(m_L, m_R)] += abs(v[i, j]) ** 2
    return dist


def pauli_basis(N):
    D = 4 ** N
    P = np.empty((D, D), dtype=complex)
    labels = []
    for idx, letters in enumerate(itertools.product(LABELS, repeat=N)):
        op = np.array([[1.0 + 0j]])
        for ch in letters:
            op = np.kron(op, PAULIS[ch])
        P[:, idx] = op.flatten() / np.sqrt(2 ** N)
        labels.append(''.join(letters))
    return P, labels


def xy_weight(lbl):
    return sum(1 for c in lbl if c in 'XY')


def yz_weight(lbl):
    return sum(1 for c in lbl if c in 'YZ')


def main():
    sigma = N * GAMMA0
    print(f"N={N}, J={J}, gamma0={GAMMA0}, sigma={sigma}, -sigma={-sigma}")

    gamma = np.full(N, GAMMA0)
    H = J * chain_H(N)
    L = build_L(H, gamma, N)
    print(f"L built ({L.shape}); diagonalizing...")
    ev, V = np.linalg.eig(L)
    re = ev.real
    im = ev.imag
    aim = np.abs(im)

    # find self-F1-paired Re=-sigma modes
    self_mask = (np.abs(re - (-sigma)) < 1e-4)
    n_self = int(self_mask.sum())
    print(f"\nself-F1-paired (Re=-sigma) modes: {n_self}")
    self_idx = np.where(self_mask)[0]
    self_aim = aim[self_idx]
    self_im = im[self_idx]
    # group by |Im|
    aim_round = np.round(self_aim, 5)
    aim_counts = Counter(aim_round.tolist())
    print(f"  |Im| distribution: {dict(sorted(aim_counts.items()))}")

    # focus on |Im|=0.15 modes (the dominant 24)
    target_aim = 0.15
    focus_mask = self_mask & (np.abs(aim - target_aim) < 1e-4)
    focus_idx = np.where(focus_mask)[0]
    print(f"\nfocus: {len(focus_idx)} modes at Re={-sigma}, |Im|={target_aim}")

    # also find slow modes: smallest non-zero |Re|
    nonzero_re = re[np.abs(re) > 1e-9]
    if nonzero_re.size:
        re_sorted_abs = np.sort(np.abs(nonzero_re))
        smallest_re = re_sorted_abs[0]
        print(f"\nsmallest non-zero |Re|: {smallest_re:.5e}")
        # expected from Admixture-Lebensader: 1.10*gamma*Q^2/N^2 = 1.10*0.05*(J/gamma)^2/36
        Q = J / GAMMA0
        predicted = 1.10 * GAMMA0 * Q ** 2 / N ** 2
        print(f"  Admixture-Lebensader prediction 1.10*gamma*Q^2/N^2 = {predicted:.5e}")
        # find modes near smallest_re
        slow_mask = (np.abs(re) > 1e-9) & (np.abs(re) < smallest_re * 1.1)
        slow_idx = np.where(slow_mask)[0]
        print(f"  modes within 10% of smallest_re: {len(slow_idx)}")
    else:
        slow_idx = np.array([], dtype=int)
        smallest_re = float('nan')

    # Build Pauli basis (vec form)
    print(f"\nbuilding Pauli basis ({4 ** N} strings)...")
    P, labels = pauli_basis(N)

    # ---------- Analyze focus (Re=-sigma, |Im|=0.15) modes ----------
    print(f"\n=== {len(focus_idx)} self-paired focus modes (Re=-sigma, |Im|={target_aim}) ===")
    # Joint popcount across all focus modes (sum weights)
    joint_total = Counter()
    for idx in focus_idx:
        v = V[:, idx]
        v_norm = v / np.linalg.norm(v)
        dist = joint_popcount_distribution(v_norm, N)
        for k, w in dist.items():
            joint_total[k] += w
    avg = {k: w / len(focus_idx) for k, w in joint_total.items()}
    print(f"  joint-popcount (m_L, m_R) distribution (avg across modes, normalized):")
    sorted_items = sorted(avg.items(), key=lambda kv: -kv[1])
    for (mL, mR), w in sorted_items[:8]:
        print(f"    (m_L={mL}, m_R={mR}): {w:.4f}")
    print(f"  total of top 8: {sum(w for _, w in sorted_items[:8]):.4f}")

    # Pauli content
    C = P.conj().T @ V[:, focus_idx]  # (D, len(focus_idx))
    total_weight = np.sum(np.abs(C) ** 2, axis=1)
    top_alpha = np.argsort(-total_weight)[:12]
    total_sum = float(total_weight.sum())
    print(f"\n  top 12 Pauli strings by total weight (of {total_sum:.2f}):")
    for a in top_alpha:
        w = float(total_weight[a])
        lbl = labels[a]
        print(f"    {lbl}: w={w:.4f}  xyW={xy_weight(lbl)}, yzW={yz_weight(lbl)}")

    # ---------- Analyze slow modes ----------
    if slow_idx.size:
        print(f"\n=== {len(slow_idx)} slow modes (smallest |Re|≈{smallest_re:.5e}) ===")
        joint_total_slow = Counter()
        for idx in slow_idx:
            v = V[:, idx]
            v_norm = v / np.linalg.norm(v)
            dist = joint_popcount_distribution(v_norm, N)
            for k, w in dist.items():
                joint_total_slow[k] += w
        avg_slow = {k: w / len(slow_idx) for k, w in joint_total_slow.items()}
        print(f"  joint-popcount distribution (avg):")
        for (mL, mR), w in sorted(avg_slow.items(), key=lambda kv: -kv[1])[:8]:
            print(f"    (m_L={mL}, m_R={mR}): {w:.4f}")

        C_slow = P.conj().T @ V[:, slow_idx]
        total_w_slow = np.sum(np.abs(C_slow) ** 2, axis=1)
        top_slow = np.argsort(-total_w_slow)[:12]
        print(f"\n  top 12 Pauli strings (of total {float(total_w_slow.sum()):.2f}):")
        for a in top_slow:
            w = float(total_w_slow[a])
            lbl = labels[a]
            print(f"    {lbl}: w={w:.4f}  xyW={xy_weight(lbl)}, yzW={yz_weight(lbl)}")


if __name__ == "__main__":
    main()
