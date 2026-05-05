"""A3 (EQ-022 Item 1' c=2 derivation) — exploration record.

Finding: σ_0 of V_inter = P_HD1† · M_H_total · P_HD3 is EXACTLY degenerate at
even N (deg=2 at N=6, N=8) and 1-dim at odd N (deg=1 at N=5, N=7).

Chain-mirror R splits the 2D top eigenspace into one R-even and one R-odd
direction with the SAME σ_0. Library-specific SVD direction inside that 2D
subspace is the obstruction blocking single-vector Tier1Derived.

This is the test bed for whoever picks up A3 next session; promotion strategy
is in `compute/RCPsiSquared.Core/F86/Item1Derivation/C2InterChannelAnalytical.cs`
PendingDerivationNote.
"""
import numpy as np


def popcount_states(N, k):
    return sorted([s for s in range(1 << N) if bin(s).count('1') == k])


def hd(p, q):
    return bin(p ^ q).count('1')


def bit_reverse(s, N):
    r = 0
    for i in range(N):
        if (s >> i) & 1:
            r |= 1 << (N - 1 - i)
    return r


def bond_flip_targets(state, N):
    out = []
    for b in range(N - 1):
        bit_hi = (state >> (N - 1 - b)) & 1
        bit_lo = (state >> (N - 2 - b)) & 1
        if bit_hi != bit_lo:
            mask = (1 << (N - 1 - b)) | (1 << (N - 2 - b))
            out.append((b, state ^ mask))
    return out


def build_inter_channel(N, n=1):
    states_p = popcount_states(N, n)
    states_q = popcount_states(N, n + 1)
    Mp = len(states_p)
    Mq = len(states_q)
    p_idx = {p: i for i, p in enumerate(states_p)}
    q_idx = {q: i for i, q in enumerate(states_q)}

    def flat(p, q):
        return p_idx[p] * Mq + q_idx[q]

    Mtot = Mp * Mq
    M_H = np.zeros((Mtot, Mtot), dtype=complex)
    for p in states_p:
        for q in states_q:
            i = flat(p, q)
            for bond, p_flip in bond_flip_targets(p, N):
                if p_flip in p_idx:
                    j = flat(p_flip, q)
                    M_H[j, i] += -1j
            for bond, q_flip in bond_flip_targets(q, N):
                if q_flip in q_idx:
                    j = flat(p, q_flip)
                    M_H[j, i] += +1j

    hd1_pairs = [(p, q) for p in states_p for q in states_q if hd(p, q) == 1]
    hd3_pairs = [(p, q) for p in states_p for q in states_q if hd(p, q) == 3]

    V_inter = np.zeros((len(hd1_pairs), len(hd3_pairs)), dtype=complex)
    for i, (p1, q1) in enumerate(hd1_pairs):
        for j, (p3, q3) in enumerate(hd3_pairs):
            V_inter[i, j] = M_H[flat(p1, q1), flat(p3, q3)]

    return hd1_pairs, hd3_pairs, V_inter


def chain_mirror_action(hd1_pairs, N):
    pmap = {(p, q): i for i, (p, q) in enumerate(hd1_pairs)}
    perm = np.zeros(len(hd1_pairs), dtype=int)
    for i, (p, q) in enumerate(hd1_pairs):
        rp, rq = bit_reverse(p, N), bit_reverse(q, N)
        perm[i] = pmap[(rp, rq)]
    return perm


print("A3 finding: top-σ degeneracy of V V† per N for c=2 (n=1):")
print()
print(f"{'N':>3s} {'σ_0':>12s} {'σ_0/(2√2)':>10s} {'deg_top':>8s}  R-parity hint")
for N in [5, 6, 7, 8]:
    hd1_pairs, hd3_pairs, V_inter = build_inter_channel(N)
    VVH = V_inter @ V_inter.conj().T
    eigvals, eigvecs = np.linalg.eigh(VVH)
    eigvals_desc = eigvals[::-1]
    top = eigvals_desc[0]
    deg = sum(1 for e in eigvals_desc if abs(e - top) / max(top, 1e-15) < 1e-8)

    # R-parity of top eigenvectors
    perm = chain_mirror_action(hd1_pairs, N)
    parity_info = ""
    for k in range(deg):
        idx = len(eigvals) - 1 - k
        v = eigvecs[:, idx]
        Rv = v[perm]
        ip = np.vdot(v, Rv).real
        parity_info += f" v_{k}<R>={ip:+.2f}"

    print(f"{N:>3d}  {np.sqrt(top):>10.6f}  {np.sqrt(top)/(2*np.sqrt(2)):>9.6f}  {deg:>6d}  {parity_info}")

print()
print("Conclusion: at even N, the top eigenspace is 2D and contains both R-even and R-odd")
print("directions with the same σ_0. SVD inside this subspace picks a library-specific")
print("direction, blocking single-vector Tier1Derived against MathNet's reference.")
print("Promotion strategy: lift test contract to projector-onto-2D-eigenspace overlap")
print("(library-independent), or perform Bogoliubov / JW reduction of the 1↔2 channel.")
