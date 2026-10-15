"""Every XY-weight block of the palindromic residual, not only the diagonal.

M = Π·L·Π⁻¹ + L + 2Σγ·I in the Pauli-string basis (Π per site: I↔X sign 1,
Y→iZ, Z→iY; Z-dephasing γ(ZρZ − ρ)), blocked by XY-weight w. The zero-immunity
theorem (docs/proofs/PROOF_ZERO_IMMUNITY.md) is about the diagonal blocks
(0,0) and (N,N); this script reads all blocks, for the numbers quoted in
experiments/V_EFFECT_BOUNDARY_LOCALIZATION.md and in the proof's section
"Where the theorem stops":

  1. Tests 3 and 5 of veffect_weight_immunity.py: the break off the diagonal.
  2. A single-site Z field and X field: the extreme blocks.
  3. YZ on every bond: no diagonal block at all.
  4. One unit XY bond: ‖M|(w,w)‖² = 32·C(N−2, w−1)·2^(N−2).
  5. The 36 two-term census of V_EFFECT_FINE_STRUCTURE at N = 3.
  6. XY − YX and XX + YY + XY on every bond, N = 3: the diagonal breaks.

The superoperator is built directly on Pauli strings (products by table), so
every Hamiltonian entry is ±2 or ±2i; the dissipator cancels from M site by site
at any rate, and γ = 0.25 (dyadic) keeps that cancellation exact in floating
point, so a zero block is compared as exactly 0.0.
"""
import itertools
import sys
from math import comb

import numpy as np

_P = [np.eye(2), np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]]), np.diag([1, -1])]
MUL = {}
for _a in range(4):
    for _b in range(4):
        _m = _P[_a] @ _P[_b]
        for _c in range(4):
            _ph = np.trace(_P[_c].conj().T @ _m) / 2
            if abs(_ph) > 1e-12:
                MUL[(_a, _b)] = (_c, complex(np.round(_ph, 12)))
PI_PERM = {0: 1, 1: 0, 2: 3, 3: 2}
PI_SIGN = {0: 1, 1: 1, 2: 1j, 3: 1j}


def _smul(s, t):
    ph, out = 1, []
    for a, b in zip(s, t):
        c, p = MUL[(a, b)]
        out.append(c)
        ph *= p
    return tuple(out), ph


def residual_blocks(N, terms, gamma=0.25):
    """terms: list of (string tuple, coefficient). Returns {(w_row, w_col): Frobenius norm}."""
    S = list(itertools.product(range(4), repeat=N))
    idx = {s: i for i, s in enumerate(S)}
    D = len(S)
    L = np.zeros((D, D), complex)
    for h, c in terms:
        for s in S:
            hs, p1 = _smul(h, s)
            _, p2 = _smul(s, h)
            v = -1j * c * (p1 - p2)
            if v != 0:
                L[idx[hs], idx[s]] += v
    for s in S:
        L[idx[s], idx[s]] += -2 * gamma * sum(1 for a in s if a in (1, 2))
    Pi = np.zeros((D, D), complex)
    for s in S:
        Pi[idx[tuple(PI_PERM[a] for a in s)], idx[s]] = np.prod([PI_SIGN[a] for a in s])
    Pi_inv = Pi.conj().T  # Π is unitary in this basis (a signed permutation)
    M = Pi @ L @ Pi_inv + L + 2 * N * gamma * np.eye(D)
    w = np.array([sum(1 for a in s if a in (1, 2)) for s in S])
    return {(a, b): float(np.linalg.norm(M[np.ix_(w == a, w == b)]))
            for a in range(N + 1) for b in range(N + 1)}


def bonds_terms(N, spec):
    """spec: list of ((i, j), 'PQ', coefficient)."""
    out = []
    for (i, j), lab, c in spec:
        s = [0] * N
        s[i] = 'IXYZ'.index(lab[0])
        s[j] = 'IXYZ'.index(lab[1])
        out.append((tuple(s), c))
    return out


def nonzero(blocks):
    return {k: round(v, 6) for k, v in blocks.items() if v != 0.0}


def main():
    print("1. Tests 3 and 5, every nonzero block")
    for N in (3, 4):
        for name, spec in [("Test 3: XX on (0,1), XY on (1,2)", [((0, 1), 'XX', 1), ((1, 2), 'XY', 1)]),
                           ("Test 5: ZZ on (0,1), XY on (1,2)", [((0, 1), 'ZZ', 1), ((1, 2), 'XY', 1)])]:
            print(f"  N={N} {name}: {nonzero(residual_blocks(N, bonds_terms(N, spec)))}")

    print("\n2. One unit field on site 0: extreme diagonal blocks (0,0), (N,N) and the largest block")
    for N in (3, 4, 5):
        for lab in ('ZI', 'XI'):
            b = residual_blocks(N, bonds_terms(N, [((0, 1), lab, 1)]))
            print(f"  N={N} {lab[0]} field: (0,0) {b[(0, 0)]:.6f}  ({N},{N}) {b[(N, N)]:.6f}  max {max(b.values()):.2e}")

    print("\n3. YZ on every bond, N=3: every nonzero block")
    print(f"  {nonzero(residual_blocks(3, bonds_terms(3, [((0, 1), 'YZ', 1), ((1, 2), 'YZ', 1)])))}")

    print("\n4. One unit XY bond on (1,2): ‖M|(w,w)‖² against 32·C(N−2, w−1)·2^(N−2)")
    for N in (3, 4, 5):
        b = residual_blocks(N, bonds_terms(N, [((1, 2), 'XY', 1)]))
        got = [round(b[(w, w)] ** 2, 9) for w in range(N + 1)]
        law = [32 * comb(N - 2, w - 1) * 2 ** (N - 2) if 1 <= w <= N - 1 else 0 for w in range(N + 1)]
        print(f"  N={N}: {got}  law {law}  {'equal' if got == law else 'DIFFERENT'}")

    print("\n5. The 36 two-term census at N=3 (terms on both bonds)")
    labs = ['XX', 'XY', 'XZ', 'YX', 'YY', 'YZ', 'ZX', 'ZY', 'ZZ']
    broken, clean, reach = 0, [], 0
    for a, b in itertools.combinations(labs, 2):
        spec = [(bd, lab, 1) for lab in (a, b) for bd in ((0, 1), (1, 2))]
        blk = residual_blocks(3, bonds_terms(3, spec))
        if max(blk.values()) == 0.0:
            continue
        broken += 1
        if max(blk[(w, w)] for w in range(4)) == 0.0:
            clean.append(f"{a}+{b}")
        if any(v != 0.0 for k, v in blk.items() if 0 in k or 3 in k):
            reach += 1
    print(f"  broken {broken}; reaching w=0 or w=3 off the diagonal {reach}; clean diagonal {len(clean)}: {clean}")

    print("\n6. Two-letter combinations on every bond, N=3: the diagonal blocks")
    for name, spec in [("XY − YX", [('XY', 1), ('YX', -1)]), ("XX + YY + XY", [('XX', 1), ('YY', 1), ('XY', 1)])]:
        blk = residual_blocks(3, bonds_terms(3, [(bd, lab, c) for lab, c in spec for bd in ((0, 1), (1, 2))]))
        print(f"  {name}: diagonal {[round(blk[(w, w)], 6) for w in range(4)]}")


if __name__ == "__main__":
    if sys.platform == "win32":
        sys.stdout.reconfigure(encoding="utf-8")
    main()
