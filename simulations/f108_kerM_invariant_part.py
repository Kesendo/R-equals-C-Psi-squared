"""F108 Part 1, the ker(M) remark: what the canonical Π still carries for YZ + ZY.

Verifier for the sentence in docs/proofs/PROOF_F108_PART1_PI2_EVEN_ALWAYS_PALINDROMIC.md ("Resolution of the
previously open ker(M) attempt"). Setting of that remark: open chain, H = sum over bonds of (YZ + ZY), uniform Z
dephasing; M = Π·L·Π⁻¹ + L + 2σ with the canonical Π. Here Π keeps L_H and reflects the dissipator, so M = 2·L_H:
ker M is H's commutant. The largest L-invariant subspace U inside it is the commutant's part that is homogeneous in
XY-weight (the dephasing is diagonal on the strings, rate 2γ per X or Y letter), and Π maps weight k to N − k inside it.

Exact throughout: L_H and the dissipator in the Pauli-string basis have Gaussian-integer entries at γ = 1 (Π is a
signed permutation with phases); every rank is a sympy rank of an integer matrix. Each check can fail.

Run: python simulations/f108_kerM_invariant_part.py
"""
from __future__ import annotations

import sys
from itertools import product
from pathlib import Path

import numpy as np
import sympy as sp

sys.path.insert(0, str(Path(__file__).parent))
from framework.symmetry import build_pi_full  # noqa: E402

RESULTS = []


def check(name, ok):
    RESULTS.append((name, bool(ok)))
    print(f"{'PASS' if ok else 'FAIL'}  {name}", flush=True)


PA = {'I': np.eye(2), 'X': np.array([[0, 1], [1, 0]]), 'Y': np.array([[0, -1j], [1j, 0]]), 'Z': np.diag([1, -1])}


def string(w):
    m = np.eye(1)
    for c in w:
        m = np.kron(m, PA[c])
    return m.astype(complex)


def words(N):
    return [''.join(w) for w in product('IXYZ', repeat=N)]         # the order of framework.symmetry.build_pi_full


def superop(N, f):
    """(F)_ab = Tr(σ_a f(σ_b)) / 2^N in the Pauli-string basis"""
    S = [string(w) for w in words(N)]
    F = np.zeros((len(S), len(S)), complex)
    for b, sb in enumerate(S):
        fb = f(sb)
        for a, sa in enumerate(S):
            F[a, b] = np.trace(sa @ fb) / 2 ** N
    return F


def exact(F):
    re, im = np.rint(F.real).astype(int), np.rint(F.imag).astype(int)
    assert np.abs(F.real - re).max() == 0 and np.abs(F.imag - im).max() == 0
    return re, im


def rank_c(F):
    """rank over C of a Gaussian-integer matrix, as half the rank over Q of its real 2x2 block form"""
    re, im = exact(F)
    return sp.Matrix(np.block([[re, -im], [im, re]])).rank() // 2


def hamiltonian(N, terms):
    return sum(string(''.join(t[0] if l == b else t[1] if l == b + 1 else 'I' for l in range(N)))
               for b in range(N - 1) for t in terms)


def weight(w):
    return sum(c in 'XY' for c in w)


def main():
    sys.stdout.reconfigure(encoding='utf-8')
    expected = {3: [1, 5, 5, 1], 4: [1, 3, 8, 3, 1], 5: [1, 3, 8, 8, 3, 1]}
    for N in (3, 4, 5):
        H = hamiltonian(N, ('YZ', 'ZY'))
        Zs = [string(''.join('Z' if l == k else 'I' for l in range(N))) for k in range(N)]
        LH = superop(N, lambda r: -1j * (H @ r - r @ H))
        D = superop(N, lambda r: sum(Z @ r @ Z - r for Z in Zs))
        Pi = build_pi_full(N)
        d = 4 ** N
        check(f'N={N} Π is unitary', np.abs(Pi @ Pi.conj().T - np.eye(d)).max() == 0)
        sig = N
        M = Pi @ (LH + D) @ Pi.conj().T + LH + D + 2 * sig * np.eye(d)
        check(f'N={N} M = 2·L_H exactly (Π keeps L_H for YZ + ZY and reflects the dephasing)',
              np.abs(M - 2 * LH).max() == 0 and np.abs(LH).max() > 0)
        # D is diagonal on the strings: the rate 2 per X or Y letter
        check(f'N={N} the dephasing is diagonal on the strings with -2·(XY-weight)',
              np.abs(D - np.diag([-2 * weight(w) for w in words(N)])).max() == 0)
        # U = the largest L-invariant subspace inside ker M = ker L_H: since D is diagonal and L_H vanishes on ker L_H,
        # U is the sum over k of (ker L_H ∩ span of the weight-k strings)
        per = []
        for k in range(N + 1):
            cols = [i for i, w in enumerate(words(N)) if weight(w) == k]
            per.append(len(cols) - rank_c(LH[:, cols]))
        check(f'N={N} U by XY-weight = {expected[N]} (sum {sum(expected[N])}), exact ranks', per == expected[N])
        # consistent with Π, which keeps L_H and maps weight k to N - k, mapping U onto itself
        check(f'N={N} the weight profile of U is symmetric under k ↔ N − k', per == per[::-1])
        if N == 3:
            # the whole kernel against its invariant part
            dimK = d - rank_c(LH)
            check(f'N=3 dim ker M = dim of the commutant of H = {dimK}, larger than dim U = {sum(per)} '
                  f'(the dephasing does not keep the whole commutant)', dimK == 24 and sum(per) == 12)
    # control: a truly Hamiltonian (XX + YY) has M = 0 at N = 3
    N = 3
    H = hamiltonian(N, ('XX', 'YY'))
    Zs = [string(''.join('Z' if l == k else 'I' for l in range(N))) for k in range(N)]
    L = superop(N, lambda r: -1j * (H @ r - r @ H)) + superop(N, lambda r: sum(Z @ r @ Z - r for Z in Zs))
    Pi = build_pi_full(N)
    check('control N=3: for XX + YY, M = 0 (truly)', np.abs(Pi @ L @ Pi.conj().T + L + 2 * N * np.eye(64)).max() == 0)
    # control: a Hamiltonian outside this case, where the residual is not 2·L_H
    H2 = hamiltonian(N, ('XY', 'YX'))
    L2 = superop(N, lambda r: -1j * (H2 @ r - r @ H2)) + superop(N, lambda r: sum(Z @ r @ Z - r for Z in Zs))
    LH2 = superop(N, lambda r: -1j * (H2 @ r - r @ H2))
    M2 = Pi @ L2 @ Pi.conj().T + L2 + 2 * N * np.eye(64)
    check('control N=3: for XY + YX the residual is not 2·L_H', np.abs(M2 - 2 * LH2).max() > 0)

    fails = [nm for nm, ok in RESULTS if not ok]
    print(f'\n{len(RESULTS) - len(fails)}/{len(RESULTS)} PASS')
    return 1 if fails else 0


if __name__ == '__main__':
    sys.exit(main())
