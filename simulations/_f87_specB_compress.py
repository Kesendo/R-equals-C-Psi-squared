#!/usr/bin/env python3
"""F87 spec-B probe 17: the omega=0 block as a COMPRESSION of the 'agreement operator' T.

NEW STRUCTURAL READING. The first-order omega=0 block is  Q := M_0 + N = P_0 T P_0, where
  - T := sum_l Ad_{Z_l}  (rho -> sum_l Z_l rho Z_l), the 'agreement operator'. On |i><j|,
        T|i><j| = (sum_l z_l(i) z_l(j)) |i><j| = (N - 2*popcount(i^j)) |i><j|.
    T is Hermitian, DIAGONAL in the computational basis, spectrum {N-2m : m=0..N}, SYMMETRIC about 0.
  - P_0 projects vec-space onto the omega=0 sector = the COMMUTANT of H (modes |E_a><E_b|, E_a=E_b),
    i.e. operators commuting with H (the centralizer algebra C(H)).
So  soft  <=>  spec(P_0 T P_0) symmetric about 0.

T's symmetry about 0 is implemented by the global-flip involution on ONE side:
  F: |i><j| -> |i><j_bar|,  j_bar = j XOR (11..1)  (right-multiply by X^{otimes N}).
  popcount(i ^ j_bar) = N - popcount(i^j)  =>  T eigenvalue N-2m -> -(N-2m).  So  F T F^-1 = -T.
Hence if the commutant subspace (range of P_0) is F-INVARIANT, then F restricts to it and
P_0 T P_0 ~ -P_0 T P_0  =>  spec symmetric about 0  =>  SOFT.

This probe TESTS the mechanism:
  (1) verify Q = P_0 T P_0 (the compression) reproduces M_0 + N bit-exact;
  (2) verify F T F^-1 = -T;
  (3) test whether  soft  <=>  [P_0, F] = 0  (commutant F-invariant)  over all 42 pairs;
  (4) if (3) holds, the converse becomes: odd cycle => H's commutant is NOT F-invariant => the
      compression loses the -symmetry => asymmetric => hard. A clean, possibly provable statement.
"""
from __future__ import annotations
import sys
from itertools import product, combinations_with_replacement
from pathlib import Path

import numpy as np
from scipy.optimize import linear_sum_assignment

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass

SCRIPT_DIR = Path(__file__).parent
sys.path.insert(0, str(SCRIPT_DIR))
import framework as fw
from framework.pauli import _build_kbody_chain, site_op

DIAG = {"I", "Z"}


def is_mixed(t):
    return any(L not in DIAG for L in t)


def build_T(N):
    """Agreement operator T = sum_l Ad_{Z_l} on vec (row-major). Diagonal: N - 2 popcount(i^j)."""
    d = 2 ** N
    diag = np.zeros(d * d)
    for i in range(d):
        for j in range(d):
            diag[i * d + j] = N - 2 * bin(i ^ j).count("1")
    return np.diag(diag).astype(complex)


def build_F(N):
    """F: |i><j| -> |i><j_bar|, j_bar = j ^ (2^N-1). Right-mult by X^{otimes N}. Permutation on vec."""
    d = 2 ** N
    full = d - 1
    F = np.zeros((d * d, d * d))
    for i in range(d):
        for j in range(d):
            F[i * d + (j ^ full), i * d + j] = 1.0
    return F


def commutant_projector(H, N, tol=1e-7):
    """Orthogonal projector P_0 onto the omega=0 sector = span of |E_a><E_b| with E_a=E_b.
    Build in vec space (row-major) from H's eigenvectors."""
    E, V = np.linalg.eigh(H)
    d = len(E)
    # omega=0 modes in the eigenbasis; map to comp-basis vec via (V (x) conj V)
    # vec of |E_a><E_b| in comp basis = (V[:,a]) kron (conj V[:,b]) (row-major: rho_{pq}=V[p,a]conjV[q,b])
    cols = []
    for a in range(d):
        for b in range(d):
            if abs(E[a] - E[b]) < tol:
                vecmode = np.kron(V[:, a], V[:, b].conj())
                cols.append(vecmode)
    Bmat = np.array(cols).T            # d^2 x n
    # orthonormalize (modes are orthonormal already since V unitary, but degenerate combos are too)
    Q_, _ = np.linalg.qr(Bmat)
    P0 = Q_ @ Q_.conj().T
    return P0


def sym_about0(M):
    s = np.linalg.eigvals(M).real
    tgt = -s
    cost = np.abs(s[:, None] - tgt[None, :])
    r, c = linear_sum_assignment(cost)
    return float(cost[r, c].mean())


def main():
    N = 4
    T = build_T(N)
    F = build_F(N)
    # (2) F T F^-1 = -T
    ftf = np.max(np.abs(F @ T @ F.T + T))
    print("=" * 84)
    print("omega=0 block = compression P_0 T P_0 of the agreement operator T")
    print("=" * 84)
    print(f"  (2) global-flip F T F^-1 = -T : residual {ftf:.1e}  (T is -symmetric via F)")
    print()

    # (1) check compression == M_0 + N for the reference pairs
    for label, pair in [("SOFT XXZ+ZXX", [('X','X','Z'),('Z','X','X')]),
                        ("HARD XXZ+XZX", [('X','X','Z'),('X','Z','X')])]:
        H = _build_kbody_chain(N, [tuple(t)+(1.0,) for t in pair])
        P0 = commutant_projector(H, N)
        Qc = P0 @ T @ P0
        # spec of compression restricted to range(P0): use eigenvalues of P0 T P0 (drop the kernel zeros)
        ev = np.sort(np.linalg.eigvals(Qc).real)
        nz = ev[np.abs(ev) > 1e-9]
        asym = sym_about0(Qc)   # includes kernel (symmetric, contributes 0)
        print(f"  {label}: rank(P0)={int(round(np.trace(P0).real))}  compression spec asym about 0 = {asym:.2e}")

    print()
    # (3) the mechanism test: soft <=> [P0, F] = 0  over all 42 pairs
    chain = fw.ChainSystem(N=N)
    terms = [t for t in product("IXYZ", repeat=3)
             if not all(L == "I" for L in t) and fw.klein_index(t) == (0, 1)]
    mixed = [t for t in terms if is_mixed(t)]
    n = agree = 0
    detail = []
    for t1, t2 in combinations_with_replacement(mixed, 2):
        if (sum(c == "Y" for c in t1) % 2) != (sum(c == "Y" for c in t2) % 2):
            continue
        n += 1
        H = _build_kbody_chain(N, [tuple(t1)+(1.0,), tuple(t2)+(1.0,)])
        P0 = commutant_projector(H, N)
        comm_PF = np.max(np.abs(P0 @ F - F @ P0))
        Finv = comm_PF < 1e-7                       # commutant F-invariant
        hard = fw.classify_pauli_pair(chain, [t1, t2], dephase_letter='Z') == 'hard'
        # mechanism predicts: F-invariant => soft
        agree += (Finv == (not hard))
        detail.append((t1, t2, hard, Finv, comm_PF))
    print(f"  (3) MECHANISM: [P_0, F]=0 (commutant F-invariant)  <=>  soft :  {agree}/{n}  "
          f"{'ALL' if agree == n else 'MISMATCH'}")
    if agree != n:
        print("      (mechanism is SUFFICIENT not equivalent; show the breakdown)")
        # count: F-invariant pairs that are soft, non-invariant that are hard
        fi_soft = sum(1 for d in detail if d[3] and not d[2])
        fi_hard = sum(1 for d in detail if d[3] and d[2])
        nfi_soft = sum(1 for d in detail if not d[3] and not d[2])
        nfi_hard = sum(1 for d in detail if not d[3] and d[2])
        print(f"      F-invariant & soft: {fi_soft}   F-invariant & hard: {fi_hard}")
        print(f"      not-F-inv  & soft: {nfi_soft}   not-F-inv  & hard: {nfi_hard}")


if __name__ == "__main__":
    main()
