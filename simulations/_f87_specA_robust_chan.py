#!/usr/bin/env python3
"""F87 spec-A robustness: the [chan]/[block]/[comm] readings at N=5 (Z) and N=4 (X,Y).

Confirms the new physical heart of the proof (−N∈spec(Φ|W_0), ω=0 block symmetric about −N,
and the {H,D}=0 2-colour nullity) are the discriminator beyond the N=4/Z anchor. For X,Y the
analysis runs in the dephase-letter eigenbasis (dephasing acts as Z there, flip is X^{⊗N}).
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
from framework.pauli import _build_kbody_chain

KLEIN = {"Z": (0, 1), "X": (1, 0), "Y": (1, 1)}
DIAG_BY_LETTER = {"Z": {"I", "Z"}, "X": {"I", "X"}, "Y": {"I", "Y"}}
H1 = np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2)


def y_to_z():
    Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
    w, v = np.linalg.eigh(Y)
    order = np.argsort(-w)
    return v[:, order].conj().T
RY = y_to_z()


def rotate(N, dl):
    if dl == "Z":
        return np.eye(2 ** N, dtype=complex)
    u1 = H1 if dl == "X" else RY
    U = u1
    for _ in range(N - 1):
        U = np.kron(U, u1)
    return U


def Zc(N):
    Z1 = np.array([[1, 0], [0, -1]], dtype=complex)
    out = []
    for l in range(N):
        ops = [np.eye(2, dtype=complex)] * N
        ops[l] = Z1
        M = ops[0]
        for o in ops[1:]:
            M = np.kron(M, o)
        out.append(M)
    return out


def readings(pair, N, dl):
    H = _build_kbody_chain(N, [tuple(t) + (1.0,) for t in pair])
    U = rotate(N, dl)
    Hr = U @ H @ U.conj().T            # dephasing is Z in this basis
    E, V = np.linalg.eigh(Hr)
    d = len(E)
    Zl = Zc(N)
    Zb = [V.conj().T @ Zl[l] @ V for l in range(N)]
    modes = [(a, b) for a in range(d) for b in range(d) if abs(E[a] - E[b]) < 1e-6]
    n = len(modes)
    # Φ|W_0 (gain) and M=Φ−N·I block
    P = np.zeros((n, n), dtype=complex)
    for j, (ap, bp) in enumerate(modes):
        for i, (a, b) in enumerate(modes):
            P[i, j] = sum(Zb[l][a, ap] * Zb[l][bp, b] for l in range(N))
    evP = np.linalg.eigvals(P).real
    chan = bool(np.any(np.abs(evP + N) < 1e-6))
    sM = np.sort(evP - N)
    tgt = -sM - 2 * N
    cost = np.abs(sM[:, None] - tgt[None, :])
    r, c = linear_sum_assignment(cost)
    block = (cost[r, c].max() < 1e-6)
    # {H,D}=0 nullity in computational(rotated) basis
    rows = []
    for i in range(d):
        for j in range(i, d):
            if abs(Hr[i, j]) > 1e-9:
                rr = np.zeros(d); rr[i] += 1.0; rr[j] += 1.0; rows.append(rr)
    if rows:
        sv = np.linalg.svd(np.array(rows), compute_uv=False)
        rank = int(np.sum(sv > 1e-7 * max(1.0, sv[0])))
        comm = (d - rank) >= 1
    else:
        comm = True
    return chan, block, comm


def main():
    print("=" * 86)
    print("F87 spec-A robustness: [chan]/[block]/[comm] ⟺ class beyond N=4/Z")
    print("=" * 86)
    for (N, k, dl) in [(4, 3, "Z"), (5, 3, "Z"), (4, 3, "X"), (4, 3, "Y")]:
        chain = fw.ChainSystem(N=N)
        diag = DIAG_BY_LETTER[dl]
        terms = [t for t in product("IXYZ", repeat=k)
                 if not all(L == "I" for L in t) and fw.klein_index(t) == KLEIN[dl]]
        mixed = [t for t in terms if any(L not in diag for L in t)]
        n = 0; ac = ab = am = 0
        for t1, t2 in combinations_with_replacement(mixed, 2):
            if (sum(c == "Y" for c in t1) % 2) != (sum(c == "Y" for c in t2) % 2):
                continue
            cls = fw.classify_pauli_pair(chain, [t1, t2], dephase_letter=dl)
            if cls == "truly":
                continue
            n += 1
            soft = (cls == "soft")
            chan, block, comm = readings([t1, t2], N, dl)
            ac += int(chan == soft); ab += int(block == soft); am += int(comm == soft)
        print(f"  N={N} {dl}-deph ({n} pairs):  "
              f"[chan] {ac}/{n}  [block] {ab}/{n}  [comm] {am}/{n}  "
              f"{'ALL' if ac==ab==am==n else 'MISMATCH'}")


if __name__ == "__main__":
    main()
