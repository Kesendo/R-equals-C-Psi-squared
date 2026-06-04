#!/usr/bin/env python3
"""F87 spec-A probe 7: the two structural facts behind the soft construction.

FACT A (the flip anti-symmetry):  FHF = −H  for EVERY Mixed+Mixed diagonal-cell pair,
  F = X^{⊗N}. Reason: X-conjugation sends Z↦−Z, Y↦−Y, X↦X, I↦I; a diagonal-cell term is
  two off-diagonal {X,Y} letters + one {I,Z}. The number of sign-flipping letters (#Y+#Z)
  per term is the parity that decides FHF=±H. Check it equals −H across all 42 pairs and
  pin the parity.

FACT B (the soft eigenmode is F·K):  given the chiral K (bipartite, KHK=−H), A:=FK obeys
  [A,H]=0 (so A∈W_0) because (FK)H(FK) = F(KHK)F = F(−H)F = −FHF = H. And A is anti-diagonal
  (F anti-diagonal, K diagonal), so Σ_l Z_l A Z_l = −N·A automatically. Hence A is the soft
  −N eigenmode. Verify [FK,H]=0 across all soft pairs (constructive soft direction).

These two facts make "bipartite ⟹ soft" at the block level a two-line algebra:
  K exists ⟹ A=FK ∈ AntiDiag∩W_0 ⟹ −N∈spec(Φ|W_0) ⟹ ω=0 block symmetric ⟹ (per-block) soft.
"""
from __future__ import annotations
import sys
from itertools import product, combinations_with_replacement
from pathlib import Path

import numpy as np
from collections import deque

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


def flip_F(N):
    Fx = np.array([[0, 1], [1, 0]], dtype=complex)
    out = Fx
    for _ in range(N - 1):
        out = np.kron(out, Fx)
    return out


def chiral_K(H, tol=1e-9):
    d = H.shape[0]
    if np.max(np.abs(np.diag(H))) > tol:
        return None
    color = {}
    for s in range(d):
        if s in color:
            continue
        color[s] = 0
        q = deque([s])
        while q:
            u = q.popleft()
            for v in np.where(np.abs(H[u]) > tol)[0]:
                if v == u:
                    continue
                if v not in color:
                    color[v] = color[u] ^ 1
                    q.append(v)
                elif color[v] == color[u]:
                    return None
    return np.diag([(-1.0) ** color.get(i, 0) for i in range(d)]).astype(complex)


def main():
    N, k = 4, 3
    chain = fw.ChainSystem(N=N)
    terms = [t for t in product("IXYZ", repeat=k)
             if not all(L == "I" for L in t) and fw.klein_index(t) == (0, 1)]
    mixed = [t for t in terms if is_mixed(t)]
    F = flip_F(N)

    factA_ok = 0
    factB_ok = 0
    n_pairs = 0
    n_soft = 0
    factA_fail = []
    parity_check = set()
    for t1, t2 in combinations_with_replacement(mixed, 2):
        if (sum(c == "Y" for c in t1) % 2) != (sum(c == "Y" for c in t2) % 2):
            continue
        cls = fw.classify_pauli_pair(chain, [t1, t2], dephase_letter="Z")
        if cls == "truly":
            continue
        n_pairs += 1
        H = _build_kbody_chain(N, [tuple(t1) + (1.0,), tuple(t2) + (1.0,)])
        # FACT A
        if np.linalg.norm(F @ H @ F + H) < 1e-9:
            factA_ok += 1
        else:
            factA_fail.append((t1, t2, np.linalg.norm(F @ H @ F + H), np.linalg.norm(F @ H @ F - H)))
        # parity per term: #Y + #Z (sign-flipping letters under X-conjugation)
        for t in (t1, t2):
            parity_check.add((sum(c == 'Y' for c in t) + sum(c == 'Z' for c in t)) % 2)
        # FACT B
        if cls == "soft":
            n_soft += 1
            K = chiral_K(H)
            if K is not None:
                A = F @ K
                if np.linalg.norm(A @ H - H @ A) < 1e-9:
                    factB_ok += 1

    print("=" * 80)
    print(f"F87 spec-A probe 7  (N={N}, k={k}, {n_pairs} pairs, {n_soft} soft)")
    print("=" * 80)
    print(f"  FACT A  (FHF = −H for every diagonal-cell pair): {factA_ok}/{n_pairs}  "
          f"{'ALL' if factA_ok==n_pairs else 'FAIL: '+str(factA_fail[:3])}")
    print(f"          per-term parity (#Y+#Z) mod 2 over all terms: {parity_check}  "
          f"(odd ⟹ each term sign-flips ⟹ FHF=−H)")
    print(f"  FACT B  ([F·K, H]=0, so soft eigenmode A=F·K ∈ W_0): {factB_ok}/{n_soft}  "
          f"{'ALL soft' if factB_ok==n_soft else 'GAP'}")
    print()
    print("  Derivation chain (block-level bipartite ⟹ soft, now fully algebraic):")
    print("    KHK=−H (chiral, bipartite) and FHF=−H (Fact A) ⟹ (FK)H(FK)=F(KHK)F=F(−H)F=H")
    print("    ⟹ [FK,H]=0 ⟹ A:=FK ∈ W_0; A anti-diagonal ⟹ Σ_l Z_l A Z_l=−N·A")
    print("    ⟹ −N ∈ spec(Φ|W_0) ⟹ ω=0 block Q symmetric about 0 ⟹ soft.")


if __name__ == "__main__":
    main()
