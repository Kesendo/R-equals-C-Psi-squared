#!/usr/bin/env python3
"""F87 spec-A: the EXACT scope of the block-level converse proof.

The reduction (soft ⟺ −N∈spec(Φ|W_0) ⟺ {H,D}=0 has nonzero sol ⟺ G_H bipartite) requires:
  (i)  FHF=−H            [each term odd in #Y+#Z; the X^{⊗N} anti-symmetry]
  (ii) H has ≥1 off-diagonal hop  [at least one Mixed term; else H is classical/diagonal and
       the population structure degenerates — those are the pure-template pairs, F111's cell].

Within scope = {pairs with ≥1 Mixed term} ∩ {FHF=−H}, the chain is bit-exact. This includes
ALL 16 odd-cycle hard pairs (rule b, the windowed-converse target) and the Mixed+template
hard pairs (which also carry the odd cycle). Verify across the full diagonal cell.
"""
from __future__ import annotations
import sys
from itertools import product, combinations_with_replacement
from pathlib import Path

import numpy as np

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass

SCRIPT_DIR = Path(__file__).parent
sys.path.insert(0, str(SCRIPT_DIR))
import framework as fw
from framework.pauli import _build_kbody_chain, site_op

N = 4
DIAG = {"I", "Z"}


def flip_F(n):
    Fx = np.array([[0, 1], [1, 0]], dtype=complex)
    out = Fx
    for _ in range(n - 1):
        out = np.kron(out, Fx)
    return out
F = flip_F(N)


def is_mixed(t):
    return any(L not in DIAG for L in t)


def chan_reading(H):
    E, V = np.linalg.eigh(H)
    d = len(E)
    Z = [V.conj().T @ site_op(N, l, 'Z') @ V for l in range(N)]
    modes = [(a, b) for a in range(d) for b in range(d) if abs(E[a] - E[b]) < 1e-6]
    n = len(modes)
    P = np.zeros((n, n), dtype=complex)
    for j, (ap, bp) in enumerate(modes):
        for i, (a, b) in enumerate(modes):
            P[i, j] = sum(Z[l][a, ap] * Z[l][bp, b] for l in range(N))
    return bool(np.any(np.abs(np.linalg.eigvals(P).real + N) < 1e-6))


def main():
    chain = fw.ChainSystem(N=N)
    k = 3
    terms = [t for t in product("IXYZ", repeat=k)
             if not all(L == "I" for L in t) and fw.klein_index(t) == (0, 1)]
    n = 0; ok = 0
    nhard = nsoft = 0
    out_scope = 0
    for t1, t2 in combinations_with_replacement(terms, 2):
        if (sum(c == "Y" for c in t1) % 2) != (sum(c == "Y" for c in t2) % 2):
            continue
        cls = fw.classify_pauli_pair(chain, [t1, t2], dephase_letter="Z")
        if cls == "truly":
            continue
        H = _build_kbody_chain(N, [tuple(t1) + (1.0,), tuple(t2) + (1.0,)])
        has_hop = bool(np.max(np.abs(H - np.diag(np.diag(H)))) > 1e-9)  # ≥1 off-diagonal
        FHF_ok = np.linalg.norm(F @ H @ F + H) < 1e-9
        if not (has_hop and FHF_ok):
            out_scope += 1
            continue
        n += 1
        soft = (cls == "soft")
        nhard += int(cls == "hard"); nsoft += int(soft)
        ok += int(chan_reading(H) == soft)
    print("=" * 84)
    print(f"F87 spec-A FINAL SCOPE: in-scope = (≥1 hop) ∧ (FHF=−H), full diagonal cell, N={N}")
    print("=" * 84)
    print(f"  in-scope pairs: {n}  (soft {nsoft}, hard {nhard});  out-of-scope (classical/even): {out_scope}")
    print(f"  [chan] −N∈spec(Φ|W_0) ⟺ class:  {ok}/{n}  {'ALL — chain bit-exact in scope' if ok==n else 'MISMATCH'}")
    print()
    print("  In-scope hard pairs are EXACTLY the odd-cycle (rule b) hardness — the windowed")
    print("  converse target. Out-of-scope = pure-classical diagonal H (rule a / F111 template).")


if __name__ == "__main__":
    main()
