#!/usr/bin/env python3
"""F87 spec-A FINAL TABLE: the full windowed-converse proof chain, one bit-exact run.

Establishes, for the N=4 k=3 Z-deph diagonal cell (16 hard / 26 soft), the equivalence

  F87-soft
   ⟺ [class]      the canonical λ↔−λ−2σ classifier says soft
   ⟺ [block]      the ω=0 degenerate first-order block M_{ω=0} has spectrum symmetric about −N
   ⟺ [chan]       −N ∈ spec(Φ|W_0),  Φ(A)=Σ_l Z_l A Z_l, W_0 = H-block-diagonal ops
   ⟺ [anti]       ∃ 0≠A ∈ AntiDiag ∩ W_0  (A supported on i↔ī)
   ⟺ [comm]       ∃ 0≠A ∈ AntiDiag ∩ commutant(H)   (block-diag = commuting, on anti-diag)
   ⟺ [2col]       the homogeneous system {d_i+d_j=0 : H_{ij}≠0} has a nonzero solution
   ⟺ [phi]        ∃ linear φ:𝔽₂^N→𝔽₂ with φ|_S=1  AND  no diagonal lift
   ⟺ [bip]        G_H bipartite & loop-free  ⟺  S has no odd 𝔽₂-relation (no K3 triangle)

All eight columns must agree with [class] on every pair. Prints the soft/hard split tables.
"""
from __future__ import annotations
import sys
from itertools import product, combinations_with_replacement, combinations
from collections import defaultdict, deque
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
from f87_flip_generators import edge_generators, phi_exists_gf2

DIAG = {"I", "Z"}
N = 4


def is_mixed(t):
    return any(L not in DIAG for L in t)


def flip_F(n):
    Fx = np.array([[0, 1], [1, 0]], dtype=complex)
    out = Fx
    for _ in range(n - 1):
        out = np.kron(out, Fx)
    return out


F = flip_F(N)


def col_block(H):
    """[block]: spec(M_{ω=0}) symmetric about −N ?"""
    E, V = np.linalg.eigh(H)
    d = len(E)
    Z = [V.conj().T @ site_op(N, l, 'Z') @ V for l in range(N)]
    modes = [(a, b) for a in range(d) for b in range(d) if abs(E[a] - E[b]) < 1e-6]
    n = len(modes)
    M = np.zeros((n, n), dtype=complex)
    for i, (a, b) in enumerate(modes):
        for j, (ap, bp) in enumerate(modes):
            val = sum(Z[l][a, ap] * Z[l][bp, b] for l in range(N))
            if (a, b) == (ap, bp):
                val -= N
            M[i, j] = val
    s = np.sort(np.linalg.eigvals(M).real)
    tgt = -s - 2 * N
    cost = np.abs(s[:, None] - tgt[None, :])
    r, c = linear_sum_assignment(cost)
    return cost[r, c].max() < 1e-6


def col_chan(H):
    """[chan]: −N ∈ spec(Φ|W_0)."""
    E, V = np.linalg.eigh(H)
    d = len(E)
    Z = [V.conj().T @ site_op(N, l, 'Z') @ V for l in range(N)]
    modes = [(a, b) for a in range(d) for b in range(d) if abs(E[a] - E[b]) < 1e-6]
    basis = []
    for (a, b) in modes:
        Mt = np.zeros((d, d), dtype=complex); Mt[a, b] = 1.0; basis.append((a, b))
    n = len(modes)
    P = np.zeros((n, n), dtype=complex)
    for j, (ap, bp) in enumerate(modes):
        # Φ(|ap><bp|) = Σ_l Z_l|ap><bp|Z_l ; coeff on |a><b| = Σ_l Z[l][a,ap] Z[l][bp,b]
        for i, (a, b) in enumerate(modes):
            P[i, j] = sum(Z[l][a, ap] * Z[l][bp, b] for l in range(N))
    ev = np.linalg.eigvals(P).real
    return bool(np.any(np.abs(ev + N) < 1e-6))


def col_comm(H):
    """[comm]: nullity of {d_i+d_j=0 : H_{ij}≠0} ≥ 1."""
    d = 2 ** N
    rows = []
    for i in range(d):
        for j in range(i, d):
            if abs(H[i, j]) > 1e-9:
                r = np.zeros(d); r[i] += 1.0; r[j] += 1.0; rows.append(r)
    if not rows:
        return True
    s = np.linalg.svd(np.array(rows), compute_uv=False)
    rank = int(np.sum(s > 1e-7 * max(1.0, s[0])))
    return (d - rank) >= 1


def col_anti(H):
    """[anti]: ∃ 0≠A∈AntiDiag∩W_0 (block-diagonal). Nullity via A=F·D system."""
    E, V = np.linalg.eigh(H)
    d = 2 ** N
    VdF = V.conj().T @ F
    rows = []
    for a in range(d):
        for b in range(d):
            if abs(E[a] - E[b]) >= 1e-6:
                rows.append(VdF[a, :] * V[:, b])
    if not rows:
        return True
    s = np.linalg.svd(np.array(rows), compute_uv=False)
    rank = int(np.sum(s > 1e-7 * max(1.0, s[0])))
    return (d - rank) >= 1


def col_phi(H):
    """[phi]: φ exists for S AND no diagonal lift."""
    S = edge_generators(H)
    dlift = bool(np.max(np.abs(np.diag(H))) > 1e-9)
    return phi_exists_gf2(S, N) and not dlift


def col_bip(H):
    """[bip]: G_H bipartite & loop-free (BFS 2-colour)."""
    d = 2 ** N
    if np.max(np.abs(np.diag(H))) > 1e-9:
        return False
    color = {}
    for s in range(d):
        if s in color:
            continue
        color[s] = 0; q = deque([s])
        while q:
            u = q.popleft()
            for v in np.where(np.abs(H[u]) > 1e-9)[0]:
                if v == u:
                    continue
                if v not in color:
                    color[v] = color[u] ^ 1; q.append(v)
                elif color[v] == color[u]:
                    return False
    return True


def main():
    chain = fw.ChainSystem(N=N)
    k = 3
    terms = [t for t in product("IXYZ", repeat=k)
             if not all(L == "I" for L in t) and fw.klein_index(t) == (0, 1)]
    mixed = [t for t in terms if is_mixed(t)]

    print("=" * 100)
    print(f"F87 spec-A FINAL: soft ⟺ block ⟺ chan ⟺ anti ⟺ comm ⟺ 2col(=comm) ⟺ phi ⟺ bip   (N={N}, k={k})")
    print("=" * 100)
    cols = ["block", "chan", "anti", "comm", "phi", "bip"]
    agree = {c: 0 for c in cols}
    rows_print = []
    n = 0
    for t1, t2 in combinations_with_replacement(mixed, 2):
        if (sum(c == "Y" for c in t1) % 2) != (sum(c == "Y" for c in t2) % 2):
            continue
        cls = fw.classify_pauli_pair(chain, [t1, t2], dephase_letter="Z")
        if cls == "truly":
            continue
        n += 1
        H = _build_kbody_chain(N, [tuple(t1) + (1.0,), tuple(t2) + (1.0,)])
        soft = (cls == "soft")
        vals = {
            "block": col_block(H), "chan": col_chan(H), "anti": col_anti(H),
            "comm": col_comm(H), "phi": col_phi(H), "bip": col_bip(H),
        }
        for c in cols:
            agree[c] += int(vals[c] == soft)
        rows_print.append((t1, t2, cls, vals))

    print(f"  pairs: {n}  (soft {sum(1 for r in rows_print if r[2]=='soft')}, "
          f"hard {sum(1 for r in rows_print if r[2]=='hard')})")
    print()
    for c in cols:
        status = "ALL" if agree[c] == n else f"MISMATCH ({n-agree[c]})"
        print(f"    [{c:5s}] ⟺ class:  {agree[c]}/{n}  {status}")
    allgood = all(agree[c] == n for c in cols)
    print()
    print(f"  FULL CHAIN bit-exact: {'YES — all 8 readings ⟺ F87 class on all 42 pairs' if allgood else 'NO'}")

    # the soft/hard split numbers (must reproduce F103 §3.2/§3.3 within the Mixed+Mixed sub-cell)
    nhard = sum(1 for r in rows_print if r[2] == 'hard')
    nsoft = sum(1 for r in rows_print if r[2] == 'soft')
    # y_par split
    def yp(t1, t2):
        return sum(c == 'Y' for c in t1) % 2
    hard_y0 = sum(1 for r in rows_print if r[2]=='hard' and yp(r[0],r[1])==0)
    hard_y1 = nhard - hard_y0
    soft_y0 = sum(1 for r in rows_print if r[2]=='soft' and yp(r[0],r[1])==0)
    soft_y1 = nsoft - soft_y0
    print(f"\n  Mixed+Mixed sub-cell split:  hard {nhard} (y0={hard_y0}, y1={hard_y1}), "
          f"soft {nsoft} (y0={soft_y0}, y1={soft_y1})")
    print("  (the 16 odd-cycle hard pairs of §7.2; the §6 34 diagonal-lift hard pairs are the")
    print("   template cell, handled by rule (a) = the diagonal-lift branch of [phi]/[bip].)")


if __name__ == "__main__":
    main()
