#!/usr/bin/env python3
"""F87 spec-A probe 3: the unpaired Perron eigenvalue as the soft/hard discriminator.

Finding from probe 2: the ω=0 block of Q always has +N as an eigenvalue (population sector
row-sum N, Perron). SOFT blocks also contain −N (so the block pairs about 0); HARD blocks
have min-eigenvalue STRICTLY ABOVE −N, leaving +N unpaired.

This probe:
  (1) Tests across ALL 16 hard + 26 soft windowed (N=4,k=3) pairs:
        soft  ⟺  −N ∈ spec(Q_{ω=0})   (min eigenvalue == −N)
        hard  ⟺  min eigenvalue of Q_{ω=0}  >  −N
  (2) Identifies the +N and (if present) −N eigenvectors. The +N eigenvector is the all-ones
      on populations. Is it pure-population? What is the −N eigenvector in the soft case?
  (3) Establishes Q_{ω=0} = G^T G structure (a Gram matrix) so spec(Q) ⊆ [?]. Actually
      Q[(a,b),(a',b')] = Σ_l Z[l]_{a,a'} Z[l]_{b',b}; on ω=0 with Z[l] Hermitian this is a
      specific contraction. Check whether ±N are the global spectral extremes (so the
      discriminator is: does the spectrum REACH its reflection-floor −N).
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

DIAG = {"I", "Z"}


def is_mixed(t):
    return any(L not in DIAG for L in t)


def omega0_Q(pair, N=4, ndig=6):
    H = _build_kbody_chain(N, [tuple(t) + (1.0,) for t in pair])
    E, V = np.linalg.eigh(H)
    d = len(E)
    Z = [V.conj().T @ site_op(N, l, 'Z') @ V for l in range(N)]
    modes0 = [(a, b) for a in range(d) for b in range(d) if abs(E[a] - E[b]) < 10**(-ndig)]
    pop = [(a, b) for (a, b) in modes0 if a == b]
    deg = [(a, b) for (a, b) in modes0 if a != b]
    order = pop + deg
    n = len(order)
    Q = np.zeros((n, n), dtype=complex)
    for i, (a, b) in enumerate(order):
        for j, (ap, bp) in enumerate(order):
            Q[i, j] = sum(Z[l][a, ap] * Z[l][bp, b] for l in range(N))
    return Q, len(pop), n


def main():
    N, k = 4, 3
    chain = fw.ChainSystem(N=N)
    terms = [t for t in product("IXYZ", repeat=k)
             if not all(L == "I" for L in t) and fw.klein_index(t) == (0, 1)]
    mixed = [t for t in terms if is_mixed(t)]

    print("=" * 90)
    print(f"F87 spec-A: −N ∈ spec(Q_{{ω=0}}) ⟺ soft   (N={N}, k={k}, Z-deph diagonal cell)")
    print("=" * 90)
    rows = []
    for t1, t2 in combinations_with_replacement(mixed, 2):
        if (sum(c == "Y" for c in t1) % 2) != (sum(c == "Y" for c in t2) % 2):
            continue
        cls = fw.classify_pauli_pair(chain, [t1, t2], dephase_letter="Z")
        if cls == "truly":
            continue
        Q, npop, n = omega0_Q([t1, t2], N)
        ev = np.sort(np.linalg.eigvals(Q).real)
        emin, emax = ev[0], ev[-1]
        has_plusN = np.any(np.abs(ev - N) < 1e-6)
        has_minusN = np.any(np.abs(ev + N) < 1e-6)
        rows.append((t1, t2, cls, emin, emax, has_plusN, has_minusN))

    # tabulate
    n_soft = sum(1 for r in rows if r[2] == "soft")
    n_hard = sum(1 for r in rows if r[2] == "hard")
    # discriminator: soft ⟺ has_minusN
    agree = sum(1 for r in rows if (r[2] == "soft") == r[6])
    all_plusN = all(r[5] for r in rows)
    print(f"  pairs: {len(rows)}  (soft {n_soft}, hard {n_hard})")
    print(f"  every block has +N eigenvalue: {all_plusN}")
    print(f"  soft ⟺ (−N ∈ spec):  {agree}/{len(rows)}  {'ALL' if agree==len(rows) else 'MISMATCH'}")
    print()
    print(f"  {'t1':>6} {'t2':>6} {'class':>5} {'min(Q)':>8} {'max(Q)':>8} {'+N?':>4} {'−N?':>4}")
    for (t1, t2, cls, emin, emax, hp, hm) in sorted(rows, key=lambda r: (r[2], r[3])):
        print(f"  {''.join(t1):>6} {''.join(t2):>6} {cls:>5} "
              f"{emin:8.4f} {emax:8.4f} {str(hp):>4} {str(hm):>4}")

    # explicit +N / −N eigenvectors for one soft and one hard
    print()
    for label, pair in [("SOFT XXZ+ZXX", [('X','X','Z'),('Z','X','X')]),
                        ("HARD XXZ+XZX", [('X','X','Z'),('X','Z','X')])]:
        Q, npop, n = omega0_Q(pair, N)
        w, U = np.linalg.eig(Q)
        w = w.real
        for tgt, name in [(N, "+N"), (-N, "−N")]:
            idx = np.where(np.abs(w - tgt) < 1e-6)[0]
            if len(idx) == 0:
                print(f"  {label}: {name} ABSENT")
                continue
            vec = U[:, idx[0]]
            pop_weight = np.linalg.norm(vec[:npop])**2
            print(f"  {label}: {name} present, |vec_pop|²={pop_weight:.4f}  "
                  f"(pure-population if ≈1)")


if __name__ == "__main__":
    main()
