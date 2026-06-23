#!/usr/bin/env python3
"""Does the §7.5/§7.6 windowed-converse derivation hold beyond k=3? Spot-check k=4, N=5 (2026-06-04).

The §7.6 premise closure (degenerate PT: L0=-i[H,.] normal; analyticity: char-poly affine in γ) is
k-agnostic. The §7.5 Perron criterion needs Fact A (FHF=-H), which holds for any k because every
diagonal-cell (Klein (0,1)) term has odd #Y+#Z (= bit_b parity, the cell's defining parity). So the
derivation should carry to windowed k=4 (k<N, i.e. N>=5). This spot-checks it: over a sample of
windowed (N=5, k=4) diagonal-cell Mixed+Mixed pairs, confirm
    Fact A: FHF = -H  (per pair)
    soft  ⟺  bipartite-in-Z-basis  ⟺  spec(L)=spec(-L-2σ) at the physical γ=0.05
the exact chain the k=3 proof rests on, now at k=4.
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
from framework.lindblad import lindbladian_pauli_dephasing
from framework.pauli import _build_kbody_chain, site_op
from f87_flip_generators import is_bipartite

DIAG = {"I", "Z"}


def fmat(N):
    X = np.array([[0, 1], [1, 0]], dtype=complex)
    F = np.array([[1]], dtype=complex)
    for _ in range(N):
        F = np.kron(F, X)
    return F


def break_at(H, gamma, N):
    L = lindbladian_pauli_dephasing(H, [gamma] * N, dephase_letter="Z")
    ev = np.linalg.eigvals(L)
    tgt = -ev - 2 * N * gamma
    cost = np.abs(ev[:, None] - tgt[None, :])
    r, c = linear_sum_assignment(cost)
    return float(cost[r, c].max())


def main():
    N, k = 5, 4
    chain = fw.ChainSystem(N=N)
    F = fmat(N)
    terms = [t for t in product("IXYZ", repeat=k)
             if not all(L == "I" for L in t) and fw.klein_index(t) == (0, 1)
             and any(L not in DIAG for L in t)]

    print("=" * 86)
    print(f"F87 windowed converse beyond k=3: spot-check N={N}, k={k} (k<N, windowed)")
    print("=" * 86)
    print(f"  diagonal-cell (0,1) Mixed+Mixed k=4 terms: {len(terms)}")
    print()

    # Fact A holds for every cell term (parity), confirm once on the term list
    factA_terms = all((sum(c == "Y" for c in t) + sum(c == "Z" for c in t)) % 2 == 1 for t in terms)
    print(f"  Fact A precondition (#Y+#Z odd for every (0,1) term, any k): {factA_terms}")
    print()

    # sample pairs until we have a few soft and a few hard
    print(f"  {'t1':>7} {'t2':>7} {'FHF=-H':>8} {'bipartite':>10} {'class':>6} "
          f"{'break@γ.05':>12} {'soft⟺bip⟺exact':>16}")
    n_soft = n_hard = 0
    all_ok = True
    for t1, t2 in combinations_with_replacement(terms, 2):
        if (sum(c == "Y" for c in t1) % 2) != (sum(c == "Y" for c in t2) % 2):
            continue
        cls = fw.classify_pauli_pair(chain, [t1, t2], dephase_letter="Z")
        want = (cls == "soft")
        if want and n_soft >= 6:
            continue
        if not want and n_hard >= 6:
            continue
        H = _build_kbody_chain(N, [tuple(t1) + (1.0,), tuple(t2) + (1.0,)])
        factA = np.max(np.abs(F @ H @ F + H)) < 1e-9
        bip = is_bipartite(H)
        brk = break_at(H, 0.05, N)
        exact = brk < 1e-6
        ok = (want == bip) and (want == exact) and factA
        all_ok = all_ok and ok
        print(f"  {''.join(t1):>7} {''.join(t2):>7} {str(factA):>8} {str(bip):>10} {cls:>6} "
              f"{brk:>12.4e} {('OK' if ok else 'FAIL'):>16}")
        if want:
            n_soft += 1
        else:
            n_hard += 1
        if n_soft >= 6 and n_hard >= 6:
            break

    print()
    print(f"  sampled {n_soft} soft + {n_hard} hard windowed k=4 pairs:  "
          f"soft ⟺ bipartite ⟺ spec-exact, all hold = {all_ok}")
    print()
    print("  => the §7.5 criterion + the §7.6 (Fact A, bipartite, spec-break) chain carry to k=4")
    print("     windowed (k<N). The derivation is not k=3-specific; only the explicit 'odd cycle =")
    print("     K3 triangle on 3 consecutive sites' characterisation is (k=4 has a larger mask set).")
    assert factA_terms and all_ok


if __name__ == "__main__":
    main()
