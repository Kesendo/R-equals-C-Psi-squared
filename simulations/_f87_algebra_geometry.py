#!/usr/bin/env python3
"""The F87 obstruction as graph homology: algebra (cycles) vs geometry (components) (2026-06-05).

Tom's reading: the boundary we mapped (what reaches the spectrum vs what stays combinatorial) is the
algebra/geometry split of the hopping graph. The windowed masks {x^w p_i} are generators of a subgroup
G = <masks> of (Z2)^N; they define the Cayley/hopping graph on the 2^N basis states. Its homology
splits by rank-nullity of the mask set:

  GEOMETRY (H0):  rank r  ->  #components = 2^(N-r)  =  dim(ker L)  =  conserved sectors  (reaches physics)
  ALGEBRA  (H1):  nullity = #masks - r  =  the relation/cycle space  ->  the obstruction (combinatorial)

Of these, only ONE crosses into the spectrum (corrected 2026-06-05; an earlier version of this header
claimed the geometry/rank does too, which is wrong):
  -N exists        <=> bipartite <=> no odd cycle (H1 PARITY) -> reaches physics (hard/soft).  [clean]
  obstruction SIZE =  the METRIC of H1 (cycle length)         -> does NOT reach physics.
  #components / rank = the GEOMETRY (H0)                       -> does NOT cleanly reach physics, because
     letter cancellations (XX+YY kills |00>-|11>) make the ACTUAL hopping graph sparser than the mask,
     so the real sector count is letter-DEPENDENT (the soft rows below show dim ker L > mask-2^(N-r)).

This script makes the split tangible and shows the geometry leg failing: for soft pairs dim ker L exceeds
the mask-based 2^(N-r), so only the algebra's parity (hard/soft) is a clean mask -> spectrum invariant.
"""
from __future__ import annotations

import sys
from itertools import combinations
from pathlib import Path

import numpy as np

if sys.platform == "win32":
    try: sys.stdout.reconfigure(encoding="utf-8")
    except Exception: pass

SCRIPT_DIR = Path(__file__).parent
sys.path.insert(0, str(SCRIPT_DIR))
import framework as fw
from framework.pauli import _build_kbody_chain
from framework.lindblad import lindbladian_pauli_dephasing

def popcount(x): return bin(x).count("1")

def gf2_rank(vectors):
    basis = []
    for v in vectors:
        for b in basis: v = min(v, v ^ b)
        if v: basis.append(v); basis.sort(reverse=True)
    return len(basis)

def windowed(mask, k, N): return [mask << w for w in range(N - k + 1)]

def min_odd_dependence(masks):
    masks = sorted(set(masks))
    for size in range(3, len(masks) + 1, 2):
        for combo in combinations(masks, size):
            x = 0
            for m in combo: x ^= m
            if x == 0: return size
    return 0

def mask_of(letters):
    m = 0
    for i, L in enumerate(letters):
        if L in ("X", "Y"): m |= 1 << i
    return m

def dim_ker_L(letter_terms, N, gamma=0.05):
    H = _build_kbody_chain(N, [tuple(t) + (1.0,) for t in letter_terms])
    L = lindbladian_pauli_dephasing(H, [gamma] * N, dephase_letter='Z')
    return int(np.sum(np.abs(np.linalg.eigvals(L)) < 1e-8))


def main():
    print("=" * 94)
    print("The F87 obstruction as graph homology: GEOMETRY (components) vs ALGEBRA (cycles)")
    print("=" * 94)

    N = 4
    # representative k=3 pairs: soft (no odd cycle) and hard (odd cycle), with letters realised.
    pairs = [
        ("soft  XXI+YYI", [('X','X','I'), ('Y','Y','I')]),
        ("hard  XYI+YIX", [('X','Y','I'), ('Y','I','X')]),
        ("hard  XIY+IXY", [('X','I','Y'), ('I','X','Y')]),
        ("soft  XYZ+YXZ", [('X','Y','Z'), ('Y','X','Z')]),
    ]
    print(f"\n  {'pair':<16} | GEOMETRY (H0)                         | ALGEBRA (H1)")
    print(f"  {'':16} | r  comps=2^(N-r)  dim ker L            | nullity  odd-cycle?  obstruction")
    print("  " + "-" * 90)
    for label, terms in pairs:
        k = len(terms[0])
        flips = []
        for t in terms: flips += windowed(mask_of(t), k, N)
        nmask = len(set(flips))
        r = gf2_rank(flips)
        comps = 1 << (N - r)
        kerd = dim_ker_L(terms, N)
        obs = min_odd_dependence(set(flips))
        nullity = nmask - r
        geo_ok = "OK" if comps == kerd else "MISMATCH"
        odd = "yes (hard)" if obs else "no (soft)"
        print(f"  {label:<16} | r={r}  comps={comps:<2}        dim ker L={kerd} [{geo_ok}]      "
              f"| {nullity:<7}  {odd:<11}  size={obs if obs else '-'}")

    print("\n  reading (corrected 2026-06-05, see PROOF_F103 §7.10-§7.11):")
    print("   The soft rows show dim ker L != #components (mask-2^(N-r)): the MASK overcounts the actual")
    print("   connectivity, because letter cancellations (XX+YY kills |00>-|11>) make H's real hopping")
    print("   graph sparser and letter-DEPENDENT. So GEOMETRY (the component/sector count) does NOT")
    print("   cleanly reach the spectrum. ALGEBRA's METRIC (the obstruction size) does not either.")
    print("   Only ONE Z2 bit reaches physics: ALGEBRA's PARITY -- whether the cycle space has an ODD")
    print("   cycle = the (1+x)-parity = bipartiteness = hard/soft = the -N mode. The combinatorics")
    print("   between algebra and geometry is rich, but only its homological parity crosses into physics.")


if __name__ == "__main__":
    main()
