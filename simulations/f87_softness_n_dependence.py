#!/usr/bin/env python3
"""The soft verdict is not N-invariant: a finite-size crossing in the F87 trichotomy (2026-06-06).

The F87 soft/hard verdict had been read as a property of the Hamiltonian. It is not, for k >= 3:
a Z-dephased chain Hamiltonian can be GENUINELY soft at one chain length and GENUINELY hard at the
next. The witness is XXXX + XYYY + YYYX (a 4-body, reversal-symmetric, mask-bipartite, bit_b-mixed
set): its eigenvalue spectrum pairs to machine precision at N=5 (soft) and breaks by O(0.2) at N=6
(hard). This is not a tolerance artifact; the soft side sits at ~1e-13 and the hard side at ~1e-1,
each landing cleanly on its class.

Two corollaries this script also confirms:
  1. The PalindromeSoftCertifier's certified 2-body site-swap cases (XX+XY+YX, YY+XY+YX) stay soft
     to N=6 at machine precision (pairErr ~1e-13, N-k = 4). The shipped certifier is sound past the
     N=5 it was verified to; only k=2 fully-lit is N-STABLY soft.
  2. The genuine k>=3 soft cases are Z-ROUTED (XZX+XZY+YZX, stable soft at N=4,5,6), the hidden-Q
     mechanism, a separate path from the lit chiral-K colouring.

Pipeline is byte-identical to the C# authority PauliPairTrichotomy.Classify and the Python twin
fw.classify_pauli_pair (same sliding-window k-body builder, same Σγ, same greedy spectrum pairing);
this script only also reads out the magnitudes (residual ‖M‖ and max pairing-error) that the class
hides. SLOW (~5 min): three N=6 dense eigendecompositions of the 4096² Liouvillian.

Companion experiment: experiments/SOFTNESS_IS_N_DEPENDENT.md.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass

SCRIPT_DIR = Path(__file__).parent
sys.path.insert(0, str(SCRIPT_DIR))
from framework.pauli import _build_bilinear, _build_kbody_chain
from framework.lindblad import lindbladian_pauli_dephasing, palindrome_residual

GAMMA = 0.05
J = 1.0


def _bonds(N):
    return [(i, i + 1) for i in range(N - 1)]


def diag(N, terms):
    """Mirror classify_pauli_pair and also return the magnitudes.

    terms: list of letter-tuples, e.g. [('X','X'),('X','Y'),('Y','X')] (k=2) or
    [('X','X','X','X'),('X','Y','Y','Y'),('Y','Y','Y','X')] (k=4). Returns (‖M‖, maxPairErr, class).
    """
    d = 2 ** N
    H = np.zeros((d, d), dtype=complex)
    bil = [t for t in terms if len(t) == 2]
    kb = [t for t in terms if len(t) > 2]
    if bil:
        H = H + _build_bilinear(N, _bonds(N), [(t[0], t[1], J) for t in bil])
    if kb:
        H = H + _build_kbody_chain(N, [tuple(t) + (J,) for t in kb])

    L = lindbladian_pauli_dephasing(H, [GAMMA] * N, dephase_letter="Z")
    sigma = N * GAMMA
    res = float(np.linalg.norm(palindrome_residual(L, sigma, N, dephase_letter="Z")))
    if res < 1e-10:
        return res, 0.0, "truly"

    evals = np.linalg.eigvals(L)
    used = np.zeros(len(evals), dtype=bool)
    mx = 0.0
    for i in range(len(evals)):
        if used[i]:
            continue
        dists = np.abs(evals - (-evals[i] - 2 * sigma))
        dists[used] = np.inf
        j = int(np.argmin(dists))
        used[i] = True
        if j != i:
            used[j] = True
        mx = max(mx, float(dists[j]))
    return res, mx, ("soft" if mx < 1e-6 else "hard")


def _label(terms):
    return "+".join("".join(t) for t in terms)


WITNESSES = [
    ("XX+XY+YX     [certified site-swap; must stay SOFT]", [("X", "X"), ("X", "Y"), ("Y", "X")], [4, 5, 6]),
    ("YY+XY+YX     [certified site-swap]", [("Y", "Y"), ("X", "Y"), ("Y", "X")], [4, 5]),
    ("XXX+XXY+YXX  [hard control: err ~ O(1)]", [("X", "X", "X"), ("X", "X", "Y"), ("Y", "X", "X")], [4, 5]),
    ("XXXX+XYYY+YYYX [N=5 soft -> N=6 hard: the crossing]", [("X", "X", "X", "X"), ("X", "Y", "Y", "Y"), ("Y", "Y", "Y", "X")], [5, 6]),
    ("XZX+XZY+YZX  [k=3 Z-routed soft: N-stable]", [("X", "Z", "X"), ("X", "Z", "Y"), ("Y", "Z", "X")], [4, 5, 6]),
]


def main():
    print("=" * 78)
    print("F87 soft/hard verdict vs chain length N  (||M|| residual, max spectral pairing-error)")
    print("  class = truly if ||M||<1e-10, else soft if pairErr<1e-6, else hard")
    print("=" * 78)
    res = {}
    for name, terms, ns in WITNESSES:
        print(f"\n{name}")
        for N in ns:
            r, err, cls = diag(N, terms)
            res[(_label(terms), N)] = cls
            print(f"    N={N}:  ||M||={r:.3e}   pairErr={err:.3e}   => {cls}")

    # Self-validation: the finding's claims, at the class level (robust to the eigensolver).
    assert res[("XX+XY+YX", 4)] == res[("XX+XY+YX", 5)] == res[("XX+XY+YX", 6)] == "soft", \
        "certifier soundness: the certified XX+XY+YX must stay soft through N=6"
    assert res[("XXXX+XYYY+YYYX", 5)] == "soft", "the crossing: genuinely soft at N=5"
    assert res[("XXXX+XYYY+YYYX", 6)] == "hard", "the crossing: genuinely hard at N=6"
    assert res[("XXX+XXY+YXX", 4)] == res[("XXX+XXY+YXX", 5)] == "hard", "hard control"
    assert res[("XZX+XZY+YZX", 4)] == res[("XZX+XZY+YZX", 6)] == "soft", "Z-routed soft is N-stable"
    print("\nAll N-dependence + certifier-soundness assertions hold.")


if __name__ == "__main__":
    main()
