#!/usr/bin/env python3
"""The soft verdict is not N-invariant: a finite-size crossing in the F87 trichotomy (2026-06-06).

The F87 soft/hard verdict had been read as a property of the Hamiltonian. It is not, for k >= 3:
a Z-dephased chain Hamiltonian can be GENUINELY soft at one chain length and GENUINELY hard at the
next. The witness is XXXX + XYYY + YYYX (a 4-body, reversal-symmetric, mask-bipartite, bit_b-mixed
set): its eigenvalue spectrum pairs to machine precision at N=5 (soft) and breaks by O(0.2) at N=6
(hard). This is not a tolerance artifact; the soft side sits at ~1e-13 and the hard side at ~1e-1,
each landing cleanly on its class.

Four corollaries this script also confirms:
  1. The PalindromeSoftCertifier's certified 2-body site-swap cases (XX+XY+YX, YY+XY+YX) stay soft
     to N=6 at machine precision (pairErr ~1e-13, N-k = 4). The shipped certifier is sound past the
     N=5 it was verified to.
  2. XZX+XZY+YZX is stably soft at N=4,5,6 for a proven reason: the period-4 golden router under
     the window-summed condition. No per-term router in the bounded candidate set reaches it, which
     makes it the HARDER of the stable cases, not the canonical one.
  3. Stability is NOT the Z's private property. The fully-lit (no-Z) k=3 sets XXY+XYX+YXX and
     XXY+YXX+YYY are soft at N=4,5,6 too, and not by grace: both are palindromized by the UNIFORM
     per-site map P4, one of the four representatives the k-body routing certifier tries
     (compute/RCPsiSquared.Diagnostics/F87/KBodyPalindromeRouting.cs), whose certificate is additive
     over windows and therefore N-independent. uniform_router_check() below builds that Pi and
     conjugates L directly, with the hard control as the discriminator. The census below
     enumerates all twelve three-term reversal-symmetric fully-lit k=3 sets: 8 hard, 2 truly,
     2 soft, the same at N=4 and N=5. So a fully-lit k>=3 soft verdict CAN be finite-size
     grace (XXXX+XYYY+YYYX is), but it is not grace by construction, and within these twelve
     what sorts the verdicts is whether a router exists.
  4. F158 reads the crossing with no spectrum: end_count_check() takes the near and far ends of
     XXXX+XYYY+YYYX as exact rational ranks, 8 and 8 at N=4,5, 8 and 4 at N=6, 8 and 0 from
     N=7, with the two lit soft sets (equal ends at every N) and the hard control XXX+XXY+YXX
     (equal ends only at N=3) as controls; Theorem 2 of the complement connection on the
     basis-state graph meets them, and square_law() checks the hand-derived holonomy of every
     square of two window flips (offset 2 closes iff the six bits under it have even parity).

Pipeline matches the C# authority PauliPairTrichotomy.Classify (same sliding-window k-body builder,
same Σγ, same greedy spectrum pairing, inlined below); this script only also reads out the
magnitudes (residual ‖M‖ and max pairing-error) that the class hides. NOTE the Python twin
fw.classify_pauli_pair no longer matches: it computes the EXACT bottleneck pairing rather than the
greedy one, so it reports smaller magnitudes for the same system (0.3 against 0.4 for IIZ+IZI at
N=4). Greedy can only overstate the error, so the soft/hard verdicts here are unaffected: against
the 1e-6 threshold the soft side sits seven orders below and the hard side five orders above, and
a factor-of-two difference in the reported magnitude cannot cross that.
SLOW: five N=6 dense eigendecompositions of the 4096² Liouvillian, plus the twelve-set census
at N=4 and N=5.

Companion experiment: experiments/SOFTNESS_IS_N_DEPENDENT.md.
"""
from __future__ import annotations

import sys
from fractions import Fraction
from itertools import combinations, product
from pathlib import Path

import numpy as np

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass

SCRIPT_DIR = Path(__file__).parent
sys.path.insert(0, str(SCRIPT_DIR))
from framework.pauli import (
    PAULI_LABELS,
    _build_bilinear,
    _build_kbody_chain,
    _k_to_indices,
    _vec_to_pauli_basis_transform,
)
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
    ("XZX+XZY+YZX  [k=3 golden-routed soft, window-summed only]", [("X", "Z", "X"), ("X", "Z", "Y"), ("Y", "Z", "X")], [4, 5, 6]),
    ("XXY+XYX+YXX  [k=3 fully lit, N-stable soft: no Z anywhere]", [("X", "X", "Y"), ("X", "Y", "X"), ("Y", "X", "X")], [4, 5, 6]),
    ("XXY+YXX+YYY  [k=3 fully lit, the second one]", [("X", "X", "Y"), ("Y", "X", "X"), ("Y", "Y", "Y")], [4, 5, 6]),
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
    assert res[("XZX+XZY+YZX", 4)] == res[("XZX+XZY+YZX", 6)] == "soft", "the golden-routed Z-middle set is soft at N=4 and N=6"
    # The lit family is not uniformly finite-size grace: a fully-lit k=3 set can
    # be soft at every length tested. Its N-independence is not read off these
    # three lengths but off the router (uniform P4, see the module docstring);
    # what the three lengths do is refute "only k=2 fully-lit is N-stably soft".
    assert res[("XXY+XYX+YXX", 4)] == res[("XXY+XYX+YXX", 5)] == res[("XXY+XYX+YXX", 6)] == "soft", \
        "fully-lit k=3 N-stable soft: the counterexample to 'only k=2 is N-stably soft'"
    assert (res[("XXY+YXX+YYY", 4)] == res[("XXY+YXX+YYY", 5)]
            == res[("XXY+YXX+YYY", 6)] == "soft"), "the second fully-lit k=3 soft set"
    print("\nAll N-dependence + certifier-soundness assertions hold.")
    lit_family_census()


def lit_family_census():
    """The twelve three-term reversal-symmetric fully-lit k=3 sets, classified.

    "Fully lit" = every letter in {X, Y}, no dephasing letter anywhere. The family
    is the 3-element subsets of the 8 length-3 XY words that are closed under
    reversal; there are exactly 12. The doc quotes the 8 hard / 2 truly / 2 soft
    split, so it is enumerated here rather than asserted from two examples.
    """
    words = ["".join(w) for w in product("XY", repeat=3)]
    sets = [s for s in combinations(words, 3) if all(w[::-1] in s for w in s)]
    print("\n" + "=" * 78)
    print(f"Census: the {len(sets)} three-term reversal-symmetric fully-lit k=3 sets")
    print("=" * 78)
    tallies = {}
    for N in (4, 5):
        tally = {"truly": 0, "soft": 0, "hard": 0}
        for s in sets:
            _, _, cls = diag(N, [tuple(w) for w in s])
            tally[cls] += 1
            if cls != "hard":
                print(f"    N={N}  {'+'.join(s):16s} {cls}")
        tallies[N] = tally
        print(f"  N={N} tally: {tally}")
    assert len(sets) == 12, f"the family has 12 members, found {len(sets)}"
    for N in (4, 5):
        assert tallies[N] == {"truly": 2, "soft": 2, "hard": 8}, \
            f"N={N} census split moved: {tallies[N]}"
    print("\nCensus split 8 hard / 2 truly / 2 soft holds at N=4 and N=5.")
    uniform_router_check()


# Three of the four per-site representatives the k-body routing certifier tries
# (KBodyPalindromeRouting.cs: P1, P4, M2, M), as maps on single-site Pauli
# letters: letter -> (phase, image letter). The fourth, the dense golden M, is
# not a signed permutation and is left to the C# certifier.
_P1 = {"I": (1, "X"), "X": (1, "I"), "Y": (1j, "Z"), "Z": (1j, "Y")}
_P4 = {"I": (1, "Y"), "X": (1j, "Z"), "Y": (1, "I"), "Z": (1j, "X")}
_M2 = {"I": (1, "Y"), "X": (-1j, "Z"), "Y": (1, "I"), "Z": (-1j, "X")}
_REPS = {"P1": _P1, "P4": _P4, "M2": _M2}


def _pi_pattern(N, pattern):
    """Pi = (x)_l rep[pattern[l mod P]] as a 4^N x 4^N matrix on Pauli coefficients.

    A period-P per-site pattern, the shape KBodyPalindromeRouting's CandidateSet
    enumerates. pattern is a tuple of representative names.
    """
    d2 = 4 ** N
    labels = ["".join(PAULI_LABELS[t] for t in _k_to_indices(k, N)) for k in range(d2)]
    index = {s: i for i, s in enumerate(labels)}
    Pi = np.zeros((d2, d2), dtype=complex)
    for k, s in enumerate(labels):
        phase = 1 + 0j
        out = []
        for site, ch in enumerate(s):
            c, o = _REPS[pattern[site % len(pattern)]][ch]
            phase *= c
            out.append(o)
        Pi[index["".join(out)], k] = phase
    return Pi


def _pauli_basis_liouvillian(N, terms):
    """L in the 4^N Pauli-string basis, plus its norm."""
    H = _build_kbody_chain(N, [tuple(t) + (J,) for t in terms])
    L_vec = lindbladian_pauli_dephasing(H, [GAMMA] * N, dephase_letter="Z")
    T = _vec_to_pauli_basis_transform(N, order="C")
    return (T.conj().T @ L_vec @ T) / (2 ** N)


def _router_residual(N, terms, pattern):
    """||Pi L Pi^-1 + L + 2*Sigma_gamma|| for the given pattern, with ||L|| for scale."""
    L = _pauli_basis_liouvillian(N, terms)
    Pi = _pi_pattern(N, pattern)
    R = Pi @ L @ np.linalg.inv(Pi) + L + 2 * (N * GAMMA) * np.eye(4 ** N)
    return float(np.linalg.norm(R)), float(np.linalg.norm(L))


def _first_router(N, terms):
    """The first period-<=2 signed-permutation pattern that palindromizes L.

    The bounded search KBodyPalindromeRouting performs, done end-to-end on the
    full Liouvillian instead of on the 4^k window. Returns its name or None;
    None means "none in this candidate set", never "no router exists".
    """
    L = _pauli_basis_liouvillian(N, terms)
    scale = float(np.linalg.norm(L))
    names = list(_REPS)
    patterns = [(r,) for r in names]
    patterns += [(a, b) for a in names for b in names]
    Id = np.eye(4 ** N)
    for pattern in patterns:
        Pi = _pi_pattern(N, pattern)
        R = Pi @ L @ np.linalg.inv(Pi) + L + 2 * (N * GAMMA) * Id
        if float(np.linalg.norm(R)) < 1e-9 * max(1.0, scale):
            return "+".join(pattern)
    return None


def uniform_router_check():
    """The two lit soft sets are palindromized by the UNIFORM map P4, exactly.

    This is what makes their softness N-independent rather than measured: the
    certificate is a per-site product and additive over windows, so it does not
    care how long the chain is. The hard control goes through the same
    construction, because a check that cannot fail is no check; there it misses
    by O(||L||).
    """
    print("\n" + "=" * 78)
    print("Uniform per-site router: ||Pi L Pi^-1 + L + 2*sigma|| for Pi = P4^(x)N")
    print("=" * 78)
    soft = [
        ("XXY+XYX+YXX", [("X", "X", "Y"), ("X", "Y", "X"), ("Y", "X", "X")]),
        ("XXY+YXX+YYY", [("X", "X", "Y"), ("Y", "X", "X"), ("Y", "Y", "Y")]),
    ]
    control = ("XXX+XXY+YXX", [("X", "X", "X"), ("X", "X", "Y"), ("Y", "X", "X")])
    for N in (4, 5):
        for name, terms in soft:
            r, nl = _router_residual(N, terms, ("P4",))
            print(f"    N={N}  {name:14s} P4  residual={r:.3e}   ||L||={nl:.3f}")
            assert r < 1e-9 * nl, f"P4 must palindromize {name} at N={N}, got {r}"
        r, nl = _router_residual(N, control[1], ("P4",))
        print(f"    N={N}  {control[0]:14s} P4  residual={r:.3e}   ||L||={nl:.3f}  (control)")
        assert r > 0.1 * nl, f"the control must NOT be routed by P4 at N={N}, got {r}"
    print("\nP4 routes both lit soft sets and not the hard control.")
    router_census()


def router_census():
    """Which of the twelve lit sets admit a period-<=2 signed-permutation router.

    The doc reads the family by routability, so the search is run here rather
    than asserted: a router exists for exactly the four non-hard sets, P4 for
    the two soft and P1 for the two truly, and for none of the eight hard ones.
    The Z-middle templates are included because the set XZX+XZY+YZX admits no
    such router while two of its three templates do on their own, which is why
    its certificate cannot be assembled template by template.
    """
    words = ["".join(w) for w in product("XY", repeat=3)]
    sets = [s for s in combinations(words, 3) if all(w[::-1] in s for w in s)]
    print("\n" + "=" * 78)
    print("Router census at N=4: first period-<=2 pattern over {P1, P4, M2}, or none")
    print("=" * 78)
    routed = {}
    for s in sets:
        terms = [tuple(w) for w in s]
        _, _, cls = diag(4, terms)
        router = _first_router(4, terms)
        routed[cls] = routed.get(cls, 0) + (1 if router else 0)
        print(f"    {'+'.join(s):16s} {cls:6s} router={router or 'none'}")
    assert routed.get("hard", 0) == 0, "no hard set may admit a router"
    assert routed.get("soft", 0) == 2 and routed.get("truly", 0) == 2, \
        f"every non-hard set must admit one: {routed}"
    print("\n  Per-template, the Z-middle set (which admits no router as a set):")
    for t in ("XZX", "XZY", "YZX"):
        router = _first_router(4, [tuple(t)])
        print(f"    {t:16s}        router={router or 'none'}")
    assert _first_router(4, [("X", "Z", "X"), ("X", "Z", "Y"), ("Y", "Z", "X")]) is None
    assert _first_router(4, [tuple("XZX")]) is None
    assert _first_router(4, [tuple("XZY")]) is not None
    print("\nRoutability splits the twelve 4 against 8, and XZX alone is the obstruction.")
    end_count_check()


# The two ends of F158, exactly. With Z jumps on every site the near end is the
# commutant of H among the {I, Z} strings and the far end among the lit {X, Y}
# strings (F158 Lemma 1); the palindrome holds iff the two have equal dimension.
# A Pauli string is an (x, z) bit pair; [T, P] = 2·T·P = ±2i·(string) when they
# anticommute, so the commutator is an integer matrix on strings and its rank is
# taken exactly over the rationals (a rank modulo a prime would only bound it).


def _string(word, offset=0):
    x = z = 0
    for i, c in enumerate(word):
        if c in "XY":
            x |= 1 << (i + offset)
        if c in "YZ":
            z |= 1 << (i + offset)
    return x, z


def _anticommute(a, b):
    return (bin(a[0] & b[1]).count("1") + bin(a[1] & b[0]).count("1")) % 2 == 1


def _commutator_sign(a, b):
    """T·P = (±i)·(a xor b) for anticommuting strings; returns the sign of the i."""
    order = {(1, 0): 0, (1, 1): 1, (0, 1): 2}   # X, Y, Z; X→Y→Z cyclic gives +i
    exponent = 0
    for s in range(max(a[0] | a[1] | b[0] | b[1], 1).bit_length()):
        p = ((a[0] >> s) & 1, (a[1] >> s) & 1)
        q = ((b[0] >> s) & 1, (b[1] >> s) & 1)
        if p != (0, 0) and q != (0, 0) and p != q:
            exponent += 1 if (order[q] - order[p]) % 3 == 1 else -1
    return 1 if exponent % 4 == 1 else -1


def _rank_q(columns):
    """Exact rank over the rationals of sparse integer columns (Gaussian elimination in Fraction)."""
    pivots = {}
    for column in columns:
        row = {k: Fraction(v) for k, v in column.items() if v}
        while row:
            k = min(row)
            if k not in pivots:
                pivots[k] = row
                break
            pivot = pivots[k]
            f = row[k] / pivot[k]
            for kk, vv in pivot.items():
                row[kk] = row.get(kk, 0) - f * vv
                if row[kk] == 0:
                    del row[kk]
    return len(pivots)


def end_counts(N, templates):
    """(near, far, colourings) for the sliding-window chain of the templates, Z jumps on every site."""
    k = len(templates[0])
    H = [_string(t, w) for w in range(N - k + 1) for t in templates]

    def kernel(strings):
        columns = []
        for s in strings:
            column = {}
            for T in H:
                if _anticommute(T, s):
                    out = (T[0] ^ s[0], T[1] ^ s[1])
                    column[out] = column.get(out, 0) + _commutator_sign(T, s)
            columns.append({key: v for key, v in column.items() if v})
        return len(strings) - _rank_q(columns)

    lit = (1 << N) - 1
    near = kernel([(0, z) for z in range(1 << N)])
    far_strings = [(lit, z) for z in range(1 << N)]
    colourings = sum(1 for s in far_strings if not any(_anticommute(T, s) for T in H))
    return near, kernel(far_strings), colourings


def _gauss_div(a, b):
    """a / b for Gaussian rationals held as (Fraction, Fraction)."""
    n = b[0] * b[0] + b[1] * b[1]
    return ((a[0] * b[0] + a[1] * b[1]) / n, (a[1] * b[0] - a[0] * b[1]) / n)


def _gauss_mul(a, b):
    return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])


def _basis_entries(N, templates):
    """H_xy of the sliding-window chain as Gaussian rationals (Fraction pairs), built bit by bit."""
    k = len(templates[0])
    H = {}
    for w in range(N - k + 1):
        for t in templates:
            flip, phase_z, n_y = _string(t, w)[0], _string(t, w)[1], t.count("Y")
            for y in range(1 << N):
                # P|y> = i^{n_Y} (−1)^{popcount(z & y)} |y xor x>, with Y = iXZ on each site
                sign = -1 if bin(phase_z & y).count("1") % 2 else 1
                re, im = [(1, 0), (0, 1), (-1, 0), (0, -1)][n_y % 4]
                key = (y ^ flip, y)
                old = H.get(key, (0, 0))
                H[key] = (old[0] + sign * re, old[1] + sign * im)
    return {key: (Fraction(v[0]), Fraction(v[1])) for key, v in H.items() if v != (0, 0)}


def complement_connection(N, templates):
    """Theorem 2 of PROOF_PALINDROME_COMPLEMENT_CONNECTION, exactly: (components of Γ_H, components of Γ,
    good components, {failed condition: count})."""
    H = _basis_entries(N, templates)
    full = (1 << N) - 1
    zero = (Fraction(0), Fraction(0))
    entry = lambda x, y: H.get((x, y), zero)

    def components(adjacent):
        label, count = {}, 0
        for s in range(1 << N):
            if s in label:
                continue
            label[s], stack = count, [s]
            while stack:
                x = stack.pop()
                for y in adjacent[x]:
                    if y not in label:
                        label[y] = count
                        stack.append(y)
            count += 1
        return label, count

    hop = {x: [y for y in range(1 << N) if y != x and entry(x, y) != zero] for x in range(1 << N)}
    _, n_hop = components(hop)
    union = {x: sorted(set(hop[x]) | {y ^ full for y in hop[x ^ full]}) for x in range(1 << N)}
    label, n_union = components(union)
    good, failed = 0, {"modulus": 0, "diagonal": 0, "holonomy": 0}
    for c in range(n_union):
        vertices = [x for x in range(1 << N) if label[x] == c]
        bad = set()
        for x in vertices:
            if entry(x, x) != entry(x ^ full, x ^ full):
                bad.add("diagonal")
            for y in union[x]:
                a, b = entry(x, y), entry(x ^ full, y ^ full)
                if a[0] ** 2 + a[1] ** 2 != b[0] ** 2 + b[1] ** 2:
                    bad.add("modulus")
        if not bad:
            section, stack = {vertices[0]: (Fraction(1), Fraction(0))}, [vertices[0]]
            while stack and not bad:
                x = stack.pop()
                for y in union[x]:
                    value = _gauss_mul(_gauss_div(entry(x, y), entry(x ^ full, y ^ full)), section[x])
                    if y not in section:
                        section[y] = value
                        stack.append(y)
                    elif section[y] != value:
                        bad.add("holonomy")
        for reason in bad:
            failed[reason] += 1
        good += not bad
    return n_hop, n_union, good, failed


def square_law(N, templates=("XXXX", "XYYY", "YYYX")):
    """The holonomy around every square of two window flips, exactly, against the law read by hand.

    For each basis state y and windows w < w', the square y → y⊕m_w → y⊕m_w⊕m_w' → y⊕m_w' → y carries
    the product of the complement connection's ratios H_xy / H_x̄ȳ. Returns the number of squares that
    break the law: offset 1 and offset ≥ 4 always flat; offset 2 flat exactly when the six bits under
    the two windows have even parity. Offset 3 is not part of the law and is skipped.
    """
    H = _basis_entries(N, templates)
    full = (1 << N) - 1
    one = (Fraction(1), Fraction(0))
    ratio = lambda x, y: _gauss_div(H[(x, y)], H[(x ^ full, y ^ full)])
    masks = [_string(templates[0], w)[0] for w in range(N - len(templates[0]) + 1)]
    broken = 0
    for y in range(1 << N):
        for i, mi in enumerate(masks):
            for j in range(i + 1, len(masks)):
                offset, mj = j - i, masks[j]
                if offset == 3:
                    continue
                hol = one
                for a, b in ((y, y ^ mi), (y ^ mi, y ^ mi ^ mj), (y ^ mi ^ mj, y ^ mj), (y ^ mj, y)):
                    hol = _gauss_mul(hol, ratio(a, b))
                flat = hol == one
                six = (y >> i) & 0b111111
                expected = offset != 2 or bin(six).count("1") % 2 == 0
                broken += flat != expected
    return broken


def end_count_check():
    """F158 reads the crossing: the far end of XXXX + XYYY + YYYX shrinks as windows overlap.

    Near 8 at every N; far 8 at N = 4, 5, then 4 at N = 6 and 0 from N = 7: soft
    exactly where the spectrum pairs, hard from N = 6. No lit string commutes with
    H at any N (XXXX and XYYY differ on three letters of one window, so a lit
    string meets them with opposite parities). The controls: the two lit soft sets
    have equal ends at every N, carried by colourings; the hard control
    XXX + XXY + YXX has equal ends only at N = 3, its own crossing.
    """
    print("\n" + "=" * 78)
    print("F158 end counts (near, far, colourings), exact over the rationals, Z jumps on every site")
    print("=" * 78)
    witness = ["XXXX", "XYYY", "YYYX"]
    expected = {4: (8, 8, 0), 5: (8, 8, 0), 6: (8, 4, 0), 7: (8, 0, 0), 8: (8, 0, 0)}
    for N, want in expected.items():
        got = end_counts(N, witness)
        print(f"    N={N}  {'+'.join(witness):16s} {got}")
        assert got == want, f"N={N}: expected {want}, got {got}"
    controls = [
        (["XXY", "XYX", "YXX"], {N: (4, 4, 4) for N in range(3, 8)}),
        (["XXY", "YXX", "YYY"], {N: (4, 4, 4) for N in range(3, 8)}),
        (["XXX", "XXY", "YXX"], {3: (4, 4, 0), 4: (4, 2, 0), 5: (4, 0, 0), 6: (4, 0, 0), 7: (4, 0, 0)}),
    ]
    for templates, table in controls:
        for N, want in table.items():
            got = end_counts(N, templates)
            print(f"    N={N}  {'+'.join(templates):16s} {got}")
            assert got == want, f"{templates} N={N}: expected {want}, got {got}"
    print("\n  Theorem 2 on the basis-state graph (components of Γ_H, of Γ, good ones, failures):")
    for N in expected:
        n_hop, n_union, good, failed = complement_connection(N, witness)
        print(f"    N={N}  {n_hop} {n_union} {good}  {failed}")
        assert (n_hop, good) == expected[N][:2], f"N={N}: Theorem 2 must meet the string route"
        assert failed["modulus"] == failed["diagonal"] == 0, "only holonomy may fail"
        assert failed["holonomy"] == n_union - good
    print("\n  The square law (offset 1 flat; offset 2 flat iff the six bits under it have even parity):")
    for N in (5, 6, 7, 8):
        broken = square_law(N)
        print(f"    N={N}  squares breaking the law: {broken}")
        assert broken == 0, f"N={N}: {broken} squares break the hand-derived law"
    print("\nThe ends agree exactly where the spectrum pairs: XXXX+XYYY+YYYX at N = 4, 5 and not")
    print("from N = 6; the soft controls at every N; the hard control only at N = 3.")


if __name__ == "__main__":
    main()
