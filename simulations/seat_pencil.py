"""The watched seat's pencil: the rank-one identity behind F157's count.

F157 counts what watching one seat j misses, blind(j) = N - dim Krylov(e_j) =
deg gcd(chi(H), chi(H struck at j)).  This script gates a third route to that
already-proved count, and the coverage claim that rides with it.

The generator is the repo's own, `f64_topology_scan.single_excitation_coherence_L`:

    L_coh(j, gamma)  =  -i H  -  2 gamma P_j ,       P_j = e_j e_j^T ,

acting on the single-excitation coherences |1_i><vac|.

WHAT IS GATED, AND WHAT IS ONLY STATED.  This distinction is the point of the
file.  Two review rounds found gates here that could not fail, both times because
a property the CONSTRUCTION forces was dressed as a measurement, so the inventory
below is deliberately short and says what each check is worth.

  GATED, and able to fail:
    G1  the pencil identity det(lambda I - H + s E_jj) = chi_H + s * chi_cut.
        The left side is an elimination determinant and the right side is
        Faddeev-LeVerrier, two routes written separately in this file, so a slip
        in either shows.  This is the one strong gate.
    G4  the zero-bond fence, on the exhaustive set the committed pricing was
        measured on, AND the containment that came out of running it: the
        Krylov blind count and the span minus one are different integers on
        1768 of the 3996 (book, profile, seat) triples, and EVERY ONE of those
        disagreements lands on a triple the two-halves face already misses.  So
        the face is right only where the two truths coincide.  That is why its
        right-count is 692 under either truth, and the coincidence is a
        consequence of the containment rather than a second finding.

  GATED, but weak, and labelled as such rather than dropped:
    G2  dim W over Q here against `blind_truth`'s GF(p) rank.  These are NOT
        independent routes: `blind_seat_span_proof.blind_basis` and
        `seat_cut_blindness.blind_truth` build the SAME Krylov matrix, line for
        line, and differ only in which field the rank is taken over.  So G2 can
        catch a bad prime pair and nothing else.  It is kept for that and is
        worth exactly that.

  STATED, NOT GATED, because a gate for them could only fail if elimination
  itself were broken:
    - P_j annihilates W.  `blind_basis` returns the nullspace of the Krylov
      matrix whose FIRST ROW is e_seat, so w[seat] = 0 is one of the equations
      solved.
    - K is L_coh-invariant, and so is W.  P_j v = e_j <e_j, v> lies in K because
      e_j does; H K is in K by construction; W = K^perp is H-invariant because H
      is Hermitian.
    - H restricted to W is self-adjoint.  An earlier version gated this as
      M^T G = G M against the basis's Gram matrix.  That test is IDENTICALLY the
      symmetry of H, since (G M)[a][b] = <b_a, H b_b> and (M^T G)[a][b] =
      <H b_a, b_b>; it holds on ANY H-invariant subspace in ANY basis, and
      `se_hamiltonian_int` writes h[a][b] and h[b][a] in one statement.  It was
      removed rather than reworded.
    - hence L_coh|_W = -i H|_W is skew-adjoint and every rate on W is exactly
      zero at every gamma, while on K no purely imaginary eigenvalue can exist,
      so C^N = W (+) K splits the generator and the undamped space is exactly W.
      gamma never multiplies anything but zero on W, which is why no gate here
      can exercise "at every gamma": that is an algebraic consequence, not a
      measurement.

Primitives are imported from `seat_cut_blindness` and `blind_seat_span_proof`
rather than restated; only the two gcd faces and the two determinant routes are
written here, so this is not a share-nothing reproduction.

Everything gated is exact: integer couplings, Fraction arithmetic, no eigensolver
and no tolerance anywhere.  Primitives are imported from `seat_cut_blindness` and
`blind_seat_span_proof` rather than restated.

Run:  python simulations/seat_pencil.py [pencil | dim | fence | all]
"""

import itertools
import sys
from fractions import Fraction

import seat_cut_blindness as S
import blind_seat_span_proof as B


# ----------------------------------------------------------------- polynomials

def _det_q(mat):
    """Exact determinant over Q by elimination.

    Deliberately NOT Faddeev-LeVerrier: G1 compares this against a
    Faddeev-LeVerrier characteristic polynomial, and two routes that shared an
    implementation would agree through a shared bug.
    """
    m = [[Fraction(x) for x in row] for row in mat]
    n = len(m)
    if n == 0:
        return Fraction(1)
    sign, det = 1, Fraction(1)
    for c in range(n):
        piv = next((r for r in range(c, n) if m[r][c] != 0), None)
        if piv is None:
            return Fraction(0)
        if piv != c:
            m[c], m[piv] = m[piv], m[c]
            sign = -sign
        det *= m[c][c]
        inv = m[c][c]
        for r in range(c + 1, n):
            f = m[r][c] / inv
            if f:
                for k in range(c, n):
                    m[r][k] -= f * m[c][k]
    return sign * det


def _charpoly(mat):
    """Monic characteristic polynomial over Q, ASCENDING coefficients."""
    n = len(mat)
    if n == 0:
        return [Fraction(1)]
    M = [[Fraction(int(i == j)) for j in range(n)] for i in range(n)]
    coeffs = [Fraction(1)]
    for k in range(1, n + 1):
        AM = [[sum(Fraction(mat[i][t]) * M[t][j] for t in range(n))
               for j in range(n)] for i in range(n)]
        c = -sum(AM[i][i] for i in range(n)) / k
        coeffs.append(c)
        M = [[AM[i][j] + (c if i == j else Fraction(0)) for j in range(n)]
             for i in range(n)]
    return coeffs[::-1]


def _pval(poly, x):
    out = Fraction(0)
    for c in reversed(poly):
        out = out * x + c
    return out


def _struck(mat, j):
    return [[mat[r][c] for c in range(len(mat)) if c != j]
            for r in range(len(mat)) if r != j]


def _sub(mat, lo, hi):
    return [row[lo:hi] for row in mat[lo:hi]]


def _gcd_deg(f, g):
    return len(B._poly_gcd(list(f), list(g))) - 1


def _bonds(js):
    return [(i, i + 1, js[i]) for i in range(len(js))]


def exhaustive_profiles(nmax=6):
    """Every profile in {0,1,2}^(N-1) for N = 3..nmax.

    At the default nmax = 6 this is the set the committed fence pricing was
    measured on (THE_BLIND_SITE.md section 11, ANALYTICAL_FORMULAS.md F4's seat
    bullets); G1 and G2 call it with nmax = 5 and make no such claim.  There, and
    G4 reproduces those numbers rather than sampling a new set: a random draw
    would give a DIFFERENT number for the same fact and leave two prices in the
    repo with no way to reconcile them.
    """
    out = []
    for n in range(3, nmax + 1):
        for js in itertools.product((0, 1, 2), repeat=n - 1):
            out.append((n, list(js)))
    return out


# ------------------------------------------------------------------- the parts

def run_pencil():
    """G1: the pencil identity."""
    print("THE PENCIL.  Watching one seat enters the generator as a rank-one")
    print("change of ONE diagonal entry, so the characteristic polynomial is")
    print("AFFINE in it:")
    print()
    print("    det(lambda I - H + s E_jj)  =  chi_H(lambda) + s * chi_cut(lambda)")
    print()
    print("  with chi_cut the polynomial of H struck at row and column j.  The")
    print("  physical case is s = 2 i gamma for the generator")
    print("  L_coh = -i H - 2 gamma P_j (f64_topology_scan.py:128), and")
    print("  s = -2 i gamma for the +i H convention of EQ-015, so the sign of s")
    print("  is convention-dependent while the identity is not.")
    print()
    print("  G1  the identity over Q, at every seat of the twenty graphs and of")
    print("      the exhaustive chain profiles, both books.  The left side is an")
    print("      elimination determinant and the right side is Faddeev-LeVerrier,")
    print("      so the two sides do not share an implementation.  Both have")
    print("      degree at most N in lambda and at most 1 in s, so agreement at")
    print("      N+1 distinct lambda and 2 distinct s PROVES the identity; three")
    print("      values of s are used, the third checking the first two.")
    print()
    tally = [0, 0]
    failures = []
    lam_nodes = [Fraction(v) for v in (0, 1, -1, 2, -2, 3, -3, 5, -5, 7, -7)]
    s_nodes = [Fraction(1), Fraction(-2), Fraction(3, 5)]
    cases = [(name, n, bonds) for name, n, bonds in S.DELETED_GRAPHS]
    cases += [(f"chain{n}:{js}", n, _bonds(js)) for n, js in exhaustive_profiles(5)]
    for zz in (True, False):
        for name, n, bonds in cases:
            h = S.se_hamiltonian_int(n, bonds, zz)
            chi = _charpoly(h)
            for seat in range(n):
                cut = _charpoly(_struck(h, seat))
                for s in s_nodes:
                    for lam in lam_nodes[:n + 1]:
                        m = [[-Fraction(h[r][c]) + (lam if r == c else 0)
                              + (s if (r == c == seat) else 0)
                              for c in range(n)] for r in range(n)]
                        ok = _det_q(m) == _pval(chi, lam) + s * _pval(cut, lam)
                        tally[1 if ok else 0] += 1
                        if not ok:
                            failures.append((name, zz, seat, str(s), str(lam)))
    print(f"  G1: pass {tally[1]:7d}   FAIL {tally[0]}   "
          f"over (book, graph, seat, s, lambda) points")
    for f in failures[:5]:
        print(f"      FAIL {f}")
    return not failures


def run_dim():
    """G2: the one thing left about W that a check can decide, and it is weak."""
    print("THE DIMENSION.  W = K^perp with K the seat's Krylov space.  The")
    print("invariance and self-adjointness steps are DEFINITIONAL and are stated")
    print("in the module docstring rather than gated; an earlier version gated")
    print("them and they could not fail.  What is left is one weak check, kept")
    print("because it is worth something and labelled so it is not read as more.")
    print()
    print("  G2  dim W over Q here, against blind_truth's GF(p) rank at two")
    print("      primes.  These build the SAME Krylov matrix, line for line, and")
    print("      differ only in the field the rank is taken over, so G2 detects a")
    print("      bad prime pair and NOTHING ELSE.  It cannot see a wrong H, a")
    print("      wrong Krylov matrix or a wrong subspace, because both sides")
    print("      would be wrong together.")
    print()
    tally = [0, 0]
    failures = []
    counted = {"triples": 0, "nonempty": 0}
    cases = [(name, n, bonds) for name, n, bonds in S.DELETED_GRAPHS]
    cases += [(f"chain{n}:{js}", n, _bonds(js)) for n, js in exhaustive_profiles(5)]
    for zz in (True, False):
        for name, n, bonds in cases:
            for seat in range(n):
                counted["triples"] += 1
                _, basis = B.blind_basis(n, bonds, seat, zz)
                ok = len(basis) == S.blind_truth(n, bonds, seat, zz)
                tally[1 if ok else 0] += 1
                if not ok:
                    failures.append((name, "ZZ" if zz else "XY", seat))
                if basis:
                    counted["nonempty"] += 1
    print(f"  G2: pass {tally[1]:6d}   FAIL {tally[0]}   "
          f"over (book, graph, seat) triples, one note each, the twenty graphs")
    print(f"      plus every profile in {{0,1,2}}^(N-1) for N = 3..5, both books;")
    print(f"      {counted['nonempty']} of the {counted['triples']} have dim W > 0")
    for f in failures[:5]:
        print(f"      FAIL {f}")
    return not failures


def run_fence():
    """G4: the committed zero-bond pricing, and the finding running it produced."""
    print("THE FENCE.  F157 carries two faces of one count: the STRUCK face")
    print("deg gcd(chi_H, chi_cut), which carries no fence, and the TWO-HALVES")
    print("face deg gcd(chi_L, chi_R), fenced to profiles with no zero bond.")
    print("The pricing is NOT new: it is committed at THE_BLIND_SITE.md section")
    print("11 and ANALYTICAL_FORMULAS.md F4's seat bullets, as 1682 zero-bond")
    print("(profile, seat) pairs PER BOOK with the Heisenberg book wrong on all")
    print("of them and XY right on 60, and 316 zero-free pairs right on both.")
    print("Counted over both books these are (book, profile, seat) TRIPLES, and")
    print("the totals below are triple counts: 1998 pairs give 3996 triples.")
    print("G4 runs the exhaustive set {0,1,2}^(N-1) for N = 3..6 that those")
    print("numbers were measured on.  This file re-implements the two gcd faces")
    print("and both determinant routes; the Hamiltonian and the two truths are")
    print("imported, so this is not a share-nothing reproduction.")
    print()
    print("  AND THE TRUTH SIDE IS NOT THE ONE SECTION 11 USED.  Section 11's")
    print("  sweep compares against the SPAN minus one; this compares against the")
    print("  KRYLOV blind count.  They differ on a large minority of the triples,")
    print("  so the tally below is computed under BOTH, and the containment of")
    print("  the disagreements is counted rather than assumed away.")
    print()
    rows = {}
    for zz in (True, False):
        for kind in ("zero-free", "zero-bond"):
            # pairs, struck right / halves right, under each of the two truths
            rows[("Heisenberg" if zz else "XY", kind)] = [0, 0, 0, 0, 0]
    disagree = 0
    leak = 0          # disagreeing triples where the halves face is RIGHT
    for n, js in exhaustive_profiles(6):
        kind = "zero-bond" if 0 in js else "zero-free"
        bonds = _bonds(js)
        for zz in (True, False):
            book = "Heisenberg" if zz else "XY"
            h = S.se_hamiltonian_int(n, bonds, zz)
            chi = _charpoly(h)
            for seat in range(n):
                krylov = S.blind_truth(n, bonds, seat, zz)
                span = S.exact_kernel_dim(n, bonds, [seat], zz) - 1
                disagree += (krylov != span)
                struck = _gcd_deg(chi, _charpoly(_struck(h, seat)))
                halves = _gcd_deg(_charpoly(_sub(h, 0, seat)),
                                  _charpoly(_sub(h, seat + 1, n)))
                r = rows[(book, kind)]
                r[0] += 1
                r[1] += (struck == krylov)
                r[2] += (halves == krylov)
                r[3] += (struck == span)
                r[4] += (halves == span)
                if krylov != span and (halves == krylov or halves == span):
                    leak += 1
    hdr = (f"{'book':12s}{'profiles':18s}{'pairs':>7}{'struck=Krylov':>15}"
           f"{'halves=Krylov':>15}{'struck=span':>13}{'halves=span':>13}")
    print(hdr)
    for book in ("Heisenberg", "XY"):
        for kind in ("zero-free", "zero-bond"):
            r = rows[(book, kind)]
            print(f"{book:12s}{kind:18s}{r[0]:>7}{r[1]:>15}{r[2]:>15}"
                  f"{r[3]:>13}{r[4]:>13}")
    tot = [sum(rows[k][i] for k in rows) for i in range(5)]
    print(f"{'total':30s}{tot[0]:>7}{tot[1]:>15}{tot[2]:>15}{tot[3]:>13}{tot[4]:>13}")
    print()
    free = sum(rows[k][0] for k in rows if k[1] == "zero-free")
    zero = tot[0] - free
    print(f"  the two-halves face is right on all {free} zero-free triples and on")
    print(f"  {tot[2] - free} of the {zero} that carry a zero bond.")
    print(f"  the two truths DISAGREE on {disagree} of the {tot[0]} triples, "
          f"{100.0 * disagree / tot[0]:.1f} %,")
    print(f"  and the two-halves face is right on {leak} of those {disagree}.")
    if leak == 0:
        print("  So the face is right ONLY where the two truths already coincide:")
        print("  every disagreement lands on a triple the face already misses.")
        print(f"  Its right-count being {tot[2]} under either truth follows from that")
        print("  containment and is not a second finding.")
    else:
        print("  So the containment does NOT hold, and the right-counts agreeing")
        print("  under the two truths would need a different explanation.")
    print()
    expected = {("Heisenberg", "zero-bond"): (1682, 0),
                ("XY", "zero-bond"): (1682, 60),
                ("Heisenberg", "zero-free"): (316, 316),
                ("XY", "zero-free"): (316, 316)}
    ok = (disagree == 1768 and leak == 0)
    if disagree != 1768:
        print(f"      the disagreement count is {disagree}, committed 1768")
    if leak != 0:
        print(f"      {leak} disagreeing triples carry a RIGHT two-halves face; "
              f"the containment is broken")
    for k, (pairs, halves_right) in expected.items():
        r = rows[k]
        # The struck face must equal the KRYLOV count everywhere, which is
        # F157's theorem.  It is NOT asserted against the span: the span
        # identity is exactly what fails off the zero-free profiles, and the
        # column above prices that failure rather than hiding it.
        if not (r[0] == pairs and r[2] == halves_right and r[4] == halves_right
                and r[1] == pairs):
            ok = False
            print(f"      {k}: got triples {r[0]}, halves {r[2]}/{r[4]}, "
                  f"struck-vs-Krylov {r[1]}; committed {pairs} and {halves_right}")
    if ok:
        print("  G4: the committed pricing is reproduced under BOTH truths; the")
        print(f"      struck face equals the Krylov count on all {tot[0]} triples")
        print(f"      and the span minus one on {tot[3]}, that gap being the span")
        print("      identity's own failure off the zero-free profiles, priced;")
        print(f"      and the {disagree} disagreements are contained in the triples")
        print("      the two-halves face already misses.")
    else:
        print("  G4: MISMATCH against the committed pricing")
    return ok


def main():
    part = sys.argv[1] if len(sys.argv) > 1 else "all"
    ok = True
    for name, fn in (("pencil", run_pencil), ("dim", run_dim),
                     ("fence", run_fence)):
        if part in (name, "all"):
            ok &= bool(fn())
            print()
    print("ALL GATES PASS" if ok else "SOME GATE FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
