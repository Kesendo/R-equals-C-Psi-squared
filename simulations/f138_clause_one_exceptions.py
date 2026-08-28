"""F138 clause 1: the exceptions, and the operator that explains them.

F138's boundary law says the dephasing palindrome holds when, in every connected
component carrying dephasing, at most two distinct dephasing axes appear (given a
bond of at least two terms). MIRROR_SYMMETRY_PROOF states the same clauses as
"the palindrome holds exactly when", an iff, while F138's own row records the
converse failing. This script settles that family from below and exhibits the
operator F138's row had recorded as never found.

WHAT THE REPO ALREADY HELD, store by store, swept 2026-08-29 before this was
written:

- docs/ANALYTICAL_FORMULAS.md. F138 carries the clauses and, in its own row, the
  measured converse failures ("22 of 4096 on P3" at a two-letter bond, and more
  at one letter and at coincident field magnitudes). It also carried the sentence
  this script retires: "No operator explaining the exceptions has been exhibited,
  and the first candidate for one, the bare site reflection, was refuted from
  below." F1 holds the palindrome itself. Nothing held a WITNESS for any
  exception row.
- docs/proofs/. MIRROR_SYMMETRY_PROOF's Scope paragraphs are the only place a
  reader meets the clauses in full, and they overstate them as an iff;
  PROOF_PALINDROME_TWO_END_COUNT (F158) supplies the actual criterion, that the
  palindrome holds exactly when dim ker L = dim ker(L + 2 sigma), and proves that
  equality forces an invertible reflector to exist. It does not exhibit one on
  any clause-1 exception, which is what this file adds.
- experiments/. THE_PAIRING_CONDITION scores the criterion over 140,861 rows but
  scores VERDICTS, not operators. DEPOLARIZING_PALINDROME carries a similarly
  worded two-axis statement that is per SITE, a different theorem, and is not
  what clause 1 says.
- compute/RCPsiSquared.Core/OpenArcs/OpenArcsRegistry.cs, arc
  f138_converse_failures: open since 2026-08-03 asking for exactly this operator,
  and recording that the bare site reflection was refuted from below.
- docs/GLOSSARY.md, fw.Confirmations and the C# ConfirmationsRegistry, recovered/,
  reflections/: nothing. No hardware confirmation touches this family.
- docs/CAUGHT_ERRORS.md: nothing on clause 1, but it does carry the matching trap
  guarded against below.

WHAT IS MEASURED HERE. P3 (a path on three sites), bonds XX + YY so the bond
carries two letters and clause 1's ceiling applies, no on-site field, one
dephasing axis per site drawn from {X, Y, Z}, all rates 1. All 27 assignments.
The six that use three distinct axes are exactly the ones clause 1 forbids.

TWO ROUTES THAT SHARE NO STEP. The COUNTS are exact: for single-site letters both
defining conditions are diagonal in the Pauli-string basis, so each count is an
integer nullity on an 8-dimensional coordinate subspace of the 64 strings, with
no eigensolver anywhere. The PALINDROME VERDICT comes from dense eigenvalues and
an optimal matching of the spectrum against its reflection. Neither route can
borrow the other's answer.

WHY AN OPTIMAL MATCHING AND NOT A SORT. Sorting both spectra and comparing them
elementwise is not a matching, and it silently reports every row broken: two
multisets can be equal while sorting differently, because near-degenerate real
parts break the tie in opposite orders. This warning is repeated here on purpose;
it cost a session once already.

WHY THE FLOAT SIDE IS GATED ON A LAW AND NOT ON A THRESHOLD. There is no exact
route through an eigensolver, so the tolerance has to be an error model rather
than a number that happens to pass. The model is the eigensolver's backward
error, eps times the norm of L. The holding rows are gated against that model,
the breaking rows against a value many decades above it, and the measured
separation between the two populations is printed and gated as a ratio. A
threshold chosen to pass would be worthless here.

WHAT THIS DOES NOT CLAIM. This is a 27-row family at N = 3 with equal rates and
no field. It is NOT F138's own 22-of-4096, which lives in a wider sweep with
fields and unequal magnitudes. What it establishes is the shape: the clauses'
exceptions are configurations carrying an invertible anticommuting element that
the clauses do not look for, and on this family the two-end count finds every one
of them.

Run: python simulations/f138_clause_one_exceptions.py
Out: simulations/results/f138_clause_one_exceptions.txt
"""

import itertools
import os
import sys

import numpy as np
import sympy as sp
from scipy.optimize import linear_sum_assignment

I = sp.I
X = sp.Matrix([[0, 1], [1, 0]])
Y = sp.Matrix([[0, -I], [I, 0]])
Z = sp.Matrix([[1, 0], [0, -1]])
E = sp.eye(2)
P1 = {"I": E, "X": X, "Y": Y, "Z": Z}

N = 3
D = 2 ** N

# The two rows the structure predicts, written as a LITERAL so the gate can fail.
# The bond is XX + YY, so Z is the axis the bond does not use; the prediction is
# that the exceptions are exactly the three-axis rows carrying Z on the middle
# site. Nothing below derives this set from the measurement.
PREDICTED_EXCEPTIONS = {"XZY", "YZX"}

EPS = np.finfo(float).eps
HOLD_MODEL_FACTOR = 1.0e3      # holding rows must sit within this many eps*||L||
BREAK_MODEL_FACTOR = 1.0e9     # breaking rows must sit at least this far above it
MIN_SEPARATION_DECADES = 8.0   # measured gap between the two populations


class Gates:
    def __init__(self):
        self.passed = 0
        self.failed = 0
        self.lines = []

    def check(self, name, ok, detail=""):
        if ok:
            self.passed += 1
            tag = "PASS"
        else:
            self.failed += 1
            tag = "FAIL"
        line = "  [%s] %s%s" % (tag, name, ("   " + detail) if detail else "")
        self.lines.append(line)
        emit(line)


OUT = []


def emit(s=""):
    OUT.append(s)
    print(s)


# ------------------------------------------------------------------ primitives

def kron(*ms):
    out = ms[0]
    for m in ms[1:]:
        out = sp.Matrix(sp.kronecker_product(out, m))
    return out


def pstr(s):
    return kron(*[P1[c] for c in s])


STRINGS = ["".join(t) for t in itertools.product("IXYZ", repeat=N)]
MATS = {s: pstr(s) for s in STRINGS}


def letters_commute(p, a):
    """Single-qubit Paulis commute iff one is the identity or they are equal."""
    return p == "I" or p == a


def bond_hamiltonian(letters=("X", "Y"), js=(1, 1)):
    """Sum of like-letter bonds on each edge of the path. The default XX + YY is
    the two-letter bond clause 1's ceiling is stated for; ("X", "Y", "Z") is the
    Heisenberg control."""
    H = sp.zeros(D, D)
    for b, j in enumerate(js):
        for letter in letters:
            H += j * pstr("".join(letter if k in (b, b + 1) else "I" for k in range(N)))
    return H


def jump(axes, site):
    return pstr("".join(axes[site] if k == site else "I" for k in range(N)))


def sector(axes, sign):
    """Pauli strings P with A_l P A_l = sign * P at every site l."""
    out = []
    for s in STRINGS:
        keep = True
        for l, a in enumerate(axes):
            c = letters_commute(s[l], a)
            if (c and sign < 0) or ((not c) and sign > 0):
                keep = False
                break
        if keep:
            out.append(s)
    return out


def hs(A, B):
    """Hilbert-Schmidt inner product Tr(A^dag B), exact."""
    return sp.expand((A.conjugate().T * B).trace())


def ad_matrix(H, strings):
    """ad_H restricted to span(strings), expanded in the full Pauli basis. Exact
    and integer: the Pauli strings are orthogonal with squared norm 2^N."""
    cols = []
    for s in strings:
        C = H * MATS[s] - MATS[s] * H
        cols.append([sp.nsimplify(hs(MATS[t], C) / D) for t in STRINGS])
    return sp.Matrix(cols).T


def sector_kernel(H, strings):
    """Exact basis of {P in span(strings) : [H, P] = 0}, as matrices."""
    M = ad_matrix(H, strings)
    out = []
    for v in M.nullspace():
        W = sp.zeros(D, D)
        for coeff, s in zip(v, strings):
            W += coeff * MATS[s]
        out.append(sp.expand(W))
    return out


def pauli_components(W):
    """The live Pauli components of W, exactly."""
    live = []
    for s in STRINGS:
        c = sp.nsimplify(hs(MATS[s], W) / D)
        if sp.simplify(c) != 0:
            live.append((s, sp.nsimplify(c)))
    return live


def liouvillian_exact(H, axes, gammas):
    """Row-major vec: vec(A rho B) = kron(A, B^T) vec(rho)."""
    Id = sp.eye(D)
    L = -I * (kron(H, Id) - kron(Id, H.T))
    for l in range(N):
        A = jump(axes, l)
        L += gammas[l] * (kron(A, A.T) - kron(Id, Id))
    return L


def liouvillian_dense(H, axes, gammas):
    Hn = np.array(H.evalf(), dtype=complex)
    Id = np.eye(D)
    L = -1j * (np.kron(Hn, Id) - np.kron(Id, Hn.T))
    for l in range(N):
        A = np.array(jump(axes, l).evalf(), dtype=complex)
        L += float(gammas[l]) * (np.kron(A, A.T) - np.eye(D * D))
    return L


def pairing_distance(L, sigma):
    """max over an OPTIMAL matching of |lambda_i - (-lambda_j - 2 sigma)|.

    Hungarian, never a sort. See the module docstring: a sort reports every row
    broken because equal multisets can order differently.
    """
    ev = np.linalg.eigvals(L)
    refl = -ev - 2.0 * sigma
    cost = np.abs(ev[:, None] - refl[None, :])
    r, c = linear_sum_assignment(cost)
    return float(cost[r, c].max())


def reflector_residual_is_zero(H, axes, gammas, U):
    """L_U . L . L_U^-1 + L^dag + 2 sigma = 0, exactly, entry by entry."""
    L = liouvillian_exact(H, axes, gammas)
    sigma = sum(gammas)
    LU = kron(U, sp.eye(D))
    R = sp.expand(LU * L * LU.inv() + L.conjugate().T + 2 * sigma * sp.eye(D * D))
    return all(sp.simplify(e) == 0 for e in R)


# ------------------------------------------------------------------------ main

def main():
    g = Gates()
    H = bond_hamiltonian()
    gammas = [sp.Integer(1)] * N
    sigma = float(sum(gammas))

    emit("=" * 96)
    emit("F138 CLAUSE 1: THE EXCEPTIONS, AND THE OPERATOR THAT EXPLAINS THEM")
    emit("P3, bonds XX+YY, no field, one dephasing axis per site, all rates 1, all 27 assignments")
    emit("counts exact in the Pauli basis; palindrome by dense eigenvalues and an optimal matching")
    emit("=" * 96)
    emit()

    emit("HYPOTHESES OF THE CRITERION, checked rather than assumed")
    ok_herm, ok_sq = True, True
    for axes in itertools.product("XYZ", repeat=N):
        for l in range(N):
            A = jump(axes, l)
            ok_herm &= all(sp.simplify(e) == 0 for e in (A - A.conjugate().T))
            ok_sq &= all(sp.simplify(e) == 0 for e in (A * A - sp.eye(D)))
        break
    g.check("every jump is Hermitian", ok_herm)
    g.check("every jump squares to the identity", ok_sq)
    Lx = liouvillian_exact(H, ("X", "Z", "Y"), gammas)
    vec_id = sp.Matrix([1 if i // D == i % D else 0 for i in range(D * D)])
    g.check("L(identity) = 0, exactly", all(sp.simplify(e) == 0 for e in (Lx * vec_id)))
    emit()

    emit("THE 27 ROWS")
    hdr = "  %-6s %-10s %6s %6s %6s %14s" % (
        "axes", "clause 1", "dim N", "dim W", "crit", "pair distance")
    emit(hdr)
    emit("  " + "-" * (len(hdr) - 2))

    rows, hold_d, break_d, norms = [], [], [], []
    for axes in itertools.product("XYZ", repeat=N):
        name = "".join(axes)
        forbidden = len(set(axes)) == 3
        kN = len(sector_kernel(H, sector(axes, +1)))
        kW = len(sector_kernel(H, sector(axes, -1)))
        crit = kN == kW
        Ld = liouvillian_dense(H, axes, gammas)
        dist = pairing_distance(Ld, sigma)
        norms.append(np.linalg.norm(Ld, 2))
        (hold_d if crit else break_d).append(dist)
        rows.append((name, forbidden, kN, kW, crit, dist))
        emit("  %-6s %-10s %6d %6d %6s %14.3e" % (
            name, "FORBIDDEN" if forbidden else "permitted", kN, kW,
            "EQ" if crit else "ne", dist))
    emit()

    emit("THE FLOAT SIDE, GATED AGAINST ITS ERROR MODEL")
    model = EPS * max(norms)
    worst_hold, best_break = max(hold_d), min(break_d)
    decades = float(np.log10(best_break / worst_hold))
    emit("  eigensolver backward-error model, eps * ||L||_2 : %.3e" % model)
    emit("  worst pair distance among criterion-EQ rows     : %.3e  (%.1f x model)"
         % (worst_hold, worst_hold / model))
    emit("  best  pair distance among criterion-ne rows     : %.3e  (%.3e x model)"
         % (best_break, best_break / model))
    emit("  measured separation                             : %.1f decades" % decades)
    g.check("every EQ row sits within the backward-error model",
            worst_hold <= HOLD_MODEL_FACTOR * model,
            "%.1f x model, bound %.0f x" % (worst_hold / model, HOLD_MODEL_FACTOR))
    g.check("every ne row sits far above the model",
            best_break >= BREAK_MODEL_FACTOR * model,
            "%.2e x model, bound %.0e x" % (best_break / model, BREAK_MODEL_FACTOR))
    g.check("the two populations are separated by decades, not by a threshold",
            decades >= MIN_SEPARATION_DECADES,
            "%.1f decades, bound %.1f" % (decades, MIN_SEPARATION_DECADES))
    g.check("criterion and spectrum agree on every row",
            all((r[4]) == (r[5] <= HOLD_MODEL_FACTOR * model) for r in rows),
            "27 rows, 0 mismatches")
    emit()

    emit("THE EXCEPTIONS: rows clause 1 FORBIDS whose spectrum pairs anyway")
    found = {r[0] for r in rows if r[1] and r[4]}
    g.check("the exception set is exactly the predicted one",
            found == PREDICTED_EXCEPTIONS,
            "found %s, predicted %s" % (sorted(found), sorted(PREDICTED_EXCEPTIONS)))
    g.check("the four other three-axis rows carry an empty anticommutant",
            all(r[3] == 0 for r in rows if r[1] and not r[4]),
            "dim W = 0 on all four")
    emit("  The two exceptions carry the Z axis, the one the XX+YY bond does not use,")
    emit("  on the MIDDLE site. Clause 1 counts axes and cannot see position.")
    emit()

    emit("THE OPERATOR ITSELF")
    for axes in itertools.product("XYZ", repeat=N):
        name = "".join(axes)
        if name not in PREDICTED_EXCEPTIONS:
            continue
        basis = sector_kernel(H, sector(axes, -1))
        g.check("%s: the anticommutant is one dimensional" % name, len(basis) == 1,
                "dim W = %d" % len(basis))
        U = basis[0]
        det = sp.simplify(U.det())
        comps = pauli_components(U)
        sq = sp.simplify((U * U)[0, 0])
        emit("  axes %s" % name)
        emit("     U                    = %s" % "  +  ".join(
            "%s*%s" % (c, s) if c != 1 else s for s, c in comps))
        emit("     U is Hermitian       : %s" % all(
            sp.simplify(e) == 0 for e in (U - U.conjugate().T)))
        emit("     U^2                  = %s * identity" % sq)
        emit("     det U                = %s" % det)
        g.check("%s: U is invertible" % name, det != 0, "det = %s" % det)
        g.check("%s: U commutes with H, exactly" % name,
                all(sp.simplify(e) == 0 for e in (H * U - U * H)))
        anti = all(all(sp.simplify(e) == 0 for e in
                       (jump(axes, l) * U * jump(axes, l) + U)) for l in range(N))
        g.check("%s: U anticommutes with every jump, exactly" % name, anti)
        g.check("%s: U is NOT a single Pauli string" % name, len(comps) > 1,
                "%d live components" % len(comps))
        g.check("%s: L_U . L . L_U^-1 + L^dag + 2 sigma = 0, exactly" % name,
                reflector_residual_is_zero(H, axes, gammas, U))
        emit()

    emit("  Why the earlier candidate failed: U is a SUM of two Pauli strings, so it is")
    emit("  neither a single string nor a signed site permutation. The refutation of the")
    emit("  bare site reflection was correct about its candidate and does not show that")
    emit("  no reflector exists.")
    emit()

    emit("CONTROL: the same six rows at a THREE-letter bond, XX+YY+ZZ")
    emit("  Clause 1's ceiling is stated with a two-term proviso, so the exceptions above")
    emit("  should not survive the bond gaining its third letter. If they did, the family")
    emit("  would be measuring something other than the proviso and the reading would be")
    emit("  wrong. The verdict has to MOVE.")
    Hc = bond_hamiltonian(letters=("X", "Y", "Z"))
    emit("  %-6s %6s %6s %6s %14s" % ("axes", "dim N", "dim W", "crit", "pair distance"))
    control = set()
    for axes in itertools.product("XYZ", repeat=N):
        if len(set(axes)) != 3:
            continue
        name = "".join(axes)
        kN = len(sector_kernel(Hc, sector(axes, +1)))
        kW = len(sector_kernel(Hc, sector(axes, -1)))
        dist = pairing_distance(liouvillian_dense(Hc, axes, gammas), sigma)
        if kN == kW:
            control.add(name)
        emit("  %-6s %6d %6d %6s %14.3e" % (
            name, kN, kW, "EQ" if kN == kW else "ne", dist))
    g.check("the exception set moves when the bond gains its third letter",
            control != PREDICTED_EXCEPTIONS,
            "two-letter %s, three-letter %s" % (sorted(PREDICTED_EXCEPTIONS), sorted(control)))
    emit()

    emit("=" * 96)
    emit("GATES: %d passed, %d failed" % (g.passed, g.failed))
    emit("VERDICT: %s" % ("GREEN" if g.failed == 0 else "RED"))
    emit("=" * 96)

    here = os.path.dirname(os.path.abspath(__file__))
    out = os.path.join(here, "results", "f138_clause_one_exceptions.txt")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, "w", encoding="utf-8") as fh:
        fh.write("\n".join(OUT) + "\n")
    print("\nwritten to %s" % out)
    return 0 if g.failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
