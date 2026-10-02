"""Exact gate for the common-carrier proposition in the complement proof.

Run: python simulations/triangle_common_carrier_gate.py
Requires SymPy. No result files are written.

For real nonzero s,t,r this proves that the common ZII-odd commutant of
H_b = s(XXI+YYI+ZZI) + t(XIX+YIY+ZIZ) + r(IXX+IYY+IZZ)
and F = IXI+IIY is a line iff (s-t)(st+sr+tr)=0, and zero otherwise.
The line's displayed generator is invertible. This does not classify carriers
that depend on a fixed field strength in H_b+hF.

Necessity uses an explicitly checked monomial minor, not a generic nullspace.
Sparse Pauli commutators are checked against independently built dense matrices.
Off-locus, unequal-field and changed-carrier inputs exercise the same equations.
"""

from itertools import product

import sympy as sp

from framework.pauli import pauli_product


LETTERS = "IXYZ"
WORDS = ["".join(w) for w in product(LETTERS, repeat=3)]
INDEX = {w: i for i, w in enumerate(WORDS)}
PAULI = {
    "I": sp.eye(2),
    "X": sp.Matrix([[0, 1], [1, 0]]),
    "Y": sp.Matrix([[0, -sp.I], [sp.I, 0]]),
    "Z": sp.diag(1, -1),
}


def dense(terms):
    return sum(
        (c * sp.kronecker_product(*(PAULI[p] for p in w)) for w, c in terms.items()),
        sp.zeros(8),
    )


def commutator(a, b):
    return a * b - b * a


def reduced_commutator(h, q):
    """Pauli coefficients of [h,q]/(2i), for real Pauli coefficients."""
    out = sp.zeros(64, 1)
    for a, ha in h.items():
        for b, qb in q.items():
            e, word = pauli_product(a, b)
            if e % 2:
                out[INDEX[word]] += (1 if e == 1 else -1) * ha * qb
    return out.applyfunc(sp.expand)


def is_zero(matrix):
    return all(sp.expand(v) == 0 for v in matrix)


def check(condition, message):
    if not condition:
        raise AssertionError(message)
    print("PASS", message)


def main():
    s, t, r = sp.symbols("s t r", real=True)
    h = {"XXI": s, "YYI": s, "ZZI": s,
         "XIX": t, "YIY": t, "ZIZ": t,
         "IXX": r, "IYY": r, "IZZ": r}
    field = {"IXI": 1, "IIY": 1}
    hb, f, jump = dense(h), dense(field), dense({"ZII": 1})
    # Independent basis of the two-free-site field commutant.
    e = [{"II": 1}, {"XI": 1}, {"IY": 1}, {"XY": 1},
         {"YX": 1, "ZZ": -1}, {"YZ": 1, "ZX": 1}]
    basis = [{p + w: c for w, c in term.items()} for p in "XY" for term in e]
    matrices = [dense(q) for q in basis]
    lit = [dense({p + a + b: 1}) for p in "XY" for a in LETTERS for b in LETTERS]
    field_map = sp.Matrix.hstack(*(commutator(f, q).reshape(64, 1) for q in lit))
    check(field_map.rank() == 20, "field condition leaves exactly 12 of the 32 lit coordinates")
    check(sp.Matrix.hstack(*(q.reshape(64, 1) for q in matrices)).rank() == 12
          and all(is_zero(commutator(f, q)) and is_zero(jump*q + q*jump) for q in matrices),
          "displayed basis spans the full field-commuting lit space")

    m = sp.Matrix.hstack(*(reduced_commutator(h, q) for q in basis))
    dense_matches = []
    for col, q in enumerate(matrices):
        rebuilt = dense({w: 2*sp.I*m[i, col] for i, w in enumerate(WORDS) if m[i, col] != 0})
        dense_matches.append(is_zero(commutator(hb, q) - rebuilt))
    check(len(dense_matches) == 12 and all(dense_matches),
          "all 12 symbolic sparse commutators equal the independent dense commutators")

    # Removing coefficient a0 leaves a minor nonzero throughout s*t*r != 0.
    selected = ["IIZ", "IXX", "IYY", "IZI", "XYZ", "YZX",
                "IYZ", "IZY", "ZIX", "ZIY", "IXZ"]
    a = m[[INDEX[w] for w in selected], :]
    minor = sp.factor(a[:, 1:].det())
    check(minor == -r**2*s**3*t**6, "necessity minor is -r^2*s^3*t^6; no exceptional nonzero-coupling branch omitted")
    coefficients = sp.Matrix([s*t**2, 0, 0, s**3, s*t**2, 0,
                              s**2*t, 0, 0, t**3, s**2*t, 0])
    check(is_zero(a*coefficients), "selected equations have the displayed line as their full kernel")
    e2 = s*t + s*r + t*r
    residual = sp.zeros(64, 1)
    for word, c in {"XIZ": s, "XZI": -s, "YIZ": -t, "YZI": t}.items():
        residual[INDEX[word]] = c*(s-t)*e2
    check(is_zero(m*coefficients - residual), "remaining equations are exactly (s-t)*(st+sr+tr) times the stated nonzero vector")

    # Independently type the operator shown in the proof, rather than derive it
    # from the sparse kernel under test.
    w = dense({"XII": s*t**2, "XXY": s**3, "XYX": s*t**2, "XZZ": -s*t**2,
               "YII": s**2*t, "YXY": t**3, "YYX": s**2*t, "YZZ": -s**2*t})
    check(is_zero(w - sum((c*q for c, q in zip(coefficients, matrices)), sp.zeros(8))),
          "printed carrier equals the unique candidate from the necessity equations")
    check(is_zero(commutator(hb, w) - 2*sp.I*(s-t)*e2*dense({"XIZ": s, "XZI": -s, "YIZ": -t, "YZI": t})),
          "dense bond commutator has the claimed factorization")
    check(is_zero(commutator(f, w)) and is_zero(jump*w + w*jump),
          "carrier commutes with the field and anticommutes with the jump")
    check(is_zero(w*w - (s*s+t*t)**3*sp.eye(8)),
          "carrier square is (s^2+t^2)^3 I, strictly positive for real nonzero s,t")
    swap = dense({"III": sp.Rational(1, 2), "IXX": sp.Rational(1, 2),
                  "IYY": sp.Rational(1, 2), "IZZ": sp.Rational(1, 2)})
    global_xy = dense({"".join(word): 1 for word in product("XY", repeat=3)})
    check(is_zero(w.subs(t, s) - s**3*swap*global_xy),
          "equal-coupling branch is the stated SWAP carrier")

    controls = [((3, 6, -2), 1), ((-1, 2, 2), 1), ((2, 2, 7), 1),
                ((2, 2, -1), 1), ((-2, -2, 1), 1),
                ((3, 6, -3), 0), ((2, -2, 3), 0), ((1, 2, 3), 0)]
    for point, expected in controls:
        sub = dict(zip((s, t, r), point))
        full = sp.Matrix.vstack(field_map, sp.Matrix.hstack(
            *(commutator(hb.subs(sub), q).reshape(64, 1) for q in lit)))
        check(32-full.rank() == expected, f"full dense common kernel at {point} has dimension {expected}")

    point = {s: 3, t: 6, r: -2}
    w0 = w.subs(point)
    check(not is_zero(commutator(hb.subs({s: 3, t: 6, r: -3}), w0)),
          "moving the coupling off both branches breaks the same commutator")
    check(not is_zero(commutator(dense({"IXI": 1, "IIY": 2}), w0)),
          "unequal fields break the candidate through the same field equation")
    check(not is_zero(commutator(hb.subs(point), w0-dense({"YYX": 3**2*6}))),
          "removing a carrier term breaks the same bond equation")
    # Quantifier control: h=0 may pair where no common carrier for all h exists.
    check(is_zero(commutator(hb.subs({s: 1, t: 2, r: 3}), dense({"XXX": 1})))
          and is_zero(jump*dense({"XXX": 1}) + dense({"XXX": 1})*jump),
          "off-locus h=0 has a carrier: common-carrier absence is not fixed-field absence")
    print("triangle common-carrier gate: ALL PASS")


if __name__ == "__main__":
    main()
