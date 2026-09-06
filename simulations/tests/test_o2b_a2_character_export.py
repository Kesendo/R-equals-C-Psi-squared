"""Exact degree-one subresultant contract for the canonical N=5 sectors."""
from pathlib import Path
import sys
import json
import copy
from fractions import Fraction
from decimal import Decimal, localcontext

import numpy as np
import pytest
import sympy as sp

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "simulations"))
import o2b_gcd_certificate as cert


@pytest.fixture(scope="module")
def exported(tmp_path_factory):
    folder = tmp_path_factory.mktemp("a2-export")
    first, second = folder / "first.json", folder / "second.json"
    cert.write_a2_json(first)
    cert.write_a2_json(second)
    assert first.read_bytes() == second.read_bytes()
    assert first.read_bytes() == (ROOT / "simulations/results/route_b_a2_n5.json").read_bytes()
    assert first.read_bytes().endswith(b"\n")
    return json.loads(first.read_text(encoding="utf-8"))


def _box(value):
    def rational(r):
        assert set(r) == {"numerator", "denominator"}
        assert all(isinstance(s, str) for s in r.values())
        q = Fraction(int(r["numerator"]), int(r["denominator"]))
        assert str(q.numerator) == r["numerator"]
        assert str(q.denominator) == r["denominator"]
        assert cert.rational_from_json(r) == q
        assert cert.rational_json(q) == r
        return q
    assert set(value) == {"reLo", "reHi", "imLo", "imHi"}
    result = tuple(tuple(rational(value[key]) for key in pair)
                   for pair in (("reLo", "reHi"), ("imLo", "imHi")))
    assert all(lo <= hi for lo, hi in result)
    return result


def _mul(a, b):
    products = [x*y for x in a for y in b]
    return min(products), max(products)


def _square(box):
    def square(a):
        return (0 if a[0] <= 0 <= a[1] else min(x*x for x in a), max(x*x for x in a))
    x, y = box
    xx, yy, xy = square(x), square(y), _mul(x, y)
    return ((xx[0]-yy[1], xx[1]-yy[0]), (2*xy[0], 2*xy[1]))


def _inside(inner, outer):
    return all(c <= a <= b <= d for (a,b),(c,d) in zip(inner,outer))


def _validate_export(data):
    assert data["schemaVersion"] == 1 and data["n"] == 5
    assert data["conventions"] == {"w": "qUnitHop^2", "qPhysicalCSharp": "qUnitHop/2"}
    assert data["exactRankCertificates"] == []
    assert all(set(sector) == {"parity", "a2Roots"} for sector in data["sectors"])
    roots = [r for sector in data["sectors"] for r in sector["a2Roots"]]
    assert len(roots) == 29
    assert [sum(r["rootKind"] == kind for r in roots) for kind in
            ("negativeReal", "positiveReal", "nonreal")] == [12, 1, 16]
    ids = []
    for root in roots:
        ids.append(root["id"])
        wbox, lbox = _box(root["wBox"]), _box(root["lambdaBox"])
        lambda_seed = tuple(Fraction(root["lambdaSeed"][axis]) for axis in ("real", "imag"))
        assert all(lo <= seed <= hi for seed,(lo,hi) in zip(lambda_seed,lbox))
        assert len(root["qLoci"]) == 2
        qboxes = []
        for locus in root["qLoci"]:
            ids.append(locus["id"])
            qb, pb = _box(locus["qUnitHopBox"]), _box(locus["qPhysicalCSharpBox"])
            assert pb == tuple(tuple(x/2 for x in a) for a in qb)
            assert _inside(_square(qb), wbox)
            assert all(b-a < Fraction(1, 2**160) for a,b in qb)
            qboxes.append(qb)
        assert qboxes[1] == tuple((-b, -a) for a,b in qboxes[0])
        for seed in [root["wSeed"], root["lambdaSeed"]] + [q[k] for q in root["qLoci"] for k in ("qUnitHopSeed", "qPhysicalCSharpSeed")]:
            assert set(seed) == {"real", "imag"}
            assert all(isinstance(x,str) for x in seed.values())
    assert len(set(ids)) == len(ids)
    for sector in data["sectors"]:
        parity = sector["parity"]
        assert [r["id"] for r in sector["a2Roots"]] == [f"N5-{parity}-A2-W-{i:03d}" for i in range(len(sector["a2Roots"]))]
        for r in sector["a2Roots"]:
            assert r["parity"] == parity
            assert [q["id"] for q in r["qLoci"]] == [r["id"]+"-Q-plus",r["id"]+"-Q-minus"]
        all_q = [(r,_box(q["qUnitHopBox"])) for r in sector["a2Roots"] for q in r["qLoci"]]
        for i,(r,qb) in enumerate(all_q):
            assert [s["id"] for s in sector["a2Roots"] if _inside(_square(qb),_box(s["wBox"]))] == [r["id"]]
            assert all(any(b<c or d<a for (a,b),(c,d) in zip(qb,other)) for _,other in all_q[i+1:])
    return roots


def test_export_schema_exact_boxes_and_inventory(exported):
    _validate_export(exported)


def test_exact_json_handles_sympy_domain_rationals_and_rejects_floats():
    assert cert.rational_json(sp.QQ(2,3)) == {"numerator": "2", "denominator": "3"}
    box = ((Fraction(1,3),Fraction(2,3)),(Fraction(0),Fraction(0)))
    encoded = cert.box_json(box)
    assert set(encoded) == {"reLo", "reHi", "imLo", "imHi"}
    assert _box(encoded) == box
    assert cert.box_contains_real(encoded,Fraction(1,2))
    with pytest.raises(TypeError):
        cert.rational_json(0.5)


def test_outward_dyadic_box_preserves_non_dyadic_endpoints():
    source = ((Fraction(1,3),Fraction(2,3)),(Fraction(-2,3),Fraction(-1,3)))
    rounded = cert._dyadic_outer_box(source,bits=8)
    assert rounded == ((Fraction(85,256),Fraction(171,256)),(Fraction(-171,256),Fraction(-85,256)))
    assert _inside(source,rounded)
    assert not _inside(source,((Fraction(86,256),Fraction(170,256)),rounded[1]))


def test_proposed_q_box_requires_exact_root_and_unambiguous_square_mapping():
    Q = sp.Poly(cert.Q2**2+1,cert.Q2)
    coarse = ((Fraction(-1),Fraction(1)),(Fraction(1,2),Fraction(3,2)))
    empty = ((Fraction(1,4),Fraction(1,3)),(Fraction(3,4),Fraction(5,4)))
    wbox = ((Fraction(-2),Fraction(0)),(Fraction(-1),Fraction(1)))
    with pytest.raises(ValueError,match="root count"):
        cert._certify_q_box(Q,empty,coarse,[wbox])
    valid = ((Fraction(-1,100),Fraction(1,100)),(Fraction(99,100),Fraction(101,100)))
    assert cert._certify_q_box(Q,valid,coarse,[wbox]) == 0
    with pytest.raises(ValueError,match="unique w"):
        cert._certify_q_box(Q,valid,coarse,[wbox,wbox])


@pytest.mark.parametrize("mutation", ["drop-lift", "half-lambda", "corrupt-endpoint"])
def test_export_mutations_are_detected(exported, mutation):
    changed = copy.deepcopy(exported)
    root = changed["sectors"][0]["a2Roots"][0]
    if mutation == "drop-lift":
        root["qLoci"].pop()
    elif mutation == "corrupt-endpoint":
        root["qLoci"][0]["qUnitHopBox"]["reLo"]["numerator"] = "99999999999999999999999999"
    else:
        with localcontext() as context:
            context.prec = 60
            root["lambdaSeed"]["real"] = str(Decimal(root["lambdaSeed"]["real"])/2)
    with pytest.raises((AssertionError, ValueError)):
        _validate_export(changed)
        _validate_seeds(changed)


def _validate_seeds(data):
    from mpmath import mp
    with mp.workdps(85):
        def z(seed):
            return mp.mpc(seed["real"], seed["imag"])
        def evaluate(poly, values):
            terms = [int(c)*mp.fprod(v**p for v,p in zip(values,m)) for m,c in poly.terms()]
            return abs(mp.fsum(terms)), mp.fsum(abs(t) for t in terms)
        for sector in data["sectors"]:
            F = cert.exact_n5_residual_polynomial(sector["parity"])
            A2 = cert.exact_n5_a2(sector["parity"])
            _, a, _ = cert.first_subresultant_linear_exact(F)
            for root in sector["a2Roots"]:
                w, l = z(root["wSeed"]), z(root["lambdaSeed"])
                assert evaluate(a,[w])[0] > 0
                for poly, vals in ((A2,[w]),(F,[l,w]),(F.diff(cert.lam),[l,w])):
                    residual, scale = evaluate(poly,vals)
                    # 50 significant decimal seed digits: Horner/monomial conditioning
                    # is bounded here by the absolute term sum, with 5 guard digits.
                    assert residual <= mp.mpf("1e-45")*max(1,scale)
                for locus in root["qLoci"]:
                    q, physical = z(locus["qUnitHopSeed"]), z(locus["qPhysicalCSharpSeed"])
                    assert abs(q*q-w) <= mp.mpf("1e-45")*max(1,abs(w))
                    assert abs(physical*2-q) <= mp.mpf("1e-45")*max(1,abs(q))


def test_export_decimal_seeds_satisfy_independent_polynomials(exported):
    _validate_seeds(exported)


def test_transformed_polynomial_exact_root_counts(exported):
    def rational(x):
        return sp.Rational(x.numerator,x.denominator)
    for sector in exported["sectors"]:
        A2 = cert.exact_n5_a2(sector["parity"])
        Q = sp.Poly(A2.as_expr().subs(cert.Q2,cert.Q2**2),cert.Q2)
        Qi = sp.Poly(Q.as_expr().subs(cert.Q2,sp.I*cert.Q2).expand(),cert.Q2)
        assert Q.degree() == 2*len(sector["a2Roots"])
        assert sp.gcd(Q,Q.diff()).degree() == 0
        for root in sector["a2Roots"]:
            for locus in root["qLoci"]:
                x,y = _box(locus["qUnitHopBox"])
                if y == (0,0):
                    count = Q.count_roots(rational(x[0]),rational(x[1]))
                elif x == (0,0):
                    count = Qi.count_roots(rational(y[0]),rational(y[1]))
                else:
                    count = Q.count_roots(rational(x[0])+sp.I*rational(y[0]),
                                         rational(x[1])+sp.I*rational(y[1]))
                assert count == 1


def test_positive_anchor_and_physical_convention(exported):
    from mpmath import mp
    root = next(r for s in exported["sectors"] for r in s["a2Roots"] if r["rootKind"] == "positiveReal")
    assert root["parity"] == "O"
    assert cert.box_contains_real(root["wBox"],"5.10083102")
    assert cert.box_contains_real(root["lambdaBox"],"-4.79196037")
    with mp.workdps(60):
        assert abs(mp.mpf(root["wSeed"]["real"])-mp.mpf("5.10083102")) < mp.mpf("5e-9")
        assert abs(mp.mpf(root["lambdaSeed"]["real"])-mp.mpf("-4.79196037")) < mp.mpf("5e-9")
        q = mp.mpf(root["qLoci"][0]["qPhysicalCSharpSeed"]["real"])
        assert abs(abs(q)-mp.sqrt(mp.mpf(root["wSeed"]["real"]))/2) < mp.mpf("1e-10")


@pytest.mark.parametrize("sector", ["E", "O"])
def test_subresultant_uses_canonical_n5_inventory(sector):
    F = cert.exact_n5_residual_polynomial(sector)
    assert cert.exact_n5_residual_polynomial(sector) is F
    inventory = cert.exact_disc_inventory_n5(sector, F)
    assert cert.exact_n5_a2(sector) is inventory["A2"]
    assert F.gens == (cert.lam, cert.Q2) and F.domain == sp.ZZ
    assert F.degree(cert.lam) == cert.FRES_DEG[5][sector]
    assert inventory["A2"].degree() == cert.DISC_TABLE[5][sector]["A2"]


def test_exact_inventory_output_keeps_character_with_its_separate_gate(capsys):
    F = cert.exact_n5_residual_polynomial("O")
    # Reuse the cached exact pencil, but execute the inventory so its output is read.
    cert.exact_disc_inventory_n5.__wrapped__("O", F)
    output = capsys.readouterr().out
    assert "order-two A2 layer; algebraic pair unique; character gated separately" in output
    assert "diabolic class" not in output


@pytest.mark.parametrize("sector", ["E", "O"])
def test_subresultant_linear_exact_matches_independent_psc1(sector):
    F = cert.exact_n5_residual_polynomial(sector)
    s1, a, b = cert.first_subresultant_linear_exact(F)
    assert s1.gens == (cert.lam, cert.Q2) and s1.domain == sp.ZZ
    assert s1.degree(cert.lam) == 1
    assert sp.expand(s1.as_expr() - a.as_expr() * cert.lam - b.as_expr()) == 0
    assert a.gens == b.gens == (cert.Q2,)
    assert a.domain == b.domain == sp.ZZ
    assert a.LC() > 0 and s1.content() == 1
    assert sp.gcd(a, cert.exact_n5_a2(sector)).degree() == 0
    common = sp.gcd(a,b)
    assert sp.gcd(common,cert.exact_n5_a2(sector)).degree() == 0
    a_reduced,b_reduced = a.exquo(common),b.exquo(common)
    assert a*b_reduced-b*a_reduced == 0
    assert a*(b_reduced+1)-b*a_reduced != 0
    good = 0
    for p in cert.make_primes(8, below=2 ** 25 - 60):
        ap = cert.sympy_poly_desc_modp(a, p)
        if int(a.LC()) % p == 0:
            continue
        psc1 = cert.psc1_poly_modp(cert.poly_to_mat(F), p)
        if cert.is_zero(psc1):
            continue
        assert cert.poly_equal_up_to_nonzero_scalar(ap, psc1, p)
        assert not cert.poly_equal_up_to_nonzero_scalar(np.zeros_like(ap), psc1, p)
        good += 1
        if good == 3:
            break
    assert good >= 3


def test_subresultant_recovers_toy_double_root_and_rejects_changed_root():
    lam, w = cert.lam, cert.Q2
    F = sp.Poly((lam - (w + 1)) ** 2 * (lam + 2), lam, w, domain="ZZ")
    _, a, b = cert.first_subresultant_linear_exact(F)
    assert sp.cancel(-b.as_expr() / a.as_expr()) == w + 1
    changed = sp.Poly((lam - (w + 2)) ** 2 * (lam + 2), lam, w, domain="ZZ")
    _, changed_a, changed_b = cert.first_subresultant_linear_exact(changed)
    changed_root = sp.cancel(-changed_b.as_expr() / changed_a.as_expr())
    assert changed_root == w + 2
    assert changed_root != w + 1


def test_subresultant_modular_comparison_rejects_zero_and_wrong_polynomials():
    p = 101
    f = sp.Poly(3 * cert.Q2 ** 2 - 2, cert.Q2, domain="ZZ")
    assert np.array_equal(cert.sympy_poly_desc_modp(f, p), [3, 0, 99])
    compare = cert.poly_equal_up_to_nonzero_scalar
    assert compare([3, 0, 99], [15, 0, 91], p)
    assert compare([0, 3, 0, 99], [15, 0, 91], p)
    assert not compare([3, 0, 99], [15, 0, 92], p)
    for zero in ([], [0], [0, 0], [p, 0]):
        assert not compare(zero, [3, 0, 99], p)
        assert not compare([3, 0, 99], zero, p)
        assert not compare(zero, zero, p)
