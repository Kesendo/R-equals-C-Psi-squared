"""Tests for the s2(C) floor and for the counted (f_leak, hop) curve.

WHY THIS FILE EXISTS. Section 7 registers the floor two ways that cannot
both hold: as the MAXIMUM over the f_leak and hop brackets, and subject
to P(s2C < floor | H1) <= 0.01. Measured 2026-08-18, the maximum breaches
its own reachability at three of the four corners (0.046-0.074 against
0.01).

WHERE THAT WENT, since these tests outlived two answers to it. The first
was to register the floor as an affine FUNCTION of the arm's own s2(U),
and the machinery for it is still exercised below. It is no longer on
the freeze path, for two reasons found after it was built: the function
is provably inert wherever it is reachable (the lemma below), and the
repo already owned the rule for this shape, namely that a LOWER cut
freezes at the LOWER envelope (RECORD_PARITY round 28, which reached it
after a cut frozen at the worst basis VOIDed 269 of 300 runs on the best
admitted device). The floor is therefore a scalar at the lower envelope,
and the function machinery is kept as measurement rather than deleted.

The second answer is the curve at the end of this file: the same Pauli
count section 10 uses to derive f_leak = 8/15 ALSO fixes hop, so the two
bracket axes are one parameter and the flown product box has corners no
channel model can produce, including the corner where the floor's
existence condition fails.

WHAT THESE TESTS PIN, and each one is a trap already walked into:

  * the affine value itself, and that b = 0 reproduces a plain scalar,
    so the scalar form is the degenerate member and the two are
    comparable;
  * the reachability criterion as a CONFIDENCE BOUND rather than a point
    estimate, because a rate of 0.01 measured on 500 draws has a standard
    error of 0.0044 and "the observed rate met the ceiling" is not the
    same statement as "the rate is below the ceiling";
  * that selection happens on the FIT sample and validation on a
    HELD-OUT one. The greedy optimum measured on the full sample reads
    0.008 in-sample and 0.012 held out: its own registered ceiling,
    missed. Section 8a registers held-out splits for exactly this and
    the harness's first version did not honour it for theta_D.
  * that the selector REFUSES rather than returning its least-bad
    candidate when nothing is admissible. A freeze that quietly degrades
    is the failure mode the whole pre-registration exists against.
"""

import importlib.util as _ilu
import math
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).parent


def _load(name, path):
    spec = _ilu.spec_from_file_location(name, path)
    mod = _ilu.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


cbr = _load("cbr_floor_test", HERE / "corner_beat_refreeze.py")


# --------------------------------------------------------- the value
def test_floor_value_is_affine():
    assert cbr.floor_value(0.0, (0.002, 0.5)) == pytest.approx(0.002)
    assert cbr.floor_value(0.001, (0.002, 0.5)) == pytest.approx(0.0025)
    assert cbr.floor_value(0.004, (0.0, 1.67)) == pytest.approx(0.00668)


def test_b_zero_reproduces_the_scalar_registration():
    """The old form is the b = 0 member: the change is a generalisation,
    not a replacement, and a reader must be able to check that."""
    for s2u in (0.0, 1e-4, 2e-3, 1.0):
        assert cbr.floor_value(s2u, (0.00389, 0.0)) == 0.00389


# ------------------------------------------- the confidence criterion
def test_upper_bound_zero_successes_matches_the_closed_form():
    """P(X = 0) = (1-p)^n = alpha has the exact root p = 1 - alpha^(1/n),
    which is the one case needing no special-function library."""
    for n in (50, 500, 1000):
        exact = 1.0 - 0.05 ** (1.0 / n)
        assert cbr.binom_upper_bound(0, n, 0.05) == pytest.approx(
            exact, rel=1e-9)


def test_upper_bound_all_successes_is_one():
    assert cbr.binom_upper_bound(7, 7, 0.05) == pytest.approx(1.0)


def test_upper_bound_is_monotone_in_the_count():
    prev = -1.0
    for k in range(0, 12):
        u = cbr.binom_upper_bound(k, 500, 0.05)
        assert u > prev
        prev = u


def test_upper_bound_inverts_the_binomial_tail():
    """From below: at the returned p, P(X <= k) must equal alpha."""
    for k, n in ((0, 500), (1, 500), (3, 500), (5, 200)):
        p = cbr.binom_upper_bound(k, n, 0.05)
        tail = sum(math.comb(n, i) * p ** i * (1 - p) ** (n - i)
                   for i in range(k + 1))
        assert tail == pytest.approx(0.05, abs=1e-9)


def test_reachability_admits_at_most_one_breach_in_five_hundred():
    """The criterion's practical content at the registered n, stated as a
    number so a change to alpha or the ceiling shows up here."""
    assert cbr.reach_admissible(0, 500)
    assert cbr.reach_admissible(1, 500)
    assert not cbr.reach_admissible(2, 500)
    assert not cbr.reach_admissible(5, 500)


def test_the_ceiling_needs_at_least_299_draws_to_be_demonstrable():
    """A DESIGN CONSTRAINT, not a detail: below this sample size even a
    flawless campaign cannot put the bound under the ceiling, so the
    selector refuses everything and the refusal is about the sample, not
    about the floor. Found by a 6-rep smoke run reporting "no admissible
    member" when the honest answer was "ask again with more draws"."""
    assert cbr.min_n_for_ceiling(0.01, 0.05) == 299
    assert not cbr.reach_admissible(0, 298)
    assert cbr.reach_admissible(0, 299)
    assert cbr.min_n_for_ceiling() == 299          # the registered pair


def test_a_point_estimate_at_the_ceiling_is_not_admissible():
    """5/500 = 0.010 meets the ceiling as a point estimate and fails as a
    bound. This is the whole reason the criterion is a bound."""
    assert 5 / 500 <= cbr.FLOOR_REACH_CEILING
    assert not cbr.reach_admissible(5, 500)


# ------------------------------------------------------- trip counts
def test_trip_counts_compare_per_rep():
    s2c = [0.005, 0.001, 0.004]
    s2u = [0.002, 0.000, 0.001]          # floors: 0.003, 0.002, 0.0025
    k, n = cbr.floor_trip_counts(s2c, s2u, (0.002, 0.5))
    assert (k, n) == (1, 3)              # only the middle rep trips


def test_trip_counts_drop_dead_reps_from_the_denominator():
    """A NaN rep is a VOID, not a non-trip: _rate's own rule, and the
    floor must not silently disagree with it."""
    s2c = [0.005, float("nan"), 0.001]
    s2u = [0.002, 0.001, 0.000]
    k, n = cbr.floor_trip_counts(s2c, s2u, (0.002, 0.5))
    assert (k, n) == (1, 2)


def test_trip_counts_reject_mismatched_columns():
    with pytest.raises(ValueError):
        cbr.floor_trip_counts([0.001, 0.002], [0.001], (0.002, 0.5))


# -------------------------------------------------- the conjunction
#
# EVERY TRIPLE BELOW SATISFIES d = s2C - s2U, and that is the point.
# The first version of these tests set d independently of the two
# variances, so it asserted rows like (d = 0.01, s2C = 0.000,
# s2U = 0.0), where d is really 0.0. They passed, they demonstrated the
# second conjunct biting, and the demonstration was of a rep the
# harness cannot produce: one_rep DEFINES d as s2C - s2U. That is what
# hid the inertness lemma below for a whole build.
def _triple(s2c, s2u):
    """Columns on the manifold the harness actually lives on."""
    return [c - u for c, u in zip(s2c, s2u)], list(s2c), list(s2u)


def test_conjunction_needs_both_conjuncts():
    """Round 9's rule: d > theta_D AND s2C >= floor, walked over the four
    truth-table rows on CONSISTENT triples. The floor must sit above
    theta_D for the second conjunct to be able to bite at all, which is
    the lemma two tests down."""
    d, s2c, s2u = _triple([0.010, 0.004, 0.003, 0.001],
                          [0.000, 0.000, 0.000, 0.000])
    #                d =   0.010  0.004  0.003  0.001
    #  conjunct 1 (>0.0035): pass   pass   fail   fail
    #  conjunct 2 (>=0.006): pass   fail   fail   fail
    r = cbr.conjunction_rate(d, s2c, s2u, 0.0035, (0.006, 0.0))
    assert r == pytest.approx(0.25)


def test_conjunction_reads_the_floor_at_the_same_rep():
    """The floor uses the rep's OWN s2(U). Two reps with identical d and
    identical s2(C) - impossible, so instead: two reps with the same
    s2(C), whose different s2(U) gives them different d AND different
    floors, and only one survives both conjuncts."""
    d, s2c, s2u = _triple([0.008, 0.008], [0.000, 0.003])
    #  d = 0.008 and 0.005; floors at b = 0.5: 0.002 and 0.0035
    r = cbr.conjunction_rate(d, s2c, s2u, 0.006, (0.002, 0.5))
    assert r == pytest.approx(0.5)   # only the first clears d > 0.006


def test_the_floor_is_inert_whenever_a_le_theta_D_and_b_le_1():
    """THE LEMMA, pinned. On the manifold d = s2C - s2U with s2U >= 0:
    passing conjunct 1 and failing conjunct 2 needs
        s2U + theta_D < s2C < a + b*s2U,
    hence a - theta_D > (1 - b)*s2U, which for b <= 1 forces a > theta_D.
    So below theta_D the floor cannot change one verdict, and the
    selector must never be allowed to report such a member as a guard.
    Checked on random consistent triples rather than by re-deriving the
    inequality in the test."""
    rng = np.random.default_rng(11)
    theta_D = 0.00291
    for a, b in ((0.0020, 0.5), (0.00291, 1.0), (0.0, 0.0),
                 (0.00289, 0.75), (0.0015, 0.25)):
        s2u = np.abs(rng.normal(0.001, 0.001, 4000))
        s2c = np.abs(rng.normal(0.004, 0.003, 4000))
        d = list(s2c - s2u)
        both = cbr.conjunction_rate(d, list(s2c), list(s2u), theta_D, (a, b))
        d_only = sum(1 for x in d if x > theta_D) / len(d)
        assert both == pytest.approx(d_only), (a, b)


def test_the_floor_can_bite_once_a_exceeds_theta_D():
    """The converse, so the lemma is a boundary and not a blanket: at
    a > theta_D the second conjunct removes reps the first one passed."""
    rng = np.random.default_rng(12)
    theta_D = 0.00291
    s2u = np.abs(rng.normal(0.001, 0.001, 4000))
    s2c = np.abs(rng.normal(0.004, 0.003, 4000))
    d = list(s2c - s2u)
    both = cbr.conjunction_rate(d, list(s2c), list(s2u), theta_D,
                                (0.0045, 0.0))
    d_only = sum(1 for x in d if x > theta_D) / len(d)
    assert both < d_only


def test_conjunction_drops_reps_dead_in_any_column():
    d, s2c, s2u = _triple([0.010, 0.010, 0.010], [0.0, 0.0, 0.0])
    s2c[1] = float("nan")
    d[2] = float("nan")
    r = cbr.conjunction_rate(d, s2c, s2u, 0.005, (0.002, 0.0))
    assert r == pytest.approx(1.0)    # one live rep, and it passes


def test_conjunction_is_none_when_nothing_is_live():
    assert cbr.conjunction_rate([float("nan")], [1.0], [1.0],
                                0.0, (0.0, 0.0)) is None


def test_conjunction_rejects_mismatched_columns():
    with pytest.raises(ValueError):
        cbr.conjunction_rate([0.1], [0.1, 0.2], [0.1], 0.0, (0.0, 0.0))


# --------------------------------------------------------- selection
def _corner(s2c, s2u):
    return {"s2C": list(s2c), "s2U": list(s2u)}


def _synthetic(n=500, sep=0.004, level=0.001, seed=7):
    """Two corners with different noise LEVELS and a common signal, the
    shape the real campaign has: H1 sits `sep` above the H0 cloud."""
    rng = np.random.default_rng(seed)
    h0, h1 = {}, {}
    for i, lv in enumerate((level, level / 8)):
        u0 = np.abs(rng.normal(lv, lv / 2, n))
        u1 = np.abs(rng.normal(lv, lv / 2, n))
        h0[f"c{i}"] = _corner(np.abs(rng.normal(lv, lv / 2, n)), u0)
        h1[f"c{i}"] = _corner(sep + np.abs(rng.normal(lv, lv / 2, n)), u1)
    return h0, h1


def test_selection_returns_an_admissible_member():
    h0, h1 = _synthetic()
    sel = cbr.select_floor_coefficients(h1, h0)
    assert sel is not None
    a, b = sel["coefficients"]
    assert a in cbr.FLOOR_A_GRID and b in cbr.FLOOR_B_GRID
    for tag in h1:
        k, n = cbr.floor_trip_counts(h1[tag]["s2C"], h1[tag]["s2U"], (a, b))
        assert cbr.reach_admissible(k, n)


def test_selection_maximises_the_weakest_h0_trip_rate():
    """The floor's registered PURPOSE (section 7 (i)): under a true null
    it should trip. Among admissible members the selector must take the
    strongest guard, not the first one it meets."""
    h0, h1 = _synthetic()
    sel = cbr.select_floor_coefficients(h1, h0)
    chosen = sel["weakest_h0_trip_fit"]
    for b in cbr.FLOOR_B_GRID:
        for a in cbr.FLOOR_A_GRID:
            ok = all(cbr.reach_admissible(
                *cbr.floor_trip_counts(h1[t]["s2C"], h1[t]["s2U"], (a, b)))
                for t in h1)
            if not ok:
                continue
            weakest = 1.0
            for t in h0:
                k, n = cbr.floor_trip_counts(
                    h0[t]["s2C"], h0[t]["s2U"], (a, b))
                weakest = min(weakest, k / n)
            assert weakest <= chosen + 1e-12


def test_selection_is_deterministic():
    h0, h1 = _synthetic()
    assert (cbr.select_floor_coefficients(h1, h0)["coefficients"]
            == cbr.select_floor_coefficients(h1, h0)["coefficients"])


def test_the_degenerate_member_is_always_admissible():
    """WRITTEN AFTER GETTING IT WRONG. The first version of the two
    tests below assumed the selector could be starved into returning
    None on realistic input. It cannot: s2(C) is a VARIANCE, so it is
    never negative, so the floor (a, b) = (0, 0) trips on nothing and
    can never breach reachability. That member is in the grid on
    purpose (it is "the floor does nothing", route (b)'s endpoint), and
    its permanent admissibility is why the operative check is the
    reported guard rather than a refusal."""
    dead = {"c0": _corner([0.0] * 500, [0.0] * 500)}
    k, n = cbr.floor_trip_counts(dead["c0"]["s2C"], dead["c0"]["s2U"],
                                 (0.0, 0.0))
    assert (k, n) == (0, 500)
    assert cbr.reach_admissible(k, n)
    assert 0.0 in cbr.FLOOR_A_GRID and 0.0 in cbr.FLOOR_B_GRID


def test_selection_refuses_when_no_grid_member_is_admissible():
    """The refusal path itself, reached with a grid that excludes the
    degenerate member: every candidate floor sits above every H1 draw."""
    h1 = {"c0": _corner([0.001] * 500, [0.0] * 500)}
    h0 = {"c0": _corner([0.001] * 500, [0.0] * 500)}
    assert cbr.select_floor_coefficients(
        h1, h0, b_grid=(0.0,), a_grid=(0.01,)) is None


def test_a_corner_that_cannot_reach_vetoes_every_working_member():
    """Registering ONE function admissible at EVERY corner is what
    replaces the max-over-bracket, and it settles section 7's own
    mismatch (the producer sentence says both brackets, the
    reachability sentence says one). A corner whose H1 draws sit at
    zero admits nothing but the degenerate member, so the selection
    must COLLAPSE to it and say so, rather than quietly keeping a
    candidate that some other corner liked."""
    h0, h1 = _synthetic()
    sel_ok = cbr.select_floor_coefficients(h1, h0)
    assert sel_ok is not None and not sel_ok["degenerate"]
    assert sel_ok["coefficients"][0] > 0.0

    poisoned_h1 = dict(h1, poison=_corner([0.0] * 500, [0.0] * 500))
    poisoned_h0 = dict(h0, poison=_corner([0.0] * 500, [0.0] * 500))
    sel_bad = cbr.select_floor_coefficients(poisoned_h1, poisoned_h0)
    assert sel_bad is not None
    assert sel_bad["coefficients"][0] == 0.0
    assert sel_bad["weakest_h0_trip_fit"] == 0.0
    assert sel_bad["degenerate"]


# -------------------------------------------------------- validation
def test_validation_reads_the_held_out_sample_and_can_fail():
    h0, h1 = _synthetic(seed=7)
    sel = cbr.select_floor_coefficients(h1, h0)
    coef = sel["coefficients"]
    good = cbr.validate_floor_coefficients(coef, h1)
    assert good["passes"]
    dead = {"c0": _corner([0.0] * 500, [0.0] * 500)}
    assert not cbr.validate_floor_coefficients(coef, dead)["passes"]


def test_validation_reports_per_corner_counts():
    h0, h1 = _synthetic()
    coef = cbr.select_floor_coefficients(h1, h0)["coefficients"]
    rep = cbr.validate_floor_coefficients(coef, h1)
    assert set(rep["corners"]) == set(h1)
    for tag, c in rep["corners"].items():
        assert c["n"] == 500
        assert 0 <= c["trips"] <= 500
        assert c["upper_bound"] >= c["rate"]


# ------------------------------------------------- the counted curve
def test_curve_endpoints_are_the_documents_own_two_counts():
    """Section 10 states BOTH numbers: 8/15 on a differing pair and
    12/15 = 0.80 on a non-differing one. The curve must reproduce them
    at its endpoints, or it is not the same object the document uses."""
    leak1, hop1 = cbr.leak_hop_at(1.0)
    leak0, hop0 = cbr.leak_hop_at(0.0)
    assert leak1 == pytest.approx(8 / 15)
    assert hop1 == pytest.approx(4 / 7)
    assert leak0 == pytest.approx(12 / 15)
    assert hop0 == pytest.approx(0.0)


def test_curve_matches_an_independent_pauli_enumeration():
    """FROM BELOW. The closed form is checked against a fresh
    enumeration written here, not against the derivation it came from:
    apply every two-qubit Pauli on every bond to every one-magnon basis
    state of a six-site chain and classify the image."""
    N = 6
    I2 = np.eye(2, dtype=complex)
    X = np.array([[0, 1], [1, 0]], dtype=complex)
    Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
    Z = np.array([[1, 0], [0, -1]], dtype=complex)
    P = {"I": I2, "X": X, "Y": Y, "Z": Z}
    names = [a + b for a in "IXYZ" for b in "IXYZ" if a + b != "II"]

    def op(name, j):
        ops = [I2] * N
        ops[j], ops[j + 1] = P[name[0]], P[name[1]]
        out = np.array([[1.0 + 0j]])
        for o in ops:
            out = np.kron(out, o)
        return out

    site_of = {1 << (N - 1 - l): l for l in range(N)}
    on = {"leave": 0, "move": 0, "stay": 0}
    off = {"leave": 0, "move": 0, "stay": 0}
    for name in names:
        for j in range(N - 1):
            M = op(name, j)
            for l in range(N):
                v = np.zeros(2 ** N, dtype=complex)
                v[1 << (N - 1 - l)] = 1.0
                w = M @ v
                idx = np.flatnonzero(np.abs(w) > 1e-12)
                assert idx.size == 1
                k = int(idx[0])
                if bin(k).count("1") != 1:
                    r = "leave"
                elif site_of[k] != l:
                    r = "move"
                else:
                    r = "stay"
                (on if l in (j, j + 1) else off)[r] += 1

    tot_on, tot_off = sum(on.values()), sum(off.values())
    for q in (0.0, 0.25, 0.5, 0.75, 1.0):
        leak = (q * on["leave"] / tot_on
                + (1 - q) * off["leave"] / tot_off)
        move = q * on["move"] / tot_on
        ins = (q * (on["move"] + on["stay"]) / tot_on
               + (1 - q) * (off["move"] + off["stay"]) / tot_off)
        expect_leak, expect_hop = cbr.leak_hop_at(q)
        assert leak == pytest.approx(expect_leak), q
        assert (move / ins if ins else 0.0) == pytest.approx(expect_hop), q


def test_the_curve_is_monotone_in_both_coordinates():
    prev_leak, prev_hop = None, None
    for i in range(11):
        q = i / 10
        leak, hop = cbr.leak_hop_at(q)
        if prev_leak is not None:
            assert leak < prev_leak      # more on-bond error -> less leaks
            assert hop > prev_hop        # ... and more of it moves
        prev_leak, prev_hop = leak, hop


def test_the_flown_boxs_anti_correlated_corner_is_off_the_curve():
    """The corner that makes the floor's existence condition fail,
    (f_leak = 8/15, hop = 0), is not reachable at ANY q: at that leak
    the count says hop = 4/7, and hop = 0 needs leak = 12/15. This is
    the whole reason the sweep set changes."""
    for i in range(1001):
        leak, hop = cbr.leak_hop_at(i / 1000)
        assert not (abs(leak - 8 / 15) < 0.02 and hop < 0.4)


def test_curve_rejects_q_outside_the_unit_interval():
    for bad in (-0.001, 1.001, 2.0):
        with pytest.raises(ValueError):
            cbr.leak_hop_at(bad)
