"""Verification for tau_max_spectral_gap.py — proven invariants + the gamma-clock result.
Run: python simulations/tau_max_spectral_gap_verify.py
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import tau_max_spectral_gap as t


def approx(a, b, tol=1e-4):
    return abs(a - b) <= tol


def test_floor_is_2gamma():
    gamma = 0.05
    for N in (2, 3, 4):
        gap = t.spectral_gap(t.build_L(t.heisenberg_chain_H(N, J=1.0), gamma, N))
        assert approx(gap, 2 * gamma), f"N={N}: gap={gap}, expected {2*gamma}"
    print("PASS test_floor_is_2gamma")


def test_palindrome_symmetry():
    gamma, N = 0.05, 3
    r = t.all_rates(t.build_L(t.heisenberg_chain_H(N, J=1.0), gamma, N))
    twoSigma = 2 * N * gamma
    for i in range(len(r)):
        assert approx(r[i] + r[-1 - i], twoSigma, tol=1e-3), \
            f"palindrome break at i={i}: {r[i]}+{r[-1-i]} != {twoSigma}"
    print("PASS test_palindrome_symmetry")


def test_gap_J_independent():
    gamma, N = 0.05, 3
    for J in (0.5, 1.0, 2.0, 4.0):
        gap = t.spectral_gap(t.build_L(t.heisenberg_chain_H(N, J=J), gamma, N))
        assert approx(gap, 2 * gamma, tol=1e-6), f"J={J}: gap={gap} != 2gamma; gap depends on J!"
    print("PASS test_gap_J_independent")


def test_gap_gamma_linear():
    N, J = 3, 1.0
    for gamma in (0.02, 0.05, 0.10, 0.20):
        gap = t.spectral_gap(t.build_L(t.heisenberg_chain_H(N, J), gamma, N))
        assert approx(gap, 2 * gamma, tol=1e-6), f"gamma={gamma}: gap={gap} != {2*gamma}"
    print("PASS test_gap_gamma_linear")


def test_gap_vanishes_at_gamma0():
    # gamma is the timekeeper: with no dephasing the gap is 0 and the clock stops.
    gap = t.spectral_gap(t.build_L(t.heisenberg_chain_H(3, 1.0), 0.0, 3))
    assert approx(gap, 0.0, tol=1e-6), f"gap at gamma=0 is {gap}, expected 0"
    assert t.clock_tau(gap) == float("inf"), "clock should be infinite at gamma=0"
    print("PASS test_gap_vanishes_at_gamma0 (clock stops)")


def test_clock_is_inverse_gap_N_and_J_independent():
    # The clock tau = 1/lambda2 = 1/(2 gamma): depends only on gamma.
    gamma = 0.05
    expect = 1.0 / (2 * gamma)
    for N in (2, 3, 4, 5):
        c = t.clock_tau(t.spectral_gap(t.build_L(t.heisenberg_chain_H(N, 1.0), gamma, N)))
        assert approx(c, expect, tol=1e-3), f"N={N}: clock={c} != {expect}"
    for J in (0.5, 1.0, 2.0, 4.0):
        c = t.clock_tau(t.spectral_gap(t.build_L(t.heisenberg_chain_H(3, J), gamma, 3)))
        assert approx(c, expect, tol=1e-3), f"J={J}: clock={c} != {expect}"
    print("PASS test_clock_is_inverse_gap_N_and_J_independent")


def test_clock_gamma_power_is_minus_one():
    gs = (0.02, 0.05, 0.10, 0.20)
    clocks = [t.clock_tau(t.spectral_gap(t.build_L(t.heisenberg_chain_H(3, 1.0), g, 3))) for g in gs]
    slope = t.loglog_slope(gs, clocks)
    assert approx(slope, -1.0, tol=0.02), f"clock gamma-power={slope}, expected -1"
    print(f"PASS test_clock_gamma_power_is_minus_one (slope={slope:+.3f})")


def test_formula_powers_are_wrong():
    # The extracted formula gives J^-1 (spurious) and gamma^-1/2 (wrong),
    # instead of the true J^0 and gamma^-1.
    gamma = 0.05
    Js = (0.5, 1.0, 2.0, 4.0)
    f_J = [t.formula_tau(t.spectral_gap(t.build_L(t.heisenberg_chain_H(3, J), gamma, 3)), J) for J in Js]
    assert approx(t.loglog_slope(Js, f_J), -1.0, tol=0.02), "formula should scale as 1/J"
    gs = (0.02, 0.05, 0.10, 0.20)
    f_g = [t.formula_tau(t.spectral_gap(t.build_L(t.heisenberg_chain_H(3, 1.0), g, 3)), 1.0) for g in gs]
    assert approx(t.loglog_slope(gs, f_g), -0.5, tol=0.02), "formula should scale as gamma^-1/2"
    print("PASS test_formula_powers_are_wrong (J^-1, gamma^-1/2)")


def test_verdict_rejected():
    a = t.run_analysis()
    assert a["verdict"] == "rejected", f"verdict={a['verdict']}, expected rejected"
    print("PASS test_verdict_rejected")


if __name__ == "__main__":
    test_floor_is_2gamma()
    test_palindrome_symmetry()
    test_gap_J_independent()
    test_gap_gamma_linear()
    test_gap_vanishes_at_gamma0()
    test_clock_is_inverse_gap_N_and_J_independent()
    test_clock_gamma_power_is_minus_one()
    test_formula_powers_are_wrong()
    test_verdict_rejected()
    print("ALL CHECKS PASSED")
