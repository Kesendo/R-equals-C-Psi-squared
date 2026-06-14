"""Resolve the XXZ handover Delta*(N) descent: does it reach the Heisenberg point Delta=1 as N->inf?

Two methods, mutually validated:
  - finite-gamma sector method (reused from xxz_delta_star.py): Delta* where the (p,p)-block
    Lebensader rate crosses the band-edge floor 2*gamma, at Q=20. Cross-check ONLY.
  - gamma->0 reduction (the lever): Delta* where gap(R)=2, R the Z-coupled classical rate matrix
    among the half-filling XXZ eigenstates, built on the SECTOR H directly (feasible to N~14-16).

The descent is fit ENTIRELY in the gamma->0 regime; the finite-gamma points are a labeled
cross-check (the Q=20 drift grows with N and would bias the limit downward).
Spec: docs/superpowers/specs/2026-06-14-xxz-delta-star-descent-design.md
"""
import sys
sys.path.insert(0, "simulations")
import numpy as np
from numpy.linalg import eigh
from scipy.optimize import brentq, curve_fit
from xxz_delta_star import xxz_Hp, lebensader_rate, delta_star, full_slowest_rate

GAMMA = 0.05  # finite-gamma cross-check regime, Q=20
PHI = 2.0 * np.cos(np.pi / 5)  # 1.61803...

# finite-gamma (Q=20) Delta* -- CROSS-CHECK ONLY, never fit input
DSTAR_Q20 = {4: 1.61789, 5: 1.52530, 6: 1.381, 7: 1.325}


def check_finite_gamma_baseline():
    """assert #2a: the reused finite-gamma sector method reproduces the canonical Q=20 Delta*."""
    for N, expected in DSTAR_Q20.items():
        got = delta_star(N, gamma=GAMMA)
        assert abs(got - expected) < 1e-3, f"finite-gamma Delta*({N})={got} != {expected}"
    print("[2a] finite-gamma sector method reproduces canonical Q=20 Delta*(4..7).  OK")


def rate_matrix_R(N, Delta, J=1.0):
    """Gamma-free classical rate matrix among the half-filling XXZ eigenstates, coupled by Z_k:
        R_ab = sum_k |<E_a|Z_k|E_b>|^2   (a != b, gain)
        R_aa = -4 sum_k Var_a(n_k)        (loss)
    Built on the SECTOR Hamiltonian (dim C(N,p)) via xxz_Hp -- NOT the full 2^N H. This is what
    makes N~14-16 feasible (a full 2^16 H is ~69 GB complex128; the sector H is 12870x12870)."""
    p = (N + 1) // 2
    H, states = xxz_Hp(N, p, Delta, J)        # sector H, dim C(N,p); states = popcount-p bitmasks
    E, V = eigh(H)                            # V[:, a] = |E_a> in the sector computational basis
    ns = len(states)
    nk = np.array([[(s >> k) & 1 for k in range(N)] for s in states])  # occupations, ns x N
    Vd = V.conj().T
    R = np.zeros((ns, ns))
    for k in range(N):
        zk = 1.0 - 2.0 * nk[:, k]             # Z_k eigenvalue +-1 per sector state (diagonal op)
        Mk = Vd @ (zk[:, None] * V)           # <E_a|Z_k|E_b>
        R += np.abs(Mk) ** 2                  # GAIN in Z_k (the load-bearing factor-4 choice)
    # the gain loop also accumulated a diagonal (the Z self-term sum_k <E_a|Z_k|E_a>^2); it is NOT
    # the generator loss and is unconditionally overwritten below by the true loss -4*Var_a(n).
    w = np.abs(V) ** 2                         # |<s|E_a>|^2, ns x ns
    mean_n = w.T @ nk                          # <n_k>_a, ns x N
    var = (mean_n - mean_n ** 2).sum(axis=1)   # sum_k Var_a(n_k)  (n^2=n so <n^2>=<n>)
    R[np.diag_indices(ns)] = -4.0 * var        # LOSS in n_k (overwrites the gain loop's diagonal)
    return R


def gapR(N, Delta):
    """Magnitude of R's slowest nonzero relaxation mode (= Lebensader rate / gamma)."""
    ev = np.sort(eigh(rate_matrix_R(N, Delta))[0])  # ascending; ev[-1] ~ 0 is the steady state
    nz = ev[ev < -1e-9]
    return -nz.max()                                 # smallest-magnitude nonzero mode


def check_R_is_generator():
    """assert #6: R real-symmetric, zero column sums, single zero mode. The Z<->n factor of 4 is
    LOAD-BEARING (gain in Z_k, loss in n_k via Var(Z)=4Var(n)); harmonizing them breaks colsum by 4x."""
    for N in (4, 5):
        R = rate_matrix_R(N, 1.5)
        assert np.max(np.abs(R - R.T)) < 1e-12, f"N={N}: R not symmetric"
        assert np.max(np.abs(R.sum(axis=0))) < 1e-10, f"N={N}: R columns do not sum to zero"
        ev = np.sort(eigh(R)[0])
        assert np.sum(np.abs(ev) < 1e-9) == 1, f"N={N}: R does not have exactly one zero mode"
    print("[6] R is a valid generator (symmetric, zero column sums, single zero mode).  OK")


DSTAR_GAMMA0 = {4: 1.61961, 5: 1.52798, 6: 1.38463, 7: 1.33007, 8: 1.27243}


def delta_star_reduction(N, lo=1.0, hi=2.5):
    """Delta* where gap(R)=2 (the gamma->0 handover condition).
    gap(R)-2 is MONOTONE-decreasing through zero on [1, 2.5] (more Delta -> more Neel order ->
    slower Lebensader -> smaller gap), so this is a ROOT-find of a monotone function. brentq is
    correct and safe -- this is NOT the multimodal Brent MINIMIZATION that failed in the
    ptf_painter arc (that was a non-convex global-min search; a bracketed monotone root-find is
    a different algorithm). f(lo)>0 (band-edge regime, gap>2), f(hi)<0 (Lebensader regime, gap<2)."""
    f = lambda D: gapR(N, D) - 2.0
    assert f(lo) > 0 > f(hi), f"N={N}: gap(R)-2 does not bracket a root on [{lo},{hi}]"
    return brentq(f, lo, hi, xtol=1e-7)


def check_reduction_matches_full():
    """assert #1: gamma*gap(R) -> the full-sector Lebensader rate as gamma->0 (ratio -> 1)."""
    for N in (4, 5):
        gap = gapR(N, 1.5)
        ratios = []
        for g in (0.05, 0.0125, 0.003125):
            full = lebensader_rate(N, 1.5, gamma=g)  # the (p,p) sector method, finite gamma
            ratios.append(full / (g * gap))
        # the ratio must approach 1 as gamma shrinks (the gamma->0 limit is exact)
        assert abs(ratios[-1] - 1.0) < 5e-3, f"N={N}: gamma*gap(R) not matching full as g->0: {ratios}"
        assert abs(ratios[-1] - 1.0) < abs(ratios[0] - 1.0) + 1e-9, f"N={N}: not converging: {ratios}"
    print("[1] gamma*gap(R) -> full Lebensader rate as gamma->0 (ratio -> 1).  OK")


def check_reduction_reproduces_gamma0_dstar():
    """assert #2b: delta_star_reduction reproduces the gamma->0 Delta* (NOT the Q=20 numbers --
    the N=4 drift 1.7e-3 exceeds a 1e-3 tolerance, so asserting vs Q=20 would falsely fail)."""
    for N, expected in DSTAR_GAMMA0.items():
        got = delta_star_reduction(N)
        assert abs(got - expected) < 1e-3, f"gamma->0 Delta*({N})={got} != {expected}"
    print("[2b] gamma->0 reduction reproduces Delta*(4..8) = "
          + ", ".join(f"{v:.5f}" for v in DSTAR_GAMMA0.values()) + ".  OK")


if __name__ == "__main__":
    check_finite_gamma_baseline()
    check_R_is_generator()
    check_reduction_matches_full()
    check_reduction_reproduces_gamma0_dstar()
