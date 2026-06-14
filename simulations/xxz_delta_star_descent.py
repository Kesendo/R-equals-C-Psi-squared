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


if __name__ == "__main__":
    check_finite_gamma_baseline()
