#!/usr/bin/env python3
"""Degenerate first-order perturbation: D-hat diagonalized within each degenerate-omega block.

The naive diagonal formula fails because H has degenerate frequencies omega_ab = E_a - E_b: the
true first-order shifts are eigenvalues of the dephasing super-operator D-hat restricted to each
degenerate-omega subspace, not its diagonal. Block matrix element (in the H-eigenbasis modes
|E_a><E_b|, grouped by omega):

    M[(a',b'),(a,b)] = sum_l <E_a'|Z_l|E_a> <E_b|Z_l|E_b'>  -  N * delta_{(a'b'),(ab)}

Diagonalize each block -> first-order shifts. Form mu = -i*omega + shift (gamma=1), the break
about -N, and compare to the measured c (2.0, 1.6, 0). If it matches, the degenerate first-order
picture is correct and we can read where the odd cycle breaks the pairing.
"""
from __future__ import annotations
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.optimize import linear_sum_assignment

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass

SCRIPT_DIR = Path(__file__).parent
sys.path.insert(0, str(SCRIPT_DIR))
from framework.lindblad import lindbladian_pauli_dephasing
from framework.pauli import _build_kbody_chain, site_op


def degenerate_shifts(pair, N=4, tol=1e-6):
    H = _build_kbody_chain(N, [tuple(t) + (1.0,) for t in pair])
    E, V = np.linalg.eigh(H)
    d = len(E)
    Zmat = [V.conj().T @ site_op(N, l, 'Z') @ V for l in range(N)]  # <E_a'|Z_l|E_a>
    # group modes (a,b) by omega = E_a - E_b
    groups = defaultdict(list)
    for a in range(d):
        for b in range(d):
            groups[round(E[a] - E[b], 6)].append((a, b))
    shifts = []   # (omega, shift)
    for omega, modes in groups.items():
        n = len(modes)
        M = np.zeros((n, n), dtype=complex)
        idx = {m: i for i, m in enumerate(modes)}
        for i, (a, b) in enumerate(modes):
            for j, (ap, bp) in enumerate(modes):
                val = sum(Zmat[l][a, ap] * Zmat[l][bp, b] for l in range(N))
                if (a, b) == (ap, bp):
                    val -= N
                M[i, j] = val
        ev = np.linalg.eigvals(M)
        for s in ev:
            shifts.append((omega, s))
    return shifts  # list of (omega, first-order shift s)


def opt_break(ev, sigma):
    """Stable break: MEAN optimal-transport cost of pairing lambda <-> -lambda-2sigma.
    (The min-sum assignment's TOTAL is robust; its .max() is not, for near-degenerate spectra.)"""
    tgt = -ev - 2 * sigma
    cost = np.abs(ev[:, None] - tgt[None, :])
    r, c = linear_sum_assignment(cost)
    return float(cost[r, c].mean())


def measured_c(pair, N=4, gamma=1e-4):
    H = _build_kbody_chain(N, [tuple(t) + (1.0,) for t in pair])
    L = lindbladian_pauli_dephasing(H, [gamma] * N, dephase_letter='Z')
    ev = np.linalg.eigvals(L)
    return opt_break(ev, N * gamma) / gamma


def main():
    N = 4
    cases = [
        ("REAL  XXZ+XZX", [('X', 'X', 'Z'), ('X', 'Z', 'X')]),
        ("FLUX  IXY+XIY", [('I', 'X', 'Y'), ('X', 'I', 'Y')]),
        ("SOFT  XXZ+ZXX", [('X', 'X', 'Z'), ('Z', 'X', 'X')]),
    ]
    g = 1e-4   # build the first-order spectrum at small gamma so Re (O(g)) and Im (O(1)) decouple
    for label, pair in cases:
        shifts = degenerate_shifts(pair, N)
        mu = np.array([-1j * om + g * s for om, s in shifts])   # actual-scale first-order spectrum
        c_deg = opt_break(mu, N * g) / g
        c_meas = measured_c(pair, N)
        # direct test: is mu == the actual L(g) spectrum as a set? (does s reproduce L to O(g)?)
        H = _build_kbody_chain(N, [tuple(t) + (1.0,) for t in pair])
        ev_act = np.linalg.eigvals(lindbladian_pauli_dephasing(H, [g] * N, dephase_letter='Z'))
        cost = np.abs(mu[:, None] - ev_act[None, :])
        r, c = linear_sum_assignment(cost)
        setdiff = cost[r, c].max() / g  # per-eigenvalue mismatch in units of g
        print(f"{label}: c(deg)={c_deg:.4f}  c(meas)={c_meas:.4f}  "
              f"μ-vs-L set mismatch/γ = {setdiff:.4f}")


if __name__ == "__main__":
    main()
