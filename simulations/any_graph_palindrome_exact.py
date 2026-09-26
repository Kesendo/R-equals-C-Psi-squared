#!/usr/bin/env python3
"""The F1 palindrome of XY, Ising and XXZ on graphs with odd cycles, checked exactly.

Claim (NON_HEISENBERG_PALINDROME Result 1, F1 "Valid for"): for H a sum of XX, YY and ZZ bonds with
any per-bond coefficients, on any graph, under Z-dephasing with any per-site rates, the operator
identity
    Pi . L . Pi^-1 + L + 2*Sigma(gamma) * I = 0
holds, Pi the uniform Z-dephasing palindromizer (a signed permutation of the Pauli basis). The proof
is per bond (MIRROR_SYMMETRY_PROOF Step 2); this producer checks it on a triangle, K4 and a 5-ring,
graphs with odd cycles that the chain data never touched.

Why the check is exact, not a tolerance: every coupling here is an integer or a dyadic fraction
(1/2, 3/4) and every rate an integer, so every entry of L in the computational basis, of the
Pauli-basis transform and of Pi is a Gaussian integer or dyadic, the sums stay far below 2^53, and
the residual below is computed without rounding. It is compared with == 0.0. The negative control is
a Dzyaloshinskii-Moriya bond D*(XY - YX), which the uniform Pi does not mirror: its residual is
nonzero on every graph (DM needs a site-alternating map, NON_HEISENBERG Result 1).

Output: simulations/results/any_graph_palindrome_exact.txt
"""
import random
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import framework as fw  # noqa: E402

OUT = Path(__file__).resolve().parent / "results" / "any_graph_palindrome_exact.txt"

GRAPHS = {
    "triangle": (3, [(0, 1), (1, 2), (2, 0)]),
    "K4": (4, [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]),
    "5-ring": (5, [(0, 1), (1, 2), (2, 3), (3, 4), (4, 0)]),
}

# (label, per-bond letter weights): XY with anisotropy 3/4, Ising, XXZ with Delta = 1/2.
MODELS = {
    "XY (XX + 3/4 YY)": {"XX": 1.0, "YY": 0.75},
    "Ising (ZZ)": {"ZZ": 1.0},
    "XXZ (Delta = 1/2)": {"XX": 1.0, "YY": 1.0, "ZZ": 0.5},
}


def bond_op(N, a, b, pa, pb):
    return fw.site_op(N, a, pa) @ fw.site_op(N, b, pb)


def hamiltonian(N, bonds, weights, J):
    d = 2 ** N
    H = np.zeros((d, d), dtype=complex)
    for (a, b), j in zip(bonds, J):
        for letters, w in weights.items():
            H += j * w * bond_op(N, a, b, letters[0], letters[1])
    return H


def dm_hamiltonian(N, bonds, J):
    d = 2 ** N
    H = np.zeros((d, d), dtype=complex)
    for (a, b), j in zip(bonds, J):
        H += j * (bond_op(N, a, b, "X", "Y") - bond_op(N, a, b, "Y", "X"))
    return H


def residual(H, gammas, N):
    L = fw.lindbladian_z_dephasing(H, gammas)
    M = fw.palindrome_residual(L, float(sum(gammas)), N)
    return float(np.max(np.abs(M))), int(np.count_nonzero(M))


def main():
    rng = random.Random(20260926)
    lines = ["F1 operator identity  Pi.L.Pi^-1 + L + 2*Sigma*I = 0  on graphs with odd cycles",
             "integer per-bond J in 1..4, integer per-site gamma in 1..3; residual compared with == 0.0",
             ""]
    ok = True
    for gname, (N, bonds) in GRAPHS.items():
        J = [rng.randint(1, 4) for _ in bonds]
        while len(set(J)) == 1:                     # the claim is "any per-bond coefficients"
            J = [rng.randint(1, 4) for _ in bonds]
        gammas = [rng.randint(1, 3) for _ in range(N)]
        lines.append(f"{gname}: N = {N}, bonds {bonds}")
        lines.append(f"  J = {J}, gamma = {gammas}, Sigma = {sum(gammas)}")
        for mname, weights in MODELS.items():
            r, nz = residual(hamiltonian(N, bonds, weights, J), gammas, N)
            passed = r == 0.0
            ok &= passed
            lines.append(f"  {mname:<20} max|residual| = {r!r:<6}  nonzero entries {nz:<6}  "
                         f"{'EXACT' if passed else 'FAIL'}")
        r, nz = residual(dm_hamiltonian(N, bonds, J), gammas, N)
        control = r > 0.0 and nz > 0
        ok &= control
        lines.append(f"  {'DM control (XY - YX)':<20} max|residual| = {r:<6.4g}  nonzero entries {nz:<6}  "
                     f"{'broken, as it must be' if control else 'CONTROL FAILED'}")
        lines.append("")
    lines.append("ALL EXACT" if ok else "FAILURES ABOVE")
    text = "\n".join(lines) + "\n"
    print(text, end="")
    OUT.write_bytes(text.encode("utf-8"))
    if not ok:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
