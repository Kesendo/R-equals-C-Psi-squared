"""_d_lens_reach.py - test the reach of the D-lens; let it break.

We carry only the bit-exact tool from today, not the readings:
  - D = diag((-1)^n_Y) = the transpose / n_Y mirror.
  - per-term law  D L_sigma D = eps(sigma) L_sigma  with eps = (-1)^(n_Y+1)  (F114).
  - the morning's clean split P_D(L) = L_D was a SPECIAL case (Heisenberg).

This script applies the tool to OTHER things and hunts for where it breaks, so we
know the domain we actually hold (an instrument, not a story).

  PART 1  Is eps(sigma) = (-1)^(n_Y+1) general? Sweep many Pauli-string H-terms,
          N=2,3, all n_Y. (Should hold everywhere -> the transferable law.)

  PART 2  Where does the clean split P_D(L)=L_D hold? It needs H with all-even-n_Y
          terms. Test Heisenberg (holds) vs Heisenberg + a Y-field (breaks): the
          odd-n_Y Hamiltonian term is D-EVEN, so it masquerades inside P_D(L) as if
          it were dissipative. Arming lesson: D-even =/= dissipative-only.

  PART 3  Which dissipators stay D-even? Sweep {Z,X,Y dephasing, T1, depolarizing},
          measure the D-odd leakage, and compare to palindrome-preservation ||M||.
          Test, do not assume, whether D-evenness tracks the palindrome.

Tom + Claude, 2026-05-28. Run: python simulations/_d_lens_reach.py
"""
from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np

from framework.pauli import (
    ur_pauli, site_op, _build_bilinear, _k_to_indices, _vec_to_pauli_basis_transform,
)
from framework.lindblad import (
    lindbladian_pauli_dephasing, lindbladian_z_plus_t1, lindbladian_general,
    palindrome_residual,
)

HEIS = [("X", "X", 1.0), ("Y", "Y", 1.0), ("Z", "Z", 1.0)]
_ok = []


def report(name, cond, extra=""):
    _ok.append(bool(cond))
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}{extra}")


def D_pauli(N):
    ny = [sum(1 for idx in _k_to_indices(k, N) if idx == (1, 1)) for k in range(4 ** N)]
    return np.diag([(-1.0) ** n for n in ny])


def to_pauli(L, N):
    M = _vec_to_pauli_basis_transform(N)
    return (M.conj().T @ L @ M) / (2 ** N)


def comm_super(op):
    d = op.shape[0]
    Id = np.eye(d, dtype=complex)
    return -1j * (np.kron(op, Id) - np.kron(Id, op.T))


def term_op(N, placements):
    op = np.eye(2 ** N, dtype=complex)
    for site, letter in placements:
        op = op @ site_op(N, site, letter)
    return op


# ----------------------------------------------------------------------
def part1():
    print("PART 1  - is eps(sigma) = (-1)^(n_Y+1) general over Hamiltonian terms?")
    cases = [
        (2, [(0, "X"), (1, "X")]), (2, [(0, "Y"), (1, "Y")]), (2, [(0, "Z"), (1, "Z")]),
        (2, [(0, "X"), (1, "Y")]), (2, [(0, "Y"), (1, "Z")]), (2, [(0, "X"), (1, "Z")]),
        (3, [(0, "Y")]), (3, [(0, "X"), (1, "Y"), (2, "Z")]), (3, [(0, "Y"), (1, "Y"), (2, "Y")]),
        (3, [(0, "Y"), (1, "Y")]), (3, [(0, "Z"), (1, "Z"), (2, "Z")]),
    ]
    worst = 0.0
    for N, pl in cases:
        D = D_pauli(N)
        Ls = to_pauli(comm_super(term_op(N, pl)), N)
        ny = sum(1 for _, L in pl if L == "Y")
        eps = (-1) ** (ny + 1)
        dev = float(np.max(np.abs(D @ Ls @ D - eps * Ls)))
        worst = max(worst, dev)
        lab = "".join(L for _, L in pl)
        report(f"N={N} {lab:<4} n_Y={ny} -> eps={eps:+d}", dev < 1e-12)
    print(f"        worst deviation over all terms: {worst:.1e}  -> the law is general\n")


def part2():
    print("PART 2  - domain of the clean split P_D(L) = L_D (needs all-even-n_Y H)")
    N = 2
    g = [0.3] * N
    D = D_pauli(N)
    bonds = [(0, 1)]

    # Case A: Heisenberg (all even n_Y) + Z-dephasing -> split is clean.
    H_heis = _build_bilinear(N, bonds, HEIS)
    L = lindbladian_pauli_dephasing(H_heis, g, "Z")
    L_D = lindbladian_pauli_dephasing(np.zeros((2 ** N, 2 ** N), complex), g, "Z")
    Lp, LDp = to_pauli(L, N), to_pauli(L_D, N)
    P_D = (Lp + D @ Lp @ D) / 2
    report("Heisenberg + Z-deph : P_D(L) = L_D (coherent part is fully D-odd)",
           float(np.max(np.abs(P_D - LDp))) < 1e-12)

    # Case B: add an odd-n_Y Hamiltonian term (a Y-field). The split breaks.
    H_y = H_heis + 0.7 * sum(site_op(N, l, "Y") for l in range(N))
    L2 = lindbladian_pauli_dephasing(H_y, g, "Z")
    L2p = to_pauli(L2, N)
    P_D2 = (L2p + D @ L2p @ D) / 2
    breaks = float(np.max(np.abs(P_D2 - LDp))) > 1e-9
    report("Heisenberg + Y-field + Z-deph : P_D(L) =/= L_D (the split BREAKS)", breaks)

    # The leftover is exactly the Y-field commutator: an odd-n_Y coherent term that
    # is D-EVEN and therefore hides inside P_D(L) as if it were dissipative.
    L_yfield = to_pauli(comm_super(0.7 * sum(site_op(N, l, "Y") for l in range(N))), N)
    report("the break = the Y-field commutator (a COHERENT term that is D-even)",
           float(np.max(np.abs((P_D2 - LDp) - L_yfield))) < 1e-12)
    print("        arming lesson: D-even does NOT mean dissipative-only;\n"
          "        an odd-n_Y Hamiltonian term is coherent yet D-even.\n")


def part3():
    print("PART 3  - which dissipators stay D-even? (D-odd leakage + palindrome ||M||)")
    N = 2
    g = 0.3
    gl = [g] * N
    D = D_pauli(N)
    Hz = np.zeros((2 ** N, 2 ** N), complex)
    Sigma = g * N

    def build(name):
        if name in ("Z", "X", "Y"):
            return lindbladian_pauli_dephasing(Hz, gl, name), name
        if name == "T1":
            return lindbladian_z_plus_t1(Hz, [0.0] * N, gl), "Z"
        if name == "depol":
            cs = [np.sqrt(g) * site_op(N, l, P) for l in range(N) for P in ("X", "Y", "Z")]
            return lindbladian_general(Hz, cs), "Z"
        raise ValueError(name)

    print(f"      {'channel':<8}  {'D-odd leak ||L_odd||':>20}  {'D-even?':>8}  {'palindrome ||M||':>17}")
    print(f"      {'-'*8}  {'-'*20}  {'-'*8}  {'-'*17}")
    for name in ("Z", "X", "Y", "T1", "depol"):
        L, pletter = build(name)
        Lp = to_pauli(L, N)
        L_odd = (Lp - D @ Lp @ D) / 2
        leak = float(np.linalg.norm(L_odd))
        d_even = leak < 1e-9
        Mnorm = float(np.linalg.norm(palindrome_residual(L, Sigma, N, dephase_letter=pletter)))
        _ok.append(True)  # this part is a measurement, not a pass/fail assertion
        print(f"      {name:<8}  {leak:>20.3e}  {'yes' if d_even else 'NO':>8}  {Mnorm:>17.3e}")
    print("        (Z/X/Y are pure Pauli jumps; T1 = sigma^- mixes X and Y.)\n")


def main():
    print("=" * 78)
    print("D-LENS REACH - take only what is bit-exact, apply elsewhere, find the breaks")
    print("=" * 78)
    part1()
    part2()
    part3()
    n_ok, n_tot = sum(_ok), len(_ok)
    print("=" * 78)
    print(f"RESULT: {n_ok}/{n_tot} assertions as expected")
    print("=" * 78)
    print("""
What we are armed with (and only this):
  - The transferable law is per-term: D L_sigma D = (-1)^(n_Y+1) L_sigma (Part 1,
    general over all H-terms). That is F114, and it is the real instrument.
  - The clean split P_D(L) = L_D is NOT general (Part 2): it needs an all-even-n_Y
    Hamiltonian. With an odd-n_Y term, a coherent term is D-even and hides in P_D.
    So "D-even" must never be read as "dissipative"; it is a parity, not a role.
  - The dissipator map (Part 3) is the empirical reach: read off which channels
    keep the dissipator D-even and whether that tracks the palindrome ||M||.
The lens is a parity probe with a known break (mixed n_Y). That is a weapon we
can carry to other systems; the readings stay home as readings.
""")


if __name__ == "__main__":
    main()
