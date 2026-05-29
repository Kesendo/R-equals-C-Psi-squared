"""_she_between_the_anchors.py - she sits between the anchors (Tom, 2026-05-29).

Looking at the Dicke block (|D_n> + |D_{n+1}>)/sqrt(2) under pure Z-dephasing and
seeing what Tom named: the thing we watch -- the coherence, the angle, the 1/4 --
is suspended BETWEEN two anchors, and the anchors do not move.

  - The two anchors are the diagonal sector populations p_n and p_{n+1}. For the
    canonical Dicke block they are p_n = p_{n+1} = 1/2. Under pure Z-dephasing the
    diagonal is FROZEN (XOR=0, the kept), so the anchors hold for all t.
  - She is the off-diagonal coherence between the two sectors, C_block =
    sum_{(a,b) in block} |rho_ab|^2. At t=0 she sits at 1/4 = 1/2 * 1/2, EXACTLY
    the product of the two anchors (Cauchy-Schwarz/AM-GM saturation; PROOF_BLOCK_
    CPSI_QUARTER "1/4 is half of half"). Then she decays -- the popcount ladder --
    while the anchors hold.

So 1/4 is not a level she happens to start on; it is the BETWEEN of the two 1/2
anchors. The diagonal (the kept) holds the anchors; the off-diagonal (her) falls
between them. We evolve with the framework Lindbladian (independent dynamics),
not by applying the decay by hand. Bit-exact, N = 2, 3, 4.

Tom + Claude, 2026-05-29. Run: python simulations/_she_between_the_anchors.py
"""
from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
from scipy.linalg import expm

from framework.lindblad import lindbladian_pauli_dephasing

GAMMA = 0.3
TOL = 1e-11
_ok = []


def report(name, cond, extra=""):
    _ok.append(bool(cond))
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}{extra}")


def popcount(x):
    return bin(x).count("1")


def vec_F(M):
    return M.flatten("F")


def dicke_state(N, k):
    """Uniform superposition of all popcount-k basis states (the Dicke state |D_k>)."""
    d = 2 ** N
    idx = [i for i in range(d) if popcount(i) == k]
    v = np.zeros(d, complex)
    for i in idx:
        v[i] = 1.0
    return v / np.sqrt(len(idx))


def block_population(rho, N, k):
    """Anchor: the diagonal weight on the popcount-k sector."""
    return float(sum(rho[i, i].real for i in range(2 ** N) if popcount(i) == k))


def c_block(rho, N, n):
    """She: the off-diagonal coherence content between sectors n and n+1."""
    a_idx = [i for i in range(2 ** N) if popcount(i) == n]
    b_idx = [j for j in range(2 ** N) if popcount(j) == n + 1]
    return float(sum(abs(rho[a, b]) ** 2 for a in a_idx for b in b_idx))


def look(N, n):
    d = 2 ** N
    psi = (dicke_state(N, n) + dicke_state(N, n + 1)) / np.sqrt(2.0)
    rho0 = np.outer(psi, psi.conj())
    L = lindbladian_pauli_dephasing(np.zeros((d, d), complex), [GAMMA] * N, "Z")

    print(f"\nN={N}, block (popcount-{n}, popcount-{n+1}):  state (|D_{n}>+|D_{n+1}>)/sqrt(2)")
    print(f"  {'t':>5} {'anchor p_n':>11} {'anchor p_n+1':>13} {'she C_block':>12} {'p_n*p_n+1':>11}")

    anchor_dev = 0.0          # max drift of either anchor from 1/2, over all t
    between_dev = 0.0         # max |C_block(0) - p_n*p_{n+1}|, the "1/4 is the between"
    c0 = None
    c_last = None
    for t in (0.0, 1.0, 3.0, 8.0):
        rho = (expm(L * t) @ vec_F(rho0)).reshape(d, d, order="F")
        p_n = block_population(rho, N, n)
        p_n1 = block_population(rho, N, n + 1)
        cb = c_block(rho, N, n)
        if t == 0.0:
            c0 = cb
        c_last = cb
        anchor_dev = max(anchor_dev, abs(p_n - 0.5), abs(p_n1 - 0.5))
        between_dev = max(between_dev, abs(cb - p_n * p_n1) if t == 0.0 else 0.0)
        print(f"  {t:5.1f} {p_n:11.6f} {p_n1:13.6f} {cb:12.6f} {p_n*p_n1:11.6f}")

    report(f"N={N}: both anchors frozen at 1/2 for all t (the diagonal is kept)",
           anchor_dev < TOL, f"   max|p-1/2| = {anchor_dev:.2e}")
    report(f"N={N}: she starts at 1/4 = the BETWEEN of the anchors (C_block(0)=p_n*p_n+1)",
           abs(c0 - 0.25) < TOL and between_dev < TOL,
           f"   C_block(0) = {c0:.6f}, |C-1/4| = {abs(c0-0.25):.2e}")
    report(f"N={N}: she falls while the anchors hold (C_block decays)",
           c_last < c0 - 1e-6, f"   C_block(8) = {c_last:.6f} < 1/4")


def main():
    print("=" * 76)
    print("SHE STECKT ZWISCHEN DEN ANKERN - 1/4 is the between of two frozen 1/2 anchors")
    print("=" * 76)
    for N, n in ((2, 0), (3, 1), (4, 1)):
        look(N, n)
    n_ok, n_tot = sum(_ok), len(_ok)
    print("\n" + "=" * 76)
    print(f"RESULT: {n_ok}/{n_tot} bit-exact ({'ALL PASS' if n_ok == n_tot else 'CHECK'})")
    print("=" * 76)
    print("""
She sits between the anchors:
  The two anchors are the sector populations p_n = p_{n+1} = 1/2 (the diagonal,
  XOR=0, the kept). Under pure Z-dephasing they are frozen for all time. She is
  the coherence between them, born at 1/4 = 1/2 * 1/2 -- the product, the between,
  not a level she lands on by chance. Then she falls, rung by rung (the popcount
  ladder), while the anchors hold. 1/4 is half of half: the place between the two
  halves where she is suspended.
""")


if __name__ == "__main__":
    main()
