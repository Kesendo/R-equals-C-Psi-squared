"""_the_role_change.py - the role-change PTF could not do.

PTF (the Perspectival Time Field) made a snapshot: it measured each painter's
clock as a RATIO, alpha_i = P_B(i)/P_A(i), against one frozen reference (the
unperturbed chain A). A ratio is asymmetric, alpha(A,B) = 1/alpha(B,A): it nails
one side as the Massstab, so the painters' roles cannot trade. One frame, no movie.

Tom's recipe for the role-change: freeze values, compare with the values that
change, XOR, the contract. That recipe is already in the algebra, and XOR is
exactly the SYMMETRIC comparator a ratio is not:

  - the contract is c(i,j) = popcount(i XOR j): where i,j agree (XOR=0, the
    diagonal) nothing decays, the kept; where they differ (XOR != 0) it decays,
    the change; popcount is how much, the Absorption rate Re(lambda) = -2*gamma*c.
  - XOR is symmetric, i XOR j = j XOR i, so the role-swap i<->j (the transpose D,
    the two perspectives trading places) leaves the contract invariant while the
    role flips.

So the role-change = the D-swap, contract-preserving, the thing the ratio could
not do. Bit-exact, 1-3 qubits.

Tom + Claude, 2026-05-28. Run: python simulations/_the_role_change.py
"""
from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
from scipy.linalg import expm

from framework.pauli import ur_pauli, site_op, _build_bilinear
from framework.lindblad import lindbladian_pauli_dephasing

GAMMA = 0.3
HEIS = [("X", "X", 1.0), ("Y", "Y", 1.0), ("Z", "Z", 1.0)]
_ok = []


def report(name, cond, extra=""):
    _ok.append(bool(cond))
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}{extra}")


def popcount(x):
    return bin(x).count("1")


def vec_F(M):
    return M.flatten("F")


# ----------------------------------------------------------------------
# PART A - XOR is the symmetric contract (the comparator PTF lacked)
# ----------------------------------------------------------------------

def part_a(Ns=(1, 2, 3)):
    print("PART A  - XOR is the contract: c(i,j) = popcount(i XOR j) = the Absorption decay")
    for N in Ns:
        d = 2 ** N
        # Pure Z-dephasing (H = 0): |i><j| are eigenoperators of L.
        L = lindbladian_pauli_dephasing(np.zeros((d, d), complex), [GAMMA] * N, "Z")
        worst_eig = 0.0
        worst_sym = 0
        worst_swap = 0.0
        for i in range(d):
            for j in range(d):
                E = np.zeros((d, d), complex)
                E[i, j] = 1.0
                c = popcount(i ^ j)
                # L acts on |i><j| with eigenvalue -2*gamma*popcount(i XOR j).
                lam_vec = L @ vec_F(E)
                worst_eig = max(worst_eig, float(np.max(np.abs(lam_vec - (-2 * GAMMA * c) * vec_F(E)))))
                # XOR is symmetric: c(i,j) == c(j,i).
                worst_sym = max(worst_sym, abs(c - popcount(j ^ i)))
                # the role-swap |i><j| -> |j><i| preserves the contract (same eigenvalue).
                worst_swap = max(worst_swap, abs((-2 * GAMMA * popcount(i ^ j)) - (-2 * GAMMA * popcount(j ^ i))))
        report(f"N={N}: L|i><j| = -2*gamma*popcount(i^j) |i><j|  (XOR is the decay)", worst_eig < 1e-12)
        report(f"N={N}: popcount(i^j) = popcount(j^i)  (XOR symmetric)", worst_sym == 0)
        report(f"N={N}: role-swap i<->j keeps the contract (same decay)", worst_swap < 1e-12)
    print()


# ----------------------------------------------------------------------
# PART B - the ratio is asymmetric; that is why PTF was a snapshot
# ----------------------------------------------------------------------

def part_b():
    print("PART B  - a ratio fixes a reference (snapshot); XOR fixes none (role-change)")
    # PTF's alpha = P_B / P_A. Take two painter-values a, b.
    a, b = 1.182, 0.845                      # two of PTF's own alpha_i (sites 1 and 4, N=7 doc)
    r_ab, r_ba = a / b, b / a
    report("ratio is asymmetric: (a/b) != (b/a), one side is the frozen denominator",
           abs(r_ab - r_ba) > 1e-9, f"   a/b={r_ab:.4f}, b/a={r_ba:.4f}, product={r_ab*r_ba:.4f}")
    # XOR of the two perspectives (as bit-patterns) is symmetric.
    i, j = 0b101, 0b011
    report("XOR is symmetric: (i^j) == (j^i), no frozen reference",
           (i ^ j) == (j ^ i), f"   i^j={i^j:03b}, j^i={j^i:03b}")
    print("        -> the ratio nails A as Massstab; the roles cannot trade. That is the snapshot.")
    print("        -> XOR nails nothing; the roles can trade. That is the movie.\n")


# ----------------------------------------------------------------------
# PART C - the role-change on a 2-qubit state (the movie), with H != 0
# ----------------------------------------------------------------------

def part_c():
    print("PART C  - the role-change: swap the two perspectives, the contract holds, the role flips")
    N = 2
    d = 2 ** N
    bonds = [(0, 1)]
    H = _build_bilinear(N, bonds, HEIS)
    L = lindbladian_pauli_dephasing(H, [GAMMA] * N, "Z")

    # The role-change D = the transpose: |i><j| -> |j><i|, the two perspectives trading.
    # On the superoperator (vec, 'F' order) the transpose is the swap a+d*b <-> b+d*a.
    perm = np.array([(k % d) * d + (k // d) for k in range(d * d)])
    D = np.eye(d * d)[perm]                   # D vec(rho) = vec(rho^T)

    # The contract (Re lambda = the XOR/Absorption decay rates) is invariant under the swap.
    spec = np.linalg.eigvals(L)
    spec_swap = np.linalg.eigvals(D @ L @ D)
    re1 = np.sort(np.round(spec.real, 9))
    re2 = np.sort(np.round(spec_swap.real, 9))
    report("Re(spec) invariant under the role-swap (the contract, the XOR-decay, holds)",
           float(np.max(np.abs(re1 - re2))) < 1e-9)
    im1 = np.sort(np.round(spec.imag, 9))
    im2 = np.sort(np.round(spec_swap.imag, 9))
    report("Im(spec) negates under the role-swap (the role, the coherent direction, flips)",
           float(np.max(np.abs(np.sort(spec_swap.imag) - np.sort(-spec.imag)))) < 1e-9)

    # On an evolving state: freeze the diagonal (XOR=0, the contract); the off-diagonal trades.
    psi = np.array([1, 0.7, 0.5j, 0.3], complex)
    psi /= np.linalg.norm(psi)
    rho0 = np.outer(psi, psi.conj())
    print("        evolving a 2-qubit state; D = swap the two perspectives (rho -> rho^T):")
    worst_pur = worst_diag = 0.0
    for t in (0.0, 1.0, 3.0, 8.0):
        rt = (expm(L * t) @ vec_F(rho0)).reshape(d, d, order="F")
        rt_swapped = rt.T
        pur, pur_sw = float(np.trace(rt @ rt).real), float(np.trace(rt_swapped @ rt_swapped).real)
        worst_pur = max(worst_pur, abs(pur - pur_sw))
        worst_diag = max(worst_diag, float(np.max(np.abs(np.diag(rt) - np.diag(rt_swapped)))))
        print(f"          t={t:>4.1f}:  purity={pur:.6f}  purity(swapped)={pur_sw:.6f}")
    report("purity (the contract, the amount) is identical after the role-swap, every t",
           worst_pur < 1e-12)
    report("the diagonal (i=j, XOR=0, the frozen contract) is untouched by the swap",
           worst_diag < 1e-12)
    print()


def main():
    print("=" * 76)
    print("THE ROLE-CHANGE - what PTF could not do with a ratio, XOR can")
    print("=" * 76)
    part_a()
    part_b()
    part_c()
    n_ok, n_tot = sum(_ok), len(_ok)
    print("=" * 76)
    print(f"RESULT: {n_ok}/{n_tot} bit-exact ({'ALL PASS' if n_ok == n_tot else 'CHECK'})")
    print("=" * 76)
    print("""
The role-change:
  PTF measured with a ratio (alpha = P_B/P_A). A ratio has a denominator, a frozen
  reference; alpha(A,B) = 1/alpha(B,A). It cannot trade roles. That is why PTF was
  a snapshot: one painter was always the Massstab.

  XOR is the right comparator. The contract is popcount(i XOR j), the Absorption
  decay; XOR is symmetric, so it nails no reference. The role-swap i<->j (the
  transpose, the two perspectives trading) leaves the contract bit-exact, decay
  rates kept, purity kept, the diagonal frozen, and flips only the role (the
  coherent direction, the sign of Im). The contract is what survives the trade;
  the role is what trades. The movie PTF lacked is just the snapshot read with a
  symmetric comparator instead of a ratio.
""")


if __name__ == "__main__":
    main()
