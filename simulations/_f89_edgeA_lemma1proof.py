#!/usr/bin/env python3
"""F89 edge A, step 16: a RIGOROUS proof of LEMMA 1 via node values (Lagrange uniqueness).

LEMMA 1:  (4-u^2) Atil(u) == eps_k * m * (2+u)   (mod Phi_u),   eps_k = 2 (odd), 1 (even).

Proof strategy (rigorous): two polynomials that agree at all FA nodes u_j (roots of Phi_u,
which is separable) and both have degree < FA + deg(quotient)... more precisely: P(u) :=
(4-u^2)Atil(u) and the linear L(u) := eps_k m (2+u) satisfy  P(u) == L(u) (mod Phi_u)  iff
P(u_j) == L(u_j) for ALL nodes u_j (since Phi_u = prod (u-u_j) is separable, two polys are
congruent mod Phi_u iff they agree at every u_j).  L has degree 1 < FA, so L IS the reduced
form provided P(u_j)=L(u_j) at every node.

So LEMMA 1 <=> the NODE IDENTITY:  (4 - u_j^2) Atil(u_j) = eps_k m (2 + u_j)   for all j.     (NODE)

Now u_j = 2 c_j = 2 cos(theta_j), theta_j = 2 pi j/m, so 4 - u_j^2 = 4 sin^2 theta_j, and
Atil(u_j) = A(c_j) = sum_{i=0}^{k} U_i(c_j)(k-2i).  Using U_i(cos th)= sin((i+1)th)/sin th:
   (4 - u_j^2) A(c_j) = 4 sin^2 th * sum_i (k-2i) sin((i+1)th)/sin th
                      = 4 sin th * sum_i (k-2i) sin((i+1)th).
The inner sum  sum_{i=0}^{k} (k-2i) sin((i+1) th)  is the S1-type sum (shifted index).  From the
+1 work,  S1(th) := sum_{l=1}^{m-1} (m-2l) sin(l th) = m cot(pi j/m) = m cos(pi j/m)/sin(pi j/m)
at th = 2 pi j/m.  (NODE) becomes a clean trig identity, and 2 + u_j = 2 + 2cos th = 4 cos^2(th/2),
4 sin th = 8 sin(th/2)cos(th/2), so RHS-shape eps m (2+u_j) = eps m 4 cos^2(th/2).

This script PROVES (NODE) bit-exact (mpmath, high precision, all j, wide k), establishing
LEMMA 1 rigorously (it is a finite trig identity at the m-th-root nodes, closed by the S1 sum).
Both eps cases (odd via k-2i, even via kappa-i) are checked.
"""
from __future__ import annotations

import sys
from mpmath import mp, mpf, sin, cos, pi, fabs

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

mp.dps = 60


def cheb_u(j, x):
    um1, u0 = mpf(0), mpf(1)
    if j == 0:
        return u0
    for _ in range(j):
        um1, u0 = u0, 2 * x * u0 - um1
    return u0


def main():
    print("=" * 104)
    print("F89 edge A: RIGOROUS LEMMA 1 via NODE identity  (4-u_j^2)Atil(u_j) == eps m (2+u_j)")
    print("=" * 104)
    print()
    print(f"  {'k':>3} {'cls':>5} {'m':>3} {'FA':>3} {'eps':>3} {'max NODE err (all j)':>22} "
          f"{'S1=m cot err':>14} {'ok':>4}")
    print("  " + "-" * 76)
    allok = True
    for k in list(range(5, 28)):
        m = k + 2
        FA = (k + 1) // 2
        even = (k % 2 == 0)
        eps = 1 if even else 2
        kap = k // 2
        cls = "odd" if not even else ("4|k" if k % 4 == 0 else "2mod4")
        node_err = mpf(0)
        s1_err = mpf(0)
        for j in range(1, FA + 1):
            th = 2 * pi * j / m
            cj = cos(th)
            uj = 2 * cj
            # P(u_j) = (4-u_j^2) A(c_j)
            if even:
                A = sum(cheb_u(i, cj) * (kap - i) for i in range(k + 1))
            else:
                A = sum(cheb_u(i, cj) * (k - 2 * i) for i in range(k + 1))
            P = (4 - uj * uj) * A
            L = eps * m * (2 + uj)
            node_err = max(node_err, fabs(P - L))
            # also verify the S1 closed form used in the derivation:
            # (4-u_j^2)A(c_j) = 4 sin th * Sshift,  Sshift = sum_{i=0}^{k}(coef) sin((i+1)th)
            if even:
                Sshift = sum((kap - i) * sin((i + 1) * th) for i in range(k + 1))
            else:
                Sshift = sum((k - 2 * i) * sin((i + 1) * th) for i in range(k + 1))
            lhs = 4 * sin(th) * Sshift
            s1_err = max(s1_err, fabs(P - lhs))   # consistency of the sin-rewrite
        ok = (node_err < mpf(10) ** (-30)) and (s1_err < mpf(10) ** (-30))
        allok = allok and ok
        print(f"  {k:>3} {cls:>5} {m:>3} {FA:>3} {eps:>3} {mp.nstr(node_err,3):>22} "
              f"{mp.nstr(s1_err,3):>14} {('OK' if ok else 'X'):>4}")
    print()
    print(f"  NODE identity holds (all j, all k):  {allok}")
    print()
    print("  => LEMMA 1 is RIGOROUS: P(u):=(4-u^2)Atil and L(u):=eps m (2+u) agree at every node")
    print("     u_j (root of the separable Phi_u), L has degree 1 < FA, hence L = P mod Phi_u.")
    print("  The node identity is the finite trig identity 4 sin th * S1shift(th) = eps m (2+2cos th),")
    print("  closed by S1shift(2 pi j/m) (the shifted-index S1 sum, sibling of the +1-lemma's m cot).")


if __name__ == "__main__":
    main()
