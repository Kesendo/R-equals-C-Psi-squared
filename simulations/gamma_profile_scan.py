#!/usr/bin/env python3
"""Scan the slowest non-kernel relaxation rate of the single-excitation flow under shaped
gamma profiles at fixed total Sum-gamma = N. Question: does the rate take closed-form values
as the dephasing shape is moved (the pure shape effect, isolated from the total)?

Rate = -max{ Re lambda : |lambda| > tol } of L'(N,Q,profile) = -iQ[H_unit, .] + sum_l p_l (Z_l . Z_l - .),
with the profile normalized to sum N. Eigenvalues are basis-independent, so this matches the
C# PostEpFlowField.SlowestRate exactly (uniform N=5 Q=20 -> 2.0; V-shape -> 1.0).
"""
import sys

import numpy as np

I2 = np.eye(2)
X = np.array([[0, 1], [1, 0]], complex)
Y = np.array([[0, -1j], [1j, 0]], complex)
Z = np.array([[1, 0], [0, -1]], complex)


def op_at(N, s, P):
    o = np.array([[1]], complex)
    for i in range(N):
        o = np.kron(o, P if i == s else I2)
    return o


def H_xy_unit(N):
    H = np.zeros((2 ** N, 2 ** N), complex)
    for b in range(N - 1):
        for P in (X, Y):
            term = np.array([[1]], complex)
            for i in range(N):
                term = np.kron(term, P if i in (b, b + 1) else I2)
            H += term
    return H


def slowest_rate(N, Q, profile):
    d = 2 ** N
    Id = np.eye(d)
    H1 = H_xy_unit(N)
    L = -1j * Q * (np.kron(Id, H1) - np.kron(H1.T, Id))
    for l in range(N):
        Zl = op_at(N, l, Z)
        L += profile[l] * (np.kron(Zl, Zl) - np.kron(Id, Id))
    w = np.linalg.eigvals(L)
    nonkernel = w[np.abs(w) > 1e-7]
    return -float(np.max(nonkernel.real))


def norm(p, N):
    p = np.array(p, float)
    return p * N / p.sum()


def main():
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    sqrt5m1 = np.sqrt(5) - 1
    profiles = {
        "uniform":          norm([1, 1, 1, 1, 1], N),
        "V mild":           norm([0.5, 0.8, 2.4, 0.8, 0.5], N),
        "V (edges .25)":    norm([0.25, 0.75, 3.0, 0.75, 0.25], N),
        "V strong":         norm([0.1, 0.5, 3.8, 0.5, 0.1], N),
        "inverse-V":        norm([2.4, 0.6, 0.2, 0.6, 2.4], N),
        "edges-quiet flat": norm([0.2, 1, 1.6, 1, 0.2], N),
        "center-quiet flat":norm([1.6, 1, 0.2, 1, 1.6], N),
    }
    Qs = [1.5, 2.0, 5.0, 20.0]
    print(f"Slowest non-kernel rate, N={N}, all profiles normalized to Sum-gamma = {N}")
    print(f"  (clean refs: 2, 1, sqrt(5)-1={sqrt5m1:.4f}, sqrt(2)={np.sqrt(2):.4f}, golden phi={(1+np.sqrt(5))/2:.4f})\n")
    head = "  " + f"{'profile':>18}  " + "  ".join(f"Q={q:<5g}" for q in Qs)
    print(head)
    for name, p in profiles.items():
        row = "  " + f"{name:>18}  " + "  ".join(f"{slowest_rate(N, q, p):>7.4f}" for q in Qs)
        print(row)

    print(f"\n  edge-depth scan at Q=20 (profile [e, x, x, x, e], 2e+3x={N}, e from 1 down):")
    print(f"  {'edge e':>7}  {'bulk x':>7}  {'rate':>8}  {'rate^2':>8}  {'2*e':>6}")
    for e in [1.0, 0.8, 0.6, 0.5, 0.4, 0.3, 0.25, 0.2, 0.1, 0.05]:
        x = (N - 2 * e) / 3.0
        p = [e, x, x, x, e]
        r = slowest_rate(N, 20.0, p)
        print(f"  {e:>7.3f}  {x:>7.3f}  {r:>8.4f}  {r*r:>8.4f}  {2*e:>6.3f}")


if __name__ == "__main__":
    main()
