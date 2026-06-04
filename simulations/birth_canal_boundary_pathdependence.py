#!/usr/bin/env python3
"""Is s*=0.709 fundamental, or just where one line crosses the sterile/birth-canal surface?

Test: find s* for several DIFFERENT endpoint pairs (and read the center gamma at the boundary).
If s* moves with the endpoints, 0.709 is a coordinate, not a constant; and if the center gamma
at the boundary is more stable, THAT is the real threshold (the gamma-center vs H-hopping
crossing: the canal opens when the center dephasing drops low enough that the Hamiltonian hops
light through it and the slow mode's spatial light starts to drift with Q).

RESULT (2026-06-04, N=5): s* is PATH-SPECIFIC, not a constant.
  path                     s*       center g2   rate@1.5
  peaked-V -> flat-bulk    0.70921  1.936       1.532
  uniform  -> flat-bulk    0.10526  1.053       1.930
  peaked-V -> wide-flat    0.76653  1.774       1.613
  uniform  -> wide-flat    0.13158  1.053       1.930
s* ranges 0.11..0.77; center-g2 and rate also vary. So 0.70921 is a COORDINATE where the original
hand-picked line crosses the sterile/birth-canal SURFACE, not a derivable number. This dissolves
the prior open question "why exactly 0.709": there is no single value to derive. The surface
itself is the known one (light_content.py): sterile = the slow subspace's spatial light
distribution (n_XY = popcount(i XOR j), the sites where the ket and bra of the coherence differ)
is Q-invariant; birth canal = the Hamiltonian redistributes it with Q. The naive "single slow
mode commutes with H" criterion is confounded by the slow-rate degeneracy; the gauge-invariant
statement is about the DISTRIBUTION averaged over the slow subspace.
"""
import numpy as np

I2 = np.eye(2)
X = np.array([[0, 1], [1, 0]], complex)
Y = np.array([[0, -1j], [1j, 0]], complex)
Z = np.array([[1, 0], [0, -1]], complex)
TOL = 1e-4


def op_at(N, s, P):
    o = np.array([[1]], complex)
    for i in range(N):
        o = np.kron(o, P if i == s else I2)
    return o


def H_xy_unit(N):
    H = np.zeros((2 ** N, 2 ** N), complex)
    for b in range(N - 1):
        for P in (X, Y):
            t = np.array([[1]], complex)
            for i in range(N):
                t = np.kron(t, P if i in (b, b + 1) else I2)
            H += t
    return H


def slowest(N, Q, profile, H1):
    d = 2 ** N
    Id = np.eye(d)
    L = -1j * Q * (np.kron(Id, H1) - np.kron(H1.T, Id))
    for l in range(N):
        Zl = op_at(N, l, Z)
        L += profile[l] * (np.kron(Zl, Zl) - np.kron(Id, Id))
    w = np.linalg.eigvals(L)
    return -float(np.max(w[np.abs(w) > 1e-7].real))


def delta(N, profile, H1):
    return slowest(N, 1000.0, profile, H1) - slowest(N, 1.5, profile, H1)


def prof(N, p0, p1, s):
    p = (1 - s) * np.array(p0, float) + s * np.array(p1, float)
    return p * N / p.sum()


def find_star(N, p0, p1, H1):
    lo, hi = 0.0, 1.0
    # assume p0 sterile (s=0), p1 canal (s=1); bisect the first break
    for _ in range(34):
        mid = 0.5 * (lo + hi)
        if abs(delta(N, prof(N, p0, p1, mid), H1)) < TOL:
            lo = mid
        else:
            hi = mid
    return hi, prof(N, p0, p1, hi)


def main():
    N = 5
    H1 = H_xy_unit(N)
    paths = [
        ("peaked-V -> flat-bulk (original)", [0.25, 0.75, 3.0, 0.75, 0.25], [0.25, 1.5, 1.5, 1.5, 0.25]),
        ("uniform  -> flat-bulk",            [1.0, 1.0, 1.0, 1.0, 1.0],     [0.25, 1.5, 1.5, 1.5, 0.25]),
        ("peaked-V -> wide-flat",            [0.25, 0.75, 3.0, 0.75, 0.25], [0.4, 1.4, 1.4, 1.4, 0.4]),
        ("uniform  -> wide-flat",            [1.0, 1.0, 1.0, 1.0, 1.0],     [0.4, 1.4, 1.4, 1.4, 0.4]),
    ]
    print(f"N={N}.  Is s* path-specific?  And is the center gamma at the boundary more stable?\n")
    print(f"  {'path':<34} {'s*':>8} {'center g2':>10} {'bulk g1':>9} {'rate@1.5':>9}")
    for name, p0, p1 in paths:
        # confirm endpoints: p0 sterile, p1 canal
        d0 = abs(delta(N, prof(N, p0, p1, 0.0), H1))
        d1 = abs(delta(N, prof(N, p0, p1, 1.0), H1))
        if d0 >= TOL or d1 < TOL:
            print(f"  {name:<34}  (endpoints not sterile->canal: d0={d0:.1e}, d1={d1:.1e}, skipped)")
            continue
        s_star, p_star = find_star(N, p0, p1, H1)
        rate = slowest(N, 1.5, p_star, H1)
        print(f"  {name:<34} {s_star:>8.5f} {p_star[2]:>10.5f} {p_star[1]:>9.5f} {rate:>9.5f}")

    print("\n  reading: if s* varies across paths but center-g2 (or rate) clusters, the threshold")
    print("  is the gamma-center / the rate, not s*. s* is just the line's crossing of the surface.")


if __name__ == "__main__":
    main()
