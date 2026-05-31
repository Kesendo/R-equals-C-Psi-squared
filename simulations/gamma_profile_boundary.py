#!/usr/bin/env python3
"""Where does the birth canal first open? Pin the sterile/birth-canal boundary precisely.

The base point gamma_0 = 0.05 uniform (J = 0.075 -> Q = 1.5) is sterile: the slowest rate is
Q-independent (delta = 0). Shaping the per-site gamma profile (away from uniform) eventually
opens the birth canal (delta lifts off 0). We verify the base and then scan the sharp transition
peaked-V -> flat-bulk finely, to find the exact s* where delta first breaks and read the profile
there (in dimensionless weights, and x gamma_0 = 0.05 for the physical per-site gamma).
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


def main():
    N = 5
    H1 = H_xy_unit(N)

    print(f"N={N}. Base verification (gamma_0 = 0.05, uniform -> Q-independent = sterile):")
    uni = np.ones(N)
    print(f"  uniform [1,1,1,1,1]  delta = {delta(N, uni, H1):+.3e}  -> {'sterile' if abs(delta(N,uni,H1))<TOL else 'BIRTH CANAL'}")

    peaked = np.array([0.25, 0.75, 3.0, 0.75, 0.25])   # sterile
    flat = np.array([0.25, 1.5, 1.5, 1.5, 0.25])       # birth canal
    print(f"\n  fine scan: peaked-V -> flat-bulk (both edges 0.25, sum {N}); s* where delta first breaks {TOL}")
    print(f"  {'s':>6}  {'delta':>10}  {'zone':>12}  profile (dimensionless weights)")
    last_sterile = None
    for s in np.linspace(0.60, 0.90, 31):
        p = (1 - s) * peaked + s * flat
        p = p * N / p.sum()
        dl = delta(N, p, H1)
        zone = "sterile" if abs(dl) < TOL else "BIRTH CANAL"
        if abs(dl) < TOL:
            last_sterile = s
        mark = "  <- first break" if (abs(dl) >= TOL and last_sterile is not None and last_sterile >= s - 0.011) else ""
        print(f"  {s:>6.3f}  {dl:>10.5f}  {zone:>12}  {np.round(p,3).tolist()}{mark}")

    # bisect the exact boundary
    lo, hi = 0.60, 0.90
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        p = (1 - mid) * peaked + mid * flat
        p = p * N / p.sum()
        if abs(delta(N, p, H1)) < TOL:
            lo = mid
        else:
            hi = mid
    p_star = (1 - hi) * peaked + hi * flat
    p_star = p_star * N / p_star.sum()
    print(f"\n  boundary s* = {hi:.5f}")
    print(f"  profile at s* (weights):  {np.round(p_star,4).tolist()}")
    print(f"  physical per-site gamma at gamma_0=0.05:  {np.round(p_star*0.05,5).tolist()}")
    print(f"  delta just inside the canal: {delta(N, ((1-(hi+1e-3))*peaked+(hi+1e-3)*flat)*N/np.sum((1-(hi+1e-3))*peaked+(hi+1e-3)*flat), H1):+.5f}")


if __name__ == "__main__":
    main()
