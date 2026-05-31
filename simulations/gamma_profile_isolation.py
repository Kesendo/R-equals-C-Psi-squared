#!/usr/bin/env python3
"""The isolation tool: where does the shortcut begin and end?

The Q-dependence of the slowest rate is the OBSERVABLE that reveals spectral isolation of the
slow mode. The DIRECT measure is the spectral gap from the slowest rate to the next: large gap =
isolated mode (J can only rotate it, the rate is Q-independent, the shortcut); small gap = a
neighbour J mixes with at second order (the rate drifts with Q, the loop). We scan a path in
profile-space (uniform -> quiet-edge flat-bulk) and read, per step:
  - delta = rate(Q=1000) - rate(Q=1.5)   the tool reading (0 = isolated/shortcut, >0 = loop)
  - gap   = (2nd slowest rate) - (slowest rate) at high Q   the isolation itself
Where delta lifts off 0 and the gap collapses: there the isolation ends, the loop begins.
"""
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
            t = np.array([[1]], complex)
            for i in range(N):
                t = np.kron(t, P if i in (b, b + 1) else I2)
            H += t
    return H


def L_dimless(N, Q, profile, H1):
    d = 2 ** N
    Id = np.eye(d)
    L = -1j * Q * (np.kron(Id, H1) - np.kron(H1.T, Id))
    for l in range(N):
        Zl = op_at(N, l, Z)
        L += profile[l] * (np.kron(Zl, Zl) - np.kron(Id, Id))
    return L


def rates(N, Q, profile, H1):
    w = np.linalg.eigvals(L_dimless(N, Q, profile, H1))
    rs = np.sort(-w[np.abs(w) > 1e-7].real)  # ascending: slowest first
    return rs


def scan(N, anchor_a, anchor_b, label):
    H1 = H_xy_unit(N)
    print(f"\n  path: {label}")
    print(f"  {'s':>5}  {'rate Q=1.5':>10}  {'rate Q=1000':>11}  {'delta':>7}  {'gap@1000':>9}  {'isolation'}")
    for s in np.linspace(0, 1, 11):
        p = (1 - s) * np.array(anchor_a, float) + s * np.array(anchor_b, float)
        p = p * N / p.sum()
        rlo = rates(N, 1.5, p, H1)[0]
        rhi = rates(N, 1000.0, p, H1)
        r0, gap = rhi[0], rhi[1] - rhi[0]
        delta = r0 - rlo
        tag = "isolated (shortcut)" if abs(delta) < 1e-4 else ("edge of loop" if abs(delta) < 0.02 else "loop")
        print(f"  {s:>5.2f}  {rlo:>10.5f}  {r0:>11.5f}  {delta:>+7.4f}  {gap:>9.5f}  {tag}")


def main():
    N = 5
    uniform = [1, 1, 1, 1, 1]
    edge_flat = [0.25, 1.5, 1.5, 1.5, 0.25]   # quiet flat-bulk edge (a loop profile)
    peaked_V = [0.25, 0.75, 3.0, 0.75, 0.25]  # quiet edges + peaked centre (an isolated profile)
    print(f"N={N}. The isolation tool: delta (Q-drift) reads it; gap is the isolation itself.")
    scan(N, uniform, edge_flat, "uniform  ->  quiet-edge flat-bulk")
    scan(N, uniform, peaked_V, "uniform  ->  quiet-edge peaked-centre")
    scan(N, peaked_V, edge_flat, "peaked-centre  ->  flat-bulk (both quiet edges 0.25)")


if __name__ == "__main__":
    main()
