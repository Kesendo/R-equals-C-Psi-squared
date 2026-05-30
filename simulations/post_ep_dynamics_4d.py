#!/usr/bin/env python3
"""The post-EP dynamics in the loop, as the 4D structure Q x bond x t.

The static F86 K-curve collapses the time axis (it is a t-integral / peak). The loop KEEPS t.
For each Q across the crossover and each site, time-evolve a single excitation (on site 0) and
track its occupation <n_site>(tau) (tau = gamma*t, the dimensionless tick; Q = J/gamma is the only
knob, so gamma drops out of L'). Below the crossover the excitation diffuses to uniform (overdamped,
forgetting); above it hops coherently and sloshes back (underdamped, remembering) before it finally
decoheres. Printed as the nested  Q -> site -> t-curve  tree , the 4D representation the Object
Manager renders (here as ASCII sparklines; each leaf would be a Curve payload drawn by --draw).
(The earlier per-bond <XX+YY> observable was Hamiltonian-invariant , Q-blind , so it showed only
the dephasing decay; the single-excitation occupation carries the rotation.)
"""
import sys

import numpy as np

I2 = np.eye(2)
X = np.array([[0, 1], [1, 0]], complex)
Y = np.array([[0, -1j], [1j, 0]], complex)
Z = np.array([[1, 0], [0, -1]], complex)
BLOCKS = "▁▂▃▄▅▆▇█"


def op_at(N, s, P):
    o = np.array([[1]], complex)
    for i in range(N):
        o = np.kron(o, P if i == s else I2)
    return o


def bond_op(N, b, P, Qop):
    o = np.array([[1]], complex)
    for i in range(N):
        o = np.kron(o, P if i == b else (Qop if i == b + 1 else I2))
    return o


def H_xy_unit(N):
    """H/J for the XY chain: sum_b (X_b X_{b+1} + Y_b Y_{b+1})."""
    H = np.zeros((2 ** N, 2 ** N), complex)
    for b in range(N - 1):
        H += bond_op(N, b, X, X) + bond_op(N, b, Y, Y)
    return H


def liouvillian_dimensionless(N, Q):
    """L' with Q = J/gamma the only knob (gamma factored out; evolve in tau = gamma*t).
    vec convention: column-stack (flatten 'F'), so [H,.] <-> I (x) H - H^T (x) I."""
    d = 2 ** N
    Id = np.eye(d)
    H1 = H_xy_unit(N)
    L = -1j * Q * (np.kron(Id, H1) - np.kron(H1.T, Id))
    for l in range(N):
        Zl = op_at(N, l, Z)
        L += np.kron(Zl, Zl) - np.kron(Id, Id)
    return L


def trajectory(N, Q, ops, taus):
    """<O>(tau) for each operator in `ops`, starting from |+>^N, via eig once + propagation."""
    d = 2 ** N
    one, zero = np.array([0, 1], complex), np.array([1, 0], complex)
    psi0 = np.array([1], complex)
    for i in range(N):
        psi0 = np.kron(psi0, one if i == 0 else zero)   # |1 0 ... 0>: a single excitation on site 0
    rho0 = np.outer(psi0, psi0.conj())
    vec0 = rho0.flatten(order="F")
    w, V = np.linalg.eig(liouvillian_dimensionless(N, Q))
    c = np.linalg.solve(V, vec0)
    out = {k: [] for k in range(len(ops))}
    for tau in taus:
        vt = V @ (np.exp(w * tau) * c)
        rho = vt.reshape(d, d, order="F")
        for k, O in enumerate(ops):
            out[k].append(float(np.trace(O @ rho).real))
    return out


def sparkline(ys):
    ys = np.asarray(ys, float)
    lo, hi = ys.min(), ys.max()
    if hi - lo < 1e-12:
        return BLOCKS[0] * len(ys)
    idx = np.clip(((ys - lo) / (hi - lo) * (len(BLOCKS) - 1)).round().astype(int), 0, len(BLOCKS) - 1)
    return "".join(BLOCKS[i] for i in idx)


def n_turns(ys):
    """Count local extrema of the detrended signal , a rough 'how much it oscillates'."""
    ys = np.asarray(ys, float)
    d = np.diff(ys)
    return int(np.sum(np.abs(np.diff(np.sign(d))) > 0))


def main():
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 4         # N is the knob: usage  script.py [N]
    taus = np.linspace(0.0, 6.0, 60)
    d = 2 ** N
    site_class = lambda b: "edge" if b in (0, N - 1) else "bulk"
    ops = [(np.eye(d) - op_at(N, b, Z)) / 2.0 for b in range(N)]       # n_b = (I - Z_b)/2 per site
    target = 1.0 / N
    print(f"post-EP dynamics  (N={N} XY chain, Z-dephasing; observable <n_site>(tau) for a single")
    print(f"  excitation on site 0; tau = gamma*t in [0,6]; Q = J/gamma the dimensionless knob)")
    print(f"  TARGET every site relaxes to = 1/N = {target:.4f}  (the fully-forgotten, equipartitioned")
    print(f"  state; it equals 1/4 only at N=4 , NOT the CΨ=1/4 fold, which is N-independent)")
    for Q in [0.5, 1.0, 1.5, 2.5]:
        regime = "overdamped (diffuses, forgets)" if Q < 1.0 else "underdamped (hops, remembers)"
        print(f"\n  Q = {Q:.2f}   ({regime})   [target 1/N = {target:.3f}]")
        traj = trajectory(N, Q, ops, taus)
        for b in range(N):
            ys = traj[b]
            turns = n_turns(ys)
            tag = "monotone" if turns <= 1 else f"oscillates ({turns} turns)"
            print(f"    site {b} ({site_class(b):4s})  {sparkline(ys)}  {ys[0]:.2f}->{ys[-1]:.2f}   {tag}")


if __name__ == "__main__":
    main()
