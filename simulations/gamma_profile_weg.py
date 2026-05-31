#!/usr/bin/env python3
"""Is the closed-form slowest rate THE path, or the ruler that shows when you leave it?

The path (the loop, the dynamics, the slosh, the memory) is Q-driven. The closed form
(10+8e)/9 for the edge family is Q-INDEPENDENT. So:
  1. Is the edge Q-independence EXACT (a genuine dynamics-free zone -> the rate is NOT the path)?
  2. Does the trajectory still slosh INSIDE the closed envelope at the edge (the path runs inside
     the ruler -> the closed form bounds the path, it is not the path)?
  3. At the center, is the rate itself Q-DEPENDENT (the path has invaded the envelope)?

If (1)+(2)+(3) hold, the closed form is the ruler/wall of the channel, and leaving it (the
Q-driven slosh / the Q-dependence) is the path.
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


def L_dimless(N, Q, profile):
    d = 2 ** N
    Id = np.eye(d)
    H1 = H_xy_unit(N)
    L = -1j * Q * (np.kron(Id, H1) - np.kron(H1.T, Id))
    for l in range(N):
        Zl = op_at(N, l, Z)
        L += profile[l] * (np.kron(Zl, Zl) - np.kron(Id, Id))
    return L


def slowest_rate(N, Q, profile):
    w = np.linalg.eigvals(L_dimless(N, Q, profile))
    return -float(np.max(w[np.abs(w) > 1e-7].real))


def n0_trajectory(N, Q, profile, taus):
    d = 2 ** N
    psi0 = np.zeros(d, complex)
    idx = 1 << (N - 1)  # site 0 excited (leftmost factor)
    psi0[idx] = 1.0
    rho0 = np.outer(psi0, psi0.conj())
    vec0 = rho0.flatten(order="F")
    n0 = (np.eye(d) - op_at(N, 0, Z)) / 2.0
    w, V = np.linalg.eig(L_dimless(N, Q, profile))
    c = np.linalg.solve(V, vec0)
    out = []
    for t in taus:
        rho = (V @ (np.exp(w * t) * c)).reshape(d, d, order="F")
        out.append(float(np.trace(n0 @ rho).real))
    return np.array(out)


def main():
    N = 5
    edge = [0.25, 1.5, 1.5, 1.5, 0.25]      # [e,x,x,x,e], e=0.25 -> closed form (10+8e)/9 = 12/9
    center = [1.6, 1.6, 0.6, 1.6, 0.6]      # quiet-ish center, sum 6 -> renormalize
    center = list(np.array(center) * 5 / sum(center))

    print(f"N={N}. Edge profile {np.round(edge,3).tolist()} (sum {sum(edge):.2f}); closed (10+8*0.25)/9 = {(10+8*0.25)/9:.6f}")
    print(f"      Center profile {np.round(center,3).tolist()} (sum {sum(center):.2f})\n")

    print("  (1) Is the rate Q-independent (dynamics-free closed zone) at the edge, but not the center?")
    print(f"  {'Q':>8}  {'edge rate':>12}  {'center rate':>12}")
    for Q in [1.5, 5.0, 20.0, 100.0, 1000.0]:
        print(f"  {Q:>8.1f}  {slowest_rate(N, Q, edge):>12.8f}  {slowest_rate(N, Q, center):>12.8f}")

    print("\n  (2) Does the trajectory slosh INSIDE the closed envelope at the edge (Q=20)?")
    r = (10 + 8 * 0.25) / 9.0
    sigma = sum(edge)  # closed rate is per unit tau already (eig of L); envelope = exp(-r*tau)
    taus = np.linspace(0, 8, 80)
    n0 = n0_trajectory(N, 20.0, edge, taus)
    target = 1.0 / N
    dev = n0 - target
    env = dev[0] * np.exp(-r * taus)  # closed-form envelope (slowest mode)
    # count sign changes of (n0 - envelope): wiggle inside the envelope = the slosh
    inside = dev - env
    turns = int(np.sum(np.abs(np.diff(np.sign(np.diff(n0)))) > 0))
    print(f"     slowest closed rate r = {r:.4f}; trajectory n0(tau) oscillation turns = {turns}")
    print(f"     {'tau':>5}  {'n0-1/N':>9}  {'envelope':>9}  {'inside(slosh)':>13}")
    for i in range(0, len(taus), 8):
        print(f"     {taus[i]:>5.2f}  {dev[i]:>9.4f}  {env[i]:>9.4f}  {inside[i]:>13.4f}")
    print("     (if 'inside' wiggles around 0 while envelope falls smoothly: the path sloshes inside the ruler)")


if __name__ == "__main__":
    main()
