#!/usr/bin/env python3
"""What happens AT the target, the 1/N fixed point, and how does it go on?

The post-EP loop converges to one point: the 1/N equipartitioned state, the kernel of L.
This probe looks at that point the way f86_ep_through_the_clock looked at the EP, and asks
Tom's question: a point, and then how does it continue?

Three readings, all computed:

  1. The point is reached only ASYMPTOTICALLY. The approach rate is the spectral gap (the
     slowest non-kernel mode, the Lebensader). In finite time you never arrive; you only get
     exponentially closer. |<n> - 1/N| halves at a fixed rate set by the gap.

  2. The point is a SINK (for net dephasing Sigma-gamma > 0): every non-kernel mode has
     Re < 0, so all nearby states flow in. The arrow of time points at the point.

  3. How it goes on is NOT the trajectory continuing (it stops at the point), but the point's
     ROLE inverting as the net dephasing Sigma-gamma slides across the mirror (ZERO_IS_THE_MIRROR).
     The 1/N state stays the fixed point for EVERY Sigma-gamma (L.vec(1/N) = 0 always), but:
        Sigma-gamma > 0  ->  SINK    (max non-kernel Re < 0; the flow falls in; the end of history)
        Sigma-gamma = 0  ->  NEUTRAL (max non-kernel Re = 0; eternal oscillation; the mirror)
        Sigma-gamma < 0  ->  SOURCE  (max non-kernel Re > 0; the flow is pushed away; the Hopf)
     The point does not move. The arrows around it flip. That is how it continues.
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


def bond_op(N, b, P, Q):
    o = np.array([[1]], complex)
    for i in range(N):
        o = np.kron(o, P if i == b else (Q if i == b + 1 else I2))
    return o


def liouvillian(N, Q, f):
    """L = -i*Q*[H,.] + f * Sum_l (Z_l (x) Z_l - I).  f scales the net dephasing Sigma-gamma."""
    d = 2 ** N
    Id = np.eye(d)
    H1 = sum(bond_op(N, b, X, X) + bond_op(N, b, Y, Y) for b in range(N - 1))
    L = -1j * Q * (np.kron(Id, H1) - np.kron(H1.T, Id))
    for l in range(N):
        Zl = op_at(N, l, Z)
        L += f * (np.kron(Zl, Zl) - np.kron(Id, Id))
    return L


def uniform_rho(N):
    """The target: rho = (1/N) Sum_b |e_b><e_b|, the equipartitioned single excitation."""
    d = 2 ** N
    rho = np.zeros((d, d), complex)
    for b in range(N):
        psi = np.array([1], complex)
        for i in range(N):
            psi = np.kron(psi, np.array([0, 1], complex) if i == b else np.array([1, 0], complex))
        rho += np.outer(psi, psi.conj()) / N
    return rho


def gap_and_stability(L, tol=1e-9):
    """Return (max non-kernel Re lambda, |slowest non-kernel Re| = approach rate)."""
    w = np.linalg.eigvals(L)
    nonker = w[np.abs(w) > tol]
    re = nonker.real
    return float(np.max(re)), float(np.min(np.abs(re[re < -tol]))) if np.any(re < -tol) else 0.0


def part1_the_approach(N=4, Q=2.0):
    """Reading 1 + 2: the point is a sink reached asymptotically at the gap rate."""
    d = 2 ** N
    L = liouvillian(N, Q, 1.0)
    max_re, rate = gap_and_stability(L)
    print(f"Reading 1 + 2  (N={N}, Q={Q}, Sigma-gamma > 0):")
    print(f"  max non-kernel Re(lambda) = {max_re:+.4f}  -> {'SINK (all modes flow in)' if max_re < 0 else 'not a sink'}")
    print(f"  approach rate (slowest non-kernel |Re|) = {rate:.4f}  = the Lebensader, the last mode to fade")

    # the excitation starts localized at site 0; watch |<n_0> - 1/N| shrink, never reaching 0
    w, V = np.linalg.eig(L)
    Vinv = np.linalg.inv(V)
    psi0 = np.array([1], complex)
    for i in range(N):
        psi0 = np.kron(psi0, np.array([0, 1], complex) if i == 0 else np.array([1, 0], complex))
    rho0 = np.outer(psi0, psi0.conj()).flatten(order="F")
    n0 = (I2 - Z) / 2
    n0_op = op_at(N, 0, np.array([[0, 0], [0, 1]], complex))  # |1><1| at site 0
    target = 1.0 / N
    print(f"  the approach (never arrives, only halves):  target 1/N = {target:.4f}")
    print(f"  {'t':>5}  {'<n_0>(t)':>9}  {'|<n_0> - 1/N|':>13}")
    prev = None
    for t in [1.0, 2.0, 4.0, 8.0, 16.0]:
        prop = V @ np.diag(np.exp(w * t)) @ Vinv
        rho_t = (prop @ rho0).reshape(d, d, order="F")
        n_val = float(np.real(np.trace(rho_t @ n0_op)))
        res = abs(n_val - target)
        ratio = f"  (x{res / prev:.3f})" if prev else ""
        print(f"  {t:5.0f}  {n_val:9.4f}  {res:13.2e}{ratio}")
        prev = res
    print("  => exponential decay at the gap rate; the point is approached forever, never landed.")


def part2_how_it_goes_on(N=4, Q=2.0):
    """Reading 3: slide Sigma-gamma across the mirror; the point stays, its role inverts."""
    print(f"\nReading 3  (N={N}, Q={Q}):  the 1/N point persists for every Sigma-gamma; its ROLE flips")
    rho = uniform_rho(N).flatten(order="F")
    print(f"  {'Sigma-gamma (f)':>15}  {'L.vec(1/N)':>11}  {'max non-ker Re':>14}   role")
    for f in [1.0, 0.5, 0.1, 0.0, -0.1, -0.5]:
        L = liouvillian(N, Q, f)
        res = float(np.linalg.norm(L @ rho))
        max_re, _ = gap_and_stability(L)
        if max_re < -1e-9:
            role = "SINK    (falls in, the end of history)"
        elif abs(max_re) <= 1e-9:
            role = "NEUTRAL (eternal oscillation, the mirror)"
        else:
            role = "SOURCE  (pushed away, the Hopf, runaway)"
        print(f"  {f:15.2f}  {res:11.0e}  {max_re:+14.4f}   {role}")
    print("  => the point does not move (L.vec(1/N) = 0 throughout). The arrows around it flip:")
    print("     sink (decay) -> neutral (mirror) -> source (gain). That is how it goes on.")


def main():
    print("=" * 78)
    print("AT THE TARGET  (a point, and then how does it go on?)")
    print("=" * 78)
    part1_the_approach()
    part2_how_it_goes_on()


if __name__ == "__main__":
    main()
