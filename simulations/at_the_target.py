#!/usr/bin/env python3
"""What happens at the sector-uniform stationary state under a formal gain continuation?

For canonical Q>0, a localized single excitation approaches the sector-uniform state. Globally,
the number-conserving Liouvillian has an (N+1)-dimensional semisimple kernel: one stationary state
per excitation-number sector. This probe follows the one-excitation representative only.

Three readings, all computed:

  1. The sector state is reached only ASYMPTOTICALLY. The printed global non-kernel spectral
     abscissa is not identified with the population-visible approach rate: its winning mode may
     have zero overlap with this preparation/readout.

  2. It attracts the transient complement at fixed one-excitation trace. Kernel-direction
     perturbations remain stationary and generally approach another sector mixture.

  3. A formal common-scalar continuation f of one positive per-site profile changes the transient
     spectral abscissa. f<0 is non-CP inverse dephasing that amplifies coherences, not an energy-gain
     bath and not a Hopf certificate. This says nothing universal about arbitrary mixed-sign profiles.
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
    """L = -i*(Q/2)*[sum(XX+YY),.] + f*sum_l(Z_l(.)Z_l-I)."""
    d = 2 ** N
    Id = np.eye(d)
    H1 = sum(bond_op(N, b, X, X) + bond_op(N, b, Y, Y) for b in range(N - 1))
    L = -1j * (Q / 2.0) * (np.kron(Id, H1) - np.kron(H1.T, Id))
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
    """Return global non-kernel spectral-abscissa data; no preparation overlap is applied."""
    w = np.linalg.eigvals(L)
    nonker = w[np.abs(w) > tol]
    re = nonker.real
    return float(np.max(re)), float(np.min(np.abs(re[re < -tol]))) if np.any(re < -tol) else 0.0


def part1_the_approach(N=4, Q=2.0):
    """Sector relaxation plus the separately labelled global spectral rate."""
    d = 2 ** N
    L = liouvillian(N, Q, 1.0)
    max_re, rate = gap_and_stability(L)
    print(f"Reading 1 + 2  (N={N}, Q={Q}, common per-site f=1):")
    print(f"  max non-kernel Re(lambda) = {max_re:+.4f}  (global transient spectral abscissa)")
    print(f"  global slowest non-kernel |Re| = {rate:.4f}; not asserted to be population-visible")

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
    print(f"  fixed-one-excitation-sector approach: target 1/N = {target:.4f}")
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
    print("  => the sampled sector residual decreases; no global-gap equality is inferred.")


def part2_how_it_goes_on(N=4, Q=2.0):
    """Common-scalar formal continuation; not arbitrary mixed-sign rates."""
    print(f"\nReading 3  (N={N}, Q={Q}): common per-site profile scalar f")
    rho = uniform_rho(N).flatten(order="F")
    print(f"  {'per-site f':>15}  {'L.vec(1/N)':>11}  {'max non-ker Re':>14}   transient role")
    for f in [1.0, 0.5, 0.1, 0.0, -0.1, -0.5]:
        L = liouvillian(N, Q, f)
        res = float(np.linalg.norm(L @ rho))
        max_re, _ = gap_and_stability(L)
        if max_re < -1e-9:
            role = "DECAY   (transient complement contracts)"
        elif abs(max_re) <= 1e-9:
            role = "NEUTRAL (eternal oscillation, the mirror)"
        else:
            role = "SOURCE  (formal inverse-dephasing amplification)"
        print(f"  {f:15.2f}  {res:11.0e}  {max_re:+14.4f}   {role}")
    print("  => the sector state remains stationary; the transient linear stability changes.")
    print("     f<0 is non-CP inverse dephasing, not an energy-gain bath or Hopf certificate.")


def main():
    print("=" * 78)
    print("AT THE TARGET  (a point, and then how does it go on?)")
    print("=" * 78)
    part1_the_approach()
    part2_how_it_goes_on()


if __name__ == "__main__":
    main()
