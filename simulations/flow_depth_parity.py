#!/usr/bin/env python3
"""The single-excitation flow as a walk across the drain-depth axis.

THE_VIEW_ONTO_THE_MEMORY sorts the Liouvillian's modes by drain depth = popcount(i XOR j).
THE_FLOW_BETWEEN_TWO_SINGULARITIES runs a single excitation into the 1/N target. This probe
ties them together and corrects an earlier conflation.

Three things, all bit-exact / state-based (no reliance on degeneracy-mixed eigenvectors):

1. drain depth IS the light: for any basis dyad |i><j|, popcount(i XOR j) == n_XY (the number
   of sites carrying X/Y rather than {I,Z}). The two axes are one.
2. the flow lives purely on EVEN depth: rho(t) of a definite-number state has odd-parity weight
   identically 0; its light <n_XY> is a transient tide (0 -> peak -> 0), and the steady state is
   100% depth-0 (the memory). The depth-1, rate-2gamma VACUUM coherence (the Liouvillian's global
   slowest non-kernel mode) is NOT used by the flow (overlap ~ 1e-16).
3. the flow's own slowest used rate is on the even ladder and is Q-DEPENDENT (even for uniform),
   distinct from the unused 2gamma vacuum floor.

Run: python simulations/flow_depth_parity.py
"""
import itertools
import sys

import numpy as np

sys.path.insert(0, "simulations")
from light_content import H_xy_unit, op_at, Z as Zm  # noqa: E402
from bound_pair_light import light_parity              # per-site, total, even-w, odd-w  # noqa: E402

PAULI_IXYZ = "IXYZ"


def depth_equals_light_identity(N=4):
    """popcount(i XOR j) == n_XY of the dyad |i><j|, bit-exact over all basis pairs."""
    worst = 0
    for i in range(2 ** N):
        for j in range(2 ** N):
            depth = bin(i ^ j).count("1")
            # n_XY of |i><j|: a site differs (X/Y) iff its bit differs between i and j
            nxy = sum(1 for b in range(N) if ((i >> b) & 1) != ((j >> b) & 1))
            worst = max(worst, abs(depth - nxy))
    return worst


def build_L(N, Q):
    d = 2 ** N
    Id = np.eye(d)
    H1 = H_xy_unit(N)
    L = -1j * Q * (np.kron(Id, H1) - np.kron(H1.T, Id))
    for l in range(N):
        Zl = op_at(N, l, Zm)
        L += np.kron(Zl, Zl) - np.kron(Id, Id)            # gamma = 1, uniform
    return L


def flow_tide_and_parity(N, Q):
    """Propagate the single-excitation flow; return (times, <n_XY>(t), odd-weight(t))."""
    d = 2 ** N
    L = build_L(N, Q)
    w, V = np.linalg.eig(L)
    Vinv = np.linalg.inv(V)
    idx = 1 << (N - 1)
    rho0 = np.zeros((d, d), complex)
    rho0[idx, idx] = 1.0
    c = Vinv @ rho0.flatten(order="F")
    ts = [0.0, 0.1, 0.3, 0.6, 1.0, 1.5, 2.5, 4.0]
    out = []
    for t in ts:
        vt = V @ (np.exp(w * t) * c)
        _, nxy, _, wo = light_parity(vt, N)
        out.append((t, nxy, wo))
    # overlap of the flow with the global slowest non-kernel mode (the vacuum coherence)
    nonk = np.abs(w) > 1e-7
    gk = int(np.argmax(np.where(nonk, w.real, -np.inf)))
    return out, -float(w[gk].real), abs(c[gk])


def flow_slowest_used_rate(N, Q):
    """Decay rate of the per-site occupation deviation from 1/N (the flow's real slowest mode)."""
    d = 2 ** N
    L = build_L(N, Q)
    w, V = np.linalg.eig(L)
    Vinv = np.linalg.inv(V)
    idx = 1 << (N - 1)
    rho0 = np.zeros((d, d), complex)
    rho0[idx, idx] = 1.0
    c = Vinv @ rho0.flatten(order="F")
    nops = [(np.eye(d) - op_at(N, k, Zm)) / 2 for k in range(N)]
    ts = np.linspace(3, 8, 11)
    dev = []
    for t in ts:
        rho = (V @ (np.exp(w * t) * c)).reshape(d, d, order="F")
        occ = np.array([np.trace(nops[k] @ rho).real for k in range(N)])
        dev.append(np.max(np.abs(occ - 1.0 / N)))
    return -np.polyfit(ts, np.log(np.array(dev)), 1)[0]


def main():
    print("1. drain depth IS the light:  max |popcount(i XOR j) - n_XY| over all dyads")
    for N in (2, 3, 4):
        print(f"     N={N}:  {depth_equals_light_identity(N)}   (0 = bit-exact identity)")

    N, Q = 5, 1.5
    print(f"\n2. the flow on the even ladder  (N={N}, Q={Q}, gamma=1 uniform):")
    tide, vac_rate, vac_overlap = flow_tide_and_parity(N, Q)
    print(f"     {'t':>6}  {'<n_XY> (tide)':>13}  {'odd-parity weight':>17}")
    for t, nxy, wo in tide:
        print(f"     {t:>6.2f}  {nxy:>13.4f}  {wo:>17.4f}")
    print(f"     -> light is a transient tide (0 -> peak -> 0); odd weight identically 0.")
    print(f"     global slowest non-kernel mode (the vacuum coherence): rate={vac_rate:.4f}, "
          f"n_XY=1, flow overlap={vac_overlap:.1e}  (UNUSED by the flow)")

    print(f"\n3. the flow's slowest USED rate is on the even ladder and Q-dependent:")
    for Q in (1.5, 1000.0):
        print(f"     N={N}, Q={Q:>7.1f}:  rate = {flow_slowest_used_rate(N, Q):.4f}"
              f"   (vs unused vacuum floor 2.0000)")


if __name__ == "__main__":
    main()
