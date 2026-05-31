#!/usr/bin/env python3
"""The new object (the two-magnon bound pair), seen in light.

The single-excitation slowest mode carries exactly 1 light quantum (n_XY=1). What does the
bound pair carry? We take the most-bound two-excitation eigenstate, form rho = |bound><bound|,
and read its light content <n_XY> (the absorption theorem's 'light mass' = exposure to the
carrier = its mortality), the per-site light, and the n_XY-parity split, across binding Delta.

Does binding change the light mass (the lifetime), and is the new object in a definite parity
sector? Compared against the free (Delta=0) two-excitation state.
"""
import itertools
import sys

import numpy as np

sys.path.insert(0, "simulations")
from light_content import PAULI                       # noqa: E402
from two_magnon_bound_state import bound_pair_vec     # noqa: E402


def light_parity(v, N):
    """Per-site <X/Y at k>, total <n_XY>, and the weight in even vs odd n_XY, of the Liouville
    vector v = vec_F(rho)."""
    d = 2 ** N
    M = v.reshape(d, d, order="F")
    norm2 = 0.0
    xy_site = np.zeros(N)
    nxy_tot = 0.0
    w_even = w_odd = 0.0
    for combo in itertools.product("IXYZ", repeat=N):
        P = np.array([[1]], complex)
        for s in combo:
            P = np.kron(P, PAULI[s])
        c = np.sum(P.conj() * M) / d
        ww = abs(c) ** 2
        norm2 += ww
        nxy = sum(1 for s in combo if s in "XY")
        nxy_tot += ww * nxy
        if nxy % 2 == 0:
            w_even += ww
        else:
            w_odd += ww
        for kk, s in enumerate(combo):
            if s in "XY":
                xy_site[kk] += ww
    return xy_site / norm2, nxy_tot / norm2, w_even / norm2, w_odd / norm2


def main():
    N = 6
    print(f"N={N}. The bound pair in light: <n_XY> (light mass), per-site light, parity, vs Delta\n")
    print(f"  {'Delta':>6}  {'<sep>note':>10}  {'<n_XY>':>7}  {'even w':>7}  {'odd w':>7}  per-site X/Y")
    for delta in [0.0, 0.5, 1.0, 2.0, 4.0, 8.0]:
        v = bound_pair_vec(N, delta)
        rho = np.outer(v, v.conj())
        xy, nxy, we, wo = light_parity(rho.flatten(order="F"), N)
        tag = "free" if delta == 0 else "bound"
        print(f"  {delta:>6.2f}  {tag:>10}  {nxy:>7.4f}  {we:>7.4f}  {wo:>7.4f}  {np.round(xy,3).tolist()}")
    print("\n  reading: <n_XY> is the bound pair's light mass (its exposure to the carrier, ~its")
    print("  decay scale). Watch whether binding (Delta up) changes the light mass, and which")
    print("  n_XY parity sector the new object lives in.")


if __name__ == "__main__":
    main()
