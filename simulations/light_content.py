#!/usr/bin/env python3
"""Is n_XY (the light content) the operator behind the sterile zone vs the birth canal?

Absorption Theorem: Re(lambda) = -2 * sum_k gamma_k * <X/Y at site k>, exactly. The rate is
ENTIRELY the dissipator's light content; the Hamiltonian only mixes the Pauli basis. So the rate
is Q-dependent only if the slow mode's light content <n_XY> changes with Q.

Test: for the slowest mode, decompose its eigenvector in the Pauli basis, read the per-site light
<X/Y at k> and the total <n_XY>, at Q=1.5 and Q=1000, for a sterile and a birth-canal profile.
  - Cross-check (validates everything): predicted rate 2*sum gamma_k <X/Y_k> == actual -Re(lambda).
  - Sterile: the per-site light DISTRIBUTION <X/Y at k> of the slow subspace is Q-independent
    (frozen). Birth canal: the distribution drifts with Q (the rate moves only through the
    gamma-weighted share 2*sum_l gamma_l*<X/Y_l>). NOTE (corrected 2026-06-10): the TOTAL
    <n_XY> is 1.00000 frozen in BOTH zones (this script's own output); the zone criterion is
    the distribution, not the total. Banked: birth_canal_boundary_pathdependence.py; the
    per-mode identity Re(lambda) = -2*gamma*light(v) is exact and test-gated in
    F8PartnerLightComplementarityTests (2026-06-10).
"""
import itertools

import numpy as np

I2 = np.eye(2)
X = np.array([[0, 1], [1, 0]], complex)
Y = np.array([[0, -1j], [1j, 0]], complex)
Z = np.array([[1, 0], [0, -1]], complex)
PAULI = {"I": I2, "X": X, "Y": Y, "Z": Z}


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


def slow_mode(N, Q, profile, H1):
    d = 2 ** N
    Id = np.eye(d)
    L = -1j * Q * (np.kron(Id, H1) - np.kron(H1.T, Id))
    for l in range(N):
        Zl = op_at(N, l, Z)
        L += profile[l] * (np.kron(Zl, Zl) - np.kron(Id, Id))
    w, V = np.linalg.eig(L)
    nonkernel = np.abs(w) > 1e-7
    k = int(np.argmax(np.where(nonkernel, w.real, -np.inf)))
    return -float(w[k].real), V[:, k]


def light(v, N):
    """Per-site <X/Y at k> and total <n_XY> of the Liouville-space eigenvector v (F-order vec)."""
    d = 2 ** N
    M = v.reshape(d, d, order="F")
    norm2 = 0.0
    xy_site = np.zeros(N)
    nxy_tot = 0.0
    for combo in itertools.product("IXYZ", repeat=N):
        P = np.array([[1]], complex)
        for s in combo:
            P = np.kron(P, PAULI[s])
        c = np.sum(P.conj() * M) / d          # Tr(P^dagger M)/d, Pauli Hermitian
        ww = abs(c) ** 2
        norm2 += ww
        nxy = sum(1 for s in combo if s in "XY")
        nxy_tot += ww * nxy
        for kk, s in enumerate(combo):
            if s in "XY":
                xy_site[kk] += ww
    return xy_site / norm2, nxy_tot / norm2


def main():
    N = 5
    H1 = H_xy_unit(N)
    profiles = {
        "peaked-V  (sterile)":   [0.25, 0.75, 3.0, 0.75, 0.25],
        "flat-bulk (birth canal)": [0.25, 1.5, 1.5, 1.5, 0.25],
    }
    print(f"N={N}. Absorption theorem: rate = 2*sum gamma_k <X/Y at k>. Is <n_XY> Q-frozen (sterile)?\n")
    for name, p in profiles.items():
        p = list(np.array(p, float) * N / np.sum(p))
        print(f"  {name}   profile {np.round(p,3).tolist()}")
        nxy_lo = nxy_hi = None
        for Q in [1.5, 1000.0]:
            rate, v = slow_mode(N, Q, p, H1)
            xy_site, nxy = light(v, N)
            pred = 2.0 * sum(p[k] * xy_site[k] for k in range(N))
            print(f"    Q={Q:>7.1f}  rate={rate:.5f}  2*Sum(g*XY)={pred:.5f}  (theorem err {abs(rate-pred):.1e})"
                  f"   <n_XY>={nxy:.5f}   per-site XY={np.round(xy_site,4).tolist()}")
            if Q == 1.5:
                nxy_lo = nxy
            else:
                nxy_hi = nxy
        drift = nxy_hi - nxy_lo
        print(f"    <n_XY> drift (Q: 1.5 -> 1000) = {drift:+.5f}   "
              f"{'FROZEN (sterile)' if abs(drift) < 1e-4 else 'DRIFTS (Hamiltonian mixes light in -> birth canal)'}\n")
    print("  reading: where <n_XY> is Q-frozen, the rate is the closed form (sterile); where the")
    print("  Hamiltonian pulls light into the slow mode (<n_XY> drifts), the rate follows -> birth canal.")


if __name__ == "__main__":
    main()
