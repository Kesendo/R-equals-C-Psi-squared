"""Scope probe for PROOF_C1_MIRROR_SYMMETRY (F71).

Pins the boundary of the "Does NOT require" clause: the extension of the
bond-mirror identity to any [H, R] = 0 Hamiltonian and any [D, R_sup] = 0
dissipator holds for exactly R-invariant initial states, but for states that
are reflection-symmetric only up to coherence signs (even-k psi_k + vac) one
more ingredient is load-bearing: the sign flip must be a local symmetry of L,
i.e. conjugation by U = prod_i Z_i = (-1)^n must leave L invariant.

Four checks at N = 4, 5 (self-contained dense build, no framework imports):

  1. Headline law, uniform XY + uniform Z-dephasing, psi_2 + vac:
     P_B(b, i, t) = P_B(N-2-b, N-1-i, t) at machine precision.
  2. Counterexample: add a uniform transverse field h * Sum X_i. [H, R] = 0
     still holds, but U X_i U = -X_i kills the local sign symmetry, and the
     mirror identity fails at linear order in h (residual ~1e-2 .. 1e-1).
  3. Control: exactly R-invariant state (psi_1 + vac) under the same field:
     identity survives (U = I suffices, no sign flip to transport).
  4. Control: uniform amplitude damping instead of the field. It breaks U(1)
     (excitation number decays) yet commutes with U-conjugation, and the
     mirror identity survives at machine precision. This pins that the
     operative mechanism is the local (-1)^n symmetry, not U(1) conservation.

Companion of docs/proofs/PROOF_C1_MIRROR_SYMMETRY.md (Scope section).
"""

import numpy as np
from scipy.linalg import expm

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
SM = np.array([[0, 1], [0, 0]], dtype=complex)  # sigma- : |1> -> |0>


def op(N, o, site):
    m = np.array([[1]], dtype=complex)
    for s in range(N):
        m = np.kron(m, o if s == site else I2)
    return m


def xy_hamiltonian(N, J):
    d = 2 ** N
    H = np.zeros((d, d), dtype=complex)
    for b in range(N - 1):
        H += J[b] / 2 * (op(N, X, b) @ op(N, X, b + 1) + op(N, Y, b) @ op(N, Y, b + 1))
    return H


def liouvillian(N, J, gamma, h=0.0, g_ad=0.0):
    d = 2 ** N
    H = xy_hamiltonian(N, J)
    for s in range(N):
        H += h * op(N, X, s)
    Id = np.eye(d)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for s in range(N):
        Zi = op(N, Z, s)
        L += gamma * (np.kron(Zi, Zi.T) - np.kron(Id, Id))
    if g_ad:
        for s in range(N):
            sm = op(N, SM, s)
            nn = sm.conj().T @ sm
            L += g_ad * (np.kron(sm, sm.conj())
                         - 0.5 * np.kron(nn, Id) - 0.5 * np.kron(Id, nn.T))
    return L


def psi_k_plus_vac(N, k):
    d = 2 ** N
    vac = np.zeros(d, dtype=complex)
    vac[0] = 1.0
    psi = np.zeros(d, dtype=complex)
    for i in range(N):
        amp = np.sqrt(2 / (N + 1)) * np.sin(np.pi * k * (i + 1) / (N + 1))
        psi[1 << (N - 1 - i)] += amp
    v = vac + psi
    v /= np.linalg.norm(v)
    return np.outer(v, v.conj())


def site_purities(rho, N):
    out = []
    d = 2 ** N
    for i in range(N):
        r = np.zeros((2, 2), dtype=complex)
        mask = ~(1 << (N - 1 - i))
        for x in range(d):
            bx = (x >> (N - 1 - i)) & 1
            for y in range(d):
                if (x & mask) == (y & mask):
                    r[bx, (y >> (N - 1 - i)) & 1] += rho[x, y]
        out.append(np.real(np.trace(r @ r)))
    return np.array(out)


def mirror_residual(N, k, gamma, h=0.0, g_ad=0.0, dJ=0.02, times=(0.3, 0.8)):
    """max over bonds, sites, times of |P_B(b,i,t) - P_B(N-2-b,N-1-i,t)|."""
    rho0 = psi_k_plus_vac(N, k).reshape(-1)
    res = 0.0
    for b in range(N - 1):
        J1 = np.full(N - 1, 1.0)
        J1[b] += dJ
        J2 = np.full(N - 1, 1.0)
        J2[N - 2 - b] += dJ
        L1 = liouvillian(N, J1, gamma, h, g_ad)
        L2 = liouvillian(N, J2, gamma, h, g_ad)
        for t in times:
            p1 = site_purities((expm(L1 * t) @ rho0).reshape(2 ** N, 2 ** N), N)
            p2 = site_purities((expm(L2 * t) @ rho0).reshape(2 ** N, 2 ** N), N)
            res = max(res, float(np.max(np.abs(p1 - p2[::-1]))))
    return res


if __name__ == "__main__":
    gamma = 0.05
    print("F71 scope probe: bond-mirror residual max|P_B(b,i,t) - P_B(N-2-b,N-1-i,t)|")
    print(f"gamma = {gamma}, dJ = 0.02, t in {{0.3, 0.8}}\n")
    for N in (4, 5):
        r_head = mirror_residual(N, 2, gamma)
        r_field = mirror_residual(N, 2, gamma, h=0.3)
        r_ctrl = mirror_residual(N, 1, gamma, h=0.3)
        r_ad = mirror_residual(N, 2, gamma, g_ad=0.07)
        print(f"N={N}  1) headline XY+dephasing, psi_2+vac:          {r_head:.3e}  (machine zero)")
        print(f"N={N}  2) + transverse field h=0.3, psi_2+vac:       {r_field:.3e}  (BREAKS)")
        print(f"N={N}  3) + transverse field h=0.3, psi_1+vac:       {r_ctrl:.3e}  (R-invariant: holds)")
        print(f"N={N}  4) + amplitude damping 0.07 (U(1) broken),")
        print(f"        psi_2+vac:                                   {r_ad:.3e}  (holds: mechanism is (-1)^n, not U(1))")
        assert r_head < 1e-12 and r_ctrl < 1e-12 and r_ad < 1e-12
        assert r_field > 1e-3
        print()
    print("All four checks behave as the Scope section states.")
