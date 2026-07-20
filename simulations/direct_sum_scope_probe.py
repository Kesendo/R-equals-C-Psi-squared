"""Scope probe for DIRECT_SUM_DECOMPOSITION and PROOF_PARITY_SELECTION_RULE.

Pins the selective breaking of the two structures the direct-sum proof rests
on, against the two perturbations its scope section names:

  structure P = n_XY parity grading: off-parity-block Frobenius norm of L
                in the Hilbert-Schmidt Pauli basis (0 iff L = L_even + L_odd).
  structure M = mirror palindrome:   ||Pi L Pi^-1 + L + 2*Sigma_gamma*I||_F.

  perturbation            parity grading P     palindrome M
  ---------------------------------------------------------
  none (Heis + Z-deph)    exactly 0            machine zero
  + amplitude damping     exactly 0 (HOLDS)    ~1.9 (BREAKS)
  + transverse field      ~5.9 (BREAKS)        machine zero (HOLDS)

The cross shows the two breaks-for items are NOT interchangeable: a
transverse field mixes the sectors because a Hamiltonian term enters the
commutator one-sided, while amplitude damping acts bilinearly (A rho Adag
puts an odd-parity operator on both sides, net parity change even, and
Adag*A = (I-Z)/2 is parity-even), so ANY jump operator whose Pauli
components share one n_XY parity preserves the grading. What T1 breaks is
the mirror (fixed arrow of time; closed forms in F82/F84), not the wall.

Heisenberg chain + uniform Z-dephasing, N = 3, dense superoperator build,
self-contained (no framework imports).
"""

import numpy as np
from itertools import product

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
SM = np.array([[0, 1], [0, 0]], dtype=complex)  # sigma- : |1> -> |0>
PAULI = {'I': I2, 'X': X, 'Y': Y, 'Z': Z}


def op(N, o, site):
    m = np.array([[1]], dtype=complex)
    for s in range(N):
        m = np.kron(m, o if s == site else I2)
    return m


def make_L(N, J=1.0, gamma=0.05, h=0.0, g_ad=0.0):
    d = 2 ** N
    H = np.zeros((d, d), dtype=complex)
    for b in range(N - 1):
        for P in (X, Y, Z):
            H += J / 4 * op(N, P, b) @ op(N, P, b + 1)
    for s in range(N):
        H += h * op(N, X, s)

    def L_of(rho):
        out = -1j * (H @ rho - rho @ H)
        for s in range(N):
            Zi = op(N, Z, s)
            out += gamma * (Zi @ rho @ Zi - rho)
        if g_ad:
            for s in range(N):
                sm = op(N, SM, s)
                nn = sm.conj().T @ sm
                out += g_ad * (sm @ rho @ sm.conj().T - 0.5 * (nn @ rho + rho @ nn))
        return out
    return L_of


def pauli_basis(N):
    labels = [''.join(p) for p in product('IXYZ', repeat=N)]
    mats = []
    for lab in labels:
        m = np.array([[1]], dtype=complex)
        for ch in lab:
            m = np.kron(m, PAULI[ch])
        mats.append(m / np.sqrt(2 ** N))
    return labels, mats


def superop_matrix(N, **kw):
    L_of = make_L(N, **kw)
    labels, mats = pauli_basis(N)
    dim = len(mats)
    M = np.zeros((dim, dim), dtype=complex)
    for j, mj in enumerate(mats):
        Lmj = L_of(mj)
        for i, mi in enumerate(mats):
            M[i, j] = np.trace(mi.conj().T @ Lmj)
    return labels, M


def pi_matrix(labels):
    # per-site Pi of MIRROR_SYMMETRY_PROOF: I->X, X->I, Y->iZ, Z->iY
    step = {'I': ('X', 1.0), 'X': ('I', 1.0), 'Y': ('Z', 1j), 'Z': ('Y', 1j)}
    dim = len(labels)
    index = {lab: i for i, lab in enumerate(labels)}
    P = np.zeros((dim, dim), dtype=complex)
    for j, lab in enumerate(labels):
        out, coef = '', 1.0 + 0j
        for ch in lab:
            nl, c = step[ch]
            out += nl
            coef *= c
        P[index[out], j] = coef
    return P


if __name__ == "__main__":
    N, gamma = 3, 0.05
    sigma = N * gamma
    print("Direct-sum scope probe (Heisenberg chain + Z-dephasing, N=3, gamma=0.05)")
    print("P = off-parity-block Frobenius of L; M = ||Pi L Pi^-1 + L + 2*sigma*I||_F\n")
    results = {}
    for name, kw in [("baseline", dict()),
                     ("+ amplitude damping 0.07", dict(g_ad=0.07)),
                     ("+ transverse field 0.3*Sum X_i", dict(h=0.3))]:
        labels, M = superop_matrix(N, **kw)
        nxy = [sum(c in 'XY' for c in lab) for lab in labels]
        off = np.sqrt(sum(abs(M[i, j]) ** 2
                          for i in range(len(labels)) for j in range(len(labels))
                          if (nxy[i] % 2) != (nxy[j] % 2)))
        P = pi_matrix(labels)
        pal = np.linalg.norm(P @ M @ np.linalg.inv(P) + M + 2 * sigma * np.eye(len(labels)))
        results[name] = (off, pal)
        print(f"{name:34s}: parity P = {off:.3e}   palindrome M = {pal:.3e}")

    p0, m0 = results["baseline"]
    p_ad, m_ad = results["+ amplitude damping 0.07"]
    p_f, m_f = results["+ transverse field 0.3*Sum X_i"]
    assert p0 < 1e-12 and m0 < 1e-12
    assert p_ad < 1e-12 and m_ad > 1e-2          # AD: grading holds, mirror breaks
    assert p_f > 1e-2 and m_f < 1e-12            # field: grading breaks, mirror holds
    print("\nCross confirmed: T1 breaks the mirror, not the wall; the field the wall, not the mirror.")
