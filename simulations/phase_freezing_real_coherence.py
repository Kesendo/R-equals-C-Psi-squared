#!/usr/bin/env python3
"""ph03 freezing (COCKPIT_SCALING OQ-096): is the central Bell pair's coherence phase EXACTLY
zero (a symmetry) or just small and shrinking with N (locality)?

The doc says: for the center Bell pair on a Heisenberg chain under uniform Z-dephasing, the phase
of the off-diagonal element |00><11| stays at "exact zero" for N>=7, but ph03 VARIES at N=5. An
exact symmetry would freeze it at every N, not only N>=7. So compute arg(rho_pair[00,11])(t)
directly for N=5,7,9 and see: machine-zero everywhere (symmetry), or nonzero at N=5 shrinking with
N (locality, the asymmetric boundary's finite-speed influence on the central pair).

Setup: |Phi+> on the central pair (c, c+1), c=(N-1)//2, tensor |+> on the rest; H = sum_bonds
(XX+YY+ZZ), uniform Z-dephasing gamma; sparse Liouvillian, expm_multiply on a t-grid.

RESULT (2026-06-04, N=5,7,9): max|Im(coherence)| ~ 1e-16 at EVERY N -> the coherence
rho_pair[00,11] is EXACTLY REAL, so ph03 = arg is in {0, pi} exactly. max|ph03| ~ 1e-14 (the
central pair stays real AND positive at all N here). PROOF: the global spin-flip P = X^(tensor N)
= Pi^2 commutes with H (Heisenberg is spin-flip symmetric) and with the Z-dephasing dissipator
(P anticommutes with each Z_l, so the dissipator is P-invariant); the initial state is P-invariant;
hence rho(t) = P rho(t) P, and P (flipping all spins) maps <00,e|rho|11,e> to <11,~e|rho|00,~e>,
giving (summed over the environment) rho_pair[00,11] = conj(rho_pair[00,11]) = real. So ph03 is in
{0, pi} for all N, exactly. The doc's "freezes only at N>=7" is a PCA-variance-gate artifact: arg
of the real coherence is a clean 0 while it stays positive, and only reads as variance when it
dips toward zero near an ESD (arg of ~0 is numerically noisy). The symmetry is universal; there is
no N-threshold. (Pi^2 = X^N and [Pi^2, L] = 0 is the framework's known Z2 parity symmetry; see
THE_OTHER_SIDE / PI_AS_TIME_REVERSAL.)
"""
from __future__ import annotations

import sys
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import expm_multiply

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

I2 = sp.identity(2, format="csr", dtype=complex)
X = sp.csr_matrix([[0, 1], [1, 0]], dtype=complex)
Y = sp.csr_matrix([[0, -1j], [1j, 0]], dtype=complex)
Z = sp.csr_matrix([[1, 0], [0, -1]], dtype=complex)


def op_at(N, site, P):
    o = None
    for i in range(N):
        f = P if i == site else I2
        o = f if o is None else sp.kron(o, f, format="csr")
    return o


def setbit(idx, qubit, val, N):
    pos = N - 1 - qubit            # qubit 0 = most significant bit
    if val:
        return idx | (1 << pos)
    return idx & ~(1 << pos)


def coherence_phase_trajectory(N, J, gamma, a, b, taus):
    d = 1 << N
    H = sp.csr_matrix((d, d), dtype=complex)
    for i in range(N - 1):
        for P in (X, Y, Z):
            H = H + J * (op_at(N, i, P) @ op_at(N, i + 1, P))
    Id = sp.identity(d, format="csr", dtype=complex)
    Idd = sp.identity(d * d, format="csr", dtype=complex)
    L = -1j * (sp.kron(Id, H, format="csr") - sp.kron(H.T, Id, format="csr"))
    for l in range(N):
        Zl = op_at(N, l, Z)
        L = L + gamma * (sp.kron(Zl, Zl, format="csr") - Idd)

    # |psi0> = (1/sqrt2)^(N-1) * [bit_a == bit_b]  (Bell+ on a,b; |+> on the rest)
    psi = np.zeros(d, dtype=complex)
    amp = (1.0 / np.sqrt(2.0)) ** (N - 1)
    for s in range(d):
        ba = (s >> (N - 1 - a)) & 1
        bb = (s >> (N - 1 - b)) & 1
        if ba == bb:
            psi[s] = amp
    rho0 = np.outer(psi, psi.conj())
    vec0 = rho0.reshape(-1, order="F")

    # env index list: ket has (a=0,b=0), bra has (a=1,b=1), env bits matched
    env_qubits = [q for q in range(N) if q not in (a, b)]
    ket_list, bra_list = [], []
    for e in range(1 << (N - 2)):
        ki = setbit(setbit(0, a, 0, N), b, 0, N)
        bi = setbit(setbit(0, a, 1, N), b, 1, N)
        for t, q in enumerate(env_qubits):
            bit = (e >> t) & 1
            ki = setbit(ki, q, bit, N)
            bi = setbit(bi, q, bit, N)
        ket_list.append(ki); bra_list.append(bi)
    ket_arr = np.array(ket_list); bra_arr = np.array(bra_list)

    evol = expm_multiply(L, vec0, start=taus[0], stop=taus[-1], num=len(taus))
    phases, ims = [], []
    for k in range(len(taus)):
        rho = evol[k].reshape(d, d, order="F")
        c = complex(np.sum(rho[ket_arr, bra_arr]))     # rho_pair[00,11], summed over env
        phases.append(np.angle(c))
        ims.append(c.imag)
    return np.array(phases), np.array(ims)


def main():
    J, gamma = 1.0, 0.05
    taus = np.linspace(0.0, 3.0, 25)
    print(f"  ph03 freezing: J={J}, gamma={gamma}, central pair (c,c+1), c=(N-1)//2\n")
    print(f"  {'N':>3} {'pair':>8} {'max|ph03|':>14} {'max|Im(coh)|':>15}   reading")
    for N in (5, 7, 9):
        c = (N - 1) // 2
        a, b = c, c + 1
        ph, im = coherence_phase_trajectory(N, J, gamma, a, b, taus)
        mph, mim = float(np.max(np.abs(ph))), float(np.max(np.abs(im)))
        reading = "EXACT zero (symmetry)" if mph < 1e-12 else ("tiny (locality?)" if mph < 1e-3 else "varies")
        print(f"  {N:>3} {f'({a},{b})':>8} {mph:>14.3e} {mim:>15.3e}   {reading}")


if __name__ == "__main__":
    main()
