"""Finite palindrome reading for XX+XY under two local channel terms.

``n_bar`` is an external jump-rate parameter, not heat created by the model.
The unique candidate centre is read from the trace; no best-centre search or
threshold-defined "orphan" population is used.
"""
import sys

import numpy as np
from framework import f1_distance_in_eps

sys.stdout.reconfigure(encoding="utf-8")

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
SM = np.array([[0, 1], [0, 0]], dtype=complex)  # |0><1|
SP = np.array([[0, 0], [1, 0]], dtype=complex)  # |1><0|
PM = {"X": X, "Y": Y, "Z": Z}


def site_op(op, k, N):
    m = np.array([[1.0 + 0j]])
    for i in range(N):
        m = np.kron(m, op if i == k else I2)
    return m


def build_H_xxxy(N, J=1.0):
    """XX + XY on every bond: a palindrome-breaking Hamiltonian (the Brecher)."""
    d = 2 ** N
    H = np.zeros((d, d), dtype=complex)
    for b in range(N - 1):
        for t in ("XX", "XY"):
            H += J * site_op(PM[t[0]], b, N) @ site_op(PM[t[1]], b + 1, N)
    return H


def build_L_thermal(H, gamma_z, gamma_amp, n_bar):
    """Thermal Lindbladian, reused verbatim from self_heating_fixpoint.py:
    -i[H,.] + Z-dephasing + amplitude/thermal channel (sigma-, sigma+ at n_bar)."""
    d = H.shape[0]
    Id = np.eye(d, dtype=complex)
    n = int(np.log2(d))
    L = -1j * (np.kron(Id, H) - np.kron(H.T, Id))
    for k in range(n):
        ops = []
        if gamma_z > 0:
            ops.append(np.sqrt(gamma_z) * site_op(Z, k, n))
        if gamma_amp > 0:
            ops.append(np.sqrt(gamma_amp * (1 + n_bar)) * site_op(SM, k, n))
        if gamma_amp > 0 and n_bar > 0:
            ops.append(np.sqrt(gamma_amp * n_bar) * site_op(SP, k, n))
        for Lop in ops:
            LdL = Lop.conj().T @ Lop
            L += np.kron(Lop.conj(), Lop)
            L -= 0.5 * np.kron(Id, LdL)
            L -= 0.5 * np.kron(LdL.T, Id)
    return L


def osc_frequencies(ev, thresh=1e-6):
    return sorted(set(round(abs(e.imag), 6) for e in ev if abs(e.imag) > thresh))


def analyse(N, n_bar, gz=0.05, ga=0.05):
    H = build_H_xxxy(N)
    L = build_L_thermal(H, gz, ga, n_bar)
    ev = np.linalg.eigvals(L)
    n_freq = len(osc_frequencies(ev))
    centre = float(np.mean(ev).real)
    distance, _ = f1_distance_in_eps(ev, -centre)
    return n_freq, centre, distance


if __name__ == "__main__":
    print("XX+XY with Z dephasing and a finite-occupation amplitude channel")
    print("n_bar is supplied externally; it is not a feedback or heat variable.\n")
    for N in (3, 4, 5):
        print(f"--- N = {N}   (L is {4 ** N}x{4 ** N}, gz = ga = 0.05) ---")
        print(f"{'n_bar':>7} {'#freq':>6} {'centre':>9} {'F1 eps*rho':>12}")
        for n_bar in (0.0, 0.5, 1.0, 2.0, 5.0):
            nf, centre, distance = analyse(N, n_bar)
            print(f"{n_bar:>7.1f} {nf:>6} {centre:>9.4f} {distance:>12.2e}")
        print()
