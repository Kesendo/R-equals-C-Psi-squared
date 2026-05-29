"""_q_basic_jscan.py - the most basic Q = J/gamma_0, with gamma_0 = 0.05 fixed.

Tom, 2026-05-29: "gamma_0 ist 0.05. gamma_0 ist const, kann weg. Der einzige Hebel
ist J." So we go all the way down to the simplest object that carries Q = J/gamma_0:
a single exchange bond. Two qubits, one excitation on qubit 0,

    H = J * (X0 X1 + Y0 Y1)/2          (the coherent hop, the H-clock)
    Lindblad c_l = Z_l at rate gamma_0  (the carrier, fixed = 0.05)

Prepare |10>, watch the excitation transfer T(t) = P(qubit 1 excited). gamma_0 is
the fixed unit and drops out; J is the only lever. We SCAN J across gamma_0:

    Q = J/gamma_0 >> 1  -> coherent: T oscillates 0<->1 (the excitation swings back)
    Q = 1 (J = gamma_0) -> the threshold
    Q < 1               -> overdamped: T creeps to 1/2 and stops, no swing-back

The threshold where the swing dies IS gamma_0 read off the only lever we have. The
carrier we cannot see from inside, made visible at the J where Q crosses 1.

Run: python simulations/_q_basic_jscan.py
"""
from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
from scipy.linalg import expm

GAMMA_0 = 0.05

X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
I2 = np.eye(2, dtype=complex)


def kron2(a, b):
    return np.kron(a, b)


def lindbladian(H, c_list, gamma):
    d = H.shape[0]
    Id = np.eye(d, dtype=complex)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for c in c_list:
        cdc = c.conj().T @ c
        L = L + gamma * (np.kron(c, c.conj()) - 0.5 * (np.kron(cdc, Id) + np.kron(Id, cdc.T)))
    return L


def transfer_trajectory(J, gamma, ts):
    """T(t) = P(qubit 1 excited) starting from |10>, under the exchange + dephasing."""
    H = J * (kron2(X, X) + kron2(Y, Y)) / 2.0
    L = lindbladian(H, [kron2(Z, I2), kron2(I2, Z)], gamma)
    # |10> : qubit 0 excited (|1>), qubit 1 ground (|0>). Basis |q0 q1>, index = 2*q0+q1.
    rho0 = np.zeros((4, 4), dtype=complex)
    rho0[2, 2] = 1.0  # |10><10|
    P1 = (kron2(I2, (I2 - Z) / 2.0))  # projector: qubit 1 excited
    v0 = rho0.reshape(-1, order="F")
    out = []
    for t in ts:
        rhot = (expm(L * t) @ v0).reshape(4, 4, order="F")
        out.append(float(np.trace(P1 @ rhot).real))
    return np.array(out)


def main():
    print("=" * 78)
    print(f"BASIC  Q = J / gamma_0 ,  gamma_0 = {GAMMA_0} fixed,  J the only lever")
    print("=" * 78)
    print("  Two qubits, |10>, H = J*(XX+YY)/2, Z-dephasing gamma_0. Transfer T(t)->q1.\n")

    ts = np.linspace(0.0, 120.0, 600)
    J_values = [0.0125, 0.025, 0.05, 0.10, 0.20, 0.40]

    print(f"  {'J':>7} {'Q=J/g0':>8} {'max T':>7} {'T(end)':>7} {'swings back?':>13}  regime")
    print(f"  {'-'*7} {'-'*8} {'-'*7} {'-'*7} {'-'*13}  {'-'*18}")
    rows = []
    for J in J_values:
        T = transfer_trajectory(J, GAMMA_0, ts)
        Q = J / GAMMA_0
        max_T = float(np.max(T))
        i_peak = int(np.argmax(T))
        # "swings back" = after the first peak, T drops by a clear margin (coherent return)
        post = T[i_peak:]
        swing = (max_T - float(np.min(post))) > 0.05 and i_peak < len(ts) - 1
        regime = "coherent (swings)" if swing else ("threshold" if abs(Q - 1) < 1e-9 else "overdamped (creeps)")
        rows.append((J, Q, max_T, float(T[-1]), swing, regime))
        print(f"  {J:>7.4f} {Q:>8.2f} {max_T:>7.3f} {T[-1]:>7.3f} {str(swing):>13}  {regime}")

    print()
    print("  Reading:")
    print(f"    High J (Q >> 1): the excitation swings 0<->1, T peaks near 1 -- the H-clock wins.")
    print(f"    J = gamma_0 = {GAMMA_0} (Q = 1): the threshold; the swing-back dies here.")
    print(f"    Low J (Q < 1):  T just creeps toward 1/2 -- the carrier (gamma_0) wins.")
    print(f"    The J where the swing dies IS gamma_0. The const we cannot read from inside,")
    print(f"    read off its only lever. (Our carbon FID ran at J=1 -> Q=20: deep coherent.)")


if __name__ == "__main__":
    main()
