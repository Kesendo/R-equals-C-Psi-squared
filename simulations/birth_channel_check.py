# -*- coding: utf-8 -*-
"""birth_channel_check.py

Is the birth channel (the depth-1 rung, rate 2*gamma) the slowest mortal mode of the WHOLE
Liouvillian, or only of the number-changing half of it?

reflections/THE_VIEW_ONTO_THE_MEMORY.md said the first, unconditionally, while its own crown
section a few paragraphs later says the longest memory is the EVEN occupation below a coupling
crossing. This builds the full 4^N Liouvillian (Heisenberg chain, uniform Z-dephasing) and reads
the smallest nonzero decay rate directly, which settles it: below the crossing the slowest mortal
mode is far under 2*gamma, and only above it does the slowest rate settle at exactly 2*gamma.

Run: python simulations/birth_channel_check.py   (a few seconds, N = 3 and 4)
"""
import numpy as np

I2 = np.eye(2)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]])
Z = np.array([[1, 0], [0, -1]], dtype=complex)


def op(n, site, m):
    out = np.array([[1]], dtype=complex)
    for s in range(n):
        out = np.kron(out, m if s == site else I2)
    return out


def liouvillian(n, J, gam):
    d = 2 ** n
    H = np.zeros((d, d), dtype=complex)
    for b in range(n - 1):
        for M in (X, Y, Z):
            H += (J / 4) * op(n, b, M) @ op(n, b + 1, M)
    L = -1j * (np.kron(H, np.eye(d)) - np.kron(np.eye(d), H.T))
    for s in range(n):
        Zs = op(n, s, Z)
        L += gam * (np.kron(Zs, Zs.conj()) - np.kron(np.eye(d), np.eye(d)))
    return L


gam = 0.05
print(f"gamma = {gam}, the birth channel (depth 1) sits at 2*gamma = {2*gam}")
for n in (3, 4):
    for J in (0.005, 0.02, 0.1, 0.5, 2.0):
        ev = np.linalg.eigvals(liouvillian(n, J, gam))
        re = -ev.real
        mortal = np.sort(re[re > 1e-9])
        slowest = mortal[0] if len(mortal) else float("nan")
        print(f"   N={n} J={J:5}: slowest MORTAL rate = {slowest:.6f}   "
              + ("BELOW 2*gamma: the birth channel is NOT the slowest here"
                 if slowest < 2 * gam - 1e-9 else "= 2*gamma: the birth channel IS the slowest"))
