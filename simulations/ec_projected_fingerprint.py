"""The projected-component fingerprint of ERROR_CORRECTION_PALINDROME Section 3.

N = 3 Heisenberg chain (J = 1, sum over XX + YY + ZZ), Z-dephasing gamma = 0.05 on
every site, L in the Pauli basis with coordinates v_s = Tr(P_s rho). P_osc is the
spectral projector onto every eigenvalue with Im(lambda) != 0 (40 of 64; the
nearest non-oscillating |Im| is at roundoff, the smallest oscillating one near 2).
A spectral projector does not depend on how the eigenvectors are normalized.

For Bell(0,1) = (|000> + |110>)/sqrt(2) and each single-qubit error E on it, the
fingerprint is the vector |(P_osc v)_s| / 2^N, one entry per Pauli string, and the
reported number is the largest entry-wise change against the error-free baseline.
On site 2 the error acts on the state as X^N does (X_2 and Y_2 up to a phase; Z_2
as the identity), and X^N x X^N commutes with L and only flips Pauli signs, so
those three changes are exactly 0; the script checks the identity on the state
exactly and prints the float reading beside it.

Output: simulations/results/ec_projected_fingerprint.txt
"""

from itertools import product as iproduct
from pathlib import Path

import numpy as np

OUT = Path(__file__).resolve().parent / "results" / "ec_projected_fingerprint.txt"

I2 = np.eye(2, dtype=complex)
SX = np.array([[0, 1], [1, 0]], dtype=complex)
SY = np.array([[0, -1j], [1j, 0]], dtype=complex)
SZ = np.array([[1, 0], [0, -1]], dtype=complex)
PAULIS = [I2, SX, SY, SZ]
NAMES = "IXYZ"


def kron_all(ops):
    m = ops[0]
    for o in ops[1:]:
        m = np.kron(m, o)
    return m


def build(n, gamma):
    d = 2 ** n
    idx = list(iproduct(range(4), repeat=n))
    pm = [kron_all([PAULIS[i] for i in s]) for s in idx]
    h = np.zeros((d, d), dtype=complex)
    for b in range(n - 1):
        for p in (SX, SY, SZ):
            ops = [I2] * n
            ops[b] = p
            ops[b + 1] = p
            h += kron_all(ops)
    num = 4 ** n
    lv = np.zeros((num, num), dtype=complex)
    for b in range(num):
        comm = -1j * (h @ pm[b] - pm[b] @ h)
        for a in range(num):
            lv[a, b] = np.trace(pm[a] @ comm) / d
    for a, s in enumerate(idx):
        lv[a, a] += -2 * gamma * sum(1 for i in s if i in (1, 2))
    return lv, pm, d


def main():
    n, gamma = 3, 0.05
    sigma = n * gamma
    lv, pm, d = build(n, gamma)
    evals, r = np.linalg.eig(lv)
    osc = np.abs(np.imag(evals + sigma)) > 1e-6
    p_osc = r @ np.diag(osc.astype(float)) @ np.linalg.inv(r)

    psi = np.zeros(d, dtype=complex)
    psi[0] = psi[6] = 1 / np.sqrt(2)

    def fingerprint(state):
        v = np.array([np.trace(p @ np.outer(state, state.conj())) for p in pm])
        return np.abs(p_osc @ v) / d

    base = fingerprint(psi)
    xn = kron_all([SX] * n)
    lines = [
        "Projected-component fingerprint, N=3 Heisenberg chain, gamma=0.05, Bell(0,1)",
        f"oscillating modes: {int(osc.sum())} of {4 ** n}; smallest oscillating |Im|: "
        f"{np.min(np.abs(np.imag(evals[osc]))):.6f}; largest non-oscillating |Im|: "
        f"{np.max(np.abs(np.imag(evals[~osc]))):.3e}",
        "",
        f"{'error':>6}  {'max change':>11}  site-2 identity E|psi> = phase * X^N|psi> (exact)",
    ]
    for site in range(n):
        for op, name in ((SX, "X"), (SY, "Y"), (SZ, "Z")):
            ops = [I2] * n
            ops[site] = op
            e = kron_all(ops)
            err = e @ psi
            change = float(np.max(np.abs(fingerprint(err) - base)))
            note = ""
            if site == n - 1:
                target = xn @ psi if name != "Z" else psi
                nz = np.flatnonzero(target)
                phase = err[nz[0]] / target[nz[0]]
                exact = bool(np.array_equal(err, phase * target)) and abs(abs(phase) - 1) == 0.0
                note = f"  {'holds' if exact else 'FAILS'} (phase {phase:.0f})"
            lines.append(f"{name}_{site:<4}  {change:>11.6f}{note}")
    text = "\n".join(lines) + "\n"
    print(text, end="")
    OUT.write_text(text, encoding="utf-8", newline="\n")


if __name__ == "__main__":
    main()
