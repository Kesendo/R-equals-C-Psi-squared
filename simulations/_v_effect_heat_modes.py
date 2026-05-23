"""Corrected first attempt: does the break-heat create paired modes or orphans?

docs/carbon master-question thread (2026-05-22, Tom + Claude).

_v_effect_iterated.py had a real error, caught by Tom: it used a pure
Z-dephasing Liouvillian, no thermal channel, so it computed the break but not
the Bruchwaerme, the heat. And it grew N by hand, so the iteration was not
heat-driven. THERMAL_BREAKING's self-heating loop: the orphan modes decay
twice as fast and that decay IS heat production; the heat itself then creates
new modes (n_bar > 0 drives 111 -> 445 frequencies at N=5).

This corrects it. Fixed N. The breaking Hamiltonian XX+XY (the Brecher: a real
break, real orphans, real Bruchwaerme). The THERMAL Liouvillian with the
sigma-/sigma+ channel and n_bar (the heat). n_bar is swept, 0 upward. The
decisive question: as the heat rises and creates new modes, are the new modes
palindromic PAIRS (new mirrors, the candidate qubit) or ORPHANS?

THERMAL_BREAKING left this open: it found the thermal channel seems to break
the pairing (Z+amplitude -> 0%) but flagged that as a possible
center-estimation artifact ("the palindromic status of thermal channels is
open"). So the pairing here is measured against the BEST-FIT reflection
centre, found per spectrum, not the assumed Z-dephasing centre.

build_L_thermal reused verbatim from simulations/self_heating_fixpoint.py;
best-centre palindrome analysis from simulations/carbon/peierls_break_structure.py.
Investigation only.
"""
import sys

import numpy as np
from scipy.optimize import minimize_scalar

sys.stdout.reconfigure(encoding="utf-8")

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
SP = np.array([[0, 1], [0, 0]], dtype=complex)
SM = np.array([[0, 0], [1, 0]], dtype=complex)
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


def residuals(ev, centre):
    """Per-eigenvalue distance to the nearest reflected partner 2*centre - lambda."""
    reflected = 2.0 * complex(centre) - ev
    return np.abs(ev[None, :] - reflected[:, None]).min(axis=1)


def analyse(N, n_bar, gz=0.05, ga=0.05):
    H = build_H_xxxy(N)
    L = build_L_thermal(H, gz, ga, n_bar)
    ev = np.linalg.eigvals(L)
    n_freq = len(osc_frequencies(ev))
    lo, hi = ev.real.min() - 1.0, ev.real.max() + 1.0
    opt = minimize_scalar(lambda c: residuals(ev, c).max(),
                          bounds=(lo, hi), method="bounded")
    centre = opt.x
    gaps = residuals(ev, centre)
    width = float(ev.real.max() - ev.real.min())
    tol = 1e-3 * max(width, 1.0)
    paired_mask = gaps < tol
    paired = int(np.sum(paired_mask))
    orphan = len(ev) - paired
    rates = -ev.real
    orph_rate = float(np.mean(rates[~paired_mask])) if orphan > 0 else 0.0
    pair_rate = float(np.mean(rates[paired_mask])) if paired > 0 else 0.0
    return n_freq, paired, orphan, centre, opt.fun, orph_rate, pair_rate


if __name__ == "__main__":
    print("Break-heat: do the heat-created modes pair (qubits) or orphan?")
    print("H = XX+XY (the Brecher), thermal channel on, n_bar swept = the heat.\n")
    for N in (3, 4, 5):
        print(f"--- N = {N}   (L is {4 ** N}x{4 ** N}, gz = ga = 0.05) ---")
        print(f"{'n_bar':>7} {'#freq':>6} {'paired':>8} {'orphan':>8} "
              f"{'centre':>9} {'max res':>10} {'orph.rate':>10} {'pair.rate':>10}")
        for n_bar in (0.0, 0.5, 1.0, 2.0, 5.0):
            nf, p, o, c, mr, orr, pr = analyse(N, n_bar)
            print(f"{n_bar:>7.1f} {nf:>6} {p:>8} {o:>8} {c:>9.4f} "
                  f"{mr:>10.4f} {orr:>10.4f} {pr:>10.4f}")
        print()
