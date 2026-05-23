"""Mega-fine scan of the bath temperature, from pure decay down into the minus.

docs/carbon master-question thread (2026-05-22, Tom + Claude).

Tom's picture: there is no cold, cold is only the absence of heat. The
wave-breaking is what creates heat; before it, the bath has no heat. So the
bath's base state is not n_bar=0, it is below it, "in the minus": a negative
temperature, which physically is not cold but a source, a pump (a negative-
temperature bath is hotter than any positive one). In the framework that is
Sum-gamma < 0, the gain side, the active vacuum.

The bosonic n_bar of _v_effect_heat_modes.py cannot represent this: n_bar in
[0, infinity) has only one side, a bosonic bath cannot be inverted. Negative
temperature needs a bounded system; our qubits are bounded, so it works, but
the bath channel must be rebuilt with two independent rates instead of n_bar:

  decay  sigma-  at rate gamma_down  (lowers |1> -> |0>)
  pump   sigma+  at rate gamma_up    (raises |0> -> |1>)

A single axis s in [+1, -1] sweeps the temperature, total rate fixed:
  gamma_down = g_tot (1+s)/2,  gamma_up = g_tot (1-s)/2
  s = +1  pure decay   (T = 0+, the old floor of _v_effect_heat_modes.py)
  s =  0  balance      (gamma_down = gamma_up, T = +-infinity)
  s = -1  pure pump    (T = 0-, deep in the minus)
  s <  0  the minus    (pump dominates, negative temperature)

Hamiltonian XX+YY, a non-breaker, so the bath is the only actor and visible.
Z-dephasing stays on (small, fixed): the framework's proven-palindrome anchor.

The question, "how does a qubit get a baby": scanning s downward, where does a
new palindromic pair first appear, a new mirror. This first scan maps the
landscape (modes, pairs, break magnitude, centre, slowest oscillation vs s) to
locate the balance and whatever happens in the minus.

build_L reuses the dissipator structure of self_heating_fixpoint.py; best-centre
palindrome analysis from peierls_break_structure.py. Investigation only.
"""
import sys

import numpy as np
from scipy.optimize import minimize_scalar

sys.stdout.reconfigure(encoding="utf-8")

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
DOWN = np.array([[0, 1], [0, 0]], dtype=complex)   # sigma- : |1> -> |0>
UP = np.array([[0, 0], [1, 0]], dtype=complex)     # sigma+ : |0> -> |1>
PM = {"X": X, "Y": Y, "Z": Z}


def site_op(op, k, N):
    m = np.array([[1.0 + 0j]])
    for i in range(N):
        m = np.kron(m, op if i == k else I2)
    return m


def chain_H(N, combo="XX+YY", J=1.0):
    """The 2-term combo on every chain bond; XX+YY is a non-breaker."""
    t1, t2 = combo.split("+")
    d = 2 ** N
    H = np.zeros((d, d), dtype=complex)
    for b in range(N - 1):
        for t in (t1, t2):
            H += J * site_op(PM[t[0]], b, N) @ site_op(PM[t[1]], b + 1, N)
    return H


def build_L_gain(H, gamma_z, gamma_down, gamma_up):
    """-i[H,.] + Z-dephasing + decay (sigma-) + pump (sigma+), independent rates.
    Dissipator structure from self_heating_fixpoint.py."""
    d = H.shape[0]
    Id = np.eye(d, dtype=complex)
    n = int(np.log2(d))
    L = -1j * (np.kron(Id, H) - np.kron(H.T, Id))
    for k in range(n):
        jumps = []
        if gamma_z > 0:
            jumps.append(np.sqrt(gamma_z) * site_op(Z, k, n))
        if gamma_down > 0:
            jumps.append(np.sqrt(gamma_down) * site_op(DOWN, k, n))
        if gamma_up > 0:
            jumps.append(np.sqrt(gamma_up) * site_op(UP, k, n))
        for Lop in jumps:
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


def analyse(N, s, g_tot=0.1, gz=0.05):
    g_down = g_tot * (1.0 + s) / 2.0
    g_up = g_tot * (1.0 - s) / 2.0
    H = chain_H(N)
    L = build_L_gain(H, gz, g_down, g_up)
    ev = np.linalg.eigvals(L)
    n_freq = len(osc_frequencies(ev))
    lo, hi = ev.real.min() - 1.0, ev.real.max() + 1.0
    opt = minimize_scalar(lambda c: residuals(ev, c).max(),
                          bounds=(lo, hi), method="bounded")
    centre = opt.x
    gaps = residuals(ev, centre)
    width = float(ev.real.max() - ev.real.min())
    tol = 1e-3 * max(width, 1.0)
    paired = int(np.sum(gaps < tol))
    rates = -ev.real
    osc_mask = np.abs(ev.imag) > 1e-6
    slow_osc = float(np.min(rates[osc_mask])) if osc_mask.any() else 0.0
    return n_freq, paired, centre, opt.fun, slow_osc


def regime(s):
    if s > 0.02:
        return "decay"
    if s < -0.02:
        return "MINUS"
    return "balance"


if __name__ == "__main__":
    print("Bath-temperature scan: pure decay -> balance -> the minus.")
    print("s=+1 pure decay (old floor), s=0 balance, s<0 the minus (pump).\n")
    s_values = [round(1.0 - 0.05 * i, 4) for i in range(41)]   # +1.0 .. -1.0
    for N in (3, 4, 5):
        print(f"--- N = {N}   (L is {4 ** N}x{4 ** N}, g_tot=0.1, gz=0.05) ---")
        print(f"{'s':>6} {'regime':>8} {'#freq':>6} {'paired':>7} "
              f"{'centre':>9} {'max res':>9} {'slow osc':>9}")
        for s in s_values:
            nf, p, c, mr, so = analyse(N, s)
            print(f"{s:>6.2f} {regime(s):>8} {nf:>6} {p:>7} "
                  f"{c:>9.4f} {mr:>9.4f} {so:>9.5f}")
        print()
