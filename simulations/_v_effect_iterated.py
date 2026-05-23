"""V-Effect iterated: do the NEW-NEW pairs survive the next step, or get wiped?

docs/carbon master-question first attempt (2026-05-22, Tom + Claude).

The reframe: to know where qubits "play carbon" we do not need a bridge to
carbon. We need the local stepping rule N -> N+1, because the rule iterated IS
the path. The V-Effect is that step: adding a bond changes the Liouvillian
frequency spectrum. pairing_structure.py found that at the first big coupling
(N=2 -> N=5, Heisenberg) every frequency is NEW relative to N=2, and the new
modes pair NEW-NEW.

The open question, never looked at: what happens to those NEW-NEW pairs at the
NEXT step? Do the frequencies of step N survive into step N+1, then a pair
carries forward, it acts as a qubit, the recursion climbs. Or are they wiped and
replaced the way the first coupling wiped the original frequencies, then each
step is a fresh start and there is no recursion.

This script iterates the incremental chain N = 2..6 and, per step, classifies
each oscillation frequency as OLD (present at N-1 within tolerance) or NEW. Two
Hamiltonians run side by side: XX+YY, a non-breaking palindrome-preserving
coupling, as a control; and XX+XY, a breaking coupling, one of the 14 that break
the palindrome at N=3 per v_effect_analysis.py. The control matters: if even the
non-breaking chain shows "all NEW every step", then matching frequencies by
value across N is confounded by the N-dependent frequency shift, and the OLD/NEW
split carries no meaning. If the control carries frequencies forward and the
breaker wipes them, the wipe is a real effect of the break.

Primitives reused from simulations/pairing_structure.py and v_effect_analysis.py.
Investigation only.
"""
import sys

import numpy as np

sys.stdout.reconfigure(encoding="utf-8")

I2 = np.eye(2, dtype=complex)
PM = {
    "I": I2,
    "X": np.array([[0, 1], [1, 0]], dtype=complex),
    "Y": np.array([[0, -1j], [1j, 0]], dtype=complex),
    "Z": np.array([[1, 0], [0, -1]], dtype=complex),
}
SZ = PM["Z"]


def site_op(op, k, N):
    m = np.array([[1.0 + 0j]])
    for i in range(N):
        m = np.kron(m, op if i == k else I2)
    return m


def chain_H(N, combo, J=1.0):
    """H = J * Sum_{bond} (term1 + term2): the 2-term combo on every chain bond.
    combo 'XX+XY' means term1 = XX, term2 = XY."""
    t1, t2 = combo.split("+")
    d = 2 ** N
    H = np.zeros((d, d), dtype=complex)
    for b in range(N - 1):
        for t in (t1, t2):
            H += J * site_op(PM[t[0]], b, N) @ site_op(PM[t[1]], b + 1, N)
    return H


def build_L(H, gamma, N):
    """Z-dephasing Lindbladian (pairing_structure.py convention)."""
    d = 2 ** N
    Id = np.eye(d, dtype=complex)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for k in range(N):
        Zk = site_op(SZ, k, N)
        L += gamma * (np.kron(Zk, Zk.conj()) - np.eye(d * d, dtype=complex))
    return L


def find_pairs(evals, Sg, tol=1e-4):
    """Palindromic partner for each eigenvalue: target = -(lambda + 2 Sg)."""
    n = len(evals)
    partners = np.full(n, -1)
    for k in range(n):
        target = -(evals[k] + 2 * Sg)
        dists = np.abs(evals - target)
        best = int(np.argmin(dists))
        if dists[best] < tol * max(1.0, abs(target)):
            partners[k] = best
    return partners


def osc_frequencies(evals, thresh=0.01):
    """Distinct oscillation frequencies |Im(lambda)| above thresh."""
    return sorted(set(round(abs(ev.imag), 6) for ev in evals
                      if abs(ev.imag) > thresh))


def classify_old_new(freqs, prev_freqs, tol=1e-3):
    """Each freq is OLD if it matches a previous-step freq within tol, else NEW."""
    if prev_freqs is None:
        return [], list(freqs)
    old, new = [], []
    for f in freqs:
        if any(abs(f - pf) < tol for pf in prev_freqs):
            old.append(f)
        else:
            new.append(f)
    return old, new


def run(combo, gamma=0.05, N_range=range(2, 7)):
    print(f"\n=== combo {combo}   (gamma = {gamma}) ===")
    print(f"{'N':>3} {'L size':>9} {'#freq':>6} {'paired':>8} {'orphan':>7} "
          f"{'OLD':>5} {'NEW':>5} {'NEW-NEW prs':>12}")
    prev_freqs = None
    for N in N_range:
        d = 2 ** N
        H = chain_H(N, combo)
        L = build_L(H, gamma, N)
        evals = np.linalg.eigvals(L)
        Sg = N * gamma
        partners = find_pairs(evals, Sg)
        n_paired = int(np.sum(partners >= 0))
        n_orphan = len(evals) - n_paired
        freqs = osc_frequencies(evals)
        old, new = classify_old_new(freqs, prev_freqs)
        new_set = set(new)
        nn = 0
        seen = set()
        for k in range(len(evals)):
            j = int(partners[k])
            if j < 0 or j == k or k in seen:
                continue
            seen.add(k)
            seen.add(j)
            fk = round(abs(evals[k].imag), 6)
            fj = round(abs(evals[j].imag), 6)
            if fk > 0.01 and fj > 0.01 and fk in new_set and fj in new_set:
                nn += 1
        print(f"{N:>3} {f'{d*d}x{d*d}':>9} {len(freqs):>6} {n_paired:>8} "
              f"{n_orphan:>7} {len(old):>5} {len(new):>5} {nn:>12}")
        prev_freqs = freqs


if __name__ == "__main__":
    print("V-Effect iterated: do the NEW-NEW pairs survive the next step?")
    print("XX+YY = non-breaking control; XX+XY = breaking coupling.")
    print("If the control also shows all-NEW each step, value-matching across N")
    print("is confounded by the frequency shift, read OLD/NEW with that caution.")
    run("XX+YY")
    run("XX+XY")
