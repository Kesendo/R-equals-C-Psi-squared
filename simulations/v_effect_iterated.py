"""V-Effect iterated: tolerance-matched frequency values across N=2..6.

The carbon master-question (2026-05-22, Tom + Claude) invites a local stepping
rule. This finite investigation compares six-decimal absolute-frequency bins,
not persistent modes: no eigenvector/projector transport is computed.
OLD means a previous-N frequency lies within 1e-3; NEW means none does.

XX+YY is the palindrome-preserving control and XX+XY the breaking comparison.
If the control also reads all-NEW, cross-N frequency matching is confounded by
the N-dependent shift. A different count between the two Hamiltonians remains
a frequency-value comparison; mode ancestry is an open question.
"""
import sys

import numpy as np

if __package__:
    from .pairing_structure import f1_partner, maximum_cardinality_f1_match
else:
    from pairing_structure import f1_partner, maximum_cardinality_f1_match

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
    print(f"{'N':>3} {'L size':>9} {'#freq':>6} {'matched':>8} {'unmatched':>9} "
          f"{'OLD':>5} {'NEW':>5} {'NEW entries':>12}")
    prev_freqs = None
    for N in N_range:
        d = 2 ** N
        H = chain_H(N, combo)
        L = build_L(H, gamma, N)
        evals = np.linalg.eigvals(L)
        Sg = N * gamma
        n_paired = maximum_cardinality_f1_match(evals, Sg, 1e-4)
        n_orphan = len(evals) - n_paired
        freqs = osc_frequencies(evals)
        old, new = classify_old_new(freqs, prev_freqs)
        new_set = set(new)
        n_new = sum(abs(ev.imag) > 0.01
                    and round(abs(ev.imag), 6) in new_set for ev in evals)
        print(f"{N:>3} {f'{d*d}x{d*d}':>9} {len(freqs):>6} {n_paired:>8} "
              f"{n_orphan:>9} {len(old):>5} {len(new):>5} {n_new:>12}")
        prev_freqs = freqs


def main():
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")
    print("V-Effect iterated: tolerance-matched frequency values across N")
    print("XX+YY = non-breaking control; XX+XY = breaking coupling.")
    print("If the control also shows all-NEW each step, value-matching across N")
    print("is confounded by the frequency shift, read OLD/NEW with that caution.")
    run("XX+YY")
    run("XX+XY")


if __name__ == "__main__":
    main()
