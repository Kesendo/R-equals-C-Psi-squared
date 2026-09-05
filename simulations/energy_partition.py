"""
Filtered-spectrum frequency and decay weight diagnostic.

For these Heisenberg chains under local Z-dephasing the full multiset pairs
by F1. Here matching excludes zero roots (|lambda| <= 1e-8), removing the
partners of roots near -2*sum_gamma. An unmatched filtered root therefore
does not demonstrate symmetry breaking. The nearest-target lookup below
tests support, not multiplicity, and is no independent F1 certificate.
Sums of |Im(lambda)| and |Re(lambda)| are spectral weights, not energies.

Script: simulations/energy_partition.py
Output: stdout (run with PYTHONIOENCODING=utf-8 on Windows)
"""
import numpy as np

# === Pauli matrices ===
I2 = np.eye(2, dtype=complex)
sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)


def site_op(op, site, N):
    ops = [I2] * N
    ops[site] = op
    r = ops[0]
    for o in ops[1:]:
        r = np.kron(r, o)
    return r


def heisenberg_H(N, J=1.0):
    d = 2 ** N
    H = np.zeros((d, d), dtype=complex)
    for i in range(N - 1):
        for P in [sx, sy, sz]:
            H += J * site_op(P, i, N) @ site_op(P, i + 1, N)
    return H


def build_liouvillian(H, gamma, N):
    d = 2 ** N
    Id = np.eye(d, dtype=complex)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for k in range(N):
        Zk = site_op(sz, k, N)
        L += gamma * (np.kron(Zk, Zk.conj()) - np.eye(d * d, dtype=complex))
    return L


def classify_modes(evals, Sg, tol=1e-6):
    """Mark whether the supplied array contains a nearby reflected target.

    Mirror formula: partner of λ is −(λ + 2Sγ)
    Targets can be reused; this support lookup is multiplicity-blind.
    Returns: paired_mask, unpaired_mask (boolean arrays)
    """
    n = len(evals)
    paired = np.zeros(n, dtype=bool)

    for k in range(n):
        target = -(evals[k] + 2 * Sg)
        dists = np.abs(evals - target)
        best = np.argmin(dists)
        if dists[best] < tol * max(1, abs(target)):
            paired[k] = True

    return paired, ~paired


def energy_partition(evals, paired_mask):
    """Sum frequency weight |Im(λ)| and decay weight |Re(λ)| by category."""
    freq = np.abs(evals.imag)
    decay = np.abs(evals.real)

    E_freq_paired = np.sum(freq[paired_mask])
    E_freq_unpaired = np.sum(freq[~paired_mask])
    E_freq_total = E_freq_paired + E_freq_unpaired

    E_decay_paired = np.sum(decay[paired_mask])
    E_decay_unpaired = np.sum(decay[~paired_mask])
    E_decay_total = E_decay_paired + E_decay_unpaired

    return {
        'freq_paired': E_freq_paired,
        'freq_unpaired': E_freq_unpaired,
        'freq_total': E_freq_total,
        'decay_paired': E_decay_paired,
        'decay_unpaired': E_decay_unpaired,
        'decay_total': E_decay_total,
    }


# =================================================================
# EXPERIMENT 1: Filtered-support scaling (N=2..5, fixed J and γ)
# =================================================================
print("=" * 70)
print("EXPERIMENT 1: Filtered-spectrum frequency weight")
print("For these Z-dephased Heisenberg chains the full multiset pairs (F1).")
print("This matching excludes zero roots (|lambda| <= 1e-8).")
print("Unmatched means no target remains after deletion; it is not F1 failure.")
print("The nearest-target lookup ignores multiplicity; weights are not energies.")
print("=" * 70)

J = 1.0
gamma = 0.1

for N in [2, 3, 4, 5]:
    d = 2 ** N
    d2 = d * d
    Sg = N * gamma  # total dephasing

    H = heisenberg_H(N, J)
    L = build_liouvillian(H, gamma, N)
    evals = np.linalg.eigvals(L)

    # Remove steady state(s) near zero
    nonzero_mask = np.abs(evals) > 1e-8
    ev = evals[nonzero_mask]

    paired, unpaired = classify_modes(ev, Sg)
    n_total = len(ev)
    n_paired = np.sum(paired)
    n_unpaired = np.sum(unpaired)
    pct = 100 * n_paired / n_total if n_total > 0 else 0

    E = energy_partition(ev, paired)

    freq_pct = 100 * E['freq_paired'] / E['freq_total'] if E['freq_total'] > 0 else 0
    decay_pct = 100 * E['decay_paired'] / E['decay_total'] if E['decay_total'] > 0 else 0

    print(f"\nN={N}  ({d2} modes, {n_total} non-steady)")
    print(f"  Filtered support matching: {n_paired}/{n_total} = {pct:.1f}%")
    print(f"  Oscillatory frequency weight (|Im(λ)|):")
    print(f"    Paired:   {E['freq_paired']:10.4f}  ({freq_pct:.1f}%)")
    print(f"    Unpaired: {E['freq_unpaired']:10.4f}  ({100-freq_pct:.1f}%)")
    print(f"    Total:    {E['freq_total']:10.4f}")
    print(f"  Decay weight (|Re(λ)|):")
    print(f"    Paired:   {E['decay_paired']:10.4f}  ({decay_pct:.1f}%)")
    print(f"    Unpaired: {E['decay_unpaired']:10.4f}  ({100-decay_pct:.1f}%)")
    print(f"    Total:    {E['decay_total']:10.4f}")


# =================================================================
# EXPERIMENT 2: γ sweep at N=3
# =================================================================
print("\n")
print("=" * 70)
print("EXPERIMENT 2: Dephasing sweep at N=3")
print("Filtered spectral weights as gamma/J changes; no state CΨ is computed.")
print("=" * 70)

N = 3
J = 1.0

header = f"{'γ':>8} {'J/γ':>8} {'pair%':>7} {'Efreq_p%':>9} {'Edecay_p%':>10} {'Efreq_tot':>10} {'Edecay_tot':>11}"
print(f"\n{header}")
print("-" * len(header))

for gamma in [0.001, 0.005, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0]:
    Sg = N * gamma
    H = heisenberg_H(N, J)
    L = build_liouvillian(H, gamma, N)
    evals = np.linalg.eigvals(L)

    ev = evals[np.abs(evals) > 1e-8]
    paired, unpaired = classify_modes(ev, Sg)

    n_total = len(ev)
    n_paired = np.sum(paired)
    pct = 100 * n_paired / n_total if n_total > 0 else 0

    E = energy_partition(ev, paired)
    fp = 100 * E['freq_paired'] / E['freq_total'] if E['freq_total'] > 0 else 0
    dp = 100 * E['decay_paired'] / E['decay_total'] if E['decay_total'] > 0 else 0

    ratio = J / gamma
    print(f"{gamma:8.3f} {ratio:8.1f} {pct:7.1f} {fp:9.1f} {dp:10.1f} {E['freq_total']:10.4f} {E['decay_total']:11.4f}")


# =================================================================
# EXPERIMENT 3: J/γ ratio sweep at N=3
# =================================================================
print("\n")
print("=" * 70)
print("EXPERIMENT 3: Coupling sweep at N=3 (fixed γ=1.0)")
print("Track filtered frequency and decay weights across coupling strengths")
print("=" * 70)

N = 3
gamma = 1.0

header = f"{'J':>8} {'J/γ':>8} {'pair%':>7} {'Efreq%':>8} {'Edecay%':>9} {'Efreq_tot':>10} {'Edecay_tot':>11} {'ratio_f/d':>10}"
print(f"\n{header}")
print("-" * len(header))

for J_val in [0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0, 20.0]:
    Sg = N * gamma
    H = heisenberg_H(N, J_val)
    L = build_liouvillian(H, gamma, N)
    evals = np.linalg.eigvals(L)

    ev = evals[np.abs(evals) > 1e-8]
    paired, unpaired = classify_modes(ev, Sg)

    n_total = len(ev)
    n_paired = np.sum(paired)
    pct = 100 * n_paired / n_total if n_total > 0 else 0

    E = energy_partition(ev, paired)
    fp = 100 * E['freq_paired'] / E['freq_total'] if E['freq_total'] > 1e-12 else 0
    dp = 100 * E['decay_paired'] / E['decay_total'] if E['decay_total'] > 1e-12 else 0
    ratio_fd = E['freq_total'] / E['decay_total'] if E['decay_total'] > 1e-12 else 0

    print(f"{J_val:8.2f} {J_val/gamma:8.2f} {pct:7.1f} {fp:8.1f} {dp:9.1f} {E['freq_total']:10.4f} {E['decay_total']:11.4f} {ratio_fd:10.4f}")

print("\n\nDone.")
