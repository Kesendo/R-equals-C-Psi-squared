"""Investigate the cross-sector pair-sum-to-6γ regularity in path-2 (3-qubit) L_super.

Empirical observation (F89c Liouvillian eigendecomposition):
  Path-2 has 8 distinct rates {0, 2, 2.556, 2.889, 3.112, 3.444, 4, 6}γ
  Cross-sector pairs sum to 6γ:
    2.556 + 3.444 = 6.000   (originating in (SE, SE) and (SE, DE) respectively)
    2.889 + 3.112 = 6.000

Hypotheses to test:
  H1: bit-flip symmetry X⊗N maps (SE, SE) ↔ (DE, DE), so they have identical
      eigenvalues. Test: compute (SE, SE) and (DE, DE) sub-block eigenvalues.
  H2: the pair-sum equals (SE, SE) sub-block trace average + (SE, DE) sub-block
      trace average.
       (SE, SE) trace = -24γ over 9 modes → mean = -8/3 γ ≈ -2.67γ
       (SE, DE) trace = -30γ over 9 modes → mean = -10/3 γ ≈ -3.33γ
       Sum of means: -6γ ✓ matches the pair-sum
  H3: at 2-qubit (path-1) the analog sum should be different.
       (SE, SE) trace = -8γ over 4 modes → mean = -2γ
       (SE, DE) trace = -8γ over 4 modes → mean = -2γ  (only overlap pairs)
       Sum: -4γ, not -6γ. So the value depends on block size + n_diff structure.
  H4: there's a one-to-one bijection between (SE, SE) and (SE, DE) eigenvalues
      such that paired sums are constant 6γ. Or the pair-sum is only mean-of-means
      (not eigenvalue-pairing).

Test by computing both sub-block spectra exactly and checking eigenvalue structure.
"""
from __future__ import annotations

import sys

import numpy as np

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

J = 1.0
GAMMA = 1.0


def pauli_at(P: np.ndarray, site: int, n_qubits: int) -> np.ndarray:
    op = np.array([[1.0]], dtype=complex)
    I2 = np.eye(2, dtype=complex)
    for q in range(n_qubits):
        op = np.kron(op, P if q == site else I2)
    return op


X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)


def build_L_super(n_qubits: int, j: float = J, gamma: float = GAMMA) -> np.ndarray:
    """L_super on vec(ρ) for path-(n-1) H_B + uniform Z-dephasing."""
    n = n_qubits
    H = np.zeros((2**n, 2**n), dtype=complex)
    for b in range(n - 1):
        H += j * (pauli_at(X, b, n) @ pauli_at(X, b + 1, n)
                  + pauli_at(Y, b, n) @ pauli_at(Y, b + 1, n))
    d = 2**n
    Id = np.eye(d, dtype=complex)
    L = -1j * (np.kron(Id, H) - np.kron(H.T, Id))
    for l in range(n):
        Zl = pauli_at(Z, l, n)
        L += gamma * (np.kron(Zl.T, Zl) - np.kron(Id, Id))
    return L


def sector_indices(n_qubits: int, p_a: int, p_b: int) -> list[int]:
    """Indices in vec(ρ) (column-major) for ρ entries (a, b) with popcount(a)=p_a, popcount(b)=p_b."""
    d = 2**n_qubits
    indices = []
    for b in range(d):
        for a in range(d):
            if bin(a).count("1") == p_a and bin(b).count("1") == p_b:
                indices.append(b * d + a)  # column-major: vec(M)[b*d + a] = M[a, b]
    return indices


def analyze(n_qubits: int) -> None:
    print(f"\n## {n_qubits}-qubit block (path-{n_qubits-1})\n")
    L = build_L_super(n_qubits)

    # Verify popcount-conservation: L should be block-diagonal in (p_a, p_b) basis
    sectors = {}
    for p_a in range(n_qubits + 1):
        for p_b in range(n_qubits + 1):
            idx = sector_indices(n_qubits, p_a, p_b)
            if not idx:
                continue
            L_sub = L[np.ix_(idx, idx)]
            sectors[(p_a, p_b)] = L_sub
            # Check that the sub-block doesn't leak: L * vec_in_sector should stay in sector
            # We just trust the U(1) structure here.

    # Pick the sectors of interest: (SE=1, SE=1), (SE=1, DE=2), (DE=2, DE=2)
    interesting = [(1, 1), (1, 2), (2, 2), (0, 1), (0, 2), (2, 1)]
    spectra = {}
    for sec in interesting:
        if sec in sectors:
            eigs = np.linalg.eigvals(sectors[sec])
            rates = -eigs.real
            freqs = eigs.imag
            spectra[sec] = (rates, freqs)
            mean_rate = np.mean(rates)
            print(f"  Sector {sec}: dim={sectors[sec].shape[0]}, trace/γ={np.trace(sectors[sec]).real:.3f}, mean rate Γ/γ={mean_rate:.4f}")
            print(f"    distinct rates Γ/γ = {sorted(set(np.round(rates, 4)))}")

    # Test H1: (SE, SE) eigenvalues == (DE, DE) eigenvalues?
    if (1, 1) in spectra and (2, 2) in spectra:
        ss_rates = sorted(np.round(spectra[(1, 1)][0], 6))
        dd_rates = sorted(np.round(spectra[(2, 2)][0], 6))
        match = ss_rates == dd_rates
        print(f"\n  H1: (SE, SE) rates == (DE, DE) rates? {match}")
        if not match:
            print(f"    (SE, SE) rates: {ss_rates}")
            print(f"    (DE, DE) rates: {dd_rates}")

    # Test H2/H4: pair-sum structure between (SE, SE) and (SE, DE)
    if (1, 1) in spectra and (1, 2) in spectra:
        ss_rates = sorted(spectra[(1, 1)][0])
        sd_rates = sorted(spectra[(1, 2)][0])
        target_sum = -np.trace(sectors[(1, 1)]).real / sectors[(1, 1)].shape[0] - np.trace(sectors[(1, 2)]).real / sectors[(1, 2)].shape[0]
        print(f"\n  Mean(SE,SE) + Mean(SE,DE) = {target_sum:.4f}γ  (the predicted 'cross-sector mean-sum' value)")
        print(f"\n  Pair-sums between (SE, SE) and (SE, DE) eigenvalues (sorted by min rate):")
        # All pairs
        pair_sums = []
        for r_ss in ss_rates:
            for r_sd in sd_rates:
                pair_sums.append((r_ss + r_sd, r_ss, r_sd))
        unique_sums = sorted(set(round(p[0], 4) for p in pair_sums))
        print(f"    distinct pair-sums (rounded to 4dp): {unique_sums[:20]} {'...' if len(unique_sums) > 20 else ''}")
        # Filter pairs that sum to target
        matching = [(s, a, b) for s, a, b in pair_sums if abs(s - target_sum) < 1e-3]
        print(f"  Pairs with sum near {target_sum:.4f}γ: {len(matching)} pairs")
        for s, a, b in matching[:10]:
            print(f"    Γ/γ = {a:.4f} + {b:.4f} = {s:.4f}")


# Path-1 (2 qubits) for baseline
analyze(2)
# Path-2 (3 qubits) for the actual question
analyze(3)
# Path-3 (4 qubits) to see how it scales
analyze(4)
