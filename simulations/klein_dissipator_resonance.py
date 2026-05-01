#!/usr/bin/env python3
"""Test the dissipator-resonance conjecture for F77 hardness.

Conjecture (from polarity-layer / Z⊗N-Brecher analysis 2026-05-01):
F77-hardness lives in the Klein cell that matches the dissipator letter.

  Z-dephasing → hard cases concentrate in Klein (0, 1)  [Z's Klein index]
  X-dephasing → hard cases concentrate in Klein (1, 0)  [X's Klein index]
  Y-dephasing → hard cases concentrate in Klein (1, 1)  [Y's Klein index]

Z-dephasing has been verified empirically (46/46 hard k=3 cases in Klein (0,1)).
This script tests the X- and Y- dephasing predictions.

Method: build L for each (Hamiltonian, dissipator) combination and apply the
spectrum-pair eigenvalue test (closest-match against λ → -λ - 2Σγ). This is
the soft/hard test independent of any specific Π operator: palindromic spectrum
is an intrinsic property of L's eigenvalue set.

The truly/soft distinction requires the right Π operator for each dissipator,
which we skip here — only soft/hard matters for the conjecture.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass

SCRIPT_DIR = Path(__file__).parent
sys.path.insert(0, str(SCRIPT_DIR))
import framework as fw
from framework.pauli import site_op, _build_kbody_chain


# ----------------------------------------------------------------------
# Spectrum-pair test (universal across dissipator letter)
# ----------------------------------------------------------------------

def spectrum_pair_max_err(L, sigma_gamma):
    """Max error in matching every eigenvalue λ of L to a partner −λ−2Σγ.

    For a "soft" (palindromic-spectrum) L this returns ~0.
    For a "hard" L it returns a finite value.
    """
    evals = np.linalg.eigvals(L)
    used = np.zeros(len(evals), dtype=bool)
    max_err = 0.0
    for i in range(len(evals)):
        if used[i]:
            continue
        target = -evals[i] - 2 * sigma_gamma
        dists = np.abs(evals - target)
        for j in range(len(evals)):
            if used[j]:
                dists[j] = np.inf
        best_j = int(np.argmin(dists))
        if best_j != i:
            used[i] = True
            used[best_j] = True
        else:
            used[i] = True
        max_err = max(max_err, float(dists[best_j]))
    return max_err


def lindbladian_pauli_dephasing(H, gamma_l, letter):
    """Single-letter dephasing Lindbladian for letter ∈ {'X', 'Y', 'Z'}.

    L = -i[H,·] + Σ_l γ_l (P_l ρ P_l - ρ),  P ∈ {X, Y, Z}.

    Each P is Hermitian with P² = I, so D[√γ·P_l]ρ = γ·(P_l ρ P_l - ρ).
    """
    if not np.allclose(H, H.conj().T):
        raise ValueError("Hamiltonian H must be Hermitian.")
    d = H.shape[0]
    N = int(round(np.log2(d)))
    Id = np.eye(d, dtype=complex)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for l, gamma in enumerate(gamma_l):
        if gamma == 0:
            continue
        Pl = site_op(N, l, letter)
        L = L + gamma * (np.kron(Pl, Pl.conj()) - np.kron(Id, Id))
    return L


def classify_pair_under_dephasing(N, terms, dephase_letter, gamma=0.1, J=1.0,
                                  spec_tol=1e-6):
    """Test whether (H from terms, dephase_letter dephasing) is F77-hard.

    Returns 'soft' if spectrum is palindromic about −Σγ, 'hard' otherwise.
    (Truly/soft distinction skipped.)
    """
    H = _build_kbody_chain(N, [tuple(t) + (J,) for t in terms])
    L = lindbladian_pauli_dephasing(H, [gamma] * N, dephase_letter)
    sigma_gamma = N * gamma
    err = spectrum_pair_max_err(L, sigma_gamma)
    return 'hard' if err > spec_tol else 'soft', err


# ----------------------------------------------------------------------
# Test pairs
# ----------------------------------------------------------------------

# Klein-cell representatives (k=3 Klein-homogeneous pairs, Y-par 0).
# Each pair has both terms in the indicated Klein cell.
TEST_PAIRS = {
    'Klein (0,0)': [('X', 'X', 'I'), ('Y', 'Y', 'I')],   # Mother sector
    'Klein (0,1)': [('I', 'I', 'Z'), ('Z', 'Z', 'Z')],   # Z-like (verified hard under Z-deph)
    'Klein (1,0)': [('I', 'I', 'X'), ('X', 'X', 'X')],   # X-like (rotated from Klein (0,1))
    'Klein (1,1)': [('I', 'I', 'Y'), ('Y', 'Y', 'Y')],   # Y-like (rotated from Klein (0,1))
}


def klein_index(letters):
    """Klein-Vierergruppe (Z₂×Z₂) index of a Pauli term."""
    bits = {'I': (0, 0), 'X': (1, 0), 'Y': (1, 1), 'Z': (0, 1)}
    a = sum(bits[L][0] for L in letters) % 2
    b = sum(bits[L][1] for L in letters) % 2
    return (a, b)


def y_parity(letters):
    return sum(1 for L in letters if L == 'Y') % 2


def enumerate_z2_homogeneous_k3(letters_alphabet=('I', 'X', 'Y', 'Z')):
    """Enumerate all unordered Klein-homogeneous + Y-par-homogeneous k=3 pairs."""
    from itertools import product
    seen = set()
    pairs = []
    for term1 in product(letters_alphabet, repeat=3):
        for term2 in product(letters_alphabet, repeat=3):
            if all(L == 'I' for L in term1) or all(L == 'I' for L in term2):
                continue
            if klein_index(term1) != klein_index(term2):
                continue
            if y_parity(term1) != y_parity(term2):
                continue
            sorted_terms = tuple(sorted([term1, term2]))
            if sorted_terms in seen:
                continue
            seen.add(sorted_terms)
            pairs.append([term1, term2])
    return pairs


def main():
    N = 4
    print(f"=== Dissipator-Resonance Test: F77-Hardness vs Klein cell ===")
    print(f"N = {N}, k = 3 (Z₂³-homogeneous pairs)")
    print(f"γ = 0.1 per site, J = 1.0")
    print()
    print("Conjecture: hardness lives in the dissipator-letter's Klein cell.")
    print("  Z's Klein index = (0, 1)  X's = (1, 0)  Y's = (1, 1)  I's = (0, 0)")
    print()

    # 1) Targeted 4×3 sanity check
    print("--- Part 1: Targeted 4×3 (one representative per Klein cell) ---")
    print(f"{'Hamiltonian':<22s} {'Z-deph':>14s} {'X-deph':>14s} {'Y-deph':>14s}")
    print("-" * 70)
    for klein_label, terms in TEST_PAIRS.items():
        results = {}
        for letter in ['Z', 'X', 'Y']:
            cls, err = classify_pair_under_dephasing(N, terms, letter)
            results[letter] = (cls, err)
        h_label = '+'.join([''.join(t) for t in terms])
        h_short = f"{klein_label} {h_label}"
        z_str = f"{results['Z'][0]} ({results['Z'][1]:.2e})"
        x_str = f"{results['X'][0]} ({results['X'][1]:.2e})"
        y_str = f"{results['Y'][0]} ({results['Y'][1]:.2e})"
        print(f"{h_short:<22s} {z_str:>14s} {x_str:>14s} {y_str:>14s}")
    print()

    # 2) Full Z₂³-homogeneous sweep × 3 dephasing letters
    print("--- Part 2: Full sweep over all Z₂³-homogeneous pairs ---")
    pairs = enumerate_z2_homogeneous_k3()
    print(f"Enumerated {len(pairs)} Klein-homogeneous + Y-par-homogeneous k=3 pairs.")
    print()

    # Count hardness by (Klein cell, dephasing letter)
    counts = {}  # (klein, letter) → {'soft': n, 'hard': n}
    for klein in [(0, 0), (0, 1), (1, 0), (1, 1)]:
        for letter in ['Z', 'X', 'Y']:
            counts[(klein, letter)] = {'soft': 0, 'hard': 0}

    import time
    start = time.time()
    for terms in pairs:
        klein = klein_index(terms[0])
        for letter in ['Z', 'X', 'Y']:
            cls, err = classify_pair_under_dephasing(N, terms, letter)
            counts[(klein, letter)][cls] += 1
    elapsed = time.time() - start

    # Print 4×3 hardness table
    print(f"Hardness counts by (Klein cell, dephasing letter)  [elapsed {elapsed:.1f}s]:")
    print()
    print(f"{'Klein cell':<14s} {'Z-deph':>14s} {'X-deph':>14s} {'Y-deph':>14s}  (#hard / #total)")
    print("-" * 75)
    for klein in [(0, 0), (0, 1), (1, 0), (1, 1)]:
        total_per_letter = {
            letter: counts[(klein, letter)]['soft'] + counts[(klein, letter)]['hard']
            for letter in ['Z', 'X', 'Y']
        }
        z_str = f"{counts[(klein, 'Z')]['hard']:>3d} / {total_per_letter['Z']:>3d}"
        x_str = f"{counts[(klein, 'X')]['hard']:>3d} / {total_per_letter['X']:>3d}"
        y_str = f"{counts[(klein, 'Y')]['hard']:>3d} / {total_per_letter['Y']:>3d}"
        klein_str = f"{klein}"
        print(f"{klein_str:<14s} {z_str:>14s} {x_str:>14s} {y_str:>14s}")

    print()
    print("Expected diagonal (conjecture):")
    print("  Z-deph: hardness only in Klein (0, 1)")
    print("  X-deph: hardness only in Klein (1, 0)")
    print("  Y-deph: hardness only in Klein (1, 1)")
    print("  Klein (0, 0): never hard (Mother sector)")


if __name__ == "__main__":
    main()
