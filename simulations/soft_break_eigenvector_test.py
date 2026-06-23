#!/usr/bin/env python3
"""Test the soft-break prediction: eigenvalue vs eigenvector pairing.

V_EFFECT_FINE_STRUCTURE found 3 truly-unbroken + 19 soft-broken + 14 hard-
broken among the 36 two-term Pauli-pair Hamiltonians. The framework's
prediction was: in soft-broken cases, eigenvalues pair (spectroscopy
sees palindromic) but eigenvectors do not pair under Π (an off-diagonal
operator-level break invisible to spectrum).

This script tests that prediction directly.

Two tests per Hamiltonian:

A. Eigenvalue pairing (spectral level): for each eigenvalue λ_i of L,
   find j minimising |λ_j − (−λ_i − 2Σγ)|. Report the worst pairing
   error across all eigenvalues.

B. Eigenvector pairing (operator level): for each eigenvalue pair (i, j),
   compute the overlap |⟨v_j | Π · v_i⟩| / (‖v_j‖ ‖Π v_i‖). For perfect
   palindrome (M = 0): overlap = 1. For soft-break: overlap < 1 even when
   eigenvalues pair.

Predicted outcome:
  3 truly-unbroken: eigenvalue pairing OK, eigenvector overlap = 1 (within tol)
  19 soft-broken:   eigenvalue pairing OK, eigenvector overlap < 1
  14 hard-broken:   eigenvalue pairing fails

If this prediction holds, the soft-break category is operationally
distinguishable from truly-unbroken via eigenvector probing, even though
spectroscopy alone treats them identically.
"""
import math
import sys
from itertools import combinations

import numpy as np

import framework as fw

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass


SINGLE_TERMS = ['XX', 'XY', 'XZ', 'YX', 'YY', 'YZ', 'ZX', 'ZY', 'ZZ']


def pair_eigenvalues(evals, Sigma_gamma, tol=1e-6):
    """Pair eigenvalues as (λ_i, λ_j) with λ_j ≈ -λ_i - 2Σγ.

    Returns:
        pairs: list of (i, j) index pairs
        pairing_errors: max |λ_j − (−λ_i − 2Σγ)| across all pairs
        unpaired: indices that have no partner within tol
    """
    n = len(evals)
    used = np.zeros(n, dtype=bool)
    pairs = []
    pairing_errors = []
    unpaired = []
    for i in range(n):
        if used[i]:
            continue
        target = -evals[i] - 2 * Sigma_gamma
        # Find best unused match
        best_j = -1
        best_dist = float('inf')
        for j in range(n):
            if used[j]:
                continue
            dist = abs(evals[j] - target)
            if dist < best_dist:
                best_dist = dist
                best_j = j
        if best_j != -1 and best_dist < tol:
            used[i] = True
            if best_j != i:
                used[best_j] = True
            pairs.append((i, best_j))
            pairing_errors.append(best_dist)
        else:
            unpaired.append(i)
    return pairs, pairing_errors, unpaired


def eigenvector_overlap(v_i, v_j, Pi):
    """Overlap |⟨v_j | Π · v_i⟩| / (‖v_j‖ · ‖Π v_i‖).

    1.0 = perfect pairing (Π takes v_i exactly to v_j up to phase)
    < 1.0 = mixing into other eigenvectors
    """
    Pi_v_i = Pi @ v_i
    norm_Pi_vi = np.linalg.norm(Pi_v_i)
    norm_vj = np.linalg.norm(v_j)
    if norm_Pi_vi < 1e-15 or norm_vj < 1e-15:
        return 0.0
    overlap = abs(np.vdot(v_j, Pi_v_i)) / (norm_Pi_vi * norm_vj)
    return float(overlap)


def subspace_overlap(v_i, partner_subspace, Pi):
    """For degenerate eigenvalues: project Π v_i onto the partner eigenspace.

    Returns the norm of the projection (1.0 if Π v_i lies entirely in the
    partner subspace, < 1.0 if it has components elsewhere).
    """
    Pi_v_i = Pi @ v_i
    norm_Pi_vi = np.linalg.norm(Pi_v_i)
    if norm_Pi_vi < 1e-15:
        return 0.0
    # Project onto partner_subspace (columns are orthonormal eigenvectors)
    if partner_subspace.shape[1] == 0:
        return 0.0
    # Gram-Schmidt orthonormalise partner_subspace columns
    Q, _ = np.linalg.qr(partner_subspace)
    proj = Q @ (Q.conj().T @ Pi_v_i)
    return float(np.linalg.norm(proj) / norm_Pi_vi)


def analyze_combo(N, gamma_l, term1, term2, Pi):
    """Run both tests on a single combo. Returns dict of metrics."""
    Sigma_gamma = sum(gamma_l)
    bonds = [(0, 1), (1, 2)]
    H = fw._build_bilinear(
        N, bonds,
        [(term1[0], term1[1], 1.0), (term2[0], term2[1], 1.0)]
    )
    L = fw.lindbladian_z_dephasing(H, gamma_l)
    M = fw._vec_to_pauli_basis_transform(N)
    L_pauli = (M.conj().T @ L @ M) / (2 ** N)

    evals, evecs = np.linalg.eig(L_pauli)

    # Test A: eigenvalue pairing
    pairs, pairing_errors, unpaired = pair_eigenvalues(evals, Sigma_gamma, tol=1e-4)
    max_pair_err = max(pairing_errors) if pairing_errors else 0.0
    n_unpaired = len(unpaired)

    # Test B: eigenvector pairing (using subspace approach for degeneracies)
    # Group eigenvalues into degenerate clusters
    eval_clusters = []
    cluster_tol = 1e-6
    used_cluster = np.zeros(len(evals), dtype=bool)
    for i in range(len(evals)):
        if used_cluster[i]:
            continue
        cluster = [i]
        used_cluster[i] = True
        for j in range(i + 1, len(evals)):
            if not used_cluster[j] and abs(evals[j] - evals[i]) < cluster_tol:
                cluster.append(j)
                used_cluster[j] = True
        eval_clusters.append(cluster)

    # For each pair (i, j), find eigenvector j's cluster and project Π v_i
    overlaps = []
    if n_unpaired == 0:
        for i, j in pairs:
            # Find cluster containing j
            cluster_with_j = next((c for c in eval_clusters if j in c), None)
            if cluster_with_j is None:
                continue
            partner_subspace = evecs[:, cluster_with_j]
            v_i = evecs[:, i]
            ov = subspace_overlap(v_i, partner_subspace, Pi)
            overlaps.append(ov)

    min_overlap = min(overlaps) if overlaps else None
    avg_overlap = sum(overlaps) / len(overlaps) if overlaps else None

    return {
        'term1': term1, 'term2': term2,
        'eigenvalue_pairing_max_err': max_pair_err,
        'unpaired_count': n_unpaired,
        'eigenvalue_paired': n_unpaired == 0 and max_pair_err < 1e-4,
        'eigenvector_overlap_min': min_overlap,
        'eigenvector_overlap_avg': avg_overlap,
        'n_pairs_tested': len(overlaps),
    }


def main():
    N = 3
    gamma = 0.1
    gamma_l = [gamma] * N
    Sigma_gamma = sum(gamma_l)
    Pi = fw.build_pi_full(N)

    op_threshold = 1e-10
    eigenvalue_pair_threshold = 1e-4
    eigenvector_overlap_threshold = 1.0 - 1e-4

    print("=" * 100)
    print(f"Soft-break test: eigenvalue vs eigenvector pairing (N={N}, γ={gamma})")
    print(f"For each of 36 combos: do eigenvalues pair? Do eigenvectors pair under Π?")
    print(f"Prediction:")
    print(f"  3 truly-unbroken: both YES")
    print(f"  19 soft-broken: eigenvalues YES, eigenvectors NO (overlap < 1)")
    print(f"  14 hard-broken: eigenvalues NO")
    print("=" * 100)

    combos = list(combinations(SINGLE_TERMS, 2))
    results = []
    for term1, term2 in combos:
        r = analyze_combo(N, gamma_l, term1, term2, Pi)
        results.append(r)

    # Categorize
    truly_unbroken = []
    soft_broken = []
    hard_broken = []
    for r in results:
        if r['eigenvalue_paired']:
            if r['eigenvector_overlap_min'] is not None and r['eigenvector_overlap_min'] >= eigenvector_overlap_threshold:
                truly_unbroken.append(r)
            else:
                soft_broken.append(r)
        else:
            hard_broken.append(r)

    print(f"\nResults:")
    print(f"  truly unbroken (both pairings intact):  {len(truly_unbroken)}")
    print(f"  soft broken (eigenvalues pair, vectors don't):  {len(soft_broken)}")
    print(f"  hard broken (eigenvalues don't pair):  {len(hard_broken)}")
    print(f"  Total: {len(truly_unbroken) + len(soft_broken) + len(hard_broken)}")

    print(f"\nExpected from V_EFFECT_FINE_STRUCTURE: 3 / 19 / 14")
    print(f"Match? truly={len(truly_unbroken)==3}, soft={len(soft_broken)==19}, hard={len(hard_broken)==14}")

    print("\n" + "=" * 100)
    print(f"TRULY UNBROKEN ({len(truly_unbroken)}): both eigenvalue and eigenvector pairing intact")
    print("=" * 100)
    print(f"{'Combo':>10s}  {'eval pair err':>16s}  {'evec overlap min':>18s}  {'evec overlap avg':>18s}")
    for r in truly_unbroken:
        ov_min = f"{r['eigenvector_overlap_min']:.6f}" if r['eigenvector_overlap_min'] is not None else "N/A"
        ov_avg = f"{r['eigenvector_overlap_avg']:.6f}" if r['eigenvector_overlap_avg'] is not None else "N/A"
        print(f"{r['term1']+'+'+r['term2']:>10s}  {r['eigenvalue_pairing_max_err']:>16.4e}  {ov_min:>18s}  {ov_avg:>18s}")

    print("\n" + "=" * 100)
    print(f"SOFT BROKEN ({len(soft_broken)}): eigenvalues pair, eigenvectors do NOT")
    print("=" * 100)
    print(f"{'Combo':>10s}  {'eval pair err':>16s}  {'evec overlap min':>18s}  {'evec overlap avg':>18s}")
    for r in soft_broken:
        ov_min = f"{r['eigenvector_overlap_min']:.6f}" if r['eigenvector_overlap_min'] is not None else "N/A"
        ov_avg = f"{r['eigenvector_overlap_avg']:.6f}" if r['eigenvector_overlap_avg'] is not None else "N/A"
        print(f"{r['term1']+'+'+r['term2']:>10s}  {r['eigenvalue_pairing_max_err']:>16.4e}  {ov_min:>18s}  {ov_avg:>18s}")

    print("\n" + "=" * 100)
    print(f"HARD BROKEN ({len(hard_broken)}): eigenvalues fail to pair")
    print("=" * 100)
    print(f"{'Combo':>10s}  {'eval pair err':>16s}  {'unpaired count':>16s}")
    for r in hard_broken:
        print(f"{r['term1']+'+'+r['term2']:>10s}  {r['eigenvalue_pairing_max_err']:>16.4e}  {r['unpaired_count']:>16d}")

    print()
    print("=" * 100)
    print("Reading: if soft-broken cases have eigenvector overlap MIN substantially")
    print("below 1.0 (e.g., 0.7-0.99) while truly-unbroken have overlap = 1.0 exactly,")
    print("the prediction is verified: soft-break is operationally distinguishable")
    print("from truly-unbroken via eigenvector probing, despite identical spectra.")
    print("=" * 100)


if __name__ == "__main__":
    main()
