#!/usr/bin/env python3
"""EQ-030 deeper Winkel: why are XY+YX and IY+YI T1-robust while YZ+ZY,
XZ+ZX (also bond-flipped) collapse?

Algebraic framing (framework-grounded):
  - Pure-Z protected observables are protected because S_λ(α) = Σ_{k∈cluster}
    V[α,k]·c[k] = 0, i.e., individual eigenmode contributions cancel within
    each *degenerate cluster* of L_Z.
  - Adding T1 = γ_T1 D_T1 perturbs the spectrum. To first order in γ_T1,
    each cluster of L_Z splits into eigenvalues offset by the eigenvalues of
    the perturbation matrix L_T1|_cluster (the projection of L_T1 onto the
    cluster's eigenspace, in L_Z's eigenbasis).
  - If L_T1|_cluster ∝ I_cluster (degenerate within the cluster), the cluster
    remains glued at γ_T1 > 0; the cancellation S_λ(α) = 0 survives, the
    observable stays protected.
  - If L_T1|_cluster has non-trivial structure (off-diagonal in V's basis,
    or distinct eigenvalues), the cluster splits, individual eigenmodes
    decouple, and the cancellation breaks for any α that relied on it.

This script computes, for each soft Hamiltonian:
  1. L_Z eigenvalue clusters
  2. L_T1 in L_Z's eigenbasis
  3. Per-cluster: off-diagonal Frobenius norm within the cluster (lifts
     the degeneracy and breaks cancellation)
  4. Cross-references against the observed drop in protected count

Hypothesis: cumulative cluster-splitting weight should correlate with the
observed drop. Robust softs (XY+YX, IY+YI) should have low cumulative
splitting; fragile ones (YZ+ZY, XZ+ZX) should have high.
"""
import math
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


def cluster_splitting_breakdown(H, gamma_l, gamma_t1_l, N, cluster_tol=1e-8):
    """For each cluster of L_Z, compute how strongly L_T1 lifts the degeneracy.

    Returns:
      clusters: list of (eigenvalue, multiplicity) pairs.
      block_norms: list of (cluster_idx, frob_norm_full, frob_offdiag) per cluster.
        - frob_norm_full = ‖L_T1|_cluster‖_F (full block weight)
        - frob_offdiag   = ‖L_T1|_cluster − diag(L_T1|_cluster)‖_F
        - frob_eigsplit  = std of eigenvalues of L_T1|_cluster (degeneracy lift)
      total_offdiag_weight: cumulative off-diagonal weight across all clusters
        of size ≥ 2.
      total_eigsplit_weight: cumulative spread of perturbation eigenvalues
        across all clusters of size ≥ 2.
    """
    L_Z = fw.lindbladian_z_dephasing(H, gamma_l)
    L_full = fw.lindbladian_z_plus_t1(H, gamma_l, gamma_t1_l)
    L_T1 = (L_full - L_Z) / max(gamma_t1_l[0], 1e-12)  # rescale to per-unit γ_T1

    # Move to Pauli basis (diagonalisation works there cleaner)
    M_basis = fw._vec_to_pauli_basis_transform(N)
    L_Z_pauli = (M_basis.conj().T @ L_Z @ M_basis) / (2 ** N)
    L_T1_pauli = (M_basis.conj().T @ L_T1 @ M_basis) / (2 ** N)

    evals, V = np.linalg.eig(L_Z_pauli)
    Vinv = np.linalg.inv(V)

    L_T1_eig = Vinv @ L_T1_pauli @ V

    n = len(evals)
    used = np.zeros(n, dtype=bool)
    clusters = []
    for i in range(n):
        if used[i]:
            continue
        cl = [i]
        used[i] = True
        for j in range(i + 1, n):
            if not used[j] and abs(evals[j] - evals[i]) < cluster_tol:
                cl.append(j)
                used[j] = True
        clusters.append(cl)

    block_norms = []
    total_offdiag_weight = 0.0
    total_eigsplit_weight = 0.0
    for ci, cl in enumerate(clusters):
        if len(cl) < 2:
            continue
        sub = L_T1_eig[np.ix_(cl, cl)]
        frob_full = float(np.linalg.norm(sub, 'fro'))
        sub_offdiag = sub - np.diag(np.diag(sub))
        frob_offdiag = float(np.linalg.norm(sub_offdiag, 'fro'))
        sub_eigs = np.linalg.eigvals(sub)
        eigsplit = float(np.std(sub_eigs.real)) + float(np.std(sub_eigs.imag))
        block_norms.append({
            'cluster_idx': ci,
            'eigenvalue': complex(evals[cl[0]]),
            'size': len(cl),
            'frob_full': frob_full,
            'frob_offdiag': frob_offdiag,
            'eigsplit': eigsplit,
        })
        total_offdiag_weight += frob_offdiag
        total_eigsplit_weight += eigsplit

    return clusters, block_norms, total_offdiag_weight, total_eigsplit_weight


def count_protected(H, gamma_l, gamma_t1_l, rho_0, N,
                    threshold=1e-9, cluster_tol=1e-8):
    if any(g != 0 for g in gamma_t1_l):
        L = fw.lindbladian_z_plus_t1(H, gamma_l, gamma_t1_l)
    else:
        L = fw.lindbladian_z_dephasing(H, gamma_l)
    M_basis = fw._vec_to_pauli_basis_transform(N)
    L_pauli = (M_basis.conj().T @ L @ M_basis) / (2 ** N)
    evals, V = np.linalg.eig(L_pauli)
    Vinv = np.linalg.inv(V)
    rho_pauli = fw.pauli_basis_vector(rho_0, N)
    c = Vinv @ rho_pauli
    n_eig = len(evals)
    used = np.zeros(n_eig, dtype=bool)
    clusters = []
    for i in range(n_eig):
        if used[i]:
            continue
        cl = [i]
        used[i] = True
        for j in range(i + 1, n_eig):
            if not used[j] and abs(evals[j] - evals[i]) < cluster_tol:
                cl.append(j)
                used[j] = True
        clusters.append(cl)
    n_protected = 0
    for alpha in range(1, 4 ** N):
        max_S = 0.0
        for cl in clusters:
            S = sum(V[alpha, k] * c[k] for k in cl)
            max_S = max(max_S, abs(S))
        if max_S < threshold:
            n_protected += 1
    return n_protected


def main():
    N = 3
    GAMMA_DEPH = 0.1
    GAMMA_T1 = 0.01  # 0.1 * γ_deph
    J = 1.0
    bonds = [(i, i + 1) for i in range(N - 1)]

    plus = np.array([1, 1], dtype=complex) / math.sqrt(2)
    minus = np.array([1, -1], dtype=complex) / math.sqrt(2)
    psi = np.kron(plus, np.kron(minus, plus))
    rho_0 = np.outer(psi, psi.conj())

    cases = [
        ('IY+YI',   [('I', 'Y', J), ('Y', 'I', J)]),  # robust, drop=0
        ('XY+YX',   [('X', 'Y', J), ('Y', 'X', J)]),  # robust, drop=1
        ('YZ+ZY',   [('Y', 'Z', J), ('Z', 'Y', J)]),  # fragile, drop=28
        ('XZ+ZX',   [('X', 'Z', J), ('Z', 'X', J)]),  # fragile, drop=29
        ('XY+XY',   [('X', 'Y', J), ('X', 'Y', J)]),  # repeat, drop=15
        ('YZ+ZZ',   [('Y', 'Z', J), ('Z', 'Z', J)]),  # asymmetric, drop=6
    ]

    print(f"T1 cluster-splitting analysis at N={N}, |+−+⟩, "
          f"γ_deph={GAMMA_DEPH}, γ_T1={GAMMA_T1}")
    print(f"For each soft Hamiltonian, computes how T1 perturbation lifts the")
    print(f"degeneracy within each L_Z eigenvalue cluster (cluster split = "
          f"protection breaks).")
    print()

    print(f"  {'case':<8s}  {'n_pure':>7s}  {'n_T1':>5s}  {'drop':>5s}  "
          f"{'#clust':>6s}  {'#multi':>6s}  {'Σ‖offdiag‖':>11s}  "
          f"{'Σ eigsplit':>11s}")
    print(f"  {'-' * 8}  {'-' * 7}  {'-' * 5}  {'-' * 5}  {'-' * 6}  "
          f"{'-' * 6}  {'-' * 11}  {'-' * 11}")

    summary = []
    for label, terms in cases:
        H = fw._build_bilinear(N, bonds, terms)
        n_pure = count_protected(H, [GAMMA_DEPH] * N, [0.0] * N, rho_0, N)
        n_t1 = count_protected(
            H, [GAMMA_DEPH] * N, [GAMMA_T1] * N, rho_0, N
        )
        drop = n_pure - n_t1

        clusters, block_norms, total_offdiag, total_eigsplit = (
            cluster_splitting_breakdown(
                H, [GAMMA_DEPH] * N, [GAMMA_T1] * N, N
            )
        )
        n_clusters = len(clusters)
        n_multi_clusters = len(block_norms)

        print(f"  {label:<8s}  {n_pure:>7d}  {n_t1:>5d}  {drop:>5d}  "
              f"{n_clusters:>6d}  {n_multi_clusters:>6d}  "
              f"{total_offdiag:>11.4e}  {total_eigsplit:>11.4e}")
        summary.append({
            'label': label,
            'drop': drop,
            'total_offdiag': total_offdiag,
            'total_eigsplit': total_eigsplit,
        })

    print()
    print("Per-case cluster breakdown (top 5 perturbed clusters by off-diag weight):")
    for label, terms in cases:
        H = fw._build_bilinear(N, bonds, terms)
        clusters, block_norms, total_offdiag, total_eigsplit = (
            cluster_splitting_breakdown(
                H, [GAMMA_DEPH] * N, [GAMMA_T1] * N, N
            )
        )
        block_norms.sort(key=lambda b: -b['frob_offdiag'])
        print(f"\n  {label}:")
        if not block_norms:
            print(f"    (no degenerate clusters of size ≥ 2)")
            continue
        print(f"    {'cluster λ':>20s}  {'size':>5s}  {'‖full‖':>10s}  "
              f"{'‖offdiag‖':>10s}  {'eigsplit':>10s}")
        for b in block_norms[:5]:
            ev = b['eigenvalue']
            ev_str = f"{ev.real:>+8.3f}{ev.imag:+8.3f}i"
            print(f"    {ev_str:>20s}  {b['size']:>5d}  "
                  f"{b['frob_full']:>10.4e}  {b['frob_offdiag']:>10.4e}  "
                  f"{b['eigsplit']:>10.4e}")

    print()
    print("Hypothesis check: does Σ‖offdiag‖ correlate with drop?")
    drops = np.array([s['drop'] for s in summary])
    offdiags = np.array([s['total_offdiag'] for s in summary])
    eigsplits = np.array([s['total_eigsplit'] for s in summary])

    if drops.std() > 0 and offdiags.std() > 0:
        rho_offdiag = float(np.corrcoef(drops, offdiags)[0, 1])
        rho_eigsplit = float(np.corrcoef(drops, eigsplits)[0, 1])
        print(f"  Pearson(drop, Σ‖offdiag‖)   = {rho_offdiag:+.4f}")
        print(f"  Pearson(drop, Σ eigsplit)   = {rho_eigsplit:+.4f}")
    print()
    print("Interpretation:")
    print("  Σ‖offdiag‖ measures total weight of T1 in cluster-coupling form")
    print("  (degeneracy-lifting). High → many clusters split → many protected")
    print("  observables leak. Low → most clusters stay glued → robustness.")


if __name__ == "__main__":
    main()
