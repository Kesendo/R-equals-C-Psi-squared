"""Π-protected observables — algebraic skeleton of which Pauli expectations stay zero.

Public API:
  pi_protected_observables(H, gamma_l, rho_0, N, threshold, cluster_tol)
"""
from __future__ import annotations

import numpy as np

from .lindblad import lindbladian_z_dephasing
from .pauli import (
    _vec_to_pauli_basis_transform, pauli_basis_vector, _pauli_label,
)


def pi_protected_observables(H, gamma_l, rho_0, N, threshold=1e-9, cluster_tol=1e-8):
    """Identify Pauli-string observables σ_α with ⟨σ_α(t)⟩ = 0 for all t.

    Under L = -i[H, ·] + Σ_l γ_l (Z_l ρ Z_l - ρ), the time-evolved expectation
    is a sum of exponentials with rates λ_k (eigenvalues of L_pauli):

        ⟨σ_α(t)⟩ = 2^N · Σ_λ S_λ(α) · exp(λ t)

    where S_λ(α) = Σ_{k: λ_k=λ} V[α, k] · c[k] sums right-eigenvector
    components within each degenerate cluster. σ_α is Π-protected iff
    S_λ(α) = 0 for every cluster.

    This is strictly weaker than "each V[α,k]·c[k] vanishes" — degenerate-
    cluster cancellations are real (e.g., ⟨X_0IZ_2⟩ = 0 for Heisenberg
    chain on |+−+⟩ via SU(2)-multiplet cancellation).

    Returns dict with:
      'protected': list of {'k', 'pauli', 'max_cluster_contribution'}
      'active':    list of same plus 'dominant_eigenvalue'
      'eigenvalues': L_pauli eigenvalues
      'n_clusters': number of distinct eigenvalue clusters

    Identity (α=0) is excluded; ⟨I⟩ = 1 trivially.
    """
    L = lindbladian_z_dephasing(H, gamma_l)
    M_basis = _vec_to_pauli_basis_transform(N)
    L_pauli = (M_basis.conj().T @ L @ M_basis) / (2 ** N)

    evals, V = np.linalg.eig(L_pauli)
    Vinv = np.linalg.inv(V)

    rho_pauli = pauli_basis_vector(rho_0, N)
    c = Vinv @ rho_pauli

    # Cluster eigenvalues
    n = len(evals)
    used = np.zeros(n, dtype=bool)
    clusters = []
    for i in range(n):
        if used[i]:
            continue
        cluster = [i]
        used[i] = True
        for j in range(i + 1, n):
            if not used[j] and abs(evals[j] - evals[i]) < cluster_tol:
                cluster.append(j)
                used[j] = True
        clusters.append(cluster)

    d2 = 4 ** N
    protected, active = [], []
    for alpha in range(1, d2):  # skip identity
        cluster_sums = []
        for cluster in clusters:
            S = sum(V[alpha, k] * c[k] for k in cluster)
            cluster_sums.append((S, cluster[0]))
        max_S = max((abs(S) for S, _ in cluster_sums), default=0.0)
        entry = {
            'k': alpha,
            'pauli': _pauli_label(alpha, N),
            'max_cluster_contribution': float(max_S),
        }
        if max_S < threshold:
            protected.append(entry)
        else:
            dom_S, dom_idx = max(cluster_sums, key=lambda x: abs(x[0]))
            entry['dominant_eigenvalue'] = complex(evals[dom_idx])
            active.append(entry)

    return {
        'protected': protected,
        'active': active,
        'eigenvalues': evals,
        'n_clusters': len(clusters),
    }
