"""Test: Is the chain (X,Y) ↔ chain (X,Z) M-spectrum universality the operator-level
shadow of a PTF-style perspectival-time-rescaling between their dynamics?

Hypothesis: chain (X,Y) and chain (X,Z) generate different L's, but their
state-trajectories are related by a per-site perspectival-time-rescaling
{α_i(t)}: site purity P_i^(XZ)(t) = P_i^(XY)(α_i · t) for some closure-respecting
α_i with Σ ln(α_i) = 0.

If true: this would be the t-parametrized unitary equivalence we couldn't find
as a static state-space rotation. The chain-universality of M would then be
a corollary of PTF's painter-equivalence between the two L's.

Run from repo root.
"""
from __future__ import annotations

import sys
sys.path.insert(0, 'simulations')
import numpy as np
from scipy.linalg import expm

import framework as fw
from framework.pauli import _build_bilinear
from framework.lindblad import lindbladian_z_dephasing


def chain_bonds(N): return [(i, i+1) for i in range(N-1)]


def build_L(N, P, Q, gamma=1.0):
    H = _build_bilinear(N, chain_bonds(N), [(P, Q, 1.0)])
    return lindbladian_z_dephasing(H, [gamma]*N)


def propagate_eigen(L, rho0_vec, t_array):
    """ρ(t) = exp(L·t)·ρ0, computed via eigendecomposition of L (works since γ>0)."""
    eigvals, eigvecs = np.linalg.eig(L)
    V_inv = np.linalg.inv(eigvecs)
    proj = V_inv @ rho0_vec
    return np.array([eigvecs @ (np.exp(eigvals * t) * proj) for t in t_array])


def vec_to_matrix(rho_vec, d):
    return rho_vec.reshape((d, d), order='F')


def site_purity(rho_full, N, site):
    """P_i(t) = Tr(ρ_i²) where ρ_i is partial trace over all sites except `site`."""
    d = 2 ** N
    rho_t = rho_full.reshape([2] * (2 * N))  # ket-then-bra indices
    # Move ket-site and bra-site axes to front
    rho_t = np.moveaxis(rho_t, [site, N + site], [0, 1])
    # Trace over remaining 2(N-1) axes pair-wise
    rho_t = rho_t.reshape(2, 2, 2 ** (N - 1), 2 ** (N - 1))
    # Partial trace: contract the second pair (env)
    rho_site = np.einsum('ijkk->ij', rho_t)
    return float(np.real(np.trace(rho_site @ rho_site)))


def main():
    N = 4
    gamma = 1.0
    L_XY = build_L(N, 'X', 'Y', gamma=gamma)
    L_XZ = build_L(N, 'X', 'Z', gamma=gamma)

    # Initial state: |+⟩^N (uniform superposition, max-coherent product state)
    plus = np.array([1, 1], dtype=complex) / np.sqrt(2)
    psi_0 = plus
    for _ in range(N - 1):
        psi_0 = np.kron(psi_0, plus)
    rho_0 = np.outer(psi_0, psi_0.conj())
    rho_0_vec = rho_0.reshape(-1, order='F')

    # Time grid
    t_array = np.linspace(0.0, 3.0, 121)

    # Propagate
    rho_XY_t = propagate_eigen(L_XY, rho_0_vec, t_array)
    rho_XZ_t = propagate_eigen(L_XZ, rho_0_vec, t_array)

    # Compute site purities
    d = 2 ** N
    P_XY = np.array([[site_purity(rho_XY_t[ti], N, s) for s in range(N)]
                     for ti in range(len(t_array))])
    P_XZ = np.array([[site_purity(rho_XZ_t[ti], N, s) for s in range(N)]
                     for ti in range(len(t_array))])

    print("=" * 78)
    print("Site purity trajectories: chain (X,Y) vs chain (X,Z)")
    print("=" * 78)
    print()
    print(f"  Initial state: |+⟩^{N}, all sites pure (P_i(0) = 1)")
    print(f"  Time grid: t ∈ [0, 3], 121 points")
    print()
    print("  Are they directly identical (P_i^XY(t) = P_i^XZ(t) for all i, t)?")
    direct_max_diff = float(np.max(np.abs(P_XY - P_XZ)))
    print(f"    Max |P_i^XY(t) - P_i^XZ(t)| over all (i, t) = {direct_max_diff:.4e}")
    print(f"    Trivially identical? {'YES' if direct_max_diff < 1e-10 else 'NO'}")
    print()

    # If not identical: try per-site time rescaling
    # P_XZ(t) = P_XY(α_i · t) for some α_i ?
    print("  Try per-site time-rescaling: find α_i so P_XZ(t) ≈ P_XY(α_i · t)")
    from scipy.interpolate import interp1d
    from scipy.optimize import minimize_scalar

    alphas = []
    fit_errors = []
    for site in range(N):
        # Interpolant for P_XY trajectory at site `site`
        P_XY_interp = interp1d(t_array, P_XY[:, site], kind='cubic',
                                bounds_error=False, fill_value=P_XY[-1, site])

        def cost(alpha):
            # Compare P_XY(α·t) to P_XZ(t)
            P_XY_rescaled = P_XY_interp(alpha * t_array)
            return float(np.sum((P_XY_rescaled - P_XZ[:, site]) ** 2))

        result = minimize_scalar(cost, bounds=(0.1, 10.0), method='bounded')
        alphas.append(result.x)
        fit_errors.append(np.sqrt(result.fun / len(t_array)))

    print(f"    Fitted α_i:        {[round(a, 4) for a in alphas]}")
    print(f"    RMS errors:        {[round(e, 4) for e in fit_errors]}")
    print(f"    Σ ln(α_i):         {sum(np.log(a) for a in alphas):.4e}")
    print(f"    All RMS errors < 0.01? {'YES' if all(e < 0.01 for e in fit_errors) else 'NO'}")
    print()

    # Also test: maybe trajectories coincide under some site-permutation
    print("  Try site-permutation: P_i^XY(t) = P_π(i)^XZ(t) for some permutation π?")
    from itertools import permutations
    best_perm = None
    best_diff = np.inf
    for perm in permutations(range(N)):
        permuted = P_XZ[:, list(perm)]
        diff = float(np.max(np.abs(P_XY - permuted)))
        if diff < best_diff:
            best_diff = diff
            best_perm = perm
    print(f"    Best permutation: {best_perm}")
    print(f"    Max |P_i^XY(t) - P_π(i)^XZ(t)| = {best_diff:.4e}")
    print(f"    Site-permutation explains it? {'YES' if best_diff < 1e-10 else 'NO'}")
    print()

    # Summary table
    print("  Per-site initial vs final purities:")
    print(f"    site | P^XY(0)  P^XY(t_max) | P^XZ(0)  P^XZ(t_max)")
    for s in range(N):
        print(f"    {s}    | {P_XY[0, s]:.4f}  {P_XY[-1, s]:.4f}     "
              f"| {P_XZ[0, s]:.4f}  {P_XZ[-1, s]:.4f}")


if __name__ == '__main__':
    main()
