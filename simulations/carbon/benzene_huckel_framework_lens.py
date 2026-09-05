"""Selected finite Hückel-graph calculations through the R=CΨ² lens.

The conventional molecule names in this script label stipulated C_N rings and
chains; they do not identify a material degree of freedom or molecular dynamics.
The calculation diagonalises a nearest-neighbour Hückel matrix. Coulson-Rushbrooke
and its K sublattice gauge concern that single-particle bipartite graph. F1 and Π
concern a Liouvillian under its stated channel. They are formal structural siblings
with different triggers and operator scopes, not one physical identity.

Conditional model vocabulary, after its coordinates and channel are selected:

  Coulson-Rushbrooke/K pairing  ↔  F1/Π pairing (different triggers)
  α (selected Hückel centre)    ↔  −Σγ analogue (Liouvillian centre)
  β (selected hopping unit)     ↔  J (selected XX+YY coupling unit)

This script builds stipulated C_N chain/ring matrices, checks their finite
pair-sum and cosine formulae, and prints selected filling-line positions. It does
not calculate aromaticity, antiaromaticity, stability, Jahn-Teller behavior, a
material bath, or a material β-to-J/gamma/Q mapping.

Run:
  PYTHONIOENCODING=utf-8 python simulations/carbon/benzene_huckel_framework_lens.py
"""
from __future__ import annotations

import sys
import numpy as np
from fractions import Fraction

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass

# Conventional Hückel labels α and β are used only as selected matrix parameters.
# We compute in β units (α = 0, β = 1), so the finite graph pair-sum centre is α = 0.

ALPHA = 0.0
BETA = 1.0


def huckel_chain(N: int) -> np.ndarray:
    """Open nearest-neighbour Hückel matrix on a selected N-site graph."""
    H = np.full((N, N), 0.0)
    np.fill_diagonal(H, ALPHA)
    for i in range(N - 1):
        H[i, i + 1] = BETA
        H[i + 1, i] = BETA
    return H


def huckel_ring(N: int) -> np.ndarray:
    """Cyclic nearest-neighbour Hückel matrix on a selected N-site graph."""
    H = huckel_chain(N)
    H[0, N - 1] = BETA
    H[N - 1, 0] = BETA
    return H


def is_bipartite_ring(N: int) -> bool:
    """A cyclic graph C_N is bipartite iff N is even."""
    return N % 2 == 0


def is_bipartite_chain(N: int) -> bool:
    """Open chains are always bipartite (alternant)."""
    return True


def coulson_rushbrooke_palindrome_check(eigvals: np.ndarray, tol: float = 1e-10) -> tuple[bool, list[tuple[float, float]]]:
    """For a selected alternant graph, pair energies as (α+x, α−x). Test:
    sort by energy and check end-to-end pair sums against 2α within ``tol``.
    Returns (is_palindromic, list of (E_pair_sum, deviation_from_2α))."""
    sorted_eigs = np.sort(eigvals)
    N = len(sorted_eigs)
    pairs = []
    is_pal = True
    for k in range(N // 2):
        low = sorted_eigs[k]
        high = sorted_eigs[N - 1 - k]
        pair_sum = low + high
        deviation = abs(pair_sum - 2 * ALPHA)
        pairs.append((pair_sum, deviation))
        if deviation > tol:
            is_pal = False
    # Middle eigenvalue if N odd
    if N % 2 == 1:
        mid = sorted_eigs[N // 2]
        deviation_mid = abs(mid - ALPHA)
        pairs.append((mid, deviation_mid))
        if deviation_mid > tol:
            is_pal = False
    return is_pal, pairs


def chain_mirror_classify(eigvecs: np.ndarray, N: int) -> list[str]:
    """For each MO column of eigvecs, classify as F71-mirror-even/odd by checking
    ψ(j) vs ψ(N-1-j). 'even' if ψ(j) = +ψ(N-1-j); 'odd' if = −ψ(N-1-j); 'mixed'
    if neither (degeneracy mixing — pick orthonormal basis that diagonalises mirror)."""
    classes = []
    for k in range(eigvecs.shape[1]):
        v = eigvecs[:, k]
        v_mirror = v[::-1]
        # Normalise sign convention by largest |entry|
        idx_max = np.argmax(np.abs(v))
        sign = np.sign(v[idx_max]) if abs(v[idx_max]) > 1e-12 else 1.0
        v = v * sign
        v_mirror = v_mirror * np.sign(v_mirror[N - 1 - idx_max]) if abs(v_mirror[N - 1 - idx_max]) > 1e-12 else v_mirror
        if np.allclose(v, v_mirror, atol=1e-8):
            classes.append("even")
        elif np.allclose(v, -v_mirror, atol=1e-8):
            classes.append("odd")
        else:
            classes.append("mixed")
    return classes


def frost_circle_predicted(N: int, topology: str = "ring") -> np.ndarray:
    """For a regular N-ring, Hückel eigenvalues = α + 2β·cos(2πk/N) for k = 0..N−1
    (Frost circle: vertices of inscribed N-gon on circle radius 2β at centre α).
    For an open chain, eigenvalues = α + 2β·cos(πk/(N+1)) for k = 1..N (the standard
    OBC sine-mode dispersion, identical to framework's XyJordanWignerModes)."""
    if topology == "ring":
        return np.array([ALPHA + 2 * BETA * np.cos(2 * np.pi * k / N) for k in range(N)])
    elif topology == "chain":
        return np.array([ALPHA + 2 * BETA * np.cos(np.pi * k / (N + 1)) for k in range(1, N + 1)])
    else:
        raise ValueError(topology)


def analyse(name: str, H: np.ndarray, topology: str, N: int, n_electrons: int):
    """One finite selected-Hückel-graph analysis pass through the framework lens."""
    eigvals, eigvecs = np.linalg.eigh(H)
    sorted_eigs = np.sort(eigvals)

    print(f"=" * 78)
    print(f"  {name} label: selected {topology} graph, N = {N} sites, {n_electrons} selected filling count")
    print(f"=" * 78)

    # 1. Frost-circle prediction match
    frost = frost_circle_predicted(N, topology)
    frost_sorted = np.sort(frost)
    max_dev = np.max(np.abs(sorted_eigs - frost_sorted))
    print(f"  Finite cosine formula:              {[f'{e:+.4f}β' for e in frost_sorted]}")
    print(f"  Hückel eigvalsh result:             {[f'{e:+.4f}β' for e in sorted_eigs]}")
    print(f"  Max deviation from Frost circle: {max_dev:.2e}  {'✓ match' if max_dev < 1e-10 else '✗ MISMATCH'}")
    print()

    # 2. Coulson-Rushbrooke palindrome around α (F1's sibling, bipartite-triggered)
    is_bipartite = is_bipartite_ring(N) if topology == "ring" else is_bipartite_chain(N)
    is_pal, pairs = coulson_rushbrooke_palindrome_check(eigvals)
    print(f"  C-R/K pair-sum check around α = {ALPHA}: " +
          f"{'✓ holds within tolerance' if is_pal else '✗ violated'}")
    print(f"  Bipartite (alternant)? {is_bipartite}  (palindrome predicted iff yes)")
    if not is_bipartite and is_pal:
        print(f"    [unexpected: palindrome holds on non-bipartite system]")
    if is_bipartite and not is_pal:
        print(f"    [unexpected: bipartite but palindrome violated — check]")
    for i, (pair_sum, deviation) in enumerate(pairs[:N // 2]):
        print(f"    pair {i}: E_lo + E_hi = {pair_sum:+.6f}β,  dev from 2α = {deviation:.2e}")
    if N % 2 == 1:
        print(f"    unpaired middle: E = {pairs[-1][0]:+.6f}β,  dev from α = {pairs[-1][1]:.2e}")
    print()

    # 3. F71 mirror classification of MOs
    classes = chain_mirror_classify(eigvecs, N)
    print(f"  F71 chain-mirror classification of MOs (j ↔ N−1−j):")
    for k, (E, c) in enumerate(zip(sorted_eigs, classes)):
        # Re-fetch the original ordering (eigvals not sorted)
        # Find which eigvec corresponds to this sorted energy
        orig_idx = np.argmin(np.abs(eigvals - E))
        original_class = classes[orig_idx]
        print(f"    MO {k}: E = {E:+.4f}β,  F71-class = {original_class}")
    print()

    # 4. Selected filling-line positions; no material orbital interpretation is assigned.
    n_doubly_filled = n_electrons // 2
    if n_electrons % 2 == 0:
        homo_E = sorted_eigs[n_doubly_filled - 1]
        lumo_E = sorted_eigs[n_doubly_filled]
        gap = lumo_E - homo_E
        print(f"  Selected filling-line structure ({n_electrons} count in {N} modes):")
        print(f"    HOMO (MO {n_doubly_filled - 1}): E = {homo_E:+.4f}β")
        print(f"    LUMO (MO {n_doubly_filled}): E = {lumo_E:+.4f}β")
        print(f"    Gap: {gap:.4f}β")
        # Selected occupied-mode sum = 2 · sum of the n_doubly_filled lowest eigenvalues.
        total_pi_energy = 2 * np.sum(sorted_eigs[:n_doubly_filled])
        # Reference: the same selected count at the diagonal value α.
        ref_energy = 2 * n_doubly_filled * ALPHA + n_electrons * 0  # n_electrons · 0 for β = 0
        delocalisation_energy = total_pi_energy - ref_energy
        print(f"    Occupied-mode sum: {total_pi_energy:+.4f}β  (= 2·Σ E_filled)")
        print(f"    Offset from diagonal reference: {delocalisation_energy:+.4f}β")
        # 4n+2 is retained as a selected filling-count label only.
        if topology == "ring":
            is_aromatic_huckel = (n_electrons - 2) % 4 == 0  # 4n+2 rule
            print(f"    Selected 4n+2 filling label: {n_electrons} count → " +
                  f"{'4n+2' if is_aromatic_huckel else '4n'} (no chemical inference)")
        print()

    # 5. Framework-lens translation
    print(f"  ─── Conditional framework-lens reading ───")
    print(f"    α (selected Hückel centre)    ↔ −Σγ analogue (Liouvillian centre)")
    print(f"    β (selected hopping unit)     ↔ J (selected XX+YY coupling unit)")
    print(f"    bipartite graph               ↔ K H K = −H, the sublattice gauge")
    print(f"    C-R/K pairing                 ↔ F1/Π sibling (different triggers; F1 topology-blind)")
    if topology == "ring" and N % 2 == 0:
        # KIntermediate Dicke anchor candidate inheritance
        m = N // 2 - 1
        print(f"    Selected even-N filling (N={N}) ↔ F86b KIntermediate candidate n ∈ {{{m}, {m+1}}}")
        print(f"      α_total(t=0) = 3/8 is F86b's stated selected-state result; material mapping unassigned")
    print()


def main():
    print()
    print("=" * 78)
    print(" Selected finite Hückel graphs through the R=CΨ² framework lens ")
    print("=" * 78)
    print()
    print(" Hückel matrix: diagonalise α·I + β·A on the stipulated finite graph.")
    print(" C-R/K: an alternant (bipartite) graph has pair energies around α.")
    print(" R=CΨ² F1 (2026): spec(L) palindromic around −Σγ under Z-deph, on any graph.")
    print()
    print(" C-R/K and F1/Π are formal sibling pairings with distinct operator scopes and")
    print(" triggers: C-R/K uses a bipartite single-particle graph; F1/Π uses its stated")
    print(" Liouvillian Hamiltonian/channel scope and is topology-blind. Their centres")
    print(" are compared only conditionally, not as a material identity.")
    print()
    print()

    # Selected conventional C6-ring label, with filling count 6.
    analyse("Benzene", huckel_ring(6), "ring", N=6, n_electrons=6)

    # Selected C10 perimeter ring; it is not a fused-ring molecular construction.
    analyse("Cyclodecapentaene (perimeter)", huckel_ring(10), "ring", N=10, n_electrons=10)

    # Selected conventional C4-chain label, with filling count 4.
    analyse("Butadiene", huckel_chain(4), "chain", N=4, n_electrons=4)

    # Selected conventional C6-chain label, with filling count 6.
    analyse("Hexatriene (open chain)", huckel_chain(6), "chain", N=6, n_electrons=6)

    # Selected C4-ring label, whose chosen filling line is degenerate at α.
    analyse("Cyclobutadiene", huckel_ring(4), "ring", N=4, n_electrons=4)

    # Selected non-bipartite C3-ring label, with filling count 2.
    analyse("Cyclopropenyl cation", huckel_ring(3), "ring", N=3, n_electrons=2)

    print("=" * 78)
    print(" Summary of framework-eye observations")
    print("=" * 78)
    print()
    print(" • The selected bipartite graph rows show C-R/K pair sums around α = 0 within")
    print("   the stated numerical tolerance. K carries that graph trigger; F1/Π has a")
    print("   different, topology-blind Liouvillian scope.")
    print()
    print(" • The selected C3 odd-ring row fails the C-R/K pair-sum check — the")
    print("   eigenvalues are {α + 2β, α − β, α − β}; pair sum (α + 2β) + (α − β) = 2α + β,")
    print("   not 2α. The bipartite trigger is gone, so K goes with it. NOT an F1 Brecher:")
    print("   F1 is topology-blind and its palindrome holds on this ring.")
    print()
    print(" • The selected 4n+2/4n filling labels sort the finite C6/C4 Hückel rows by")
    print("   their filling-line position, not by C-R/K pair symmetry. No aromaticity,")
    print("   antiaromaticity, stability, Jahn-Teller, or material mechanism follows.")
    print("   A Klein-4 reading is an open model proposal requiring a selected Liouvillian.")
    print()
    print(" • The selected even-N C6 graph has a candidate comparison to F86b's stated")
    print("   KIntermediate window n ∈ {2, 3}; F86b's γ = 1/2 row gives α_total = 3/8.")
    print("   This does not identify a benzene state. Any F98 comparison needs an explicitly")
    print("   selected open-system operator, channel, state, and producer.")


if __name__ == "__main__":
    main()
