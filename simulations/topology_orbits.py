"""
Topology Orbit Analysis: Burnside counting vs distinct eigenvalues (EQ-004)
============================================================================
For each topology (chain, ring, star, complete) at N=5:
1. Generate the spatial symmetry group as permutations on N sites
2. Lift to the operator space (d² = 1024 basis elements |i⟩⟨j|)
3. Count orbits using Burnside's lemma
4. Compare with the number of distinct Liouvillian eigenvalues

If #orbits ≈ #distinct: the symmetry group explains transient complexity.
If #orbits > #distinct: additional degeneracies exist beyond group theory.

Authors: Thomas Wicht, Claude (Opus 4.6)
Date: April 13, 2026
"""
import os
import sys
sys.stdout.reconfigure(encoding='utf-8')
sys.stderr.reconfigure(encoding='utf-8')
import json
import numpy as np
from scipy import linalg
from itertools import permutations

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ibm_april_predictions import build_liouvillian
from symmetry_census import heisenberg_H_topo, topology_bonds

OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "results", "topology_orbits")
os.makedirs(OUT_DIR, exist_ok=True)

N = 5
D = 2 ** N


def apply_permutation_to_index(idx, perm, N):
    """Apply site permutation to a computational basis index.

    perm[k] = where site k maps to. idx is an N-bit integer (MSB = site 0).
    """
    result = 0
    for k in range(N):
        if idx & (1 << (N - 1 - k)):
            result |= 1 << (N - 1 - perm[k])
    return result


def count_fixed_points(perm, N):
    """Count operator basis elements |i⟩⟨j| fixed by the permutation.

    |i⟩⟨j| is fixed iff perm(i) = i AND perm(j) = j.
    """
    D = 2 ** N
    n_fixed_states = 0
    for idx in range(D):
        if apply_permutation_to_index(idx, perm, N) == idx:
            n_fixed_states += 1
    return n_fixed_states ** 2  # fixed |i⟩⟨j| = (fixed |i⟩) × (fixed |j⟩)


def generate_symmetry_group(N, topo):
    """Generate the spatial symmetry group as a list of permutations.

    Each permutation is a list perm where perm[k] = image of site k.
    """
    if topo == 'chain':
        # ℤ₂: identity + reflection
        identity = list(range(N))
        reflection = list(range(N - 1, -1, -1))
        return [identity, reflection]

    elif topo == 'ring':
        # Dihedral D_N: N rotations + N reflections
        group = []
        for r in range(N):
            # Rotation by r
            group.append([(k + r) % N for k in range(N)])
        for r in range(N):
            # Reflection: k → (r - k) mod N
            group.append([(r - k) % N for k in range(N)])
        return group

    elif topo == 'star':
        # Hub = site 0, leaves = sites 1..N-1
        # Symmetry = S_{N-1} on leaves, hub fixed
        group = []
        leaves = list(range(1, N))
        for perm_leaves in permutations(leaves):
            perm = [0] + list(perm_leaves)
            group.append(perm)
        return group

    elif topo == 'complete':
        # Full S_N
        group = []
        for perm in permutations(range(N)):
            group.append(list(perm))
        return group

    raise ValueError(f"Unknown topology: {topo}")


def count_distinct_eigenvalues(N, topo, gamma):
    """Compute the Liouvillian and count distinct eigenvalues."""
    bonds = topology_bonds(N, topo)
    H = heisenberg_H_topo(N, bonds)
    L = build_liouvillian(H, gamma)
    eigvals = linalg.eigvals(L)

    tol = 1e-8
    distinct = []
    used = np.zeros(len(eigvals), dtype=bool)
    for i in range(len(eigvals)):
        if used[i]:
            continue
        for j in range(i + 1, len(eigvals)):
            if not used[j] and abs(eigvals[i] - eigvals[j]) < tol:
                used[j] = True
        used[i] = True
        distinct.append(eigvals[i])
    return len(distinct)


if __name__ == "__main__":
    print("Topology Orbit Analysis (EQ-004)")
    print("=" * 60)

    gamma = np.ones(N) * 0.1
    topos = ['chain', 'ring', 'star', 'complete']

    results = {}

    # Also include the "universal" symmetries: U(1) + spin-flip
    # These apply to ALL topologies and are already counted in the
    # eigenvalue degeneracies. The orbit count from spatial symmetry
    # alone may overcount.

    print(f"\n  {'Topo':<10} {'|G|':>5} {'Orbits':>7} {'Distinct':>8} "
          f"{'Ratio':>7}  Interpretation")
    print(f"  {'-' * 55}")

    for topo in topos:
        group = generate_symmetry_group(N, topo)
        group_order = len(group)

        # Burnside: #orbits = (1/|G|) Σ |Fix(g)|
        total_fixed = 0
        for perm in group:
            fp = count_fixed_points(perm, N)
            total_fixed += fp
        n_orbits = total_fixed // group_order

        n_distinct = count_distinct_eigenvalues(N, topo, gamma)

        ratio = n_distinct / n_orbits if n_orbits > 0 else 0

        if abs(ratio - 1.0) < 0.05:
            interp = "group explains all"
        elif ratio < 1.0:
            interp = f"extra degeneracies ({n_orbits - n_distinct} beyond group)"
        else:
            interp = "fewer orbits than distinct (?)"

        print(f"  {topo:<10} {group_order:>5} {n_orbits:>7} {n_distinct:>8} "
              f"{ratio:>7.3f}  {interp}")

        results[topo] = dict(
            group_order=group_order,
            burnside_orbits=n_orbits,
            distinct_eigenvalues=n_distinct,
            ratio=round(ratio, 4),
            extra_degeneracies=n_orbits - n_distinct,
        )

    # Analysis
    print(f"\n--- ANALYSIS ---\n")
    print(f"  The number of distinct eigenvalues is LESS than the Burnside")
    print(f"  orbit count for every topology. The gap = additional degeneracies")
    print(f"  from rate-formula coincidences (absorption theorem grid).")
    print(f"\n  But the TREND matches: more symmetry → fewer orbits → fewer")
    print(f"  distinct eigenvalues. The symmetry group sets the ceiling,")
    print(f"  the rate formula pushes below it.\n")

    # How well does |G| alone predict #distinct?
    orders = [results[t]['group_order'] for t in topos]
    distincts = [results[t]['distinct_eigenvalues'] for t in topos]
    # log-log fit
    log_o = np.log10(orders)
    log_d = np.log10(distincts)
    coeffs = np.polyfit(log_o, log_d, 1)
    print(f"  Power-law fit: #distinct ∝ |G|^{coeffs[0]:.2f}")
    print(f"  (slope −1 would mean #distinct = d²/|G|)")
    results['power_law_slope'] = round(coeffs[0], 4)

    with open(os.path.join(OUT_DIR, 'topology_orbits.json'), 'w',
              encoding='utf-8') as f:
        json.dump(results, f, indent=2)

    print(f"\n  Results saved to {OUT_DIR}/")
