"""Investigate the structural origin of M's SVD for 2-body bilinears.

Three structural questions to test:

H1 (per-bond additivity): For non-overlapping bonds (star, disjoint pairs),
   does M = Σ_bonds M_bond ⊗ I_others (additive like the single-body case)?

H2 (max-uniform XX+XY): For star topology, H = X_hub ⊗ Σ_leaves (X+Y),
   factorizing as (single hub op) ⊗ (sum on leaves). What does M look like?
   At chain we observed: N=3 all SVs = 2√2; N=4 two clusters at 2√5 and 2.
   Predict the general N pattern.

H3 (YZ ≡ XY+YX SVD-equivalence): What identifies these two as same M?
   Common Π-orbit? Same bilinear-class? Test the equivalence at N=5.

Run from repo root.
"""
from __future__ import annotations

import sys
sys.path.insert(0, 'simulations')
import numpy as np

import framework as fw
from framework.pauli import _build_bilinear, site_op
from framework.lindblad import lindbladian_z_dephasing, palindrome_residual


def chain_bonds(N): return [(i, i+1) for i in range(N-1)]
def star_bonds(N): return [(0, i) for i in range(1, N)]
def disjoint_bonds(N):
    """Non-overlapping bond pairs: (0,1), (2,3), (4,5), ..."""
    return [(2*i, 2*i+1) for i in range(N // 2)]


def build_M_from_terms(N, terms, bonds):
    bilinear = [(t[0], t[1], 1.0) for t in terms]
    H = _build_bilinear(N, bonds, bilinear)
    L = lindbladian_z_dephasing(H, [1.0]*N)
    return palindrome_residual(L, N*1.0, N)


def cluster_svs(M, tol=1e-6):
    svs = np.linalg.svd(M, compute_uv=False)
    out = []
    for s in svs:
        placed = False
        for i, (v, c) in enumerate(out):
            if abs(s - v) < tol:
                out[i] = (v, c+1); placed = True; break
        if not placed:
            out.append((s, 1))
    return sorted(out, key=lambda x: -x[0])


# ----------------------------------------------------------------------
# H1: Per-bond additivity test (star + disjoint bonds)
# ----------------------------------------------------------------------

def test_H1_per_bond_additivity():
    print("=" * 78)
    print("H1: Per-bond additivity (star, disjoint)")
    print("=" * 78)
    print()
    print("If M = Σ_bonds M_bond ⊗ I_others, the SVDs of M should match")
    print("the additive prediction (sums of single-bond eigenvalues).")
    print()

    # Even N for disjoint, all N for star
    for terms_label, terms in [("XX+XY (hard)", [('X','X'), ('X','Y')]),
                               ("YZ (soft)", [('Y','Z')]),
                               ("YZ+ZY (soft)", [('Y','Z'), ('Z','Y')])]:
        print(f"  Terms: {terms_label}")
        for N, topology, bonds_fn in [(4, "disjoint(2 bonds)", disjoint_bonds),
                                      (4, "star(3 bonds)", star_bonds),
                                      (5, "star(4 bonds)", star_bonds)]:
            bonds = bonds_fn(N)

            # Direct M
            M_full = build_M_from_terms(N, terms, bonds)
            direct = cluster_svs(M_full)

            # Per-bond additive prediction:
            # Compute M for each bond alone (treating as single-bond H), then
            # verify if eigenvalues sum like single-body case.
            per_bond_evs = []
            for b in bonds:
                M_b = build_M_from_terms(N, terms, [b])
                evs = np.linalg.eigvals(M_b)
                per_bond_evs.append(evs)

            # If additive: eigenvalues of M = sum of one ev from each bond's M
            # But each M_b is on the full 4^N space, so sum of eigenvalues works
            # only if they share eigenbasis — usually not the case.
            # Instead, check if ‖M_full‖² = Σ_b ‖M_b‖² (Frobenius additivity,
            # which holds iff bond-supports are disjoint AND M_b's are F-orthogonal)
            sum_sq = sum(float(np.linalg.norm(M_b)**2)
                         for b in bonds
                         for M_b in [build_M_from_terms(N, terms, [b])])
            actual_sq = float(np.linalg.norm(M_full)**2)
            ratio = actual_sq / sum_sq if sum_sq > 0 else float('nan')
            print(f"    N={N} {topology}: ‖M_full‖²={actual_sq:.0f}, "
                  f"Σ ‖M_b‖²={sum_sq:.0f}, ratio={ratio:.4f}")
        print()


# ----------------------------------------------------------------------
# H2: Max-uniform XX+XY pattern across N
# ----------------------------------------------------------------------

def test_H2_max_uniform_xxxy():
    print("=" * 78)
    print("H2: XX+XY max-uniform pattern (chain)")
    print("=" * 78)
    print()

    for N in [2, 3, 4, 5]:
        bonds = chain_bonds(N) if N >= 2 else []
        if not bonds: continue
        M = build_M_from_terms(N, [('X','X'), ('X','Y')], bonds)
        clusters = cluster_svs(M)
        norm_sq = float(np.linalg.norm(M)**2)
        print(f"  N={N} chain ({len(bonds)} bonds): ‖M‖²={norm_sq:.0f}, "
              f"clusters={[(round(v,4), m) for v, m in clusters]}")

    print()
    print("  Pattern hypothesis:")
    print("  N=2 chain: 1 bond, B=1")
    print("  N=3 chain: 2 bonds, B=2: all SVs 2√2 ≈ 2.828 (uniform!)")
    print("  N=4 chain: 3 bonds, B=3: SVs at 2√5 and 2")
    print("  N=5 chain: 4 bonds, B=4: ?")


# ----------------------------------------------------------------------
# H3: YZ ≡ XY+YX equivalence
# ----------------------------------------------------------------------

def test_H3_yz_xyyx_equivalence():
    print("=" * 78)
    print("H3: YZ vs XY+YX SVD equivalence (chain)")
    print("=" * 78)
    print()

    print("  At N=3, 4 we observed identical SVD clusters. Test N=5:")
    for N in [3, 4, 5]:
        bonds = chain_bonds(N)
        M_yz = build_M_from_terms(N, [('Y','Z')], bonds)
        M_xy = build_M_from_terms(N, [('X','Y'), ('Y','X')], bonds)
        cl_yz = cluster_svs(M_yz)
        cl_xy = cluster_svs(M_xy)
        match = (len(cl_yz) == len(cl_xy)
                 and all(abs(a[0]-b[0]) < 1e-5 and a[1] == b[1]
                         for a, b in zip(cl_yz, cl_xy)))
        norm_yz = float(np.linalg.norm(M_yz)**2)
        norm_xy = float(np.linalg.norm(M_xy)**2)
        print(f"    N={N}: YZ ‖M‖²={norm_yz:.0f}, XY+YX ‖M‖²={norm_xy:.0f}, "
              f"clusters_match={'YES' if match else 'NO'}")
        if N == 5:
            print(f"      YZ clusters:    {[(round(v,3), m) for v, m in cl_yz[:6]]}")
            print(f"      XY+YX clusters: {[(round(v,3), m) for v, m in cl_xy[:6]]}")

    print()
    print("  Test more pairs of related Hamiltonians:")
    print("  (a) YZ vs ZY single-term (a chain-mirror image)")
    print("  (b) XY single-term vs YX single-term")
    print()
    pairs = [
        ("YZ", [('Y','Z')], "ZY", [('Z','Y')]),
        ("XY", [('X','Y')], "YX", [('Y','X')]),
        ("XX+XY", [('X','X'), ('X','Y')], "XX+YX", [('X','X'), ('Y','X')]),
    ]
    for label_a, terms_a, label_b, terms_b in pairs:
        for N in [3, 4]:
            bonds = chain_bonds(N)
            M_a = build_M_from_terms(N, terms_a, bonds)
            M_b = build_M_from_terms(N, terms_b, bonds)
            cl_a = cluster_svs(M_a)
            cl_b = cluster_svs(M_b)
            match = (len(cl_a) == len(cl_b)
                     and all(abs(x[0]-y[0]) < 1e-5 and x[1] == y[1]
                             for x, y in zip(cl_a, cl_b)))
            print(f"    N={N}: {label_a} vs {label_b}: clusters_match={'YES' if match else 'NO'}")


def main():
    test_H1_per_bond_additivity()
    test_H2_max_uniform_xxxy()
    test_H3_yz_xyyx_equivalence()


if __name__ == '__main__':
    main()
