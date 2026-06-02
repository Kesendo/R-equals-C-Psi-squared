#!/usr/bin/env python3
"""Star spectral compactness (Reading 2): does S_{N-1} representation theory
account for the star's distinct-eigenvalue count?

STAR_SPECTRUM_COMPACTNESS observes that the Heisenberg star (hub + N-1 leaves)
under uniform Z-dephasing has a strikingly compact Liouvillian spectrum: at N=8
only 2275 of 65536 eigenvalues are distinct (3.5%), ~30x more degenerate than
the chain. The hub-and-leaves geometry has automorphism group S_{N-1} (permuting
the N-1 leaves); Reading 2 asks whether the 2275 follows from S_{N-1} irrep
multiplicities.

Mechanism under test. L commutes with the leaf-permutation S_{N-1} (H_star and
the uniform dephasing are both leaf-symmetric), so L block-diagonalises into
S_{N-1}-isotypic components. The Liouville (operator) space is
    C^4 (hub, S_{N-1}-trivial)  ⊗  (C^4)^{⊗(N-1)} (leaves, S_{N-1} permutes factors).
By Schur-Weyl for the symmetric-group permutation action, the multiplicity of the
S_{N-1}-irrep λ ⊢ (N-1) in (C^4)^{⊗(N-1)} is dim S^λ(C^4) (the number of SSYT of
shape λ with entries in {1,2,3,4}); the hub multiplies it by 4. So L restricts to
an M_λ × M_λ block on the multiplicity space, M_λ = 4·dim S^λ(C^4), and each of its
eigenvalues appears with degeneracy d_λ = dim(S_{N-1}-irrep λ). Hence

    distinct eigenvalues of L  ≤  Σ_{λ ⊢ (N-1)} M_λ   (equality iff every block
    is non-degenerate and no two blocks share a value).

This probe computes the prediction Σ_λ M_λ (exact combinatorics, any N) and the
actual star distinct count (dense eig, N=3..6), and compares, to see whether
S_{N-1} alone is the dominant accounting or whether the common U(1)/palindrome
structure reduces it further.

Finding (2026-06-02). S_{N-1} gives a clean upper bound, distinct ≤ Σ M_λ (= 4400
at N=8), but the exact 2275 lies below it (≈0.52×). The gap is cross-block value
coincidence from the non-normal dephasing spectrum, NOT a representation-theoretic
effect: the U(1) (Sz_left, Sz_right) refinement does not lower Σ M_λ, and the
distinct = |Re|·|Im| factorisation is refuted (Im(λ_L) ≠ H-gaps since L is
non-normal). So S_{N-1} is necessary and dominant but not sufficient for the exact
count; the derivable object is the Σ M_λ bound. See
hypotheses/STAR_SPECTRUM_COMPACTNESS.md (Reading 2).
"""
from __future__ import annotations

import sys
from math import factorial

import numpy as np

sys.path.insert(0, 'simulations')
import framework as fw  # noqa: E402

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass


# ---------------------------------------------------------------------------
# Partition / Young-diagram combinatorics
# ---------------------------------------------------------------------------
def partitions(n, max_part=None):
    """All integer partitions of n (as non-increasing tuples)."""
    if max_part is None:
        max_part = n
    if n == 0:
        yield ()
        return
    for first in range(min(n, max_part), 0, -1):
        for rest in partitions(n - first, first):
            yield (first,) + rest


def hook_lengths(lam):
    """Hook length of every cell of the Young diagram λ. Returns list of ints."""
    cols = [0] * (lam[0] if lam else 0)
    for r in lam:
        for j in range(r):
            cols[j] += 1
    hooks = []
    for i, row in enumerate(lam):
        for j in range(row):
            arm = row - (j + 1)            # cells to the right in the row
            leg = cols[j] - (i + 1)        # cells below in the column
            hooks.append(arm + leg + 1)
    return hooks


def dim_irrep(lam):
    """dim of the S_n irrep V_λ (hook-length formula)."""
    n = sum(lam)
    prod = 1
    for h in hook_lengths(lam):
        prod *= h
    return factorial(n) // prod


def dim_schur(lam, d):
    """dim S^λ(C^d) = # SSYT of shape λ with entries in [d] (hook-content formula).
    Zero if λ has more than d rows."""
    if len(lam) > d:
        return 0
    # content (j - i) with 1-indexed (i,j); arrange as exact rational via product.
    num = 1
    den = 1
    hooks = hook_lengths(lam)
    k = 0
    for i, row in enumerate(lam):       # i = row index (0-based)
        for j in range(row):            # j = col index (0-based)
            num *= (d + (j - i))        # content = (j+1) - (i+1) = j - i
            den *= hooks[k]
            k += 1
    assert num % den == 0, f"non-integer dim_schur for {lam}"
    return num // den


def repn_prediction(N, d=4):
    """The S_{N-1} accounting: blocks M_λ = d_hub·dim S^λ(C^d), degeneracy d_λ.
    Returns (sum_M, rows) where rows = [(λ, d_λ, dim_schur, M_λ)]."""
    m = N - 1
    rows = []
    total_dim = 0
    sum_M = 0
    for lam in partitions(m):
        ds = dim_schur(lam, d)
        if ds == 0:
            continue
        dlam = dim_irrep(lam)
        M = d * ds                      # hub factor d=4
        rows.append((lam, dlam, ds, M))
        sum_M += M
        total_dim += M * dlam
    return sum_M, total_dim, rows


# ---------------------------------------------------------------------------
# Actual star spectrum
# ---------------------------------------------------------------------------
def star_distinct(N, gamma_0=0.5, J=1.0, tol=1e-6):
    """Distinct Liouvillian eigenvalues of the Heisenberg star (cluster within tol).
    Returns (n_total, n_distinct, degeneracy_multiset)."""
    chain = fw.ChainSystem(N, gamma_0=gamma_0, J=J, topology='star', H_type='heisenberg')
    ev = np.linalg.eigvals(chain.L)
    keys = np.round(ev / tol).astype(np.complex128)   # bin to a tol grid
    uniq, counts = np.unique(keys, return_counts=True)
    from collections import Counter
    degmult = Counter(int(c) for c in counts)
    return len(ev), len(uniq), degmult


def analyze_structure(N, gamma_0=0.5, J=1.0, tol=1e-6):
    """Decompose the star spectrum into Re (decay ladder) x Im (H-gaps).

    Im(λ_L) = ω_α − ω_β are differences of H_star eigenvalues; Re(λ_L) is the
    dephasing decay. Tests whether distinct(L) factorises as |distinct Re| x
    |distinct Im|, and whether |distinct Im(L)| equals |distinct H-gaps|."""
    chain = fw.ChainSystem(N, gamma_0=gamma_0, J=J, topology='star', H_type='heisenberg')
    ev = np.linalg.eigvals(chain.L)

    def nbin(arr):
        return len(np.unique(np.round(np.asarray(arr) / tol).astype(np.int64)))

    n_re = nbin(ev.real)
    n_im = nbin(ev.imag)
    n_dist = len(np.unique(np.round(ev / tol).astype(np.complex128)))

    # H_star spectrum: distinct energies (SU(2) multiplets) and distinct gaps.
    hev = np.linalg.eigvalsh(chain.H)              # H is Hermitian
    n_h = nbin(hev)
    gaps = (hev[:, None] - hev[None, :]).ravel()
    n_gap = nbin(gaps)

    return {
        'n_dist': n_dist, 'n_re': n_re, 'n_im': n_im,
        'prod': n_re * n_im, 'n_h': n_h, 'n_gap': n_gap,
    }


def main():
    print("Star spectral compactness (Reading 2): S_{N-1} representation theory")
    print("  Heisenberg star (hub + N-1 leaves), J=1, γ=0.5, uniform Z-dephasing")
    print("  Prediction: distinct ≤ Σ_λ M_λ,  M_λ = 4·dim S^λ(C^4),  λ ⊢ (N-1)")
    print()

    print(f"  {'N':>2s} {'4^N':>7s} {'Σ M_λ (pred)':>13s} {'actual distinct':>16s} "
          f"{'ratio':>7s} {'#irreps':>8s}")
    print(f"  {'-'*2} {'-'*7} {'-'*13} {'-'*16} {'-'*7} {'-'*8}")

    for N in (3, 4, 5, 6):
        sum_M, total_dim, rows = repn_prediction(N)
        assert total_dim == 4 ** N, f"dim check failed N={N}: {total_dim} != {4**N}"
        n_tot, n_dist, degmult = star_distinct(N)
        assert n_dist <= sum_M, f"upper bound violated at N={N}: {n_dist} > {sum_M}"
        ratio = n_dist / sum_M if sum_M else float('nan')
        print(f"  {N:>2d} {4**N:>7d} {sum_M:>13d} {n_dist:>16d} {ratio:>7.3f} {len(rows):>8d}")

    # N=8 prediction (the headline anchor; actual = 2275 from the doc)
    print()
    sum_M8, total_dim8, rows8 = repn_prediction(8)
    assert total_dim8 == 4 ** 8
    print(f"  N=8 prediction Σ M_λ = {sum_M8}  (λ ⊢ 7, {len(rows8)} irreps); "
          f"doc actual distinct = 2275")
    print(f"  N=8 ratio actual/pred = {2275 / sum_M8:.3f}")
    print()

    # Re x Im factorisation test: is distinct(L) = |distinct Re| x |distinct Im|?
    print("  Re (decay ladder) x Im (H-gaps) factorisation:")
    print(f"  {'N':>2s} {'distinct':>8s} {'|Re|':>5s} {'|Im|':>5s} {'|Re|·|Im|':>9s} "
          f"{'dist/prod':>9s} {'|H eig|':>7s} {'|H gaps|':>8s} {'Im==gaps?':>9s}")
    for N in (3, 4, 5, 6):
        s = analyze_structure(N)
        eq = 'yes' if s['n_im'] == s['n_gap'] else 'NO'
        print(f"  {N:>2d} {s['n_dist']:>8d} {s['n_re']:>5d} {s['n_im']:>5d} "
              f"{s['prod']:>9d} {s['n_dist']/s['prod']:>9.3f} {s['n_h']:>7d} "
              f"{s['n_gap']:>8d} {eq:>9s}")
    print()

    # Detail for N=4: the per-irrep block table + the actual degeneracy multiset.
    print("  --- N=4 detail (S_3 on 3 leaves) ---")
    _, _, rows4 = repn_prediction(4)
    print(f"    {'λ':>10s} {'d_λ':>4s} {'dim S^λ(C^4)':>13s} {'M_λ=4·dim':>10s}")
    for lam, dlam, ds, M in rows4:
        print(f"    {str(lam):>10s} {dlam:>4d} {ds:>13d} {M:>10d}")
    n_tot, n_dist, degmult = star_distinct(4)
    print(f"    actual: {n_tot} total, {n_dist} distinct; degeneracy multiset "
          f"(deg: count) = {dict(sorted(degmult.items()))}")
    print(f"    S_3 irrep dims present: {sorted({dlam for _, dlam, _, _ in rows4})}")
    print()
    print("Done.")


if __name__ == "__main__":
    main()
