using System;
using System.Collections.Generic;

namespace RCPsiSquared.Core.F89PathK;

/// <summary>The joint-popcount block lattice {(p,w)} of the N-site chain Liouvillian and the fold-lattice
/// D₄ machinery of docs/proofs/PROOF_CODIM1_BY_ADDITIVITY.md §7 (the fold-lattice lemma, gate
/// BlockLatticeFoldGroupTests): G = ⟨t, f_P, f_Q⟩ ≅ D₄ acts by transpose (p,w) ↦ (w,p), Klein
/// (p,w) ↦ (N−p,N−w) (same λ), and the two folds (p,w) ↦ (p,N−w), (N−p,w) (λ ↦ −λ−2N); membership of the
/// s-invariant pair {λ_A, μ = −λ_A−2N} is constant on orbits, so exclusion work lives on the quotient
/// region {p ≤ w, p + w ≤ N}. The rate window is the Bendixson interval Re spec ⊆ [−2·n_max, −2·n_min]
/// (the window-shell lemma, gate WindowShellLemmaTests). Block label convention: (p,w) are the two
/// popcount weights in the argument order of <see cref="WeightCoherenceBlock.Build(int,int,int,System.Numerics.Complex)"/>,
/// i.e. Build(n, p, w) IS the (p,w) block — the naming WeightCoherenceBlockTests uses for L(1,2). A
/// consistent swap of the two labels is unitary-silent for spectra, σ_min, windows, and R-parity
/// (transpose ∘ bipartite gauge), so the convention is cosmetic; it is pinned here once.</summary>
public static class BlockLattice
{
    /// <summary>All (N+1)² blocks (p,w) ∈ [0,N]².</summary>
    public static IEnumerable<(int P, int W)> AllBlocks(int n)
    {
        for (int p = 0; p <= n; p++)
            for (int w = 0; w <= n; w++)
                yield return (p, w);
    }

    /// <summary>The quotient region {p ≤ w, p + w ≤ N} of the D₄ block-lattice action. NOTE: an
    /// over-covering transversal, not a strict fundamental domain — some orbits intersect it more than
    /// once (N=5: 12 blocks for 6 orbits; e.g. (1,2) and its fold image (1,3) both lie inside). Safe for
    /// the census: over-covering yields redundant probes, never a missed member, and the two-shift probe
    /// {λ_A, μ} distinguishes the band representative from the fold representative.</summary>
    public static bool InFundamentalDomain(int n, int p, int w) => p <= w && p + w <= n;

    /// <summary>Quotient-region blocks in lexicographic order.</summary>
    public static IEnumerable<(int P, int W)> FundamentalDomain(int n)
    {
        foreach (var (p, w) in AllBlocks(n))
            if (InFundamentalDomain(n, p, w)) yield return (p, w);
    }

    /// <summary>The 8 D₄ images of (p,w) with the spectral cocycle: FoldParity 0 = same spectrum
    /// (identity, transpose, Klein, transpose·Klein), FoldParity 1 = folded spectrum λ ↦ −λ−2N
    /// (bra fold, ket fold, and their transpose partners). Duplicates occur on the fixed lines
    /// (p = w, p + w = N, and the even-N mid-lines) — callers dedup if needed.</summary>
    public static IReadOnlyList<(int P, int W, int FoldParity)> OrbitImages(int n, int p, int w) =>
        new (int, int, int)[]
        {
            (p, w, 0), (w, p, 0), (n - p, n - w, 0), (n - w, n - p, 0),
            (p, n - w, 1), (n - p, w, 1), (n - w, p, 1), (w, n - p, 1),
        };

    /// <summary>n_diff on block (p,w) ranges over [|p−w|, min(p+w, 2N−p−w)] in steps of 2
    /// (the window-combinatorics shell lemma; pinned from below in BlockLatticeTests).</summary>
    public static (int NMin, int NMax) NdiffRange(int n, int p, int w) =>
        (Math.Abs(p - w), Math.Min(p + w, 2 * n - p - w));

    /// <summary>The Bendixson rate window at real q: Re spec(p,w) ⊆ [−2·n_max, −2·n_min].</summary>
    public static (double Lo, double Hi) RateWindow(int n, int p, int w)
    {
        var (nMin, nMax) = NdiffRange(n, p, w);
        return (-2.0 * nMax, -2.0 * nMin);
    }

    /// <summary>Distance of a real shift to the rate window (0 inside): the analytic exclusion margin.
    /// At real q, σ_min(L−s) ≥ dist(s, W(L)) ≥ this margin (Bendixson; the Hermitian part of the block
    /// is exactly the −2·n_diff diagonal at real q, the hopping being anti-Hermitian).</summary>
    public static double WindowDistance(int n, int p, int w, double reShift)
    {
        var (lo, hi) = RateWindow(n, p, w);
        if (reShift < lo) return lo - reShift;
        if (reShift > hi) return reShift - hi;
        return 0.0;
    }

    /// <summary>dim(p,w) = C(N,p)·C(N,w).</summary>
    public static long Dim(int n, int p, int w) => Binomial(n, p) * Binomial(n, w);

    static long Binomial(int n, int k)
    {
        if (k < 0 || k > n) return 0;
        long r = 1;
        for (int i = 1; i <= k; i++) r = r * (n - k + i) / i;   // r stays C(n−k+i, i), each division exact
        return r;
    }
}
