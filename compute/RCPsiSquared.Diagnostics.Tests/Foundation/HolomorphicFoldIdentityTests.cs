using System;
using System.Linq;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.F89PathK;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;
using Xunit.Abstractions;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>The holomorphic-fold identity (premise 1 of the complex-q resultant route, remainder R1), now
/// DERIVED and shown to be FULL-SPECTRUM, not residual-specific. The diagonal core (p_c,p_c) = (3,3) carries the
/// HOLOMORPHIC fold μ = −λ − 2N of every (2,3)-band eigenvalue λ, exactly and at every q, real AND complex:
///
/// <code>  spec(3,3)(q)  =  { −z − 2N : z ∈ spec(2,3)(q) }   as multisets, ∀q  </code>
///
/// <para>This is a COROLLARY of three exact Tier-1 similarities, the §7 diamond maps composed
/// (<c>docs/proofs/PROOF_CODIM1_BY_ADDITIVITY.md</c> §7):
/// <list type="number">
/// <item>the climbing W-step (§5): spec(1,2)(q) ⊆ spec(2,3)(q) (injective intertwiner, σ_min = √2);</item>
/// <item>the transpose (§7.2): spec(3,2)(q) = conj(spec(2,3)(q̄));</item>
/// <item>F89d (§7.4, <see cref="RCPsiSquared.Core.Symmetry.F89CrossFoldSimilarityClaim"/>):
///       L(3,3)(q̄) = −P·conj(L(3,2)(q))·Pᵀ − 2N, hence spec(3,3)(q̄) = −conj(spec(3,2)(q)) − 2N.</item>
/// </list>
/// Composing 2 then 3 (substitute q→q̄), the TWO conjugations cancel HOLOMORPHICALLY at every q:
/// spec(3,3)(q) = −conj(spec(3,2)(q̄)) − 2N = −conj(conj(spec(2,3)(q))) − 2N = −spec(2,3)(q) − 2N. With link 1,
/// every (1,2) eigenvalue λ (residual or AT) folds: −λ − 2N ∈ spec(3,3)(q). So the R1 chain's former "one
/// underived entry" is derived, and it is stronger than the resultant needs (the full band, not just the F_18
/// residual, folds).</para>
///
/// <para>CORRECTION (2026-07-03): the earlier reading "the fold is RESIDUAL-SPECIFIC; the full-spectrum identity
/// spec(3,3) = −spec(band) − 2N is FALSE (≈ 44)" was a SORTED-ZIP ARTIFACT. A sort-by-(Re,Im)-then-zip
/// comparison is unreliable here: the fold negates Im, and the band spectrum has real-part rungs of multiplicity
/// up to 10 (self-conjugate at real q), so within a rung the two sides order oppositely and the zip pairs +Im
/// against −Im, manufacturing an O(2·max|Im|) ≈ 44 gap. Under a rigorous multiset metric (optimal matching, or a
/// rounded sort key) the distance is ~10⁻¹³ at every q. The full-spectrum identity HOLDS; the fold is not
/// residual-specific. Verified independently (numpy, 40-digit mpmath at the branch locus, negative controls) and
/// in the real C# builder here. <see cref="SortedZip_IsAnArtifact_Guard"/> pins this so it is not re-introduced.</para>
///
/// <para>R3 (the gap byte-identity) falls out of the same identity: the near-defective pair's eigenvalue
/// separation is preserved by the fold (|(−λ₁−2N)−(−λ₂−2N)| = |λ₁−λ₂|), and the (1,2) seed pair sits in
/// spec(3,3) by links 1+fold, so the (3,3) core gap equals the (1,2) seed gap identically
/// (<see cref="R3_GapByteIdentity_SeedPairGap_AppearsIn_Core33"/>).</para>
///
/// <para>Run: <c>dotnet test "compute/RCPsiSquared.Diagnostics.Tests" --filter "Category=FOLDIDENTITY"
/// --logger "console;verbosity=detailed"</c>.</para></summary>
public sealed class HolomorphicFoldIdentityTests
{
    private readonly ITestOutputHelper _out;
    public HolomorphicFoldIdentityTests(ITestOutputHelper o) => _out = o;

    private static Complex[] Spec(int n, int wk, int wb, Complex q) =>
        Matrix<Complex>.Build.DenseOfArray(WeightCoherenceBlock.Build(n, wk, wb, q)).Evd().EigenValues.ToArray();

    // Rigorous multiset distance: greedy nearest-with-removal. For equal multisets it returns ~machine zero and
    // is a valid witness of set-equality; it is NOT fooled by the plane-reflecting fold the way sort-zip is.
    private static double MultisetDist(Complex[] a, Complex[] b)
    {
        if (a.Length != b.Length) return double.NaN;
        var used = new bool[b.Length];
        double worst = 0;
        foreach (var z in a)
        {
            int best = -1; double bd = double.PositiveInfinity;
            for (int j = 0; j < b.Length; j++)
                if (!used[j]) { double d = (z - b[j]).Magnitude; if (d < bd) { bd = d; best = j; } }
            used[best] = true;
            worst = Math.Max(worst, bd);
        }
        return worst;
    }

    // Each element of `sub` has a match in `super` (multiset subset), reporting the worst nearest distance.
    private static double SubsetDist(Complex[] sub, Complex[] super) =>
        sub.Max(z => super.Min(w => (z - w).Magnitude));

    private static readonly Complex[] Qs =
    {
        new(3, 0), new(2, 0), new(1.8, 0.4), new(0.9938, 0.1825),
        new(2, 1), new(0.37, -1.4), new(0, 1),
    };

    // ── The headline: the FULL-SPECTRUM holomorphic identity, real AND complex q (replaces the artifact test). ──
    [Fact]
    [Trait("Category", "FOLDIDENTITY")]
    public void FullSpectrumHolomorphicFold_Band23_To_Core33_RealAndComplexQ()
    {
        const int N = 5, twoN = 2 * N;
        foreach (var q in Qs)
        {
            var band23 = Spec(N, 2, 3, q);
            var core33 = Spec(N, 3, 3, q);
            var fold = band23.Select(z => -z - twoN).ToArray();
            double d = MultisetDist(core33, fold);
            _out.WriteLine($"q={q}: multiset |spec(3,3) − (−spec(2,3)−2N)| = {d:E2}");
            Assert.True(d < 1e-9, $"the FULL-SPECTRUM holomorphic fold must hold at q={q}; multiset dist {d:E2}");
        }
    }

    // ── Every (1,2) eigenvalue folds (subsumes the F_18 residual roots the resultant composes against). ──
    [Fact]
    [Trait("Category", "FOLDIDENTITY")]
    public void AllOf_Block12_FoldsInto_Core33_AndSoDoTheExactResidualRoots()
    {
        const int N = 5, twoN = 2 * N;
        foreach (var q in Qs)
        {
            var block12 = Spec(N, 1, 2, q);            // dim 50: AT strands ⊎ residual
            var core33 = Spec(N, 3, 3, q);             // dim 100
            double d = SubsetDist(block12.Select(z => -z - twoN).ToArray(), core33);
            Assert.True(d < 1e-9, $"all of spec(1,2) must fold into spec(3,3) at q={q}; worst {d:E2}");
        }
        // The exact-arithmetic F_18 residual roots (mirror basis, ×2-cleared) are the resultant-relevant subset;
        // they fold too (this was the original passing gate, now a corollary of the full identity above).
        foreach (int q0 in new[] { 2, 3, 5, 7 })
        {
            var q = new Complex(q0, 0);
            var residual = PathKMonodromyScout.ResidualRootsExact(N - 1, q);
            var core33 = Spec(N, 3, 3, q);
            double d = SubsetDist(residual.Select(lam => -lam - twoN).ToArray(), core33);
            _out.WriteLine($"q0={q0}: {residual.Length} residual roots; worst min|spec(3,3) − (−λ−2N)| = {d:E2}");
            Assert.True(d < 1e-8, $"the F_18 residual roots must fold at q0={q0}; worst {d:E2}");
        }
    }

    // ── The three exact links that DERIVE the fold (so it is a corollary, not an observation). ──
    [Fact]
    [Trait("Category", "FOLDIDENTITY")]
    public void DerivationLinks_W_Transpose_F89d_ComposeToTheFold()
    {
        const int N = 5, twoN = 2 * N;
        foreach (var q in Qs)
        {
            var qb = Complex.Conjugate(q);
            // Link 1 (W): spec(1,2)(q) ⊆ spec(2,3)(q).
            double dW = SubsetDist(Spec(N, 1, 2, q), Spec(N, 2, 3, q));
            Assert.True(dW < 1e-9, $"W link spec(1,2)⊆spec(2,3) at q={q}; worst {dW:E2}");

            // Link 2 (transpose): spec(3,2)(q) = conj(spec(2,3)(q̄)).
            double dT = MultisetDist(Spec(N, 3, 2, q), Spec(N, 2, 3, qb).Select(Complex.Conjugate).ToArray());
            Assert.True(dT < 1e-9, $"transpose link at q={q}; {dT:E2}");

            // Link 3 (F89d): the bit-exact matrix identity L(3,3)(q̄) = −P·conj(L(3,2)(q))·Pᵀ − 2N.
            var l32 = WeightCoherenceBlock.Build(N, 3, 2, q);
            var l33qb = WeightCoherenceBlock.Build(N, 3, 3, qb);
            int[] p = WeightCoherenceBlock.BraComplementPermutation(N, 3, 2);   // (3,2) index → (3,3) index
            int d32 = l32.GetLength(0);
            double f89d = 0;
            for (int t = 0; t < d32; t++)
                for (int u = 0; u < d32; u++)
                {
                    // (−P conj(L32) Pᵀ − 2N)[p[t], p[u]] = −conj(L32[t,u]) − 2N·δ_tu
                    var rhs = -Complex.Conjugate(l32[t, u]) - (t == u ? new Complex(twoN, 0) : Complex.Zero);
                    f89d = Math.Max(f89d, (l33qb[p[t], p[u]] - rhs).Magnitude);
                }
            Assert.True(f89d < 1e-9, $"F89d matrix identity at q={q}; {f89d:E2}");

            _out.WriteLine($"q={q}: W {dW:E2}, transpose {dT:E2}, F89d(matrix) {f89d:E2}");
        }
    }

    // ── R3: the gap byte-identity is a corollary (the fold preserves the seed pair's eigenvalue separation). ──
    [Fact]
    [Trait("Category", "FOLDIDENTITY")]
    public void R3_GapByteIdentity_SeedPairGap_AppearsIn_Core33()
    {
        const int N = 5, twoN = 2 * N;
        // The two real N=5 defective loci and the seed/core eigenvalue centres.
        foreach (var (q0, lamA) in new[] { (0.620878, -4.6189), (1.077615, -3.7917) })
        {
            var q = new Complex(q0, 0);
            var lamB = -lamA - twoN;                                 // holomorphic fold (real at a real locus)
            var seed = Spec(N, 1, 2, q).OrderBy(z => (z - lamA).Magnitude).Take(2).ToArray();
            var core = Spec(N, 3, 3, q).OrderBy(z => (z - lamB).Magnitude).Take(2).ToArray();
            double seedGap = (seed[0] - seed[1]).Magnitude;
            double coreGap = (core[0] - core[1]).Magnitude;
            // the folded seed pair must be the core pair (so the gaps are equal by |·|-preservation)
            double foldMatch = SubsetDist(seed.Select(z => -z - twoN).ToArray(), core);
            _out.WriteLine($"q*={q0}: seed gap {seedGap:E4}, core gap {coreGap:E4}, |Δ| {Math.Abs(seedGap - coreGap):E2}, foldMatch {foldMatch:E2}");
            Assert.True(Math.Abs(seedGap - coreGap) < 1e-6, $"R3 gap byte-identity at q*={q0}: {seedGap:E4} vs {coreGap:E4}");
            Assert.True(foldMatch < 1e-6, $"the folded seed pair must be the core pair at q*={q0}; {foldMatch:E2}");
        }
    }

    // ── Guard: the sorted-zip comparison is an ARTIFACT (spurious ≈44 at real q); do not re-introduce it. ──
    [Fact]
    [Trait("Category", "FOLDIDENTITY")]
    public void SortedZip_IsAnArtifact_Guard()
    {
        const int N = 5, twoN = 2 * N;
        var q = new Complex(3, 0);
        var core33 = Spec(N, 3, 3, q);
        var fold = Spec(N, 2, 3, q).Select(z => -z - twoN).ToArray();

        double sortedZip = core33.OrderBy(z => z.Real).ThenBy(z => z.Imaginary)
            .Zip(fold.OrderBy(z => z.Real).ThenBy(z => z.Imaginary), (a, b) => (a - b).Magnitude).Max();
        double rounded = core33.OrderBy(z => Math.Round(z.Real, 6)).ThenBy(z => Math.Round(z.Imaginary, 6))
            .Zip(fold.OrderBy(z => Math.Round(z.Real, 6)).ThenBy(z => Math.Round(z.Imaginary, 6)),
                 (a, b) => (a - b).Magnitude).Max();
        double multiset = MultisetDist(core33, fold);

        _out.WriteLine($"raw sorted-zip = {sortedZip:E2} (ARTIFACT), rounded-key sorted-zip = {rounded:E2}, multiset = {multiset:E2}");
        Assert.True(sortedZip > 1.0, "raw sorted-zip is spuriously large at real q (the artifact this guards against)");
        Assert.True(rounded < 1e-9 && multiset < 1e-9,
            "the same data matches to machine zero under a rounded key / multiset metric: the fold IS full-spectrum");
    }
}
