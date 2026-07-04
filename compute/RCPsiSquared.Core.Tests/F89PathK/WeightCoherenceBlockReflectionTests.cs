using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.F89PathK;
using RCPsiSquared.Core.Numerics;
using Xunit;

namespace RCPsiSquared.Core.Tests.F89PathK;

/// <summary>The general-(p,w) site-reflection parity split behind the step-3 census's ~¼ LU-cost lever:
/// R = the site reflection i ↦ N−1−i acting as |a⟩⟨b| ↦ |rev(a)⟩⟨rev(b)|, an involution commuting
/// ENTRY-WISE with the uniform chain block; the orthonormal even/odd orbit bases split the spectrum
/// exactly (spec(full) = spec(even) ⊎ spec(odd)) and σ_min(full−s) = min over sectors. The sector
/// builder assembles DIRECTLY in sector coordinates (never materializes the full block — the property
/// that fits the N=11 blocks under the LP64 dim ≤ 46340 wall); pinned here against the full-block
/// transform at N=5/6, including a block with reflection fixed points ((3,3): palindromic masks).</summary>
public class WeightCoherenceBlockReflectionTests
{
    [Theory]
    [InlineData(5, 2, 3)]
    [InlineData(5, 1, 3)]
    [InlineData(6, 2, 2)]
    public void ReflectionPermutation_IsInvolution_AndCommutesExactly(int n, int wKet, int wBra)
    {
        var perm = WeightCoherenceBlock.ReflectionPermutation(n, wKet, wBra);
        for (int t = 0; t < perm.Length; t++) Assert.Equal(t, perm[perm[t]]);

        var q = new Complex(0.83, 0.11); // generic complex q: commutation is exact entry-wise
        var l = WeightCoherenceBlock.Build(n, wKet, wBra, q);
        int d = perm.Length;
        for (int i = 0; i < d; i++)
            for (int j = 0; j < d; j++)
                Assert.True((l[perm[i], perm[j]] - l[i, j]).Magnitude < 1e-15);
    }

    [Theory]
    [InlineData(5, 2, 3)]
    [InlineData(5, 3, 3)]   // has reflection fixed points (palindromic masks) — the even/odd skew case
    public void Sectors_SplitTheSpectrumExactly(int n, int wKet, int wBra)
    {
        var q = new Complex(0.620878, 0);
        var full = Matrix<Complex>.Build.DenseOfArray(WeightCoherenceBlock.Build(n, wKet, wBra, q));
        var (even, dEven) = WeightCoherenceBlock.BuildReflectionSectorColumnMajor(n, wKet, wBra, q, odd: false);
        var (odd, dOdd) = WeightCoherenceBlock.BuildReflectionSectorColumnMajor(n, wKet, wBra, q, odd: true);
        Assert.Equal(full.RowCount, dEven + dOdd);

        var fullEv = full.Evd().EigenValues.Enumerate().ToList();
        var sectorEv = ToMatrix(even, dEven).Evd().EigenValues.Enumerate()
            .Concat(ToMatrix(odd, dOdd).Evd().EigenValues.Enumerate()).ToList();
        Assert.True(MultisetDistance(fullEv, sectorEv) < 1e-9);
    }

    [Fact]
    public void SigmaMin_FullEqualsSectorMinimum()
    {
        int n = 5, wKet = 2, wBra = 3;
        var q = new Complex(0.620878, 0);
        var shift = new Complex(-4.618886, 0);
        var full = WeightCoherenceBlock.Build(n, wKet, wBra, q);
        double sFull = ShiftedSigmaMin.Estimate(full, shift, maxIter: 500).SigmaMin;

        double sSector = new[] { false, true }.Min(o =>
        {
            var (a, d) = WeightCoherenceBlock.BuildReflectionSectorColumnMajor(n, wKet, wBra, q, o);
            for (int i = 0; i < d; i++) a[(long)i * d + i] -= shift;
            return ShiftedSigmaMin.EstimateColumnMajor(a, d, maxIter: 500).SigmaMin;
        });
        Assert.Equal(sFull, sSector, Math.Max(sFull * 1e-2, 1e-8));
    }

    static Matrix<Complex> ToMatrix(Complex[] colMajor, int d)
    {
        var m = Matrix<Complex>.Build.Dense(d, d);
        for (int c = 0; c < d; c++)
            for (int r = 0; r < d; r++)
                m[r, c] = colMajor[(long)c * d + r];
        return m;
    }

    // greedy nearest-with-removal (the sorted-zip-artifact-safe multiset metric)
    static double MultisetDistance(List<Complex> a, List<Complex> b)
    {
        double worst = 0;
        var pool = b.ToList();
        foreach (var x in a)
        {
            int best = Enumerable.Range(0, pool.Count).OrderBy(i => (pool[i] - x).Magnitude).First();
            worst = Math.Max(worst, (pool[best] - x).Magnitude);
            pool.RemoveAt(best);
        }
        return worst;
    }
}
