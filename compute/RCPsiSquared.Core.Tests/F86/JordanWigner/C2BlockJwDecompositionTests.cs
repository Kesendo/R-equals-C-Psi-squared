using System;
using System.Linq;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.F86.JordanWigner;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Resonance;
using Xunit;

namespace RCPsiSquared.Core.Tests.F86.JordanWigner;

public class C2BlockJwDecompositionTests
{
    [Theory]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    [InlineData(8)]
    public void F73SumRuleResidual_BelowTolerance(int N)
    {
        var block = new CoherenceBlock(N: N, n: 1, gammaZero: 0.05);
        var decomp = C2BlockJwDecomposition.Build(block);
        Assert.True(decomp.F73SumRuleResidual < C2BlockJwDecomposition.F73Tolerance,
            $"N={N}: F73 sum-rule residual {decomp.F73SumRuleResidual:G3} ≥ tolerance");
    }

    [Theory]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    [InlineData(8)]
    public void EndpointBond_HasMoreUniformKDistribution_ThanInnermostInteriorBond(int N)
    {
        // Physical claim: Endpoint bonds (b ∈ {0, N−2}) sit at the OBC boundary, where the
        // sine modes that confine away from edges (low-k) are weakest; the absolute mass of
        // the bond's bilinear coefficients C_b(k1, k2) is therefore distributed more
        // uniformly across all k. The innermost Interior bond (deepest into the chain,
        // farthest from any boundary) is most localized to the smoothest sine modes and
        // carries the smallest L1 mass.
        //
        // The ‖C_b‖_F is *constant* (= √2 by sine-mode orthonormality), so we discriminate by
        // L1 norm Σ_{k1,k2} |C_b[k1, k2]|. Verified pattern (Python cross-check at unit J):
        //   N=5: Endpoint = 5.887, innermost Interior(b=2) = 5.774 → Endpoint > Interior
        //   N=6: Endpoint = 7.054, innermost Interior(b=2) = 5.484 → Endpoint > Interior
        //   N=7: Endpoint = 8.200, innermost Interior(b=3) = 7.983 → Endpoint > Interior
        //   N=8: Endpoint = 9.330, innermost Interior(b=3) = 7.147 → Endpoint > Interior
        //
        // The spec's "high-k mass fraction over rows k1 ≥ N/2" metric is symmetry-protected
        // (exactly 0.5 for even N, near-0.5 for odd N) and cannot discriminate; switched to
        // the L1-norm-vs-innermost-Interior comparison that captures the same physical claim.
        var block = new CoherenceBlock(N: N, n: 1, gammaZero: 0.05);
        var decomp = C2BlockJwDecomposition.Build(block);
        var endpoint = decomp.Bonds.First(b => b.BondClass == BondClass.Endpoint);
        var innermostInterior = InnermostInterior(decomp.Bonds, N);
        if (innermostInterior is null) return;  // N ≤ 4: no Interior bonds
        double endpointL1 = L1Mass(endpoint.Coefficients);
        double interiorL1 = L1Mass(innermostInterior.Coefficients);
        Assert.True(endpointL1 > interiorL1,
            $"N={N}: Endpoint bond should have larger L1 mass than innermost Interior; endpoint={endpointL1:F3}, interior={interiorL1:F3}");
    }

    private static BondJwCoefficients? InnermostInterior(IReadOnlyList<BondJwCoefficients> bonds, int N)
    {
        // Innermost = bond closest to the chain center (N−2)/2; pick by minimal absolute
        // distance from center among Interior bonds.
        double center = (N - 2) / 2.0;
        BondJwCoefficients? best = null;
        double bestDist = double.PositiveInfinity;
        foreach (var bond in bonds)
        {
            if (bond.BondClass != BondClass.Interior) continue;
            double d = Math.Abs(bond.Bond - center);
            if (d < bestDist) { bestDist = d; best = bond; }
        }
        return best;
    }

    private static double L1Mass(Matrix<double> coeffs)
    {
        double total = 0;
        for (int i = 0; i < coeffs.RowCount; i++)
            for (int j = 0; j < coeffs.ColumnCount; j++)
                total += Math.Abs(coeffs[i, j]);
        return total;
    }

    [Theory]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    [InlineData(8)]
    public void NumBonds_MatchesBlockNumBonds(int N)
    {
        var block = new CoherenceBlock(N: N, n: 1, gammaZero: 0.05);
        var decomp = C2BlockJwDecomposition.Build(block);
        Assert.Equal(block.NumBonds, decomp.Bonds.Count);
    }

    [Theory]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    [InlineData(8)]
    public void BondClassTagging_FirstAndLastEndpoint_OthersInterior(int N)
    {
        var block = new CoherenceBlock(N: N, n: 1, gammaZero: 0.05);
        var decomp = C2BlockJwDecomposition.Build(block);
        int last = decomp.Bonds.Count - 1;
        Assert.Equal(BondClass.Endpoint, decomp.Bonds[0].BondClass);
        Assert.Equal(BondClass.Endpoint, decomp.Bonds[last].BondClass);
        for (int b = 1; b < last; b++)
            Assert.Equal(BondClass.Interior, decomp.Bonds[b].BondClass);
    }

    [Theory]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    [InlineData(8)]
    public void Coefficients_AreSymmetric_InK1K2(int N)
    {
        var block = new CoherenceBlock(N: N, n: 1, gammaZero: 0.05);
        var decomp = C2BlockJwDecomposition.Build(block);
        foreach (var bond in decomp.Bonds)
        {
            for (int k1 = 0; k1 < N; k1++)
                for (int k2 = 0; k2 < N; k2++)
                {
                    double a = bond.Coefficients[k1, k2];
                    double b = bond.Coefficients[k2, k1];
                    Assert.Equal(a, b, precision: 12);
                }
        }
    }

    [Fact]
    public void Tier_IsTier2Verified()
    {
        var block = new CoherenceBlock(N: 5, n: 1, gammaZero: 0.05);
        var decomp = C2BlockJwDecomposition.Build(block);
        Assert.Equal(Tier.Tier2Verified, decomp.Tier);
    }

    [Fact]
    public void Anchor_References_PROOF_F86_QPEAK_JwTrack()
    {
        var block = new CoherenceBlock(N: 5, n: 1, gammaZero: 0.05);
        var decomp = C2BlockJwDecomposition.Build(block);
        Assert.Contains("PROOF_F86_QPEAK", decomp.Anchor);
        Assert.Contains("JW track", decomp.Anchor);
    }

    [Fact]
    public void Build_RejectsNonC2Block()
    {
        // c=2 stratum requires n=1 for small N; pick a block with c != 2.
        // For N=4, n=2: c = min(2, 4-1-2) + 1 = min(2, 1) + 1 = 2. Still c=2.
        // For N=5, n=2: c = min(2, 5-1-2) + 1 = min(2, 2) + 1 = 3. c=3 ≠ 2.
        var block = new CoherenceBlock(N: 5, n: 2, gammaZero: 0.05);
        Assert.Equal(3, block.C);
        Assert.Throws<ArgumentException>(() => C2BlockJwDecomposition.Build(block));
    }

    [Theory]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    [InlineData(8)]
    public void Modes_AreXyJordanWignerModes_WithMatchingN(int N)
    {
        var block = new CoherenceBlock(N: N, n: 1, gammaZero: 0.05);
        var decomp = C2BlockJwDecomposition.Build(block);
        Assert.NotNull(decomp.Modes);
        Assert.Equal(N, decomp.Modes.N);
    }
}
