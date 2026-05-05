using System;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.Decomposition;
using RCPsiSquared.Core.F86.Item1Derivation;
using RCPsiSquared.Core.Knowledge;
using Xunit;

namespace RCPsiSquared.Core.Tests.F86.Item1Derivation;

public class C2InterChannelAnalyticalTests
{
    [Theory]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    [InlineData(8)]
    public void U0V0_OverlapWithNumericalSvd_AtMachinePrecision(int N)
    {
        var block = new CoherenceBlock(N, n: 1, gammaZero: 0.05);
        var analytical = C2InterChannelAnalytical.Build(block);
        var numerical = InterChannelSvd.Build(block, hd1: 1, hd2: 3);

        var u0Overlap = (analytical.U0.Conjugate() * numerical.U0InFullBlock).Magnitude;
        var v0Overlap = (analytical.V0.Conjugate() * numerical.V0InFullBlock).Magnitude;

        // BOTH paths (Tier1Derived analytical or Tier2Verified numerical) must pass this:
        // Tier1Derived: ansatz overlaps numerical at machine precision.
        // Tier2Verified: U0/V0 ARE the numerical vectors, overlap is exactly 1 to ~1e-15.
        Assert.InRange(u0Overlap, 1.0 - 1e-10, 1.0 + 1e-10);
        Assert.InRange(v0Overlap, 1.0 - 1e-10, 1.0 + 1e-10);
    }

    [Theory]
    [InlineData(5)]
    [InlineData(8)]
    public void Sigma0_MatchesNumericalSvd(int N)
    {
        var block = new CoherenceBlock(N, n: 1, gammaZero: 0.05);
        var analytical = C2InterChannelAnalytical.Build(block);
        var numerical = InterChannelSvd.Build(block, hd1: 1, hd2: 3);

        Assert.InRange(analytical.Sigma0, numerical.Sigma0 - 1e-10, numerical.Sigma0 + 1e-10);
    }

    [Fact]
    public void Tier_IsTier1DerivedOrTier2Verified_NotOther()
    {
        var block = new CoherenceBlock(N: 7, n: 1, gammaZero: 0.05);
        var analytical = C2InterChannelAnalytical.Build(block);
        Assert.True(analytical.Tier == Tier.Tier1Derived || analytical.Tier == Tier.Tier2Verified,
            $"Expected Tier1Derived or Tier2Verified; got {analytical.Tier}");
    }

    [Fact]
    public void IsAnalyticallyDerived_ConsistentWithTier()
    {
        var block = new CoherenceBlock(N: 7, n: 1, gammaZero: 0.05);
        var analytical = C2InterChannelAnalytical.Build(block);
        if (analytical.IsAnalyticallyDerived)
            Assert.Equal(Tier.Tier1Derived, analytical.Tier);
        else
            Assert.Equal(Tier.Tier2Verified, analytical.Tier);
    }

    [Fact]
    public void Constructor_ThrowsIfNotC2()
    {
        var block = new CoherenceBlock(N: 5, n: 2, gammaZero: 0.05);  // c=3
        Assert.Throws<ArgumentException>(() => C2InterChannelAnalytical.Build(block));
    }
}
