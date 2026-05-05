using System;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.Decomposition;
using RCPsiSquared.Core.F86.Item1Derivation;
using RCPsiSquared.Core.Knowledge;
using Xunit;

namespace RCPsiSquared.Core.Tests.F86.Item1Derivation;

public class C2ChannelUniformAnalyticalTests
{
    [Theory]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    [InlineData(8)]
    public void Vectors_OverlapWithFourModeBasis_AtMachinePrecision(int N)
    {
        var block = new CoherenceBlock(N, n: 1, gammaZero: 0.05);
        var analytical = new C2ChannelUniformAnalytical(block);
        var fourMode = FourModeBasis.Build(block);

        // FourModeBasis columns: 0 = |c_1⟩, 1 = |c_3⟩, 2 = |u_0⟩, 3 = |v_0⟩.
        var c1Numerical = fourMode.BasisMatrix.Column(0);
        var c3Numerical = fourMode.BasisMatrix.Column(1);

        // Channel-uniform vectors are real positive uniform, so overlap should be
        // 1.0 to machine precision (no phase ambiguity).
        var overlap1 = (analytical.C1Vector.Conjugate() * c1Numerical).Magnitude;
        var overlap3 = (analytical.C3Vector.Conjugate() * c3Numerical).Magnitude;

        Assert.InRange(overlap1, 1.0 - 1e-12, 1.0 + 1e-12);
        Assert.InRange(overlap3, 1.0 - 1e-12, 1.0 + 1e-12);
    }

    [Theory]
    [InlineData(5)]
    [InlineData(8)]
    public void Vectors_AreUnitNormalized(int N)
    {
        var block = new CoherenceBlock(N, n: 1, gammaZero: 0.05);
        var analytical = new C2ChannelUniformAnalytical(block);
        Assert.InRange(analytical.C1Vector.L2Norm(), 1.0 - 1e-12, 1.0 + 1e-12);
        Assert.InRange(analytical.C3Vector.L2Norm(), 1.0 - 1e-12, 1.0 + 1e-12);
    }

    [Fact]
    public void Tier_IsTier1Derived()
    {
        var block = new CoherenceBlock(N: 5, n: 1, gammaZero: 0.05);
        var analytical = new C2ChannelUniformAnalytical(block);
        Assert.Equal(Tier.Tier1Derived, analytical.Tier);
    }

    [Fact]
    public void Constructor_ThrowsIfNotC2()
    {
        var block = new CoherenceBlock(N: 5, n: 2, gammaZero: 0.05);  // c=3
        Assert.Throws<ArgumentException>(() => new C2ChannelUniformAnalytical(block));
    }
}
