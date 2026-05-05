using System;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.Decomposition;
using RCPsiSquared.Core.F86.Item1Derivation;
using RCPsiSquared.Core.Knowledge;
using Xunit;

namespace RCPsiSquared.Core.Tests.F86.Item1Derivation;

public class C2BondCouplingTests
{
    [Theory]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    [InlineData(8)]
    public void ProbeBlock_AllBondsAllEntries_MatchFourModeEffective(int N)
    {
        var block = new CoherenceBlock(N, n: 1, gammaZero: 0.05);
        var coupling = C2BondCoupling.Build(block);
        var effective = FourModeEffective.Build(block);

        for (int b = 0; b < block.NumBonds; b++)
        {
            var expected = effective.MhPerBondEff[b];
            for (int alpha = 0; alpha < 2; alpha++)
                for (int beta = 0; beta < 2; beta++)
                {
                    var actual = coupling.ProbeBlockEntry(b, alpha, beta);
                    var diff = (actual - expected[alpha, beta]).Magnitude;
                    Assert.True(diff < 1e-12,
                        $"V_b[{alpha},{beta}] at N={N}, b={b}: expected {expected[alpha, beta]}, got {actual}, diff={diff:E2}");
                }
        }
    }

    [Theory]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    [InlineData(8)]
    public void ProbeBlock_OffDiagonalSumRule_F73Generalisation(int N)
    {
        // Σ_b V_b[α, β] = 0 for α ≠ β (M_H_total is diagonal in channel-uniform basis)
        var block = new CoherenceBlock(N, n: 1, gammaZero: 0.05);
        var coupling = C2BondCoupling.Build(block);

        Complex sum01 = Complex.Zero;
        Complex sum10 = Complex.Zero;
        for (int b = 0; b < block.NumBonds; b++)
        {
            sum01 += coupling.ProbeBlockEntry(b, 0, 1);
            sum10 += coupling.ProbeBlockEntry(b, 1, 0);
        }
        Assert.True(sum01.Magnitude < 1e-12, $"Σ_b V_b[0,1] should vanish (F73), got {sum01}");
        Assert.True(sum10.Magnitude < 1e-12, $"Σ_b V_b[1,0] should vanish (F73), got {sum10}");
    }

    [Fact]
    public void Tier_IsTier1Derived()
    {
        var block = new CoherenceBlock(N: 5, n: 1, gammaZero: 0.05);
        var coupling = C2BondCoupling.Build(block);
        Assert.Equal(Tier.Tier1Derived, coupling.Tier);
    }

    [Fact]
    public void Build_ThrowsIfNotC2()
    {
        var block = new CoherenceBlock(N: 5, n: 2, gammaZero: 0.05);  // c=3
        Assert.Throws<ArgumentException>(() => C2BondCoupling.Build(block));
    }
}
