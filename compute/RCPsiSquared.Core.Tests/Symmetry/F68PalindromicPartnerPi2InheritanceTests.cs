using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class F68PalindromicPartnerPi2InheritanceTests
{
    private static F68PalindromicPartnerPi2Inheritance BuildClaim()
    {
        var ladder = new Pi2DyadicLadderClaim();
        var memoryLoop = new Pi2I4MemoryLoopClaim();
        var f1 = new F1Pi2Inheritance(new RCPsiSquared.Core.F1.F1PalindromeIdentity(), ladder, memoryLoop);
        var qubitAnchor = new QubitDimensionalAnchorClaim();
        var f66 = new F66PoleModesPi2Inheritance(ladder, qubitAnchor);
        var f65 = new F65XxChainSpectrumPi2Inheritance(ladder, f66);
        var f67 = new F67BondingBellPairPi2Inheritance(ladder, f65);
        return new F68PalindromicPartnerPi2Inheritance(ladder, f1, f67);
    }

    [Fact]
    public void Tier_IsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived, BuildClaim().Tier);
    }

    [Fact]
    public void SumCoefficient_IsExactlyTwo()
    {
        // 2 = a_0 on the dyadic ladder = polynomial root d.
        Assert.Equal(2.0, BuildClaim().SumCoefficient, precision: 14);
    }

    [Theory]
    // F68 closed form: α_p = 2γ₀ − α_b. F67 gives α_b = (4γ₀/(N+1))·sin²(π/(N+1)).
    // N=3, γ=0.05: α_b = 0.025, α_p = 0.1 − 0.025 = 0.075
    // N=5, γ=0.05: α_b = 0.05/6 ≈ 0.008333, α_p = 0.1 − 0.008333 ≈ 0.091667
    [InlineData(3, 0.05, 0.075)]
    [InlineData(5, 0.05, 0.1 - 0.05 / 6.0)]
    [InlineData(3, 1.0, 1.5)]      // α_p = 2 − 0.5 = 1.5
    [InlineData(5, 1.0, 2.0 - 1.0 / 6.0)]
    public void PartnerRate_MatchesClosedForm(int N, double gammaZero, double expected)
    {
        Assert.Equal(expected, BuildClaim().PartnerRate(N, gammaZero), precision: 12);
    }

    [Theory]
    // The palindromic sum α_b + α_p = 2γ₀ should hold exactly (machine precision)
    // for all N, all γ₀.
    [InlineData(3, 0.05)]
    [InlineData(4, 0.05)]
    [InlineData(5, 0.05)]
    [InlineData(7, 0.1)]
    [InlineData(13, 0.025)]
    public void PalindromicSumHolds_AcrossAllN(int N, double gammaZero)
    {
        Assert.True(BuildClaim().PalindromicSumHolds(N, gammaZero));
    }

    [Theory]
    [InlineData(3, false)]   // N=3 is rank-2 fourfold-degenerate
    [InlineData(4, true)]
    [InlineData(5, true)]
    [InlineData(7, true)]
    public void IsRankOneOperational_HoldsForNGreaterEqual4(int N, bool expected)
    {
        Assert.Equal(expected, BuildClaim().IsRankOneOperational(N));
    }

    [Fact]
    public void PartnerRate_PlusBondingRate_EqualsTwoGamma()
    {
        // The structural relationship: α_b + α_p = 2γ₀ exactly.
        var f = BuildClaim();
        for (int N = 3; N <= 7; N++)
        {
            double gamma = 0.05;
            double sum = f.PalindromicSum(N, gamma);
            Assert.Equal(2.0 * gamma, sum, precision: 14);
        }
    }

    [Fact]
    public void PartnerRate_NLessThanTwo_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            BuildClaim().PartnerRate(N: 1, gammaZero: 0.05));
    }

    [Fact]
    public void PartnerRate_NegativeGamma_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            BuildClaim().PartnerRate(N: 3, gammaZero: -0.05));
    }

    [Fact]
    public void IsRankOneOperational_NLessThanTwo_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            BuildClaim().IsRankOneOperational(N: 1));
    }

    [Fact]
    public void Constructor_NullLadder_Throws()
    {
        var ladder = new Pi2DyadicLadderClaim();
        var memoryLoop = new Pi2I4MemoryLoopClaim();
        var f1 = new F1Pi2Inheritance(new RCPsiSquared.Core.F1.F1PalindromeIdentity(), ladder, memoryLoop);
        var f66 = new F66PoleModesPi2Inheritance(ladder, new QubitDimensionalAnchorClaim());
        var f65 = new F65XxChainSpectrumPi2Inheritance(ladder, f66);
        var f67 = new F67BondingBellPairPi2Inheritance(ladder, f65);
        Assert.Throws<ArgumentNullException>(() =>
            new F68PalindromicPartnerPi2Inheritance(null!, f1, f67));
    }

    [Fact]
    public void Constructor_NullF1_Throws()
    {
        var ladder = new Pi2DyadicLadderClaim();
        var f66 = new F66PoleModesPi2Inheritance(ladder, new QubitDimensionalAnchorClaim());
        var f65 = new F65XxChainSpectrumPi2Inheritance(ladder, f66);
        var f67 = new F67BondingBellPairPi2Inheritance(ladder, f65);
        Assert.Throws<ArgumentNullException>(() =>
            new F68PalindromicPartnerPi2Inheritance(ladder, null!, f67));
    }

    [Fact]
    public void Constructor_NullF67_Throws()
    {
        var ladder = new Pi2DyadicLadderClaim();
        var memoryLoop = new Pi2I4MemoryLoopClaim();
        var f1 = new F1Pi2Inheritance(new RCPsiSquared.Core.F1.F1PalindromeIdentity(), ladder, memoryLoop);
        Assert.Throws<ArgumentNullException>(() =>
            new F68PalindromicPartnerPi2Inheritance(ladder, f1, null!));
    }
}
