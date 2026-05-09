using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class F65XxChainSpectrumPi2InheritanceTests
{
    private static F65XxChainSpectrumPi2Inheritance BuildClaim()
    {
        var ladder = new Pi2DyadicLadderClaim();
        var qubitAnchor = new QubitDimensionalAnchorClaim();
        var f66 = new F66PoleModesPi2Inheritance(ladder, qubitAnchor);
        return new F65XxChainSpectrumPi2Inheritance(ladder, f66);
    }

    [Fact]
    public void Tier_IsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived, BuildClaim().Tier);
    }

    [Fact]
    public void NumeratorCoefficient_IsExactlyFour()
    {
        Assert.Equal(4.0, BuildClaim().NumeratorCoefficient, precision: 14);
    }

    [Fact]
    public void UpperBoundCoefficient_IsExactlyTwo()
    {
        Assert.Equal(2.0, BuildClaim().UpperBoundCoefficient, precision: 14);
    }

    [Fact]
    public void LowerBoundCoefficient_IsExactlyZero()
    {
        Assert.Equal(0.0, BuildClaim().LowerBoundCoefficient, precision: 14);
    }

    [Theory]
    // F65 verified table from ANALYTICAL_FORMULAS:
    // N=3: α/γ₀ ∈ {1/2, 1, 1/2}
    // N=5: α/γ₀ ∈ {1/6, 1/2, 2/3, 1/2, 1/6}
    [InlineData(3, 1, 1.0, 0.5)]
    [InlineData(3, 2, 1.0, 1.0)]
    [InlineData(3, 3, 1.0, 0.5)]
    [InlineData(5, 1, 1.0, 1.0 / 6.0)]
    [InlineData(5, 2, 1.0, 0.5)]
    [InlineData(5, 3, 1.0, 2.0 / 3.0)]
    [InlineData(5, 4, 1.0, 0.5)]
    [InlineData(5, 5, 1.0, 1.0 / 6.0)]
    public void SingleExcitationRate_MatchesF65VerifiedTable(int N, int k, double gammaZero, double expected)
    {
        Assert.Equal(expected, BuildClaim().SingleExcitationRate(N, k, gammaZero), precision: 8);
    }

    [Theory]
    [InlineData(3)]
    [InlineData(5)]
    [InlineData(7)]
    [InlineData(15)]
    public void RatesLieInF66Interval_AcrossN(int N)
    {
        Assert.True(BuildClaim().RatesLieInF66Interval(N, 0.05));
    }

    [Theory]
    [InlineData(3, 1)]
    [InlineData(5, 2)]
    [InlineData(7, 3)]
    [InlineData(11, 5)]
    public void MirrorSymmetryHolds_AcrossNK(int N, int k)
    {
        Assert.True(BuildClaim().MirrorSymmetryHolds(N, k, 0.05));
    }

    [Theory]
    [InlineData(3, 4.0 / 4.0)]    // 4/(3+1) = 1
    [InlineData(5, 4.0 / 6.0)]    // 4/(5+1) = 2/3
    [InlineData(7, 4.0 / 8.0)]    // 4/(7+1) = 1/2
    public void MaxRateCoefficient_Is4OverNPlus1(int N, double expected)
    {
        Assert.Equal(expected, BuildClaim().MaxRateCoefficient(N), precision: 12);
    }

    [Theory]
    // F65's bonding-mode population matches F75's exactly
    [InlineData(5, 2, 0)]
    [InlineData(5, 2, 2)]
    [InlineData(7, 3, 1)]
    public void BondingModePopulation_MatchesF75ClosedForm(int N, int k, int site)
    {
        var f = BuildClaim();
        double expected = (2.0 / (N + 1)) * Math.Pow(Math.Sin(Math.PI * k * (site + 1) / (N + 1)), 2);
        Assert.Equal(expected, f.BondingModePopulation(N, k, site), precision: 12);
    }

    [Fact]
    public void SingleExcitationRate_KOutOfRange_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().SingleExcitationRate(5, 0, 1.0));
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().SingleExcitationRate(5, 6, 1.0));
    }

    [Fact]
    public void Constructor_NullLadder_Throws()
    {
        var ladder = new Pi2DyadicLadderClaim();
        var qubitAnchor = new QubitDimensionalAnchorClaim();
        var f66 = new F66PoleModesPi2Inheritance(ladder, qubitAnchor);
        Assert.Throws<ArgumentNullException>(() =>
            new F65XxChainSpectrumPi2Inheritance(null!, f66));
    }

    [Fact]
    public void Constructor_NullF66_Throws()
    {
        Assert.Throws<ArgumentNullException>(() =>
            new F65XxChainSpectrumPi2Inheritance(new Pi2DyadicLadderClaim(), null!));
    }
}
