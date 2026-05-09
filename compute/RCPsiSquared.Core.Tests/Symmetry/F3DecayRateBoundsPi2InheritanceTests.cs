using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class F3DecayRateBoundsPi2InheritanceTests
{
    private static F3DecayRateBoundsPi2Inheritance BuildClaim()
    {
        var ladder = new Pi2DyadicLadderClaim();
        var f50 = new F50WeightOneDegeneracyPi2Inheritance(ladder);
        return new F3DecayRateBoundsPi2Inheritance(ladder, f50);
    }

    [Fact]
    public void Tier_IsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived, BuildClaim().Tier);
    }

    [Fact]
    public void RateCoefficient_IsExactlyTwo()
    {
        Assert.Equal(2.0, BuildClaim().RateCoefficient, precision: 14);
    }

    [Theory]
    [InlineData(0.05, 0.1)]
    [InlineData(0.5, 1.0)]
    [InlineData(1.0, 2.0)]
    public void MinRate_Equals2Gamma(double gammaZero, double expected)
    {
        Assert.Equal(expected, BuildClaim().MinRate(gammaZero), precision: 14);
    }

    [Theory]
    [InlineData(2, 0.05, 0.1)]    // 2·1·0.05 = 0.1
    [InlineData(3, 0.05, 0.2)]    // 2·2·0.05 = 0.2
    [InlineData(5, 1.0, 8.0)]     // 2·4·1 = 8
    [InlineData(7, 0.5, 6.0)]     // 2·6·0.5 = 6
    public void MaxRate_Equals2NMinus1Gamma(int N, double gammaZero, double expected)
    {
        Assert.Equal(expected, BuildClaim().MaxRate(N, gammaZero), precision: 14);
    }

    [Theory]
    [InlineData(2, 0.05, 0.0)]    // 2·0·0.05 = 0
    [InlineData(3, 0.05, 0.1)]    // 2·1·0.05 = 0.1
    [InlineData(5, 1.0, 6.0)]     // 2·3·1 = 6
    public void Bandwidth_Equals2NMinus2Gamma(int N, double gammaZero, double expected)
    {
        Assert.Equal(expected, BuildClaim().Bandwidth(N, gammaZero), precision: 14);
    }

    [Theory]
    [InlineData(3, 0.05, 0.3)]    // 2·3·0.05 = 0.3
    [InlineData(5, 0.1, 1.0)]     // 2·5·0.1 = 1.0
    [InlineData(1, 0.5, 1.0)]     // 2·1·0.5 = 1.0
    public void XorBoundary_Equals2NGamma(int N, double gammaZero, double expected)
    {
        Assert.Equal(expected, BuildClaim().XorBoundary(N, gammaZero), precision: 14);
    }

    [Fact]
    public void MinRateMatchesF50_AcrossAllGamma()
    {
        var f = BuildClaim();
        Assert.True(f.MinRateMatchesF50(0.05));
        Assert.True(f.MinRateMatchesF50(0.5));
        Assert.True(f.MinRateMatchesF50(1.0));
    }

    [Theory]
    [InlineData(2, 0.05)]
    [InlineData(5, 0.1)]
    [InlineData(7, 1.0)]
    public void BandwidthIsMaxMinusMin(int N, double gammaZero)
    {
        Assert.True(BuildClaim().BandwidthIsMaxMinusMin(N, gammaZero));
    }

    [Theory]
    [InlineData(2, 0.05)]
    [InlineData(3, 0.05)]
    [InlineData(5, 0.1)]
    public void XorBoundaryAboveMaxRate(int N, double gammaZero)
    {
        Assert.True(BuildClaim().XorBoundaryAboveMaxRate(N, gammaZero));
    }

    [Fact]
    public void HybridMinimumRates_TableHasN4AndN5()
    {
        var table = BuildClaim().HybridMinimumRates;
        Assert.Equal(2, table.Count);
        Assert.Equal((4, 0.98), table[0]);
        Assert.Equal((5, 0.62), table[1]);
    }

    [Fact]
    public void MinRate_NegativeGamma_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().MinRate(-0.05));
    }

    [Fact]
    public void MaxRate_NLessThanTwo_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().MaxRate(N: 1, gammaZero: 0.05));
    }

    [Fact]
    public void Constructor_NullLadder_Throws()
    {
        var ladder = new Pi2DyadicLadderClaim();
        var f50 = new F50WeightOneDegeneracyPi2Inheritance(ladder);
        Assert.Throws<ArgumentNullException>(() =>
            new F3DecayRateBoundsPi2Inheritance(null!, f50));
    }

    [Fact]
    public void Constructor_NullF50_Throws()
    {
        Assert.Throws<ArgumentNullException>(() =>
            new F3DecayRateBoundsPi2Inheritance(new Pi2DyadicLadderClaim(), null!));
    }
}
