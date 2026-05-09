using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class F66PoleModesPi2InheritanceTests
{
    private static F66PoleModesPi2Inheritance BuildClaim() =>
        new F66PoleModesPi2Inheritance(
            new Pi2DyadicLadderClaim(),
            new QubitDimensionalAnchorClaim());

    [Fact]
    public void Tier_IsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived, BuildClaim().Tier);
    }

    [Fact]
    public void LowerPoleAlpha_IsZero()
    {
        Assert.Equal(0.0, BuildClaim().LowerPoleAlpha, precision: 14);
    }

    [Fact]
    public void UpperPoleCoefficient_IsExactlyTwo()
    {
        Assert.Equal(2.0, BuildClaim().UpperPoleCoefficient, precision: 14);
    }

    [Fact]
    public void UpperPoleCoefficientMatchesQubitDimension_HoldsExactly()
    {
        Assert.True(BuildClaim().UpperPoleCoefficientMatchesQubitDimension());
    }

    [Theory]
    [InlineData(0.0, 0.0)]
    [InlineData(0.05, 0.1)]
    [InlineData(0.5, 1.0)]
    [InlineData(1.0, 2.0)]
    [InlineData(2.5, 5.0)]
    public void UpperPoleAlpha_IsTwoTimesGammaZero(double gammaZero, double expected)
    {
        Assert.Equal(expected, BuildClaim().UpperPoleAlpha(gammaZero), precision: 12);
    }

    [Theory]
    [InlineData(0.05, 0.1)]
    [InlineData(1.0, 2.0)]
    public void DissipationIntervalWidth_EqualsUpperPoleAlpha(double gammaZero, double expected)
    {
        var f = BuildClaim();
        Assert.Equal(expected, f.DissipationIntervalWidth(gammaZero), precision: 12);
        Assert.Equal(f.UpperPoleAlpha(gammaZero), f.DissipationIntervalWidth(gammaZero), precision: 14);
    }

    [Theory]
    [InlineData(3, 4)]
    [InlineData(4, 5)]
    [InlineData(5, 6)]
    [InlineData(6, 7)]
    [InlineData(7, 8)]
    public void EndpointMultiplicity_IsNPlus1(int N, int expected)
    {
        Assert.Equal(expected, BuildClaim().EndpointMultiplicity(N));
    }

    [Theory]
    [InlineData(3, 0, 3)]
    [InlineData(4, 0, 4)]
    [InlineData(5, 0, 5)]
    [InlineData(7, 0, 7)]
    public void PalindromicWeightPairOfPoles_IsZeroAndN(int N, int expectedLow, int expectedHigh)
    {
        var (l, u) = BuildClaim().PalindromicWeightPairOfPoles(N);
        Assert.Equal(expectedLow, l);
        Assert.Equal(expectedHigh, u);
    }

    [Theory]
    [InlineData(3, 0, 3, true)]
    [InlineData(5, 0, 5, true)]
    [InlineData(5, 0, 4, false)]
    [InlineData(5, 1, 5, false)]
    public void PalindromePartnershipHolds_MatchesExpected(int N, int lower, int upper, bool expected)
    {
        Assert.Equal(expected, BuildClaim().PalindromePartnershipHolds(N, lower, upper));
    }

    [Fact]
    public void UpperPoleAlpha_NegativeGammaZero_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().UpperPoleAlpha(-0.1));
    }

    [Fact]
    public void EndpointMultiplicity_NLessThanThree_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().EndpointMultiplicity(2));
    }

    [Fact]
    public void PalindromicWeightPairOfPoles_NLessThanOne_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().PalindromicWeightPairOfPoles(0));
    }

    [Fact]
    public void Constructor_NullLadder_Throws()
    {
        Assert.Throws<ArgumentNullException>(() =>
            new F66PoleModesPi2Inheritance(
                ladder: null!,
                qubitAnchor: new QubitDimensionalAnchorClaim()));
    }

    [Fact]
    public void Constructor_NullQubitAnchor_Throws()
    {
        Assert.Throws<ArgumentNullException>(() =>
            new F66PoleModesPi2Inheritance(
                ladder: new Pi2DyadicLadderClaim(),
                qubitAnchor: null!));
    }
}
