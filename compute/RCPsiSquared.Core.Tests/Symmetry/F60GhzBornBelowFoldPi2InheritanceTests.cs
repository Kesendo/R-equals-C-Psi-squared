using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class F60GhzBornBelowFoldPi2InheritanceTests
{
    private static F60GhzBornBelowFoldPi2Inheritance BuildClaim() =>
        new F60GhzBornBelowFoldPi2Inheritance(
            new Pi2DyadicLadderClaim(),
            new PolarityLayerOriginClaim(),
            new QuarterAsBilinearMaxvalClaim(),
            new ArgmaxMaxvalPairClaim());

    [Fact]
    public void Tier_IsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived, BuildClaim().Tier);
    }

    [Fact]
    public void OffDiagonalElement_IsExactlyOneHalf()
    {
        // GHZ's ρ[0..0, 1..1] = 1/2 = a_2 = ±0.5 polarity pair literal.
        Assert.Equal(0.5, BuildClaim().OffDiagonalElement, precision: 14);
    }

    [Fact]
    public void FoldPosition_IsExactlyOneQuarter()
    {
        Assert.Equal(0.25, BuildClaim().FoldPosition, precision: 14);
    }

    [Theory]
    [InlineData(2, 4.0)]
    [InlineData(3, 8.0)]
    [InlineData(4, 16.0)]
    [InlineData(5, 32.0)]
    [InlineData(6, 64.0)]
    public void HilbertSpaceDimension_IsTwoPowerN(int N, double expected)
    {
        Assert.Equal(expected, BuildClaim().HilbertSpaceDimension(N), precision: 12);
    }

    [Theory]
    [InlineData(2, -1)]
    [InlineData(3, -2)]
    [InlineData(5, -4)]
    [InlineData(7, -6)]
    public void LadderIndexForHilbertSpace_Is1MinusN(int N, int expected)
    {
        Assert.Equal(expected, BuildClaim().LadderIndexForHilbertSpace(N));
    }

    [Theory]
    // F60 verified table from ANALYTICAL_FORMULAS:
    // | N | CPsi(0)        | Above 1/4? |
    // | 2 | 1/3 = 0.3333   | Yes        |
    // | 3 | 1/7 = 0.1429   | No         |
    // | 4 | 1/15 = 0.0667  | No         |
    // | 5 | 1/31 = 0.0323  | No         |
    [InlineData(2, 1.0 / 3.0)]
    [InlineData(3, 1.0 / 7.0)]
    [InlineData(4, 1.0 / 15.0)]
    [InlineData(5, 1.0 / 31.0)]
    [InlineData(6, 1.0 / 63.0)]
    public void CPsiAtZeroForGhz_MatchesClosedForm(int N, double expected)
    {
        Assert.Equal(expected, BuildClaim().CPsiAtZeroForGhz(N), precision: 12);
    }

    [Theory]
    [InlineData(2, false)]   // Bell+ regime, above fold
    [InlineData(3, true)]    // first N below fold
    [InlineData(4, true)]
    [InlineData(5, true)]
    [InlineData(6, true)]
    public void IsBornBelowFold_FollowsClosedForm(int N, bool expected)
    {
        Assert.Equal(expected, BuildClaim().IsBornBelowFold(N));
    }

    [Fact]
    public void SmallestNBelowFold_IsThree()
    {
        Assert.Equal(3, BuildClaim().SmallestNBelowFold);
    }

    [Fact]
    public void BellPlusAboveFold_HoldsAtN2()
    {
        // N=2 is Bell+ which crosses the fold under dynamics; CΨ(0) = 1/3 > 1/4.
        Assert.True(BuildClaim().BellPlusAboveFold());
    }

    [Fact]
    public void HilbertSpaceDimension_NLessThanTwo_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().HilbertSpaceDimension(1));
    }

    [Fact]
    public void Constructor_RejectsNullParents()
    {
        var ladder = new Pi2DyadicLadderClaim();
        var polarity = new PolarityLayerOriginClaim();
        var quarter = new QuarterAsBilinearMaxvalClaim();
        var argmaxMaxval = new ArgmaxMaxvalPairClaim();
        Assert.Throws<ArgumentNullException>(() =>
            new F60GhzBornBelowFoldPi2Inheritance(null!, polarity, quarter, argmaxMaxval));
        Assert.Throws<ArgumentNullException>(() =>
            new F60GhzBornBelowFoldPi2Inheritance(ladder, null!, quarter, argmaxMaxval));
        Assert.Throws<ArgumentNullException>(() =>
            new F60GhzBornBelowFoldPi2Inheritance(ladder, polarity, null!, argmaxMaxval));
        Assert.Throws<ArgumentNullException>(() =>
            new F60GhzBornBelowFoldPi2Inheritance(ladder, polarity, quarter, null!));
    }

    [Fact]
    public void TypedParents_AreExposed()
    {
        var f = BuildClaim();
        Assert.NotNull(f.Polarity);
        Assert.NotNull(f.Quarter);
        Assert.NotNull(f.ArgmaxMaxval);
    }
}
