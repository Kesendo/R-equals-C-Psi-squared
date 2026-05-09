using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class F56CriticalSlowingPi2InheritanceTests
{
    private static F56CriticalSlowingPi2Inheritance BuildClaim()
    {
        var ladder = new Pi2DyadicLadderClaim();
        var quarter = new QuarterAsBilinearMaxvalClaim();
        var half = new HalfAsStructuralFixedPointClaim();
        return new F56CriticalSlowingPi2Inheritance(ladder, quarter, half);
    }

    [Fact]
    public void Tier_IsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived, BuildClaim().Tier);
    }

    [Fact]
    public void HalfPrefactor_IsExactlyOneHalf()
    {
        Assert.Equal(0.5, BuildClaim().HalfPrefactor, precision: 14);
    }

    [Fact]
    public void FourFactor_IsExactlyFour()
    {
        Assert.Equal(4.0, BuildClaim().FourFactor, precision: 14);
    }

    [Fact]
    public void SixteenFactor_IsExactlySixteen()
    {
        Assert.Equal(16.0, BuildClaim().SixteenFactor, precision: 14);
    }

    [Fact]
    public void NegFourTransient_IsExactlyMinusFour()
    {
        Assert.Equal(-4.0, BuildClaim().NegFourTransient, precision: 14);
    }

    [Fact]
    public void CardioidCuspPosition_IsExactlyOneQuarter()
    {
        Assert.Equal(0.25, BuildClaim().CardioidCuspPosition, precision: 14);
    }

    [Theory]
    // α(tol) = -4 + (1/2)·ln(16·tol)
    // tol = 1: α = -4 + 0.5·ln(16) = -4 + 0.5·2.7726 = -4 + 1.3863 ≈ -2.6137
    // tol = 0.1: α = -4 + 0.5·ln(1.6) = -4 + 0.5·0.4700 = -3.7650
    // tol = 1e-10: α = -4 + 0.5·ln(1.6e-9) = -4 + 0.5·(-20.255) = -14.128
    [InlineData(1.0, -2.6137056388801093)]
    [InlineData(0.1, -3.7649975058299856)]
    public void Alpha_MatchesClosedForm(double tol, double expected)
    {
        Assert.Equal(expected, BuildClaim().Alpha(tol), precision: 6);
    }

    [Theory]
    // K(ε, tol) = (1/2)·ln(4ε/tol) + α(tol)·√ε
    // For ε = 0.01, tol = 0.001:
    //   log term = 0.5·ln(0.04/0.001) = 0.5·ln(40) = 0.5·3.689 = 1.844
    //   α(0.001) = -4 + 0.5·ln(0.016) = -4 + 0.5·(-4.135) = -6.068
    //   sqrt term = -6.068·0.1 = -0.6068
    //   K = 1.844 - 0.607 = 1.237
    [InlineData(0.01, 0.001, 1.23768)]
    public void IterationCount_MatchesClosedForm(double epsilon, double tol, double expected)
    {
        Assert.Equal(expected, BuildClaim().IterationCount(epsilon, tol), precision: 4);
    }

    [Fact]
    public void CardioidDistance_EqualsQuarterMinusEpsilon()
    {
        // c = 1/4 − ε
        var f = BuildClaim();
        Assert.Equal(0.24, f.CardioidDistance(0.01), precision: 14);
        Assert.Equal(0.0, f.CardioidDistance(0.25), precision: 14);
    }

    [Fact]
    public void SixteenIsFourSquared()
    {
        Assert.True(BuildClaim().SixteenIsFourSquared());
    }

    [Fact]
    public void CardioidCuspMatchesQuarter()
    {
        Assert.True(BuildClaim().CardioidCuspMatchesQuarter());
    }

    [Fact]
    public void IterationCount_NonPositiveEpsilon_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().IterationCount(epsilon: 0.0, tol: 0.001));
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().IterationCount(epsilon: -0.01, tol: 0.001));
    }

    [Fact]
    public void IterationCount_NonPositiveTol_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().IterationCount(epsilon: 0.01, tol: 0.0));
    }

    [Fact]
    public void Alpha_NonPositiveTol_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().Alpha(tol: 0.0));
    }

    [Fact]
    public void CardioidDistance_NonPositiveEpsilon_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().CardioidDistance(0.0));
    }

    [Fact]
    public void Constructor_NullLadder_Throws()
    {
        Assert.Throws<ArgumentNullException>(() =>
            new F56CriticalSlowingPi2Inheritance(null!, new QuarterAsBilinearMaxvalClaim(), new HalfAsStructuralFixedPointClaim()));
    }

    [Fact]
    public void Constructor_NullQuarter_Throws()
    {
        Assert.Throws<ArgumentNullException>(() =>
            new F56CriticalSlowingPi2Inheritance(new Pi2DyadicLadderClaim(), null!, new HalfAsStructuralFixedPointClaim()));
    }

    [Fact]
    public void Constructor_NullHalf_Throws()
    {
        Assert.Throws<ArgumentNullException>(() =>
            new F56CriticalSlowingPi2Inheritance(new Pi2DyadicLadderClaim(), new QuarterAsBilinearMaxvalClaim(), null!));
    }
}
