using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class F73SpatialSumPurityClosurePi2InheritanceTests
{
    private static F73SpatialSumPurityClosurePi2Inheritance BuildClaim()
    {
        var ladder = new Pi2DyadicLadderClaim();
        var f70 = new F70DeltaNSelectionRulePi2Inheritance(ladder);
        var f72 = new F72BlockDiagonalPurityPi2Inheritance(ladder, f70);
        return new F73SpatialSumPurityClosurePi2Inheritance(ladder, f70, f72);
    }

    [Fact]
    public void Tier_IsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived, BuildClaim().Tier);
    }

    [Fact]
    public void ClosureBaseline_IsExactlyOneHalf()
    {
        Assert.Equal(0.5, BuildClaim().ClosureBaseline, precision: 14);
    }

    [Fact]
    public void DecayRateCoefficient_IsExactlyFour()
    {
        Assert.Equal(4.0, BuildClaim().DecayRateCoefficient, precision: 14);
    }

    [Fact]
    public void PerSitePrefactor_IsExactlyTwo()
    {
        Assert.Equal(2.0, BuildClaim().PerSitePrefactor, precision: 14);
    }

    [Fact]
    public void SpatialSumClosure_AtTZero_IsExactlyOneHalf()
    {
        Assert.Equal(0.5, BuildClaim().SpatialSumClosure(0.05, 0.0), precision: 14);
    }

    [Theory]
    [InlineData(0.05, 0.0, 0.5)]
    [InlineData(0.05, 20.0, 9.157819e-3)]   // verified value
    [InlineData(0.10, 5.0, 0.5 * 0.135335283)]  // 0.5·e^{-2}
    public void SpatialSumClosure_MatchesClosedForm(double gammaZero, double t, double expected)
    {
        Assert.Equal(expected, BuildClaim().SpatialSumClosure(gammaZero, t), precision: 6);
    }

    [Fact]
    public void SpatialSumClosure_LargeT_GoesToZero()
    {
        Assert.Equal(0.0, BuildClaim().SpatialSumClosure(0.05, 1000.0), precision: 6);
    }

    [Fact]
    public void ClosureAtTZeroIsHalf_HoldsExactly()
    {
        Assert.True(BuildClaim().ClosureAtTZeroIsHalf());
    }

    [Fact]
    public void DecayRateMatchesF25AndF76_HoldsExactly()
    {
        Assert.True(BuildClaim().DecayRateMatchesF25AndF76());
    }

    [Fact]
    public void VerifiedValueAtN5Gamma0p05T20_MatchesAnalyticalFormulas()
    {
        // (1/2) · exp(-4 · 0.05 · 20) = (1/2) · exp(-4) = 0.5 · 0.018316 = 0.009158
        Assert.Equal(9.157819e-3, BuildClaim().VerifiedValueAtN5Gamma0p05T20(), precision: 6);
    }

    [Fact]
    public void SpatialSumClosure_NegativeT_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().SpatialSumClosure(0.05, -0.1));
    }

    [Fact]
    public void SpatialSumClosure_NegativeGamma_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().SpatialSumClosure(-0.05, 0.1));
    }

    [Fact]
    public void Constructor_NullLadder_Throws()
    {
        var ladder = new Pi2DyadicLadderClaim();
        var f70 = new F70DeltaNSelectionRulePi2Inheritance(ladder);
        var f72 = new F72BlockDiagonalPurityPi2Inheritance(ladder, f70);
        Assert.Throws<ArgumentNullException>(() =>
            new F73SpatialSumPurityClosurePi2Inheritance(null!, f70, f72));
    }

    [Fact]
    public void Constructor_NullF70_Throws()
    {
        var ladder = new Pi2DyadicLadderClaim();
        var f70 = new F70DeltaNSelectionRulePi2Inheritance(ladder);
        var f72 = new F72BlockDiagonalPurityPi2Inheritance(ladder, f70);
        Assert.Throws<ArgumentNullException>(() =>
            new F73SpatialSumPurityClosurePi2Inheritance(ladder, null!, f72));
    }

    [Fact]
    public void Constructor_NullF72_Throws()
    {
        var ladder = new Pi2DyadicLadderClaim();
        var f70 = new F70DeltaNSelectionRulePi2Inheritance(ladder);
        Assert.Throws<ArgumentNullException>(() =>
            new F73SpatialSumPurityClosurePi2Inheritance(ladder, f70, null!));
    }
}
