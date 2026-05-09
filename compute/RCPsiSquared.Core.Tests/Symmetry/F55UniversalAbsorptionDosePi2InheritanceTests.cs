using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class F55UniversalAbsorptionDosePi2InheritanceTests
{
    private static F55UniversalAbsorptionDosePi2Inheritance BuildClaim()
    {
        var ladder = new Pi2DyadicLadderClaim();
        var f50 = new F50WeightOneDegeneracyPi2Inheritance(ladder);
        return new F55UniversalAbsorptionDosePi2Inheritance(ladder, f50);
    }

    [Fact]
    public void Tier_IsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived, BuildClaim().Tier);
    }

    [Fact]
    public void RateMinCoefficient_IsExactlyTwo()
    {
        Assert.Equal(2.0, BuildClaim().RateMinCoefficient, precision: 14);
    }

    [Fact]
    public void KDeath_EqualsLn10()
    {
        Assert.Equal(Math.Log(10.0), F55UniversalAbsorptionDosePi2Inheritance.KDeath, precision: 14);
    }

    [Fact]
    public void KDeath_ApproximatelyTwoPoint303()
    {
        Assert.Equal(2.303, F55UniversalAbsorptionDosePi2Inheritance.KDeath, precision: 3);
    }

    [Fact]
    public void AbsorptionFraction_IsExactlyZeroPoint99()
    {
        Assert.Equal(0.99, F55UniversalAbsorptionDosePi2Inheritance.AbsorptionFraction, precision: 14);
    }

    [Fact]
    public void KDeathOverKFold_IsApproximately2Point3()
    {
        Assert.Equal(2.3, F55UniversalAbsorptionDosePi2Inheritance.KDeathOverKFoldApprox, precision: 1);
    }

    [Theory]
    [InlineData(1, 2)]
    [InlineData(3, 4)]
    [InlineData(5, 6)]
    [InlineData(20, 21)]
    public void ImmortalModeCount_EqualsNPlusOne(int N, int expected)
    {
        Assert.Equal(expected, BuildClaim().ImmortalModeCount(N));
    }

    [Theory]
    [InlineData(1.0, 2.30258509299)]   // ln(10) at γ=1
    [InlineData(0.05, 46.05170186)]    // ln(10)/0.05
    public void AbsorptionTime_EqualsLn10OverGamma(double gammaZero, double expected)
    {
        Assert.Equal(expected, BuildClaim().AbsorptionTime(gammaZero), precision: 6);
    }

    [Fact]
    public void DerivationConsistencyHolds_AcrossAllGamma()
    {
        var f = BuildClaim();
        // 2γ · t_death = ln(100) for any γ.
        Assert.True(f.DerivationConsistencyHolds(0.05));
        Assert.True(f.DerivationConsistencyHolds(0.5));
        Assert.True(f.DerivationConsistencyHolds(1.0));
        Assert.True(f.DerivationConsistencyHolds(0.001));
    }

    [Fact]
    public void RateMinMatchesF50_AcrossAllGamma()
    {
        // F55 rate_min = F50 |EigenvaluePosition| at any γ.
        var f = BuildClaim();
        Assert.True(f.RateMinMatchesF50(0.05));
        Assert.True(f.RateMinMatchesF50(0.5));
        Assert.True(f.RateMinMatchesF50(1.0));
    }

    [Fact]
    public void AbsorptionTime_NonPositiveGamma_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().AbsorptionTime(0.0));
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().AbsorptionTime(-0.1));
    }

    [Fact]
    public void ImmortalModeCount_NLessThanOne_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().ImmortalModeCount(N: 0));
    }

    [Fact]
    public void Constructor_NullLadder_Throws()
    {
        var ladder = new Pi2DyadicLadderClaim();
        var f50 = new F50WeightOneDegeneracyPi2Inheritance(ladder);
        Assert.Throws<ArgumentNullException>(() =>
            new F55UniversalAbsorptionDosePi2Inheritance(null!, f50));
    }

    [Fact]
    public void Constructor_NullF50_Throws()
    {
        Assert.Throws<ArgumentNullException>(() =>
            new F55UniversalAbsorptionDosePi2Inheritance(new Pi2DyadicLadderClaim(), null!));
    }
}
