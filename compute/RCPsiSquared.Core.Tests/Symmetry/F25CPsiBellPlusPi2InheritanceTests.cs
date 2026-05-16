using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class F25CPsiBellPlusPi2InheritanceTests
{
    private static F25CPsiBellPlusPi2Inheritance BuildClaim() =>
        new F25CPsiBellPlusPi2Inheritance(new Pi2DyadicLadderClaim(), new QuarterAsBilinearMaxvalClaim());

    [Fact]
    public void Tier_IsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived, BuildClaim().Tier);
    }

    [Fact]
    public void DecayRateCoefficient_IsExactlyFour()
    {
        // 4 in e^{-4γt} = a_{-1} on dyadic ladder.
        Assert.Equal(4.0, BuildClaim().DecayRateCoefficient, precision: 14);
    }

    [Fact]
    public void Coefficient2_IsExactlyTwo()
    {
        Assert.Equal(2.0, BuildClaim().Coefficient2, precision: 14);
    }

    [Fact]
    public void CrossingThreshold_IsExactlyOneQuarter()
    {
        Assert.Equal(0.25, BuildClaim().CrossingThreshold, precision: 14);
    }

    [Fact]
    public void BellPlusFCross_Is08612()
    {
        Assert.Equal(0.8612, BuildClaim().BellPlusFCross, precision: 4);
    }

    [Fact]
    public void BellPlusKInvariant_Is00374()
    {
        Assert.Equal(0.0374, BuildClaim().BellPlusKInvariant, precision: 4);
    }

    [Fact]
    public void BellPlusF57Prefactor_Is1080088()
    {
        // F57's Bell+ prefactor 1.080088 = 2 / 1.851701; the "2" IS Coefficient2 = a_0.
        Assert.Equal(1.080088, BuildClaim().BellPlusF57Prefactor, precision: 6);
    }

    [Fact]
    public void CPsiAtTime_AtZero_IsOneThird()
    {
        // CΨ(0) = 1·(1+1)/6 = 1/3 for Bell+ initial state.
        Assert.Equal(1.0 / 3.0, BuildClaim().CPsiAtTime(0.05, 0.0), precision: 12);
    }

    [Fact]
    public void CPsiAtTime_LargeT_GoesToZero()
    {
        // f → 0 as t → ∞; CΨ → 0.
        Assert.Equal(0.0, BuildClaim().CPsiAtTime(0.05, 1000.0), precision: 6);
    }

    [Fact]
    public void CPsiAtTime_DecreasesMonotonically()
    {
        var f = BuildClaim();
        double cpsi0 = f.CPsiAtTime(0.05, 0.0);
        double cpsi1 = f.CPsiAtTime(0.05, 0.5);
        double cpsi2 = f.CPsiAtTime(0.05, 5.0);

        Assert.True(cpsi1 < cpsi0);
        Assert.True(cpsi2 < cpsi1);
    }

    [Fact]
    public void DCPsiDtAtTime_IsNegativeForPositiveT()
    {
        var f = BuildClaim();
        Assert.True(f.DCPsiDtAtTime(0.05, 0.1) < 0.0);
        Assert.True(f.DCPsiDtAtTime(0.05, 1.0) < 0.0);
    }

    [Fact]
    public void BellPlusInitialIsOneThird_HoldsExactly()
    {
        Assert.True(BuildClaim().BellPlusInitialIsOneThird());
    }

    [Fact]
    public void CrossingFConsistency_HoldsAtBellPlusFCross()
    {
        // At f = 0.8612, CΨ ≈ 0.25 (the fold)
        Assert.True(BuildClaim().CrossingFConsistency());
    }

    [Fact]
    public void CPsiAtTime_NegativeT_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().CPsiAtTime(0.05, -0.1));
    }

    [Fact]
    public void CPsiAtTime_NegativeGamma_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().CPsiAtTime(-0.05, 0.1));
    }

    [Fact]
    public void Constructor_RejectsNullParents()
    {
        var ladder = new Pi2DyadicLadderClaim();
        var quarter = new QuarterAsBilinearMaxvalClaim();
        Assert.Throws<ArgumentNullException>(() => new F25CPsiBellPlusPi2Inheritance(null!, quarter));
        Assert.Throws<ArgumentNullException>(() => new F25CPsiBellPlusPi2Inheritance(ladder, null!));
    }

    [Fact]
    public void TypedParents_AreExposed()
    {
        var f = BuildClaim();
        Assert.NotNull(f.Quarter);
    }
}
