using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class F25CPsiBellPlusPi2InheritanceTests
{
    /// <summary>Double-precision machine epsilon: the unit the exactness gates
    /// below are stated in, since the only deviation left in F25's Bell+ chain
    /// is float rounding of ∛ and ln.</summary>
    private const double MachineEps = 2.220446049250313e-16;

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
    public void BellPlusFCross_SolvesTheCubic()
    {
        // The claim is not "f* is 0.8612" but "f* solves f(1+f²) = 3/2". Gate the
        // equation, not the tabulation: the residual is float rounding of ∛ only.
        double f = BuildClaim().BellPlusFCross;
        Assert.True(Math.Abs(f * (1.0 + f * f) - 1.5) < 8.0 * MachineEps,
            $"cubic residual {Math.Abs(f * (1.0 + f * f) - 1.5):E3} exceeds 8 eps");
        Assert.Equal(0.8612240997395737, f, precision: 14);
    }

    [Fact]
    public void BellPlusFCross_IsTheOnlyRealRootInTheUnitInterval()
    {
        // f³ + f − 3/2 is strictly increasing, so the root is unique and the
        // solver cannot have landed on a sibling: it brackets and nothing else does.
        double f = BuildClaim().BellPlusFCross;
        Assert.True(f > 0.0 && f < 1.0);
        Assert.True((f - 1e-6) * (1.0 + (f - 1e-6) * (f - 1e-6)) < 1.5);
        Assert.True((f + 1e-6) * (1.0 + (f + 1e-6) * (f + 1e-6)) > 1.5);
    }

    [Fact]
    public void BellPlusKInvariant_IsMinusLogFStarOverFour()
    {
        // K = γ·t_cross from f* = e^{−4γt}: 0.037350132494447214, which is the
        // 0.03735 ANALYTICAL_FORMULAS calls unrounded and prints as 0.0374.
        Assert.Equal(0.037350132494447214, BuildClaim().BellPlusKInvariant, precision: 14);
    }

    [Fact]
    public void BellPlusF57Prefactor_IsTwoOverTheCrossingDerivative()
    {
        // F57's Bell+ prefactor = 2 / |dCΨ/dt|_{t_cross}/γ; the "2" IS Coefficient2 = a_0.
        var f = BuildClaim();
        Assert.Equal(1.851701200347236, f.BellPlusDCPsiDtMagnitudePerGamma, precision: 14);
        Assert.Equal(1.0800878671056402, f.BellPlusF57Prefactor, precision: 14);
        Assert.Equal(f.Coefficient2 / f.BellPlusDCPsiDtMagnitudePerGamma, f.BellPlusF57Prefactor);
    }

    [Fact]
    public void BellPlusDCPsiDtMagnitude_IsTheLiveDerivativeAtTheCrossingTime()
    {
        // The per-γ magnitude is not a second formula: it is DCPsiDtAtTime read at
        // t = K/γ, divided by γ. Checked at three γ, exactly the same double each time.
        var f = BuildClaim();
        foreach (double gamma in new[] { 0.01, 0.05, 10.0 })
        {
            double live = Math.Abs(f.DCPsiDtAtTime(gamma, f.BellPlusKInvariant / gamma)) / gamma;
            Assert.Equal(f.BellPlusDCPsiDtMagnitudePerGamma, live, precision: 14);
        }
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
        // CΨ(f*) = 1/4, the fold.
        Assert.True(BuildClaim().CrossingFConsistency());
    }

    [Fact]
    public void CrossingFConsistency_FailsOnTheFourDigitTabulation()
    {
        // The gate has to be able to fail: at 8 eps the old 4-decimal f* = 0.8612
        // misses the fold by 3.7e-6, sixteen thousand million eps.
        double tabulated = 0.8612;
        double cpsi = tabulated * (1.0 + tabulated * tabulated) / 6.0;
        Assert.True(Math.Abs(cpsi - 0.25) > 8.0 * MachineEps);
    }

    [Theory]
    [InlineData(0.01)]
    [InlineData(0.05)]
    [InlineData(100.0)]
    public void CrossingTimeConsistency_HoldsAtEveryGamma(double gamma)
    {
        // K is the γ-invariant, so t = K/γ lands on the fold whatever γ is.
        Assert.True(BuildClaim().CrossingTimeConsistency(gamma));
    }

    [Fact]
    public void CrossingTimeConsistency_RejectsNonPositiveGamma()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().CrossingTimeConsistency(0.0));
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().CrossingTimeConsistency(-1.0));
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
