using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class F57DwellTimeQuarterPi2InheritanceTests
{
    /// <summary>Double-precision machine epsilon: the unit F57's exactness gates
    /// are stated in, the F58 ↔ F57 identity being algebraic at the crossing.</summary>
    private const double MachineEps = 2.220446049250313e-16;

    private static F57DwellTimeQuarterPi2Inheritance BuildClaim()
    {
        var ladder = new Pi2DyadicLadderClaim();
        var quarter = new QuarterAsBilinearMaxvalClaim();
        var f25 = new F25CPsiBellPlusPi2Inheritance(ladder, quarter);
        return new F57DwellTimeQuarterPi2Inheritance(ladder, quarter, f25, new ArgmaxMaxvalPairClaim());
    }

    [Fact]
    public void Tier_IsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived, BuildClaim().Tier);
    }

    [Fact]
    public void CrossingThreshold_IsExactlyOneQuarter()
    {
        Assert.Equal(0.25, BuildClaim().CrossingThreshold, precision: 14);
    }

    [Fact]
    public void WindowDoublingFactor_IsExactlyTwo()
    {
        Assert.Equal(2.0, BuildClaim().WindowDoublingFactor, precision: 14);
    }

    [Fact]
    public void ThresholdMatchesQuarterAnchor_HoldsExactly()
    {
        Assert.True(BuildClaim().ThresholdMatchesQuarterAnchor());
    }

    [Fact]
    public void BellPlusKDwellPrefactor_IsTheMotherClaimsValue()
    {
        // Not re-tabulated here: F57 reads F25's 2 / |dCΨ/dt|_{t_cross}/γ.
        var f = BuildClaim();
        Assert.Equal(f.F25.BellPlusF57Prefactor, f.BellPlusKDwellPrefactor);
        Assert.Equal(1.0800878671056402, f.BellPlusKDwellPrefactor, precision: 14);
    }

    [Fact]
    public void BellPlusW2AtCrossing_IsHalfFStarSquared()
    {
        // W₂ = f*²/2 from ρ = ¼(II + ZZ + f·XX − f·YY) with weights Σ c_P²/4.
        var f = BuildClaim();
        double fStar = f.F25.BellPlusFCross;
        Assert.Equal(fStar * fStar / 2.0, f.BellPlusW2AtCrossing);
        Assert.Equal(0.3708534749861196, f.BellPlusW2AtCrossing, precision: 14);
    }

    [Fact]
    public void BellPlusW0_IsOneHalf()
    {
        // The frozen diagonal: (1² + 1²)/4, and Z-dephasing never touches it.
        Assert.Equal(0.5, BuildClaim().BellPlusW0);
    }

    [Theory]
    [InlineData(0.01, 1.080088, 0.01080088)]
    [InlineData(0.05, 1.080088, 0.0540044)]
    [InlineData(0.1, 1.0, 0.1)]
    [InlineData(0.0, 1.5, 0.0)]
    public void KDwell_IsPrefactorTimesDelta(double delta, double prefactor, double expected)
    {
        Assert.Equal(expected, BuildClaim().KDwell(delta, prefactor), precision: 8);
    }

    [Theory]
    [InlineData(0.01, 0.05, 1.080088, 0.21601760)]   // 1.080088·0.01/0.05 = 0.2160176
    [InlineData(0.05, 1.0, 1.080088, 0.0540044)]
    [InlineData(0.1, 0.5, 1.0, 0.2)]
    public void TDwell_IsKDwellOverGamma(double delta, double gamma, double prefactor, double expected)
    {
        Assert.Equal(expected, BuildClaim().TDwell(delta, gamma, prefactor), precision: 8);
    }

    [Fact]
    public void TDwell_IsGammaInvariantForKDwell()
    {
        // K_dwell = γ · t_dwell is γ-independent: scale γ by ×10, t_dwell scales by /10,
        // K_dwell stays put.
        var f = BuildClaim();
        double delta = 0.05;
        double prefactor = 1.080088;
        double k1 = 0.1 * f.TDwell(delta, 0.1, prefactor);
        double k10 = 10.0 * f.TDwell(delta, 10.0, prefactor);
        Assert.Equal(k1, k10, precision: 14);
        Assert.Equal(f.KDwell(delta, prefactor), k1, precision: 14);
    }

    [Fact]
    public void EvenWeightPrefactor_AtW2_0_IsTwo()
    {
        // W₂ = 0 (no light-face content) → prefactor = 2/1 = 2 (degenerate case)
        Assert.Equal(2.0, BuildClaim().EvenWeightPrefactor(0.0), precision: 14);
    }

    [Fact]
    public void EvenWeightPrefactor_AtBellPlusW2_IsTheSameNumberNotANearOne()
    {
        // F58 at W₂ = f*²/2 IS F57's prefactor: substituting turns (2+4W₂)/(1+6W₂)
        // into (2+2f²)/(1+3f²), equal to 3/(f(1+3f²)) exactly when f(1+f²) = 3/2.
        // So this gate tests the F25 crossing equation; only rounding is left.
        var f = BuildClaim();
        Assert.True(f.EvenWeightPrefactorReducesToF57());
        Assert.Equal(f.BellPlusKDwellPrefactor, f.EvenWeightPrefactor(f.BellPlusW2AtCrossing), precision: 14);
    }

    [Fact]
    public void EvenWeightPrefactor_AtTheRoundedW2_IsTenOrdersWorse()
    {
        // The gate above can fail, and this is what makes it fail: the four-digit 0.3709
        // round-trip costs 3.578e-5 against 4.441e-16 for the exact weight, a factor 8.06e10.
        // The 1e9 asked for below leaves an order of magnitude of headroom under that, so a
        // change of libm does not turn the control red while a real regression still does.
        var f = BuildClaim();
        double exact = Math.Abs(f.EvenWeightPrefactor(f.BellPlusW2AtCrossing) - f.BellPlusKDwellPrefactor);
        double rounded = Math.Abs(f.EvenWeightPrefactor(0.3709) - f.BellPlusKDwellPrefactor);
        Assert.True(exact < 8.0 * MachineEps);
        Assert.True(rounded > 1e-5);
        Assert.True(rounded / exact > 1e9);
    }

    [Fact]
    public void TwoSectorPrefactor_AtBellPlus_IsTheSameIdentity()
    {
        // Bell+ as F59 special case: k = 2, W₀ = 1/2, W_k = f*²/2.
        var f = BuildClaim();
        double prefactor = f.TwoSectorPrefactor(k: 2, w0: f.BellPlusW0, wk: f.BellPlusW2AtCrossing);
        Assert.Equal(f.BellPlusKDwellPrefactor, prefactor, precision: 14);
    }

    [Fact]
    public void EvenWeightPrefactor_AndTwoSectorPrefactor_AgreeForBellPlus()
    {
        // F58's even-weight reading and F59 at k = 2, W₀ = 1/2 are the same
        // expression term for term, so they agree bit for bit, not to six digits.
        var f = BuildClaim();
        double w2 = f.BellPlusW2AtCrossing;
        Assert.Equal(f.EvenWeightPrefactor(w2), f.TwoSectorPrefactor(k: 2, w0: f.BellPlusW0, wk: w2));
    }

    [Fact]
    public void KDwell_NegativeDelta_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().KDwell(delta: -0.01, prefactor: 1.0));
    }

    [Fact]
    public void TDwell_ZeroGamma_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().TDwell(delta: 0.05, gamma: 0.0, prefactor: 1.0));
    }

    [Fact]
    public void TDwell_NegativeGamma_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().TDwell(delta: 0.05, gamma: -0.1, prefactor: 1.0));
    }

    [Fact]
    public void EvenWeightPrefactor_OutOfRangeW2_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().EvenWeightPrefactor(w2: -0.1));
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().EvenWeightPrefactor(w2: 1.1));
    }

    [Fact]
    public void TwoSectorPrefactor_KZero_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().TwoSectorPrefactor(k: 0, w0: 0.5, wk: 0.5));
    }

    [Fact]
    public void TwoSectorPrefactor_NegativeWeights_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().TwoSectorPrefactor(k: 2, w0: -0.1, wk: 0.5));
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().TwoSectorPrefactor(k: 2, w0: 0.5, wk: -0.1));
    }

    [Fact]
    public void TwoSectorPrefactor_DegenerateDenominator_Throws()
    {
        // W₀ + 3·W_k = 0 only when both are 0
        Assert.Throws<ArgumentException>(() => BuildClaim().TwoSectorPrefactor(k: 2, w0: 0.0, wk: 0.0));
    }

    [Fact]
    public void Constructor_RejectsNullParents()
    {
        var ladder = new Pi2DyadicLadderClaim();
        var quarter = new QuarterAsBilinearMaxvalClaim();
        var f25 = new F25CPsiBellPlusPi2Inheritance(ladder, quarter);
        var argmaxMaxval = new ArgmaxMaxvalPairClaim();
        Assert.Throws<ArgumentNullException>(() =>
            new F57DwellTimeQuarterPi2Inheritance(null!, quarter, f25, argmaxMaxval));
        Assert.Throws<ArgumentNullException>(() =>
            new F57DwellTimeQuarterPi2Inheritance(ladder, null!, f25, argmaxMaxval));
        Assert.Throws<ArgumentNullException>(() =>
            new F57DwellTimeQuarterPi2Inheritance(ladder, quarter, null!, argmaxMaxval));
        Assert.Throws<ArgumentNullException>(() =>
            new F57DwellTimeQuarterPi2Inheritance(ladder, quarter, f25, null!));
    }

    [Fact]
    public void TypedParents_AreExposed()
    {
        var f = BuildClaim();
        Assert.NotNull(f.F25);
        Assert.NotNull(f.ArgmaxMaxval);
    }
}
