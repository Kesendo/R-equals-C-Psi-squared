using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class F57DwellTimeQuarterPi2InheritanceTests
{
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
    public void BellPlusKDwellPrefactor_IsClosedFormValue()
    {
        Assert.Equal(1.080088, BuildClaim().BellPlusKDwellPrefactor, precision: 6);
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
    public void EvenWeightPrefactor_AtBellPlusW2_RecoversBellPlusPrefactor()
    {
        // Bell+ has W₂ = 0.3709 per F58 → prefactor ≈ 1.080088
        var f = BuildClaim();
        double w2 = 0.3709;
        double prefactor = f.EvenWeightPrefactor(w2);
        Assert.Equal(1.080088, prefactor, precision: 4);
    }

    [Fact]
    public void TwoSectorPrefactor_AtBellPlus_RecoversBellPlusPrefactor()
    {
        // Bell+ as F59 special case: k=2, W₀=1/2, W₂=0.3709
        var f = BuildClaim();
        double prefactor = f.TwoSectorPrefactor(k: 2, w0: 0.5, wk: 0.3709);
        Assert.Equal(1.080088, prefactor, precision: 4);
    }

    [Fact]
    public void EvenWeightPrefactor_AndTwoSectorPrefactor_AgreeForBellPlus()
    {
        // F58 even-weight reading = F59 special case (k=2, W₀=1/2)
        var f = BuildClaim();
        double w2 = 0.3709;
        double f58 = f.EvenWeightPrefactor(w2);
        double f59 = f.TwoSectorPrefactor(k: 2, w0: 0.5, wk: w2);
        Assert.Equal(f58, f59, precision: 6);
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
