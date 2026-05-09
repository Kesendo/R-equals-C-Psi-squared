using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class F57DwellTimeQuarterPi2InheritanceRegistrationTests
{
    private static ClaimRegistryBuilder BuildBaseRegistry() =>
        new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterPi2DyadicLadder()
            .RegisterF25CPsiBellPlusPi2Inheritance();

    [Fact]
    public void RegisterF57_AddsClaim()
    {
        var registry = BuildBaseRegistry()
            .RegisterF57DwellTimeQuarterPi2Inheritance()
            .Build();

        Assert.True(registry.Contains<F57DwellTimeQuarterPi2Inheritance>());
    }

    [Fact]
    public void RegisterF57_TierIsTier1Derived()
    {
        var registry = BuildBaseRegistry()
            .RegisterF57DwellTimeQuarterPi2Inheritance()
            .Build();

        Assert.Equal(Tier.Tier1Derived,
            registry.Get<F57DwellTimeQuarterPi2Inheritance>().Tier);
    }

    [Fact]
    public void RegisterF57_AncestorsContainBothPi2Anchors()
    {
        var registry = BuildBaseRegistry()
            .RegisterF57DwellTimeQuarterPi2Inheritance()
            .Build();

        var ancestors = registry.AncestorsOf<F57DwellTimeQuarterPi2Inheritance>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(Pi2DyadicLadderClaim), ancestors);
        Assert.Contains(typeof(QuarterAsBilinearMaxvalClaim), ancestors);
    }

    [Fact]
    public void RegisterF57_CrossingThresholdIsExactlyOneQuarter()
    {
        var registry = BuildBaseRegistry()
            .RegisterF57DwellTimeQuarterPi2Inheritance()
            .Build();

        Assert.Equal(0.25, registry.Get<F57DwellTimeQuarterPi2Inheritance>().CrossingThreshold, precision: 14);
    }

    [Fact]
    public void RegisterF57_WindowDoublingFactorIsExactlyTwo()
    {
        var registry = BuildBaseRegistry()
            .RegisterF57DwellTimeQuarterPi2Inheritance()
            .Build();

        Assert.Equal(2.0, registry.Get<F57DwellTimeQuarterPi2Inheritance>().WindowDoublingFactor, precision: 14);
    }

    [Fact]
    public void RegisterF57_ThresholdMatchesQuarterAnchor()
    {
        var registry = BuildBaseRegistry()
            .RegisterF57DwellTimeQuarterPi2Inheritance()
            .Build();

        Assert.True(registry.Get<F57DwellTimeQuarterPi2Inheritance>().ThresholdMatchesQuarterAnchor());
    }

    [Theory]
    [InlineData(0.01, 0.05, 1.080088, 0.21601760)]
    [InlineData(0.05, 1.0, 1.080088, 0.0540044)]
    public void RegisterF57_TDwellAcrossRegistry(double delta, double gamma, double prefactor, double expected)
    {
        var registry = BuildBaseRegistry()
            .RegisterF57DwellTimeQuarterPi2Inheritance()
            .Build();

        Assert.Equal(expected, registry.Get<F57DwellTimeQuarterPi2Inheritance>().TDwell(delta, gamma, prefactor), precision: 8);
    }

    [Fact]
    public void RegisterF57_BellPlusF58F59AgreeAcrossRegistry()
    {
        var registry = BuildBaseRegistry()
            .RegisterF57DwellTimeQuarterPi2Inheritance()
            .Build();

        var f = registry.Get<F57DwellTimeQuarterPi2Inheritance>();
        double w2 = 0.3709;
        Assert.Equal(f.EvenWeightPrefactor(w2), f.TwoSectorPrefactor(k: 2, w0: 0.5, wk: w2), precision: 6);
    }

    [Fact]
    public void RegisterF57_AncestorsTransitivelyReachPolynomialFoundation()
    {
        // CrossingThreshold (1/4) and WindowDoublingFactor (2) both trace back to
        // d²−2d=0: a_3 = (1/d)² and a_0 = d via the dyadic ladder.
        var registry = BuildBaseRegistry()
            .RegisterF57DwellTimeQuarterPi2Inheritance()
            .Build();

        var ancestors = registry.AncestorsOf<F57DwellTimeQuarterPi2Inheritance>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(PolynomialFoundationClaim), ancestors);
        Assert.Contains(typeof(QubitDimensionalAnchorClaim), ancestors);
    }
}
