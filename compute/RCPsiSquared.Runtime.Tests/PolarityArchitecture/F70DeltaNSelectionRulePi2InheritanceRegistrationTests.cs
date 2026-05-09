using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class F70DeltaNSelectionRulePi2InheritanceRegistrationTests
{
    private static ClaimRegistryBuilder BuildBaseRegistry() =>
        new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterPi2DyadicLadder();

    [Fact]
    public void RegisterF70_AddsClaim()
    {
        var registry = BuildBaseRegistry()
            .RegisterF70DeltaNSelectionRulePi2Inheritance()
            .Build();

        Assert.True(registry.Contains<F70DeltaNSelectionRulePi2Inheritance>());
    }

    [Fact]
    public void RegisterF70_TierIsTier1Derived()
    {
        var registry = BuildBaseRegistry()
            .RegisterF70DeltaNSelectionRulePi2Inheritance()
            .Build();

        Assert.Equal(Tier.Tier1Derived,
            registry.Get<F70DeltaNSelectionRulePi2Inheritance>().Tier);
    }

    [Fact]
    public void RegisterF70_AncestorsContainLadder()
    {
        var registry = BuildBaseRegistry()
            .RegisterF70DeltaNSelectionRulePi2Inheritance()
            .Build();

        var ancestors = registry.AncestorsOf<F70DeltaNSelectionRulePi2Inheritance>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(Pi2DyadicLadderClaim), ancestors);
    }

    [Fact]
    public void RegisterF70_AncestorsTransitivelyReachPolynomialFoundation()
    {
        var registry = BuildBaseRegistry()
            .RegisterF70DeltaNSelectionRulePi2Inheritance()
            .Build();

        var ancestors = registry.AncestorsOf<F70DeltaNSelectionRulePi2Inheritance>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(PolynomialFoundationClaim), ancestors);
    }

    [Fact]
    public void RegisterF70_SingleSiteMaxDeltaNIsOne()
    {
        var registry = BuildBaseRegistry()
            .RegisterF70DeltaNSelectionRulePi2Inheritance()
            .Build();

        Assert.Equal(1.0, registry.Get<F70DeltaNSelectionRulePi2Inheritance>().SingleSiteMaxDeltaN, precision: 14);
    }

    [Fact]
    public void RegisterF70_PairMaxDeltaNIsTwo()
    {
        var registry = BuildBaseRegistry()
            .RegisterF70DeltaNSelectionRulePi2Inheritance()
            .Build();

        Assert.Equal(2.0, registry.Get<F70DeltaNSelectionRulePi2Inheritance>().PairMaxDeltaN, precision: 14);
    }

    [Theory]
    [InlineData(1)]
    [InlineData(2)]
    [InlineData(3)]
    public void RegisterF70_PartialTraceMaxDeltaNAcrossRegistry(int k)
    {
        var registry = BuildBaseRegistry()
            .RegisterF70DeltaNSelectionRulePi2Inheritance()
            .Build();

        Assert.Equal(k, registry.Get<F70DeltaNSelectionRulePi2Inheritance>().PartialTraceMaxDeltaN(k));
    }

    [Fact]
    public void RegisterF70_BothLadderAnchorsHoldAcrossRegistry()
    {
        var registry = BuildBaseRegistry()
            .RegisterF70DeltaNSelectionRulePi2Inheritance()
            .Build();

        var f70 = registry.Get<F70DeltaNSelectionRulePi2Inheritance>();
        Assert.True(f70.SingleSiteThresholdMatchesSelfMirrorPivot());
        Assert.True(f70.PairThresholdMatchesPolynomialRoot());
    }
}
