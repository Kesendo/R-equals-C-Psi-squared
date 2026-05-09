using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class F26CPsiPauliChannelsPi2InheritanceRegistrationTests
{
    private static ClaimRegistryBuilder BuildBaseRegistry() =>
        new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterPi2DyadicLadder()
            .RegisterF25CPsiBellPlusPi2Inheritance();

    [Fact]
    public void RegisterF26_AddsClaim()
    {
        var registry = BuildBaseRegistry()
            .RegisterF26CPsiPauliChannelsPi2Inheritance()
            .Build();

        Assert.True(registry.Contains<F26CPsiPauliChannelsPi2Inheritance>());
    }

    [Fact]
    public void RegisterF26_TierIsTier1Derived()
    {
        var registry = BuildBaseRegistry()
            .RegisterF26CPsiPauliChannelsPi2Inheritance()
            .Build();

        Assert.Equal(Tier.Tier1Derived,
            registry.Get<F26CPsiPauliChannelsPi2Inheritance>().Tier);
    }

    [Fact]
    public void RegisterF26_AncestorsContainF25()
    {
        var registry = BuildBaseRegistry()
            .RegisterF26CPsiPauliChannelsPi2Inheritance()
            .Build();

        var ancestors = registry.AncestorsOf<F26CPsiPauliChannelsPi2Inheritance>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(F25CPsiBellPlusPi2Inheritance), ancestors);
    }

    [Fact]
    public void RegisterF26_AncestorsTransitivelyReachPolynomialFoundation()
    {
        var registry = BuildBaseRegistry()
            .RegisterF26CPsiPauliChannelsPi2Inheritance()
            .Build();

        var ancestors = registry.AncestorsOf<F26CPsiPauliChannelsPi2Inheritance>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(PolynomialFoundationClaim), ancestors);
    }

    [Fact]
    public void RegisterF26_DecayRateCoefficientIsFour()
    {
        var registry = BuildBaseRegistry()
            .RegisterF26CPsiPauliChannelsPi2Inheritance()
            .Build();

        Assert.Equal(4.0, registry.Get<F26CPsiPauliChannelsPi2Inheritance>().DecayRateCoefficient, precision: 14);
    }

    [Fact]
    public void RegisterF26_RecoversF25AtSingleChannelLimit()
    {
        var registry = BuildBaseRegistry()
            .RegisterF26CPsiPauliChannelsPi2Inheritance()
            .Build();

        Assert.True(registry.Get<F26CPsiPauliChannelsPi2Inheritance>().RecoversF25AtSingleChannelLimit(0.05, 0.5));
    }

    [Fact]
    public void RegisterF26_DepolarizingChannelMonotonic()
    {
        var registry = BuildBaseRegistry()
            .RegisterF26CPsiPauliChannelsPi2Inheritance()
            .Build();

        var f26 = registry.Get<F26CPsiPauliChannelsPi2Inheritance>();
        Assert.True(f26.DepolarizingCPsi(0.05, 1.0) < f26.DepolarizingCPsi(0.05, 0.0));
    }
}
