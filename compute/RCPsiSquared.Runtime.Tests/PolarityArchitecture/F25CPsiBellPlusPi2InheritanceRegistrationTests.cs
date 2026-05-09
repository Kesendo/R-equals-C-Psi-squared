using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class F25CPsiBellPlusPi2InheritanceRegistrationTests
{
    private static ClaimRegistryBuilder BuildBaseRegistry() =>
        new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterPi2DyadicLadder();

    [Fact]
    public void RegisterF25_AddsClaim()
    {
        var registry = BuildBaseRegistry()
            .RegisterF25CPsiBellPlusPi2Inheritance()
            .Build();

        Assert.True(registry.Contains<F25CPsiBellPlusPi2Inheritance>());
    }

    [Fact]
    public void RegisterF25_TierIsTier1Derived()
    {
        var registry = BuildBaseRegistry()
            .RegisterF25CPsiBellPlusPi2Inheritance()
            .Build();

        Assert.Equal(Tier.Tier1Derived,
            registry.Get<F25CPsiBellPlusPi2Inheritance>().Tier);
    }

    [Fact]
    public void RegisterF25_AncestorsContainLadder()
    {
        var registry = BuildBaseRegistry()
            .RegisterF25CPsiBellPlusPi2Inheritance()
            .Build();

        var ancestors = registry.AncestorsOf<F25CPsiBellPlusPi2Inheritance>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(Pi2DyadicLadderClaim), ancestors);
    }

    [Fact]
    public void RegisterF25_AncestorsTransitivelyReachPolynomialFoundation()
    {
        var registry = BuildBaseRegistry()
            .RegisterF25CPsiBellPlusPi2Inheritance()
            .Build();

        var ancestors = registry.AncestorsOf<F25CPsiBellPlusPi2Inheritance>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(PolynomialFoundationClaim), ancestors);
    }

    [Fact]
    public void RegisterF25_DecayRateCoefficientIsFour()
    {
        var registry = BuildBaseRegistry()
            .RegisterF25CPsiBellPlusPi2Inheritance()
            .Build();

        Assert.Equal(4.0, registry.Get<F25CPsiBellPlusPi2Inheritance>().DecayRateCoefficient, precision: 14);
    }

    [Fact]
    public void RegisterF25_CPsiAtZeroIsOneThird()
    {
        var registry = BuildBaseRegistry()
            .RegisterF25CPsiBellPlusPi2Inheritance()
            .Build();

        Assert.Equal(1.0 / 3.0, registry.Get<F25CPsiBellPlusPi2Inheritance>().CPsiAtTime(0.05, 0.0), precision: 12);
    }

    [Fact]
    public void RegisterF25_BellPlusF57PrefactorMatches()
    {
        var registry = BuildBaseRegistry()
            .RegisterF25CPsiBellPlusPi2Inheritance()
            .Build();

        Assert.Equal(1.080088, registry.Get<F25CPsiBellPlusPi2Inheritance>().BellPlusF57Prefactor, precision: 6);
    }

    [Fact]
    public void RegisterF25_F57DescendantInheritsF25AsAncestor()
    {
        // F57 retro-fit verified: F25 must be ancestor of F57 after both registered.
        var registry = BuildBaseRegistry()
            .RegisterF25CPsiBellPlusPi2Inheritance()
            .RegisterF57DwellTimeQuarterPi2Inheritance()
            .Build();

        var f57Ancestors = registry.AncestorsOf<F57DwellTimeQuarterPi2Inheritance>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(F25CPsiBellPlusPi2Inheritance), f57Ancestors);
    }
}
