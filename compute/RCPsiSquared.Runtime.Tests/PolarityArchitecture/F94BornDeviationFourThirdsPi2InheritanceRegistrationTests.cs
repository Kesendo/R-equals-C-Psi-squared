using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class F94BornDeviationFourThirdsPi2InheritanceRegistrationTests
{
    private static ClaimRegistryBuilder BuildBaseRegistry() =>
        new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterPi2DyadicLadder();

    [Fact]
    public void RegisterF94_AddsClaim()
    {
        var registry = BuildBaseRegistry()
            .RegisterF94BornDeviationFourThirdsPi2Inheritance()
            .Build();

        Assert.True(registry.Contains<F94BornDeviationFourThirdsPi2Inheritance>());
    }

    [Fact]
    public void RegisterF94_TierIsTier1Derived()
    {
        var registry = BuildBaseRegistry()
            .RegisterF94BornDeviationFourThirdsPi2Inheritance()
            .Build();

        Assert.Equal(Tier.Tier1Derived,
            registry.Get<F94BornDeviationFourThirdsPi2Inheritance>().Tier);
    }

    [Fact]
    public void RegisterF94_AncestorsContainPi2DyadicLadder()
    {
        // The "4" in 4/3 is a_{-1} = 4 on the Pi2 dyadic ladder — typed parent.
        var registry = BuildBaseRegistry()
            .RegisterF94BornDeviationFourThirdsPi2Inheritance()
            .Build();

        var ancestors = registry.AncestorsOf<F94BornDeviationFourThirdsPi2Inheritance>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(Pi2DyadicLadderClaim), ancestors);
    }

    [Fact]
    public void RegisterF94_CoefficientIsFourThirdsAcrossRegistry()
    {
        var registry = BuildBaseRegistry()
            .RegisterF94BornDeviationFourThirdsPi2Inheritance()
            .Build();

        var f = registry.Get<F94BornDeviationFourThirdsPi2Inheritance>();
        Assert.Equal(4.0 / 3.0, f.Coefficient, precision: 15);
    }

    [Fact]
    public void RegisterF94_DeltaDominantSampleAcrossRegistry()
    {
        var registry = BuildBaseRegistry()
            .RegisterF94BornDeviationFourThirdsPi2Inheritance()
            .Build();

        var f = registry.Get<F94BornDeviationFourThirdsPi2Inheritance>();
        // Carrier-sweep anchor point (Q=20, K=0.0143):
        double Q = 20.0, K = 0.0143;
        double expected = (4.0 / 3.0) * Q * Q * K * K * K;
        Assert.Equal(expected, f.DeltaDominant(Q, K), precision: 14);
    }
}
