using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class F3DecayRateBoundsPi2InheritanceRegistrationTests
{
    private static ClaimRegistryBuilder BuildBaseRegistry() =>
        new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterPi2DyadicLadder()
            .RegisterF50WeightOneDegeneracyPi2Inheritance();

    [Fact]
    public void RegisterF3_AddsClaim()
    {
        var registry = BuildBaseRegistry()
            .RegisterF3DecayRateBoundsPi2Inheritance()
            .Build();

        Assert.True(registry.Contains<F3DecayRateBoundsPi2Inheritance>());
    }

    [Fact]
    public void RegisterF3_TierIsTier1Derived()
    {
        var registry = BuildBaseRegistry()
            .RegisterF3DecayRateBoundsPi2Inheritance()
            .Build();

        Assert.Equal(Tier.Tier1Derived,
            registry.Get<F3DecayRateBoundsPi2Inheritance>().Tier);
    }

    [Fact]
    public void RegisterF3_AncestorsContainF50AndPi2DyadicLadder()
    {
        var registry = BuildBaseRegistry()
            .RegisterF3DecayRateBoundsPi2Inheritance()
            .Build();

        var ancestors = registry.AncestorsOf<F3DecayRateBoundsPi2Inheritance>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(F50WeightOneDegeneracyPi2Inheritance), ancestors);
        Assert.Contains(typeof(Pi2DyadicLadderClaim), ancestors);
    }

    [Theory]
    [InlineData(5, 0.05, 0.4)]   // max = 2·4·0.05 = 0.4
    [InlineData(7, 1.0, 12.0)]   // max = 2·6·1 = 12
    public void RegisterF3_MaxRateAcrossRegistry(int N, double gammaZero, double expected)
    {
        var registry = BuildBaseRegistry()
            .RegisterF3DecayRateBoundsPi2Inheritance()
            .Build();

        Assert.Equal(expected, registry.Get<F3DecayRateBoundsPi2Inheritance>().MaxRate(N, gammaZero), precision: 14);
    }

    [Fact]
    public void RegisterF3_MinRateMatchesF50AcrossRegistry()
    {
        var registry = BuildBaseRegistry()
            .RegisterF3DecayRateBoundsPi2Inheritance()
            .Build();

        Assert.True(registry.Get<F3DecayRateBoundsPi2Inheritance>().MinRateMatchesF50(0.05));
    }

    [Fact]
    public void RegisterF3_BandwidthIsMaxMinusMinAcrossRegistry()
    {
        var registry = BuildBaseRegistry()
            .RegisterF3DecayRateBoundsPi2Inheritance()
            .Build();

        var f3 = registry.Get<F3DecayRateBoundsPi2Inheritance>();
        for (int N = 2; N <= 10; N++)
            Assert.True(f3.BandwidthIsMaxMinusMin(N, 0.05));
    }

    [Fact]
    public void RegisterF3_WithoutF50_Throws()
    {
        Assert.Throws<InvariantViolationException>(() =>
            new ClaimRegistryBuilder()
                .RegisterPi2Family()
                .RegisterPi2DyadicLadder()
                // Missing: RegisterF50WeightOneDegeneracyPi2Inheritance
                .RegisterF3DecayRateBoundsPi2Inheritance()
                .Build());
    }
}
