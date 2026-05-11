using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class F65XxChainSpectrumPi2InheritanceRegistrationTests
{
    private static ClaimRegistryBuilder BuildBaseRegistry() =>
        new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterPi2DyadicLadder()
            .RegisterAbsorptionTheoremClaim()
            .RegisterF66PoleModesPi2Inheritance();

    [Fact]
    public void RegisterF65_AddsClaim()
    {
        var registry = BuildBaseRegistry()
            .RegisterF65XxChainSpectrumPi2Inheritance()
            .Build();

        Assert.True(registry.Contains<F65XxChainSpectrumPi2Inheritance>());
    }

    [Fact]
    public void RegisterF65_TierIsTier1Derived()
    {
        var registry = BuildBaseRegistry()
            .RegisterF65XxChainSpectrumPi2Inheritance()
            .Build();

        Assert.Equal(Tier.Tier1Derived,
            registry.Get<F65XxChainSpectrumPi2Inheritance>().Tier);
    }

    [Fact]
    public void RegisterF65_AncestorsContainF66()
    {
        var registry = BuildBaseRegistry()
            .RegisterF65XxChainSpectrumPi2Inheritance()
            .Build();

        var ancestors = registry.AncestorsOf<F65XxChainSpectrumPi2Inheritance>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(F66PoleModesPi2Inheritance), ancestors);
        Assert.Contains(typeof(Pi2DyadicLadderClaim), ancestors);
    }

    [Theory]
    [InlineData(3, 2, 1.0, 1.0)]    // peak at k=2 for N=3
    [InlineData(5, 3, 1.0, 2.0 / 3.0)]
    public void RegisterF65_SingleExcitationRateAcrossRegistry(int N, int k, double gammaZero, double expected)
    {
        var registry = BuildBaseRegistry()
            .RegisterF65XxChainSpectrumPi2Inheritance()
            .Build();

        Assert.Equal(expected, registry.Get<F65XxChainSpectrumPi2Inheritance>().SingleExcitationRate(N, k, gammaZero), precision: 8);
    }

    [Fact]
    public void RegisterF65_RatesLieInF66Interval()
    {
        var registry = BuildBaseRegistry()
            .RegisterF65XxChainSpectrumPi2Inheritance()
            .Build();

        Assert.True(registry.Get<F65XxChainSpectrumPi2Inheritance>().RatesLieInF66Interval(7, 0.05));
    }
}
