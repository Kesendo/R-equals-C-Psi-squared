using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Spectrum;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;
using RCPsiSquared.Runtime.Spectrum;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class F2W1DispersionPi2InheritanceRegistrationTests
{
    private static ClaimRegistryBuilder BuildBaseRegistry() =>
        new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterPi2DyadicLadder()
            .RegisterW1Dispersion(N: 5, J: 1.0, gammaZero: 0.05);

    [Fact]
    public void RegisterF2_AddsClaim()
    {
        var registry = BuildBaseRegistry()
            .RegisterF2W1DispersionPi2Inheritance()
            .Build();

        Assert.True(registry.Contains<F2W1DispersionPi2Inheritance>());
    }

    [Fact]
    public void RegisterF2_TierIsTier1Derived()
    {
        var registry = BuildBaseRegistry()
            .RegisterF2W1DispersionPi2Inheritance()
            .Build();

        Assert.Equal(Tier.Tier1Derived,
            registry.Get<F2W1DispersionPi2Inheritance>().Tier);
    }

    [Fact]
    public void RegisterF2_AncestorsContainW1AndPi2DyadicLadder()
    {
        var registry = BuildBaseRegistry()
            .RegisterF2W1DispersionPi2Inheritance()
            .Build();

        var ancestors = registry.AncestorsOf<F2W1DispersionPi2Inheritance>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(W1Dispersion), ancestors);
        Assert.Contains(typeof(Pi2DyadicLadderClaim), ancestors);
    }

    [Theory]
    // N=3, k=1: 4·(1−cos(π/3)) = 4·(1−0.5) = 2
    // N=4, k=2: 4·(1−cos(π/2)) = 4·(1−0) = 4
    [InlineData(3, 1.0, 1, 2.0)]
    [InlineData(4, 1.0, 2, 4.0)]
    public void RegisterF2_FrequencyAcrossRegistry(int N, double J, int k, double expected)
    {
        var registry = BuildBaseRegistry()
            .RegisterF2W1DispersionPi2Inheritance()
            .Build();

        Assert.Equal(expected, registry.Get<F2W1DispersionPi2Inheritance>().Frequency(N, J, k), precision: 12);
    }

    [Fact]
    public void RegisterF2_MatchesW1DispersionParentAcrossRegistry()
    {
        var registry = BuildBaseRegistry()
            .RegisterF2W1DispersionPi2Inheritance()
            .Build();

        Assert.True(registry.Get<F2W1DispersionPi2Inheritance>().MatchesW1DispersionParent());
    }

    [Fact]
    public void RegisterF2_WithoutW1Dispersion_Throws()
    {
        Assert.Throws<InvariantViolationException>(() =>
            new ClaimRegistryBuilder()
                .RegisterPi2Family()
                .RegisterPi2DyadicLadder()
                // Missing: RegisterW1Dispersion
                .RegisterF2W1DispersionPi2Inheritance()
                .Build());
    }
}
