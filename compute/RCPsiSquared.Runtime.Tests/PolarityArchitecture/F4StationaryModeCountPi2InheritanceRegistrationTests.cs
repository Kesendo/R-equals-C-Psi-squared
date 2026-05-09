using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class F4StationaryModeCountPi2InheritanceRegistrationTests
{
    private static ClaimRegistryBuilder BuildBaseRegistry() =>
        new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterPi2DyadicLadder();

    [Fact]
    public void RegisterF4_AddsClaim()
    {
        var registry = BuildBaseRegistry()
            .RegisterF4StationaryModeCountPi2Inheritance()
            .Build();

        Assert.True(registry.Contains<F4StationaryModeCountPi2Inheritance>());
    }

    [Fact]
    public void RegisterF4_TierIsTier1Derived()
    {
        var registry = BuildBaseRegistry()
            .RegisterF4StationaryModeCountPi2Inheritance()
            .Build();

        Assert.Equal(Tier.Tier1Derived,
            registry.Get<F4StationaryModeCountPi2Inheritance>().Tier);
    }

    [Fact]
    public void RegisterF4_AncestorsContainPi2DyadicLadder()
    {
        var registry = BuildBaseRegistry()
            .RegisterF4StationaryModeCountPi2Inheritance()
            .Build();

        var ancestors = registry.AncestorsOf<F4StationaryModeCountPi2Inheritance>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(Pi2DyadicLadderClaim), ancestors);
    }

    [Theory]
    [InlineData(2, 10)]
    [InlineData(3, 24)]
    [InlineData(4, 54)]
    [InlineData(5, 120)]
    public void RegisterF4_StationaryModeCountAcrossRegistry(int N, int expected)
    {
        var registry = BuildBaseRegistry()
            .RegisterF4StationaryModeCountPi2Inheritance()
            .Build();

        Assert.Equal(expected, registry.Get<F4StationaryModeCountPi2Inheritance>().StationaryModeCount(N));
    }

    [Fact]
    public void RegisterF4_SchurWeylIdentityHoldsAcrossRegistry()
    {
        var registry = BuildBaseRegistry()
            .RegisterF4StationaryModeCountPi2Inheritance()
            .Build();

        var f4 = registry.Get<F4StationaryModeCountPi2Inheritance>();
        for (int N = 1; N <= 10; N++)
            Assert.True(f4.SchurWeylDimensionIdentityHolds(N));
    }

    [Fact]
    public void RegisterF4_WithoutPi2DyadicLadder_Throws()
    {
        Assert.Throws<InvariantViolationException>(() =>
            new ClaimRegistryBuilder()
                .RegisterPi2Family()
                // Missing: RegisterPi2DyadicLadder
                .RegisterF4StationaryModeCountPi2Inheritance()
                .Build());
    }
}
