using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class F56CriticalSlowingPi2InheritanceRegistrationTests
{
    private static ClaimRegistryBuilder BuildBaseRegistry() =>
        new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterPi2DyadicLadder();

    [Fact]
    public void RegisterF56_AddsClaim()
    {
        var registry = BuildBaseRegistry()
            .RegisterF56CriticalSlowingPi2Inheritance()
            .Build();
        Assert.True(registry.Contains<F56CriticalSlowingPi2Inheritance>());
    }

    [Fact]
    public void RegisterF56_TierIsTier1Derived()
    {
        var registry = BuildBaseRegistry()
            .RegisterF56CriticalSlowingPi2Inheritance()
            .Build();
        Assert.Equal(Tier.Tier1Derived, registry.Get<F56CriticalSlowingPi2Inheritance>().Tier);
    }

    [Fact]
    public void RegisterF56_AncestorsContainAllThreeParents()
    {
        var registry = BuildBaseRegistry()
            .RegisterF56CriticalSlowingPi2Inheritance()
            .Build();
        var ancestors = registry.AncestorsOf<F56CriticalSlowingPi2Inheritance>()
            .Select(c => c.GetType()).ToHashSet();
        Assert.Contains(typeof(Pi2DyadicLadderClaim), ancestors);
        Assert.Contains(typeof(QuarterAsBilinearMaxvalClaim), ancestors);
        Assert.Contains(typeof(HalfAsStructuralFixedPointClaim), ancestors);
    }

    [Theory]
    [InlineData(0.01, 0.001, 1.23768)]
    public void RegisterF56_IterationCountAcrossRegistry(double epsilon, double tol, double expected)
    {
        var registry = BuildBaseRegistry()
            .RegisterF56CriticalSlowingPi2Inheritance()
            .Build();
        Assert.Equal(expected, registry.Get<F56CriticalSlowingPi2Inheritance>().IterationCount(epsilon, tol), precision: 4);
    }

    [Fact]
    public void RegisterF56_SixteenIsFourSquaredAcrossRegistry()
    {
        var registry = BuildBaseRegistry()
            .RegisterF56CriticalSlowingPi2Inheritance()
            .Build();
        Assert.True(registry.Get<F56CriticalSlowingPi2Inheritance>().SixteenIsFourSquared());
    }

    [Fact]
    public void RegisterF56_CardioidCuspMatchesQuarterAcrossRegistry()
    {
        var registry = BuildBaseRegistry()
            .RegisterF56CriticalSlowingPi2Inheritance()
            .Build();
        Assert.True(registry.Get<F56CriticalSlowingPi2Inheritance>().CardioidCuspMatchesQuarter());
    }

    [Fact]
    public void RegisterF56_WithoutPi2DyadicLadder_Throws()
    {
        Assert.Throws<InvariantViolationException>(() =>
            new ClaimRegistryBuilder()
                .RegisterPi2Family()
                .RegisterF56CriticalSlowingPi2Inheritance()
                .Build());
    }
}
