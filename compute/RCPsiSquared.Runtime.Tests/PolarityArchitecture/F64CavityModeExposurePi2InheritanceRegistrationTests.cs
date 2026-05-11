using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class F64CavityModeExposurePi2InheritanceRegistrationTests
{
    private static ClaimRegistryBuilder BuildBaseRegistry() =>
        new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterPi2DyadicLadder()
            .RegisterAbsorptionTheoremClaim();

    [Fact]
    public void RegisterF64_AddsClaim()
    {
        var registry = BuildBaseRegistry()
            .RegisterF64CavityModeExposurePi2Inheritance()
            .Build();

        Assert.True(registry.Contains<F64CavityModeExposurePi2Inheritance>());
    }

    [Fact]
    public void RegisterF64_TierIsTier1Derived()
    {
        var registry = BuildBaseRegistry()
            .RegisterF64CavityModeExposurePi2Inheritance()
            .Build();

        Assert.Equal(Tier.Tier1Derived,
            registry.Get<F64CavityModeExposurePi2Inheritance>().Tier);
    }

    [Fact]
    public void RegisterF64_AncestorsContainBothParents()
    {
        var registry = BuildBaseRegistry()
            .RegisterF64CavityModeExposurePi2Inheritance()
            .Build();

        var ancestors = registry.AncestorsOf<F64CavityModeExposurePi2Inheritance>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(Pi2DyadicLadderClaim), ancestors);
        Assert.Contains(typeof(QuarterAsBilinearMaxvalClaim), ancestors);
    }

    [Theory]
    [InlineData(0.5, 0.2)]
    [InlineData(1.0, 0.25)]
    [InlineData(2.0, 0.1)]
    public void RegisterF64_N3GValueAcrossRegistry(double r, double expected)
    {
        var registry = BuildBaseRegistry()
            .RegisterF64CavityModeExposurePi2Inheritance()
            .Build();

        Assert.Equal(expected, registry.Get<F64CavityModeExposurePi2Inheritance>().N3GValue(r), precision: 12);
    }

    [Fact]
    public void RegisterF64_UniformJMatchesQuarterAcrossRegistry()
    {
        var registry = BuildBaseRegistry()
            .RegisterF64CavityModeExposurePi2Inheritance()
            .Build();

        Assert.True(registry.Get<F64CavityModeExposurePi2Inheritance>().UniformJSpecialValueMatchesQuarter());
    }

    [Fact]
    public void RegisterF64_N3CrossoverContinuityHoldsAcrossRegistry()
    {
        var registry = BuildBaseRegistry()
            .RegisterF64CavityModeExposurePi2Inheritance()
            .Build();

        Assert.True(registry.Get<F64CavityModeExposurePi2Inheritance>().N3CrossoverContinuityHolds());
    }

    [Fact]
    public void RegisterF64_WithoutPi2DyadicLadder_Throws()
    {
        Assert.Throws<InvariantViolationException>(() =>
            new ClaimRegistryBuilder()
                .RegisterPi2Family()
                // Missing: RegisterPi2DyadicLadder
                .RegisterF64CavityModeExposurePi2Inheritance()
                .Build());
    }
}
