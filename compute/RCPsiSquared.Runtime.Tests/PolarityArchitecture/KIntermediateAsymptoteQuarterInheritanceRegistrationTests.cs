using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class KIntermediateAsymptoteQuarterInheritanceRegistrationTests
{
    private static ClaimRegistryBuilder BuildBaseRegistry() =>
        new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterPi2DyadicLadder()
            .RegisterDickeSuperpositionQuarterPi2Inheritance();

    [Fact]
    public void RegisterKIntermediateAsymptoteQuarterInheritance_AddsClaim()
    {
        var registry = BuildBaseRegistry()
            .RegisterKIntermediateAsymptoteQuarterInheritance()
            .Build();

        Assert.True(registry.Contains<KIntermediateAsymptoteQuarterInheritance>());
    }

    [Fact]
    public void RegisterKIntermediateAsymptoteQuarterInheritance_TierIsTier1Derived()
    {
        var registry = BuildBaseRegistry()
            .RegisterKIntermediateAsymptoteQuarterInheritance()
            .Build();

        Assert.Equal(Tier.Tier1Derived,
            registry.Get<KIntermediateAsymptoteQuarterInheritance>().Tier);
    }

    [Fact]
    public void RegisterKIntermediateAsymptoteQuarterInheritance_AncestorsContainAllThreeParents()
    {
        var registry = BuildBaseRegistry()
            .RegisterKIntermediateAsymptoteQuarterInheritance()
            .Build();

        var ancestors = registry.AncestorsOf<KIntermediateAsymptoteQuarterInheritance>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(QuarterAsBilinearMaxvalClaim), ancestors);
        Assert.Contains(typeof(HalfAsStructuralFixedPointClaim), ancestors);
        Assert.Contains(typeof(DickeSuperpositionQuarterPi2Inheritance), ancestors);
    }

    [Fact]
    public void RegisterKIntermediateAsymptoteQuarterInheritance_WithoutDickeSuperposition_Throws()
    {
        Assert.Throws<InvariantViolationException>(() =>
            new ClaimRegistryBuilder()
                .RegisterPi2Family()
                .RegisterPi2DyadicLadder()
                // Missing: RegisterDickeSuperpositionQuarterPi2Inheritance
                .RegisterKIntermediateAsymptoteQuarterInheritance()
                .Build());
    }
}
