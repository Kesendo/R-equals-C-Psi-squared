using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class CanonicalTrigAnchorPi2InheritanceRegistrationTests
{
    private static ClaimRegistryBuilder BuildBaseRegistry() =>
        new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterPi2DyadicLadder()
            .RegisterDickeSuperpositionQuarterPi2Inheritance()
            .RegisterKIntermediateAsymptoteQuarterInheritance();

    [Fact]
    public void RegisterCanonicalTrigAnchorPi2Inheritance_AddsClaim()
    {
        var registry = BuildBaseRegistry()
            .RegisterCanonicalTrigAnchorPi2Inheritance()
            .Build();

        Assert.True(registry.Contains<CanonicalTrigAnchorPi2Inheritance>());
    }

    [Fact]
    public void RegisterCanonicalTrigAnchorPi2Inheritance_TierIsTier1Derived()
    {
        var registry = BuildBaseRegistry()
            .RegisterCanonicalTrigAnchorPi2Inheritance()
            .Build();

        Assert.Equal(Tier.Tier1Derived,
            registry.Get<CanonicalTrigAnchorPi2Inheritance>().Tier);
    }

    [Fact]
    public void RegisterCanonicalTrigAnchorPi2Inheritance_AncestorsContainAllThreeParents()
    {
        var registry = BuildBaseRegistry()
            .RegisterCanonicalTrigAnchorPi2Inheritance()
            .Build();

        var ancestors = registry.AncestorsOf<CanonicalTrigAnchorPi2Inheritance>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(HalfAsStructuralFixedPointClaim), ancestors);
        Assert.Contains(typeof(QuarterAsBilinearMaxvalClaim), ancestors);
        Assert.Contains(typeof(KIntermediateAsymptoteQuarterInheritance), ancestors);
    }

    [Fact]
    public void RegisterCanonicalTrigAnchorPi2Inheritance_WithoutKIntermediate_Throws()
    {
        Assert.Throws<InvariantViolationException>(() =>
            new ClaimRegistryBuilder()
                .RegisterPi2Family()
                .RegisterPi2DyadicLadder()
                .RegisterDickeSuperpositionQuarterPi2Inheritance()
                // Missing: RegisterKIntermediateAsymptoteQuarterInheritance
                .RegisterCanonicalTrigAnchorPi2Inheritance()
                .Build());
    }
}
