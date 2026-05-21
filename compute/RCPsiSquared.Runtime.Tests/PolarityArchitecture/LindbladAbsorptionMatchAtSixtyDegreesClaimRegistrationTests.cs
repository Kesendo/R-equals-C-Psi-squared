using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class LindbladAbsorptionMatchAtSixtyDegreesClaimRegistrationTests
{
    private static ClaimRegistryBuilder BuildBaseRegistry() =>
        new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterPi2DyadicLadder()
            .RegisterDickeSuperpositionQuarterPi2Inheritance()
            .RegisterKIntermediateAsymptoteQuarterInheritance()
            .RegisterCanonicalTrigAnchorPi2Inheritance()
            .RegisterAbsorptionTheoremClaim()
            .RegisterF95AngleAtQuadraticZeroPi2Inheritance();

    [Fact]
    public void RegisterLindbladAbsorptionMatchAtSixtyDegreesClaim_AddsClaim()
    {
        var registry = BuildBaseRegistry()
            .RegisterLindbladAbsorptionMatchAtSixtyDegreesClaim()
            .Build();

        Assert.True(registry.Contains<LindbladAbsorptionMatchAtSixtyDegreesClaim>());
    }

    [Fact]
    public void RegisterLindbladAbsorptionMatchAtSixtyDegreesClaim_TierIsTier1Derived()
    {
        var registry = BuildBaseRegistry()
            .RegisterLindbladAbsorptionMatchAtSixtyDegreesClaim()
            .Build();

        Assert.Equal(Tier.Tier1Derived,
            registry.Get<LindbladAbsorptionMatchAtSixtyDegreesClaim>().Tier);
    }

    [Fact]
    public void RegisterLindbladAbsorptionMatchAtSixtyDegreesClaim_AncestorsContainAllThreeParents()
    {
        var registry = BuildBaseRegistry()
            .RegisterLindbladAbsorptionMatchAtSixtyDegreesClaim()
            .Build();

        var ancestors = registry.AncestorsOf<LindbladAbsorptionMatchAtSixtyDegreesClaim>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(F95AngleAtQuadraticZeroPi2Inheritance), ancestors);
        Assert.Contains(typeof(AbsorptionTheoremClaim), ancestors);
        Assert.Contains(typeof(CanonicalTrigAnchorPi2Inheritance), ancestors);
    }

    [Fact]
    public void RegisterLindbladAbsorptionMatchAtSixtyDegreesClaim_WithoutCanonicalTrig_Throws()
    {
        Assert.Throws<InvariantViolationException>(() =>
            new ClaimRegistryBuilder()
                .RegisterPi2Family()
                .RegisterPi2DyadicLadder()
                .RegisterDickeSuperpositionQuarterPi2Inheritance()
                .RegisterKIntermediateAsymptoteQuarterInheritance()
                .RegisterAbsorptionTheoremClaim()
                .RegisterF95AngleAtQuadraticZeroPi2Inheritance()
                // Missing: RegisterCanonicalTrigAnchorPi2Inheritance
                .RegisterLindbladAbsorptionMatchAtSixtyDegreesClaim()
                .Build());
    }
}
