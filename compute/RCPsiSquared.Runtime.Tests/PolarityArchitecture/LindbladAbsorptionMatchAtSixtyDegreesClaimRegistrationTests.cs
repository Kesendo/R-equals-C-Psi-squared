using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class LindbladAbsorptionMatchAtSixtyDegreesClaimRegistrationTests
{
    private static ClaimRegistryBuilder BuildBaseRegistry() =>
        new ClaimRegistryBuilder()
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
    public void RegisterLindbladAbsorptionMatchAtSixtyDegreesClaim_HasOnlyF95Ancestor()
    {
        var registry = BuildBaseRegistry()
            .RegisterLindbladAbsorptionMatchAtSixtyDegreesClaim()
            .Build();

        var ancestors = registry.AncestorsOf<LindbladAbsorptionMatchAtSixtyDegreesClaim>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(F95AngleAtQuadraticZeroPi2Inheritance), ancestors);
        Assert.DoesNotContain(typeof(AbsorptionTheoremClaim), ancestors);
        Assert.DoesNotContain(typeof(CanonicalTrigAnchorPi2Inheritance), ancestors);
        Assert.Single(ancestors);
    }

    [Fact]
    public void RegisterLindbladAbsorptionMatchAtSixtyDegreesClaim_WithoutF95_Throws()
    {
        Assert.Throws<InvariantViolationException>(() =>
            new ClaimRegistryBuilder()
                .RegisterLindbladAbsorptionMatchAtSixtyDegreesClaim()
                .Build());
    }
}
