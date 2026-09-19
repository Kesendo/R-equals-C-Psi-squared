using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Knowledge;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class ApproachFamilyCarrierClaimRegistrationTests
{
    [Fact]
    public void DefaultRegistry_ContainsApproachFamilyCarrierClaim()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        Assert.True(registry.Contains<ApproachFamilyCarrierClaim>());
    }

    [Fact]
    public void ApproachFamilyCarrierClaim_HasExactlyTwoTypedParents()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();

        var directParents = registry.EdgesInto<ApproachFamilyCarrierClaim>()
            .Select(edge => edge.Parent).ToHashSet();

        Assert.Equal(
            new[]
            {
                typeof(AbsorptionTheoremClaim),
                typeof(F25CPsiBellPlusPi2Inheritance),
            }.OrderBy(type => type.FullName),
            directParents.OrderBy(type => type.FullName));
        Assert.DoesNotContain(typeof(TwoReadingsClaim), directParents);
        Assert.Equal(2, directParents.Count);

        var ancestors = registry.AncestorsOf<ApproachFamilyCarrierClaim>()
            .Select(claim => claim.GetType()).ToHashSet();
        Assert.DoesNotContain(typeof(UniversalCarrierClaim), ancestors);
        Assert.DoesNotContain(typeof(PolynomialDiscriminantAnchorClaim), ancestors);
    }
}
