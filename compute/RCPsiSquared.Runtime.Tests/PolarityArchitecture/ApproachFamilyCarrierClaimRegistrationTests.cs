using RCPsiSquared.Core.F86.Item1Derivation;
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
    public void ApproachFamilyCarrierClaim_HasFourTypedParents()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();

        var ancestors = registry.AncestorsOf<ApproachFamilyCarrierClaim>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(UniversalCarrierClaim), ancestors);
        Assert.Contains(typeof(C2BareDoubledPtfClosedForm), ancestors);
        Assert.Contains(typeof(TwoReadingsClaim), ancestors);
        Assert.Contains(typeof(F25CPsiBellPlusPi2Inheritance), ancestors);
    }
}
