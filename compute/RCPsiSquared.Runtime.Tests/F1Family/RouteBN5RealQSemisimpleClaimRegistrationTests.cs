using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Knowledge;

namespace RCPsiSquared.Runtime.Tests.F1Family;

/// <summary>Wiring tests for the parameterless F164 claim. Its prior-work anchors are provenance,
/// not typed claim dependencies, so it has no direct parents or ancestors.</summary>
public class RouteBN5RealQSemisimpleClaimRegistrationTests
{
    [Fact]
    public void TheClaimIsRegisteredAndTier1Derived()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        Assert.True(registry.Contains<RouteBN5RealQSemisimpleClaim>());
        var claim = registry.Get<RouteBN5RealQSemisimpleClaim>();
        Assert.Equal(Tier.Tier1Derived, claim.Tier);
        Assert.Contains("F164 certifies all 26 R-odd loci", claim.Summary);
        Assert.Contains("38 distinct exact-certified loci", claim.Summary);
        Assert.Contains("20 nonreal-q R-even loci numerical-only", claim.Summary);
    }

    [Fact]
    public void TheClaimHasNoTypedParentsOrAncestors()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var claim = registry.Get<RouteBN5RealQSemisimpleClaim>();
        Assert.Empty(registry.AncestorsOf<RouteBN5RealQSemisimpleClaim>());
        Assert.DoesNotContain(claim.Children, child => child is Claim);
    }

    [Fact]
    public void TheClaimConstructsWithoutParameters()
        => Assert.Equal(Tier.Tier1Derived, new RouteBN5RealQSemisimpleClaim().Tier);
}
