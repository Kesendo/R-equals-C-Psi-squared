using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Knowledge;

namespace RCPsiSquared.Runtime.Tests.F1Family;

/// <summary>Wiring tests for <see cref="RouteBN5RealQSemisimpleClaim"/> (F164), whose one typed parent is
/// <see cref="F89Path3OcticEpClaim"/>: the same (1,2) coherence block one chain shorter, the same character
/// question, the same verdict reached numerically rather than by exact arithmetic.</summary>
public class RouteBN5RealQSemisimpleClaimRegistrationTests
{
    [Fact]
    public void TheClaimIsRegisteredAndTier1Derived()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        Assert.True(registry.Contains<RouteBN5RealQSemisimpleClaim>());
        Assert.Equal(Tier.Tier1Derived, registry.Get<RouteBN5RealQSemisimpleClaim>().Tier);
    }

    [Fact]
    public void TheTypedParentIsTheOcticSiblingAndIsShared()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var claim = registry.Get<RouteBN5RealQSemisimpleClaim>();
        var ancestors = registry.AncestorsOf<RouteBN5RealQSemisimpleClaim>().Select(c => c.GetType()).ToHashSet();
        Assert.Contains(typeof(F89Path3OcticEpClaim), ancestors);
        Assert.Same(registry.Get<F89Path3OcticEpClaim>(), claim.OcticSibling);
    }

    /// <summary>The tier rule: a Tier1Derived claim may not rest on a weaker parent.</summary>
    [Fact]
    public void TheParentIsNotWeakerThanTheChild()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        Assert.Equal(Tier.Tier1Derived, registry.Get<F89Path3OcticEpClaim>().Tier);
    }

    [Fact]
    public void TheClaimRejectsANullParent()
        => Assert.Throws<ArgumentNullException>(() => new RouteBN5RealQSemisimpleClaim(null!));
}
