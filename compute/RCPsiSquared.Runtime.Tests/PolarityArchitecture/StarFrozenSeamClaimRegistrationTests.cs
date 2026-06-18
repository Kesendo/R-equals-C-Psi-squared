using System.Linq;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Knowledge;
using Xunit;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class StarFrozenSeamClaimRegistrationTests
{
    [Fact]
    public void DefaultRegistry_ContainsStarFrozenSeamClaim()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        Assert.True(registry.Contains<StarFrozenSeamClaim>());
    }

    [Fact]
    public void StarFrozenSeamClaim_ResolvesFromDefaultRegistry_WithSharedTypedParents()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();

        var claim = registry.Get<StarFrozenSeamClaim>();
        Assert.NotNull(claim);

        // the registration must wire the SHARED registry instances, not fresh ones
        Assert.Same(registry.Get<StructuralCeilingClaim>(), claim.Ceiling);
        Assert.Same(registry.Get<SecondClockRegimeClaim>(), claim.Regime);
    }

    [Fact]
    public void StarFrozenSeamClaim_HasBothParents()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();

        var ancestors = registry.AncestorsOf<StarFrozenSeamClaim>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(StructuralCeilingClaim), ancestors);
        Assert.Contains(typeof(SecondClockRegimeClaim), ancestors);
    }
}
