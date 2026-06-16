using System.Linq;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Knowledge;
using Xunit;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class SecondClockRegimeClaimRegistrationTests
{
    [Fact]
    public void DefaultRegistry_ContainsSecondClockRegimeClaim()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        Assert.True(registry.Contains<SecondClockRegimeClaim>());
    }

    [Fact]
    public void SecondClockRegimeClaim_ResolvesFromDefaultRegistry_WithSharedTypedParents()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();

        var claim = registry.Get<SecondClockRegimeClaim>();
        Assert.NotNull(claim);

        // the registration must wire the SHARED registry instances, not fresh ones
        Assert.Same(registry.Get<CoherenceHorizonClaim>(), claim.Horizon);
        Assert.Same(registry.Get<StructuralCeilingClaim>(), claim.Ceiling);
    }

    [Fact]
    public void SecondClockRegimeClaim_HasBothRegimeParents()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();

        var ancestors = registry.AncestorsOf<SecondClockRegimeClaim>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(CoherenceHorizonClaim), ancestors);
        Assert.Contains(typeof(StructuralCeilingClaim), ancestors);
    }
}
