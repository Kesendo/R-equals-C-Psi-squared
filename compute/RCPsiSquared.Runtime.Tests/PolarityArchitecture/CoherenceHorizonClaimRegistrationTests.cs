using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Knowledge;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class CoherenceHorizonClaimRegistrationTests
{
    [Fact]
    public void DefaultRegistry_ContainsCoherenceHorizonClaim()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        Assert.True(registry.Contains<CoherenceHorizonClaim>());
    }

    [Fact]
    public void CoherenceHorizonClaim_ResolvesFromDefaultRegistry_WithSharedTypedParents()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();

        var claim = registry.Get<CoherenceHorizonClaim>();
        Assert.NotNull(claim);

        // The registration must wire the SHARED registry instances, not fresh ones.
        Assert.Same(registry.Get<ClockHandLadderClaim>(), claim.Horizon);
        Assert.Same(registry.Get<F2bXyChainSpectrumPi2Inheritance>(), claim.BandEdge);
    }

    [Fact]
    public void CoherenceHorizonClaim_HasTwoTypedParents()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();

        var ancestors = registry.AncestorsOf<CoherenceHorizonClaim>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(ClockHandLadderClaim), ancestors);
        Assert.Contains(typeof(F2bXyChainSpectrumPi2Inheritance), ancestors);
    }
}
