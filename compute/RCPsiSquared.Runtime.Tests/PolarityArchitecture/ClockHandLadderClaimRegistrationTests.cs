using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Knowledge;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class ClockHandLadderClaimRegistrationTests
{
    [Fact]
    public void DefaultRegistry_ContainsClockHandLadderClaim()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        Assert.True(registry.Contains<ClockHandLadderClaim>());
    }

    [Fact]
    public void ClockHandLadderClaim_ResolvesFromDefaultRegistry_WithSharedTypedParents()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();

        var claim = registry.Get<ClockHandLadderClaim>();
        Assert.NotNull(claim);

        // The registration must wire the SHARED registry instances, not fresh ones.
        Assert.Same(registry.Get<F2bXyChainSpectrumPi2Inheritance>(), claim.BandEdge);
        Assert.Same(registry.Get<AbsorptionTheoremClaim>(), claim.Absorption);
        Assert.Same(registry.Get<UniversalCarrierClaim>(), claim.Carrier);
    }

    [Fact]
    public void ClockHandLadderClaim_HasThreeTypedParents()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();

        var ancestors = registry.AncestorsOf<ClockHandLadderClaim>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(F2bXyChainSpectrumPi2Inheritance), ancestors);
        Assert.Contains(typeof(AbsorptionTheoremClaim), ancestors);
        Assert.Contains(typeof(UniversalCarrierClaim), ancestors);
    }
}
