using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Foundation;
using RCPsiSquared.Diagnostics.Knowledge;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class HandoverFloorClaimRegistrationTests
{
    [Fact]
    public void DefaultRegistry_ContainsHandoverFloorClaim()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        Assert.True(registry.Contains<HandoverFloorClaim>());
    }

    [Fact]
    public void HandoverFloorClaim_ResolvesFromDefaultRegistry_WithSharedTypedParents()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();

        var claim = registry.Get<HandoverFloorClaim>();
        Assert.NotNull(claim);

        // The registration must wire the SHARED registry instances, not fresh ones.
        Assert.Same(registry.Get<AbsorptionTheoremClaim>(), claim.Survival);
        Assert.Same(registry.Get<F50WeightOneDegeneracyPi2Inheritance>(), claim.Floor);
        Assert.Same(registry.Get<CoherenceHorizonClaim>(), claim.ChainSolution);
    }

    [Fact]
    public void HandoverFloorClaim_HasThreeTypedParents()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();

        var ancestors = registry.AncestorsOf<HandoverFloorClaim>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(AbsorptionTheoremClaim), ancestors);
        Assert.Contains(typeof(F50WeightOneDegeneracyPi2Inheritance), ancestors);
        Assert.Contains(typeof(CoherenceHorizonClaim), ancestors);
    }

    [Fact]
    public void HandoverFloorClaim_BatteryAllPass()
    {
        // The 5-case battery (the F50 floor, chain = Q*(N) at N=4 and the trace-dressing gap at N=6,
        // the distinct ring (2,2) seam, the ring growth) is computed live via the witness.
        var claim = HandoverFloorClaim.Shared;
        Assert.Equal(claim.Cases.Count, claim.PassCount);
    }
}
