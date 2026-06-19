using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Foundation;
using RCPsiSquared.Diagnostics.Knowledge;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

/// <summary>(D) the closure functional, typed as <see cref="SurvivorDiffusionGradientClaim"/>
/// (Tier1Candidate; the eigenvalue-level dual of the stone). Wiring + battery guard: the claim resolves
/// from the default registry with its two SHARED typed parents (AbsorptionTheoremClaim = the -2gamma rate,
/// SurvivalIncompletenessMirrorClaim = the (A) survivor), and its live battery (the diffusion-Rayleigh
/// dRe ~ grad^2 law at N=4) all-passes.</summary>
public class SurvivorDiffusionGradientClaimRegistrationTests
{
    [Fact]
    public void DefaultRegistry_ContainsSurvivorDiffusionGradientClaim()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        Assert.True(registry.Contains<SurvivorDiffusionGradientClaim>());
    }

    [Fact]
    public void Claim_ResolvesWithSharedTypedParents()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var claim = registry.Get<SurvivorDiffusionGradientClaim>();
        Assert.NotNull(claim);
        Assert.Same(registry.Get<AbsorptionTheoremClaim>(), claim.RateLaw);
        Assert.Same(registry.Get<SurvivalIncompletenessMirrorClaim>(), claim.Survivor);
    }

    [Fact]
    public void Claim_HasBothTypedParents()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var ancestors = registry.AncestorsOf<SurvivorDiffusionGradientClaim>()
            .Select(c => c.GetType()).ToHashSet();
        Assert.Contains(typeof(AbsorptionTheoremClaim), ancestors);
        Assert.Contains(typeof(SurvivalIncompletenessMirrorClaim), ancestors);
    }

    [Fact]
    public void Claim_BatteryAllPass()
    {
        // the live battery: the diffusion-Rayleigh law at N=4 -- dRe ~ grad^2 (slope ~2), the ratio
        // bond-independent, and the gradient quiet at the no-flux chain ends.
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var claim = registry.Get<SurvivorDiffusionGradientClaim>();
        Assert.Equal(claim.Cases.Count, claim.PassCount);
    }
}
