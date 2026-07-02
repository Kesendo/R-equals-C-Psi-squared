using System;
using System.Linq;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Knowledge;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>Wiring audit for <see cref="SpectatorIntertwinerClaim"/> (Theorem B of
/// PROOF_CODIM1_BY_ADDITIVITY): the site-summed spectator W(ρ) = Σ_l c_l†ρc_l is an exact part-by-part
/// intertwiner, block-shifting (p,q̃)→(p+1,q̃+1) and carrying the Jordan block up the diamond. Two typed
/// parents, both Tier1Derived: <see cref="F89CrossFoldSimilarityClaim"/> (F89d, the antiunitary orbit leg)
/// and <see cref="AbsorptionTheoremClaim"/> (the rate law behind A = −2·diag(n_diff)).</summary>
public class SpectatorIntertwinerClaimRegistrationTests
{
    [Fact]
    public void BuildDefault_ContainsClaim()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        Assert.True(registry.Contains<SpectatorIntertwinerClaim>());
    }

    [Fact]
    public void Claim_IsTier1Derived()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        Assert.Equal(Tier.Tier1Derived, registry.Get<SpectatorIntertwinerClaim>().Tier);
    }

    [Fact]
    public void Claim_Ancestors_ContainBothTypedParents()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var ancestors = registry.AncestorsOf<SpectatorIntertwinerClaim>()
            .Select(c => c.GetType()).ToHashSet();
        Assert.Contains(typeof(F89CrossFoldSimilarityClaim), ancestors);
        Assert.Contains(typeof(AbsorptionTheoremClaim), ancestors);
    }

    [Fact]
    public void SigmaMinClimbingRung_IsSqrtTwo()
    {
        // The gate-pinned injectivity of the (1,2)→(2,3) rung at N=5: σ_min(W) = √2 exactly.
        Assert.Equal(Math.Sqrt(2.0), SpectatorIntertwinerClaim.SigmaMinClimbingRungN5, 15);
    }
}
