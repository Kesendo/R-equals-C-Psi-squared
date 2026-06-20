using System.Linq;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Knowledge;
using RCPsiSquared.Diagnostics.Ptf;

namespace RCPsiSquared.Diagnostics.Tests.Ptf;

/// <summary>F124, the band-edge transition invariant (the reading-grammar arc's frame-theoretic capstone):
/// ‖M‖_F² + λ_min(MMᵀ) = z = 2 for the open chain's band-edge carrier, with λ_min = E the Dirichlet-edge
/// coupling. Covers: tier (Tier1Derived), the summary content, the anchor, the two typed parent edges
/// (KPartnerSelectionRuleClaim + ClockHandLadderClaim), the live gate-first battery (all cases pass), and
/// the registry wiring + ancestry.</summary>
public class BandEdgeTransitionInvariantClaimTests
{
    private static BandEdgeTransitionInvariantClaim Build() => BandEdgeTransitionInvariantClaim.Build();

    [Fact]
    public void Claim_TierIsTier1Derived()
    {
        // Both typed parents are Tier1Derived, so the tier-inheritance invariant leaves the child Tier1Derived.
        Assert.Equal(Tier.Tier1Derived, Build().Tier);
    }

    [Fact]
    public void Build_WiresBothParents_BothTier1Derived()
    {
        var claim = Build();
        Assert.NotNull(claim.KPartner);
        Assert.NotNull(claim.Clock);
        Assert.Equal(Tier.Tier1Derived, claim.KPartner.Tier);
        Assert.Equal(Tier.Tier1Derived, claim.Clock.Tier);

        var kids = ((IInspectable)claim).Children.Select(c => c.GetType()).ToHashSet();
        Assert.Contains(typeof(KPartnerSelectionRuleClaim), kids);
        Assert.Contains(typeof(ClockHandLadderClaim), kids);
    }

    [Fact]
    public void Summary_NamesIdentityAndLambdaMinAndFrame()
    {
        var summary = Build().Summary;
        Assert.Contains("‖M‖_F²", summary);
        Assert.Contains("λ_min", summary);
        Assert.Contains("K-partner", summary);
    }

    [Fact]
    public void Anchor_CitesProofAndFormulaAndWitness()
    {
        var anchor = Build().Anchor;
        Assert.Contains("PROOF_HANDSHAKE_TRANSITION_INVARIANT.md", anchor);
        Assert.Contains("F124", anchor);
        Assert.Contains("BandEdgeTransitionInvariantWitness.cs", anchor);
    }

    [Fact]
    public void Battery_AllPass()
    {
        var claim = Build();
        Assert.NotEmpty(claim.Cases);
        Assert.Equal(claim.Cases.Count, claim.PassCount);
        foreach (var c in claim.Cases)
            Assert.True(c.Passes, $"battery case '{c.Name}': expected {c.Expected}, got {c.Actual}");
    }

    // --- Registration ---

    [Fact]
    public void BuildDefault_ContainsBandEdgeTransitionInvariantClaim()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        Assert.True(registry.Contains<BandEdgeTransitionInvariantClaim>());
    }

    [Fact]
    public void Claim_Ancestors_ContainBothParents()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var ancestors = registry.AncestorsOf<BandEdgeTransitionInvariantClaim>()
            .Select(c => c.GetType()).ToHashSet();
        Assert.Contains(typeof(KPartnerSelectionRuleClaim), ancestors);
        Assert.Contains(typeof(ClockHandLadderClaim), ancestors);
    }
}
