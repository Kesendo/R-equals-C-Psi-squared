using System.Linq;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Diagnostics.Knowledge;
using RCPsiSquared.Diagnostics.Ptf;

namespace RCPsiSquared.Diagnostics.Tests.Ptf;

/// <summary>The K-partner selection rule (the reading-grammar arc's first DERIVED result):
/// the carrier ψ_1 never couples to its K-partner ψ_N through any bond defect, ⟨ψ_N|V_b|ψ_1⟩ = 0,
/// a two-line corollary of <see cref="ChiralMirrorTrajectoryClaim"/>; the consequence is that the
/// defect-decoder's location dictionary has rank N−2, which IS its sign-location ambiguity (the
/// K-partner channel is the dictionary's null direction). Covers: tier, summary content, the
/// double-anchor, the typed parent edge, and the live battery (all cases pass).</summary>
public class KPartnerSelectionRuleClaimTests
{
    private static KPartnerSelectionRuleClaim Build() =>
        new(chiralMirror: new ChiralMirrorTrajectoryClaim());

    [Fact]
    public void Claim_TierIsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived, Build().Tier);
    }

    [Fact]
    public void Summary_NamesSelectionRuleAndRankAndKPartnership()
    {
        var summary = Build().Summary;
        Assert.Contains("selection rule", summary, System.StringComparison.OrdinalIgnoreCase);
        Assert.Contains("rank", summary, System.StringComparison.OrdinalIgnoreCase);
        Assert.Contains("K-partner", summary, System.StringComparison.OrdinalIgnoreCase);
    }

    [Fact]
    public void Anchor_CitesBothDirections()
    {
        var anchor = Build().Anchor;
        Assert.Contains("HANDSHAKE_GEOMETRY.md", anchor);
        Assert.Contains("DefectDecoder.cs", anchor);
        Assert.Contains("k_partner_selection_rule.py", anchor);
    }

    [Fact]
    public void Children_ContainChiralMirrorTrajectoryClaim()
    {
        var claim = Build();
        var children = ((RCPsiSquared.Core.Inspection.IInspectable)claim).Children
            .Select(c => c.GetType()).ToHashSet();
        Assert.Contains(typeof(ChiralMirrorTrajectoryClaim), children);
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
    public void BuildDefault_ContainsKPartnerSelectionRuleClaim()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        Assert.True(registry.Contains<KPartnerSelectionRuleClaim>());
    }

    [Fact]
    public void Claim_Ancestors_ContainChiralMirrorTrajectoryClaim()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var ancestors = registry.AncestorsOf<KPartnerSelectionRuleClaim>()
            .Select(c => c.GetType()).ToHashSet();
        Assert.Contains(typeof(ChiralMirrorTrajectoryClaim), ancestors);
    }
}
