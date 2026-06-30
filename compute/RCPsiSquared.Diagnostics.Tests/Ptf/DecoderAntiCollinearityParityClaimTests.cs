using System.Linq;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Diagnostics.Knowledge;
using RCPsiSquared.Diagnostics.Ptf;

namespace RCPsiSquared.Diagnostics.Tests.Ptf;

/// <summary>The empirical Q-instance refinement of the defect-reading equivariance: the decoder's painted
/// α-dictionary worst anti-collinearity follows odd-N parity (a Q=20 painting effect; the confuser is the
/// distance-2 bond pair, not the mirror pair). Covers tier (Tier2Empirical, weaker than the Tier1Derived
/// parent), summary content, the parent edge (DefectReadingEquivarianceClaim), and registration into
/// BuildDefault. The numeric grounding is the verifier DictionaryParityInvestigationTests.</summary>
public class DecoderAntiCollinearityParityClaimTests
{
    private static DecoderAntiCollinearityParityClaim Build() =>
        new(equivariance: new DefectReadingEquivarianceClaim(
            new KPartnerSelectionRuleClaim(new ChiralMirrorTrajectoryClaim())));

    [Fact]
    public void Claim_TierIsTier2Empirical()
    {
        Assert.Equal(Tier.Tier2Empirical, Build().Tier);
    }

    [Fact]
    public void Summary_NamesParityAndPainting()
    {
        var summary = Build().Summary;
        Assert.Contains("parity", summary, System.StringComparison.OrdinalIgnoreCase);
        Assert.Contains("painting", summary, System.StringComparison.OrdinalIgnoreCase);
        Assert.Contains("distance-2", summary, System.StringComparison.OrdinalIgnoreCase);
    }

    [Fact]
    public void Children_ContainDefectReadingEquivarianceParent()
    {
        var children = ((RCPsiSquared.Core.Inspection.IInspectable)Build()).Children
            .Select(c => c.GetType()).ToHashSet();
        Assert.Contains(typeof(DefectReadingEquivarianceClaim), children);
    }

    [Fact]
    public void BuildDefault_ContainsDecoderAntiCollinearityParityClaim()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        Assert.True(registry.Contains<DecoderAntiCollinearityParityClaim>());
    }

    [Fact]
    public void Claim_Ancestors_ContainDefectReadingEquivarianceClaim()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var ancestors = registry.AncestorsOf<DecoderAntiCollinearityParityClaim>()
            .Select(c => c.GetType()).ToHashSet();
        Assert.Contains(typeof(DefectReadingEquivarianceClaim), ancestors);
    }
}
