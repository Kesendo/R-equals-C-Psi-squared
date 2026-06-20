using System.Linq;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Diagnostics.Knowledge;
using RCPsiSquared.Diagnostics.Ptf;

namespace RCPsiSquared.Diagnostics.Tests.Ptf;

/// <summary>M3: the defect-reading map is spatial-reflection equivariant,
/// M[N−2−b,k] = (−1)^{k−1}M[b,k], the sign-location confusability a closed-form parity-weighted
/// mode sum. Covers tier, summary content, the parent edge (KPartnerSelectionRuleClaim), the live
/// battery (all cases pass), and registration into BuildDefault.</summary>
public class DefectReadingEquivarianceClaimTests
{
    private static DefectReadingEquivarianceClaim Build() =>
        new(kPartner: new KPartnerSelectionRuleClaim(new ChiralMirrorTrajectoryClaim()));

    [Fact]
    public void Claim_TierIsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived, Build().Tier);
    }

    [Fact]
    public void Summary_NamesEquivarianceAndCosine()
    {
        var summary = Build().Summary;
        Assert.Contains("equivariance", summary, System.StringComparison.OrdinalIgnoreCase);
        Assert.Contains("cos", summary, System.StringComparison.OrdinalIgnoreCase);
    }

    [Fact]
    public void Children_ContainKPartnerSelectionRuleClaim()
    {
        var children = ((RCPsiSquared.Core.Inspection.IInspectable)Build()).Children
            .Select(c => c.GetType()).ToHashSet();
        Assert.Contains(typeof(KPartnerSelectionRuleClaim), children);
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

    [Fact]
    public void BuildDefault_ContainsDefectReadingEquivarianceClaim()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        Assert.True(registry.Contains<DefectReadingEquivarianceClaim>());
    }

    [Fact]
    public void Claim_Ancestors_ContainKPartnerSelectionRuleClaim()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var ancestors = registry.AncestorsOf<DefectReadingEquivarianceClaim>()
            .Select(c => c.GetType()).ToHashSet();
        Assert.Contains(typeof(KPartnerSelectionRuleClaim), ancestors);
    }
}
