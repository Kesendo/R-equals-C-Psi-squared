using System.Linq;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Diagnostics.Knowledge;
using RCPsiSquared.Diagnostics.Ptf;

namespace RCPsiSquared.Diagnostics.Tests.Ptf;

/// <summary>The empirical Q-instance refinement of the defect-reading equivariance: the decoder's painted
/// α-dictionary is anti-collinear at every N read, the seesaw's (cos = 1 − 2w on mirror pairs). Covers tier
/// (Tier2Empirical, weaker than the Tier1Derived parent), summary content, the parent edge
/// (DefectReadingEquivarianceClaim), and registration into BuildDefault. The numeric grounding is
/// SingleExcitationPainterDictionaryTests.</summary>
public class DecoderAntiCollinearitySeesawClaimTests
{
    private static DecoderAntiCollinearitySeesawClaim Build() =>
        new(equivariance: new DefectReadingEquivarianceClaim(
            new KPartnerSelectionRuleClaim(new ChiralMirrorTrajectoryClaim())));

    [Fact]
    public void Claim_TierIsTier2Empirical()
    {
        Assert.Equal(Tier.Tier2Empirical, Build().Tier);
    }

    [Fact]
    public void Summary_NamesTheIdentityAndNotAParityLaw()
    {
        var summary = Build().Summary;
        Assert.Contains("1 − 2w", summary);
        Assert.Contains("mirror-odd", summary);
        Assert.DoesNotContain("distance-2", summary);
        Assert.DoesNotContain("odd-N parity", summary);
        Assert.Contains("carrier channels", summary);
        Assert.DoesNotContain("seesaw is open", summary);
    }

    [Fact]
    public void Children_ContainDefectReadingEquivarianceParent()
    {
        var children = ((RCPsiSquared.Core.Inspection.IInspectable)Build()).Children
            .Select(c => c.GetType()).ToHashSet();
        Assert.Contains(typeof(DefectReadingEquivarianceClaim), children);
    }

    [Fact]
    public void BuildDefault_ContainsDecoderAntiCollinearitySeesawClaim()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        Assert.True(registry.Contains<DecoderAntiCollinearitySeesawClaim>());
    }

    [Fact]
    public void Claim_Ancestors_ContainDefectReadingEquivarianceClaim()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var ancestors = registry.AncestorsOf<DecoderAntiCollinearitySeesawClaim>()
            .Select(c => c.GetType()).ToHashSet();
        Assert.Contains(typeof(DefectReadingEquivarianceClaim), ancestors);
    }
}
