using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Diagnostics.F87;

namespace RCPsiSquared.Diagnostics.Tests.F87;

public class DissipatorAxisSelectsPolarityClaimTests
{
    [Fact]
    public void Claim_IsTier1Derived()
    {
        var c = new DissipatorAxisSelectsPolarityClaim();
        Assert.Equal(Tier.Tier1Derived, c.Tier);
    }

    [Fact]
    public void Claim_AnchorReferences_AllThreeBridgeSources()
    {
        var c = new DissipatorAxisSelectsPolarityClaim();
        Assert.Contains("THE_POLARITY_LAYER.md", c.Anchor);
        Assert.Contains("DissipatorResonanceLaw", c.Anchor);
        Assert.Contains("PolarityLayerOriginClaim", c.Anchor);
    }

    [Fact]
    public void Claim_HasNamedCascadeChildren_FromWhatToProvenance()
    {
        var c = new DissipatorAxisSelectsPolarityClaim();
        IInspectable node = c;
        var labels = node.Children.Select(ch => ch.DisplayName).ToList();
        // tier + anchor metadata (from Claim base) + the cascade
        Assert.Contains("what is being differentiated", labels);
        Assert.Contains("where F87-hardness lives", labels);
        Assert.Contains("the selector (Z → bit_b axis)", labels);
        Assert.Contains("the selector (X → bit_a axis)", labels);
        Assert.Contains("the selector (Y → both axes)", labels);
        Assert.Contains("two readings unified (Brecher ↔ Hardness)", labels);
        Assert.Contains("γ-as-light bridge", labels);
        Assert.Contains("operational anchor", labels);
    }

    [Fact]
    public void Claim_IsRegisteredInF87KnowledgeBase_AsTier1Property()
    {
        var chain = new ChainSystem(N: 3, J: 1.0, GammaZero: 0.05);
        var kb = new F87KnowledgeBase(chain);
        Assert.NotNull(kb.DissipatorAxisSelectsPolarity);
        Assert.Equal(Tier.Tier1Derived, kb.DissipatorAxisSelectsPolarity.Tier);
    }

    [Fact]
    public void F87KnowledgeBase_TierOneGroup_HasThreeClaims()
    {
        var chain = new ChainSystem(N: 3, J: 1.0, GammaZero: 0.05);
        var kb = new F87KnowledgeBase(chain);
        IInspectable root = kb;
        var tier1Group = root.Children.FirstOrDefault(c => c.DisplayName.StartsWith("Tier 1"));
        Assert.NotNull(tier1Group);
        // 3 children: Trichotomy, DissipatorResonance, DissipatorAxisSelectsPolarity
        Assert.Equal(3, tier1Group!.Children.Count());
    }
}
