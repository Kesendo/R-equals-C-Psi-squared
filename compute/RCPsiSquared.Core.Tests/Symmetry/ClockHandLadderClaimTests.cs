using System.Linq;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using Xunit;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class ClockHandLadderClaimTests
{
    [Fact]
    public void Claim_IsTier1Derived_AndNamesTheTwoClocks()
    {
        // Tier1Derived since 2026-06-16: the gap-dominance lemma is proven (PROOF_CHAIN_GAP_DOMINANCE,
        // free fermions); all three parents (F2b, Absorption, UniversalCarrier) are Tier1Derived.
        var c = ClockHandLadderClaim.Shared;
        Assert.Equal(Tier.Tier1Derived, c.Tier);
        Assert.Contains("ω_mem", c.Summary);
    }

    [Fact]
    public void Claim_CitesItsMarkdownAndCSharpHomes_BothDirections()
    {
        var c = ClockHandLadderClaim.Shared;
        Assert.Contains("ANALYTICAL_FORMULAS.md", c.Anchor);
        Assert.Contains("ON_HOW_GAMMA_BECAME_THE_TICK", c.Anchor);
        Assert.Contains("ClockHandLadderWitness", c.Anchor); // C# witness named for discoverability
    }

    [Fact]
    public void Claim_HasTheThreeTypedParents_F2b_Absorption_Carrier()
    {
        var c = ClockHandLadderClaim.Shared;
        var children = ((IInspectable)c).Children.ToList();
        Assert.Contains(children, ch => ch is F2bXyChainSpectrumPi2Inheritance);
        Assert.Contains(children, ch => ch is AbsorptionTheoremClaim);
        Assert.Contains(children, ch => ch is UniversalCarrierClaim);
    }
}
