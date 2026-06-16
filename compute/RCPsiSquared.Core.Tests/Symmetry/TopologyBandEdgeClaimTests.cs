using System.Linq;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using Xunit;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class TopologyBandEdgeClaimTests
{
    [Fact]
    public void Build_WiresParents_AndIsTier1Candidate()
    {
        var claim = TopologyBandEdgeClaim.Build();
        // Tier1Candidate, NOT Derived: the typed parent edge to ClockHandLadderClaim (Tier1Candidate)
        // forces the tier down by the registry's tier-inheritance invariant (a child can be no stronger
        // than its weakest parent). The band-edge LAW is exact and the floor is the Tier1Derived
        // Absorption Theorem, but the gap-dominance map reuses the chain's open gap-dominance proof.
        // The sibling CoherenceHorizonClaim (the other ClockHandLadder child) is Tier1Candidate too.
        Assert.Equal(Tier.Tier1Candidate, claim.Tier);
        Assert.NotNull(claim.ChainInstance);   // ClockHandLadderClaim (the chain Im instance)
        Assert.NotNull(claim.Absorption);      // AbsorptionTheoremClaim (the Re=-2γ floor)
        var kids = ((IInspectable)claim).Children.ToList();
        Assert.Contains(claim.ChainInstance, kids);
        Assert.Contains(claim.Absorption, kids);
    }
}
