using System.Linq;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using Xunit;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class TopologyBandEdgeClaimTests
{
    [Fact]
    public void Build_WiresParents_AndIsTier1Derived()
    {
        var claim = TopologyBandEdgeClaim.Build();
        // Tier1Derived since 2026-06-16: both typed parents are now Tier1Derived. The band-edge LAW is exact,
        // the Re=−2γ floor is the Absorption Theorem, and the chain gap-dominance that capped the parent
        // ClockHandLadderClaim is proven (PROOF_CHAIN_GAP_DOMINANCE), so ClockHandLadder graduated and the
        // inherited cap is lifted. (The sibling CoherenceHorizonClaim graduated 2026-07-19: its own open
        // piece — the ring 2-excitation (2,2)/(N−2,N−2) doublet V-Effect seam — was resolved by the
        // reviewed PROOF_RING_HANDOVER_SLOPE.)
        Assert.Equal(Tier.Tier1Derived, claim.Tier);
        Assert.NotNull(claim.ChainInstance);   // ClockHandLadderClaim (the chain Im instance)
        Assert.NotNull(claim.Absorption);      // AbsorptionTheoremClaim (the Re=-2γ floor)
        var kids = ((IInspectable)claim).Children.ToList();
        Assert.Contains(claim.ChainInstance, kids);
        Assert.Contains(claim.Absorption, kids);
    }
}
