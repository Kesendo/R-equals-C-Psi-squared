using System;
using System.Linq;
using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.F86;
using RCPsiSquared.Core.Knowledge;
using Xunit;

namespace RCPsiSquared.Core.Tests.F86;

/// <summary>F86 ↔ FRAGILE_BRIDGE meta-claim tests. Pins the empirical witness table
/// (c=2 N=5..8 Petermann-K sweep, 2026-05-06), the parity-asymmetry prediction from A3,
/// the Tier2Verified honesty (analytic continuation gap), and KB integration.</summary>
public class LocalGlobalEpLinkTests
{
    [Fact]
    public void LocalGlobalEpLink_IsTier2Verified()
    {
        var link = LocalGlobalEpLink.Build();
        Assert.Equal(Tier.Tier2Verified, link.Tier);
    }

    [Fact]
    public void LocalGlobalEpLink_HasFourWitnesses_C2N5To8()
    {
        var link = LocalGlobalEpLink.Build();
        Assert.Equal(4, link.Witnesses.Count);
        Assert.Equal(5, link.Witnesses[0].N);
        Assert.Equal(8, link.Witnesses[3].N);
    }

    [Fact]
    public void LocalGlobalEpLink_ParityAsymmetry_OddDominatesEven()
    {
        var link = LocalGlobalEpLink.Build();
        var odd = link.Witnesses.Where(w => w.Parity == "odd").Select(w => w.MaxKGlobal);
        var even = link.Witnesses.Where(w => w.Parity == "even").Select(w => w.MaxKGlobal);
        // A3 prediction: σ_0 R-even/R-odd degeneracy at even N → parity asymmetry
        // Empirical: odd-N peaks dominate even-N peaks by 2-4× at every step
        Assert.True(odd.Min() > even.Max(), "Odd-N min K should exceed even-N max K");
    }

    [Fact]
    public void LocalGlobalEpLink_AboveFragileBridgeKBallpark_ByN7()
    {
        var link = LocalGlobalEpLink.Build();
        var n7 = link.Witnesses.Single(w => w.N == 7);
        // FRAGILE_BRIDGE Petermann K=403 is the global-EP near-singularity ballpark.
        // By N=7 the c=2 local EP K is ~6× above this on the real Q axis,
        // empirically confirming "real-axis hit" (not "off-axis siblings").
        Assert.True(n7.MaxKGlobal > 403.0 * 5.0,
            $"By N=7 c=2 K should exceed 5×K_FRAGILE; got {n7.MaxKGlobal}");
    }

    [Fact]
    public void LocalGlobalEpLink_PendingDerivationNote_NamesAnalyticContinuation()
    {
        var link = LocalGlobalEpLink.Build();
        Assert.NotNull(link.PendingDerivationNote);
        Assert.Contains("complex-γ", link.PendingDerivationNote!, StringComparison.OrdinalIgnoreCase);
        Assert.Contains("analytic continuation", link.PendingDerivationNote!, StringComparison.OrdinalIgnoreCase);
    }

    [Fact]
    public void F86KnowledgeBase_ExposesLocalGlobalEpLink()
    {
        // The link is a meta-claim (block-independent); should be exposed at the KB root
        // for any block, not just c=2.
        var block = new CoherenceBlock(N: 5, n: 1, gammaZero: 0.05);
        var kb = new F86KnowledgeBase(block);
        Assert.NotNull(kb.LocalGlobalEpLink);
        Assert.Equal(Tier.Tier2Verified, kb.LocalGlobalEpLink.Tier);

        // Also verify availability for a c=1 block (block-independent meta-claim).
        var c1Block = new CoherenceBlock(N: 3, n: 0, gammaZero: 0.05);
        var c1Kb = new F86KnowledgeBase(c1Block);
        Assert.NotNull(c1Kb.LocalGlobalEpLink);
    }
}
