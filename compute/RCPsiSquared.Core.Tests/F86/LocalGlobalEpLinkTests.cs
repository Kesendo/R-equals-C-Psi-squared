using System;
using System.Linq;
using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.F86;
using RCPsiSquared.Core.Knowledge;
using Xunit;

namespace RCPsiSquared.Core.Tests.F86;

/// <summary>F86 ↔ FRAGILE_BRIDGE meta-claim tests, post the 2026-06-21 F86a-retraction.
/// The claim is now an <see cref="Tier.OpenQuestion"/>: the full Σγ=N·γ₀ block is genuinely
/// non-normal on the real Q axis but has NO eigenvalue coalescence there (eigenvalues simple),
/// so the prior real-axis defective-EP reading is dropped along with the grid-sensitive peak
/// magnitudes ("6×", K=2384.7, the within-parity growth law, the parity asymmetry). The four
/// witness rows are retained ONLY as a cautionary non-normality record. These tests pin the
/// corrected tier, the retained-rows count, the correction note, and KB integration; they do
/// NOT pin any magnitude (those numbers are grid artifacts).</summary>
public class LocalGlobalEpLinkTests
{
    [Fact]
    public void LocalGlobalEpLink_IsOpenQuestion()
    {
        var link = LocalGlobalEpLink.Build();
        Assert.Equal(Tier.OpenQuestion, link.Tier);
    }

    [Fact]
    public void LocalGlobalEpLink_HasFourCautionaryWitnesses_C2N5To8()
    {
        // The four rows are RETAINED as a cautionary record of genuine non-normality near
        // each Q_peak — NOT as EP evidence and NOT as a magnitude law. We pin only their
        // presence and identity (N values), never their magnitudes (grid-sensitive, dropped).
        var link = LocalGlobalEpLink.Build();
        Assert.Equal(4, link.Witnesses.Count);
        Assert.Equal(5, link.Witnesses[0].N);
        Assert.Equal(8, link.Witnesses[3].N);
    }

    [Fact]
    public void LocalGlobalEpLink_PendingDerivationNote_RecordsCorrectionAndOpenQuestion()
    {
        var link = LocalGlobalEpLink.Build();
        Assert.NotNull(link.PendingDerivationNote);
        // The correction: genuine non-normality, NOT an eig artifact; eigenvalues simple on
        // the real axis (no coalescence); magnitudes retracted; off-axis EP is OPEN.
        Assert.Contains("non-normal", link.PendingDerivationNote!, StringComparison.OrdinalIgnoreCase);
        Assert.Contains("artifact-free", link.PendingDerivationNote!, StringComparison.OrdinalIgnoreCase);
        Assert.Contains("RETRACTED", link.PendingDerivationNote!, StringComparison.OrdinalIgnoreCase);
        Assert.Contains("OPEN QUESTION", link.PendingDerivationNote!, StringComparison.OrdinalIgnoreCase);
    }

    [Fact]
    public void F86KnowledgeBase_ExposesLocalGlobalEpLink_AsOpenQuestion()
    {
        // The link is a meta-claim (block-independent); should be exposed at the KB root
        // for any block, not just c=2.
        var block = new CoherenceBlock(N: 5, n: 1, gammaZero: 0.05);
        var kb = new F86KnowledgeBase(block);
        Assert.NotNull(kb.LocalGlobalEpLink);
        Assert.Equal(Tier.OpenQuestion, kb.LocalGlobalEpLink.Tier);

        // Also verify availability for a c=1 block (block-independent meta-claim).
        var c1Block = new CoherenceBlock(N: 3, n: 0, gammaZero: 0.05);
        var c1Kb = new F86KnowledgeBase(c1Block);
        Assert.NotNull(c1Kb.LocalGlobalEpLink);
    }
}
