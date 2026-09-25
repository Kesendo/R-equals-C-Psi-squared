using System;
using System.Linq;
using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.F86;
using RCPsiSquared.Core.Knowledge;
using Xunit;

namespace RCPsiSquared.Core.Tests.F86;

/// <summary>F86 ↔ FRAGILE_BRIDGE meta-claim tests, post the 2026-06-21 F86a-retraction.
/// The claim is now an <see cref="Tier.OpenQuestion"/>: the full Σγ=N·γ₀ block is genuinely
/// non-normal at the sampled real-Q peaks, whose eigenvalues are simple. Separately certified
/// real-axis F89 seeds are outside that coarse sweep; the grid-sensitive peak
/// magnitudes ("6×", K=2384.7, the within-parity growth law, the parity asymmetry) are dropped.
/// The four witness rows are retained ONLY as a cautionary non-normality record. These tests pin
/// the corrected tier, the retained-rows count, the correction note, the scope of the
/// FRAGILE_BRIDGE threshold statement, and KB integration; they do NOT pin any Petermann magnitude
/// (those numbers are grid artifacts). The FRAGILE_BRIDGE threshold itself is recomputed in
/// <see cref="FragileBridgeThresholdTests"/>.</summary>
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
        // The current boundary: the sampled peak magnitudes are not EP evidence; certified
        // real-axis seeds and the distinct open off-axis question remain separate.
        Assert.Contains("non-normal", link.PendingDerivationNote!, StringComparison.OrdinalIgnoreCase);
        Assert.Contains("grid-sensitive", link.PendingDerivationNote!, StringComparison.OrdinalIgnoreCase);
        Assert.Contains("remains open", link.PendingDerivationNote!, StringComparison.OrdinalIgnoreCase);
        Assert.Contains("DISTINCT off-real-axis", link.PendingDerivationNote!, StringComparison.OrdinalIgnoreCase);
    }

    [Fact]
    public void FragileBridgeSurface_CarriesTheCertifiedCouplingsAndTheChainLength()
    {
        // The threshold's character is certified at J_bridge = 1.0 and 1.9
        // (simulations/fragile_bridge_ep_signature.py, section 1) and recomputed from the C# builders in FragileBridgeThresholdTests, which also
        // checks the located γ* and λ* against the constants this claim carries. Here only the scope
        // the claim states it in: the chain length, and the couplings where the full-L Jordan count
        // was read.
        var link = LocalGlobalEpLink.Build();
        string surface = string.Join("\n", new[]
        {
            link.Name,
            link.DisplayName,
            link.Summary,
            link.PendingDerivationNote,
            link.GlobalInstanceAnchor
        }.Concat(link.Children.Select(child => $"{child.DisplayName}\n{child.Summary}")));

        Assert.Contains("two qubits per chain", link.GlobalInstanceAnchor, StringComparison.Ordinal);
        Assert.Contains("at J_bridge = 1.0 and 1.9", link.GlobalInstanceAnchor, StringComparison.Ordinal);
        Assert.Contains("two qubits per chain", link.PendingDerivationNote, StringComparison.Ordinal);
        Assert.Contains("toy 2x2 EP", surface, StringComparison.OrdinalIgnoreCase);
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
