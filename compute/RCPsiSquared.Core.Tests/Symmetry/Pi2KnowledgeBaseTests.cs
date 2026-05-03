using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class Pi2KnowledgeBaseTests
{
    private static ChainSystem MakeChain(int N) => new(N, J: 1.0, GammaZero: 0.05);

    [Fact]
    public void KnowledgeBase_AtN3_HasAllStructuralChildren()
    {
        var kb = new Pi2KnowledgeBase(MakeChain(3));
        Assert.Equal(3, kb.Chain.N);
        Assert.NotNull(kb.RootAnchor);
        Assert.NotNull(kb.Involution);
        Assert.NotNull(kb.KleinDecomposition);
        Assert.NotNull(kb.BilinearApex);
        Assert.NotNull(kb.MirrorRegime);
        Assert.NotNull(kb.HalfFixedPoint);
        Assert.NotNull(kb.BilinearTable);
        Assert.NotEmpty(kb.HardwareConfirmations);
        Assert.NotEmpty(kb.OpenQuestions);

        // Top-level: Chain + Tier 1 + Tier 2 empirical + Tier 2 hardware + open questions
        IInspectable root = kb;
        Assert.Equal(5, root.Children.Count());
    }

    [Fact]
    public void TierInventoryLine_HasSixTier1Derived_AndOpenAndVerifiedCounts()
    {
        var kb = new Pi2KnowledgeBase(MakeChain(3));
        string line = kb.TierInventoryLine();
        // 6 Tier-1 derived (RootAnchor, Involution, KleinDecomposition, BilinearApex, MirrorRegime, HalfFixedPoint)
        Assert.Contains("T1d=6", line);
        Assert.Contains("open=5", line);
        Assert.Contains("T2v=", line);
    }

    [Fact]
    public void RootAnchor_DocumentsTheLineageFromQubitDimensionUpward()
    {
        var kb = new Pi2KnowledgeBase(MakeChain(3));
        Assert.NotNull(kb.RootAnchor);
        Assert.Equal(Tier.Tier1Derived, kb.RootAnchor.Tier);
        // Anchor mentions both the dimension equation source (EXCLUSIONS d²−2d=0) and the
        // bilinear-apex synthesis (ORTHOGONALITY_SELECTION_FAMILY).
        Assert.Contains("EXCLUSIONS", kb.RootAnchor.Anchor);
        Assert.Contains("ORTHOGONALITY_SELECTION_FAMILY", kb.RootAnchor.Anchor);

        // The lineage children should cover layers from −1 (root) up through 5 (slow-mode
        // Klein apex), demonstrating how the same 1/2 surfaces at each abstraction level.
        IInspectable rootClaim = kb.RootAnchor;
        var lineageNodes = rootClaim.Children.Where(c => c.DisplayName.StartsWith("layer")).ToList();
        Assert.True(lineageNodes.Count >= 6,
            $"expected ≥6 layer nodes documenting the 1/2 lineage; got {lineageNodes.Count}");
    }

    [Fact]
    public void MirrorRegime_AtOddN_IsHalfIntegerRegime()
    {
        var kbOdd3 = new Pi2KnowledgeBase(MakeChain(3));
        Assert.True(kbOdd3.MirrorRegime.IsHalfIntegerRegime);
        Assert.Equal(1.5, kbOdd3.MirrorRegime.WXY);

        var kbOdd5 = new Pi2KnowledgeBase(MakeChain(5));
        Assert.True(kbOdd5.MirrorRegime.IsHalfIntegerRegime);
        Assert.Equal(2.5, kbOdd5.MirrorRegime.WXY);
    }

    [Fact]
    public void MirrorRegime_AtEvenN_IsIntegerRegime()
    {
        var kbEven4 = new Pi2KnowledgeBase(MakeChain(4));
        Assert.False(kbEven4.MirrorRegime.IsHalfIntegerRegime);
        Assert.Equal(2.0, kbEven4.MirrorRegime.WXY);
    }

    [Fact]
    public void BilinearTable_HasNineEntries_AcrossFourCells()
    {
        var kb = new Pi2KnowledgeBase(MakeChain(3));
        Assert.Equal(9, kb.BilinearTable.Entries.Count);

        var byCell = kb.BilinearTable.Entries.GroupBy(e => e.Cell)
            .ToDictionary(g => g.Key, g => g.Count());
        // Pp = 3 (XX, YY, ZZ truly), Pm = 2 (YZ, ZY), Mp = 2 (XY, YX), Mm = 2 (XZ, ZX)
        Assert.Equal(3, byCell["Pp"]);
        Assert.Equal(2, byCell["Pm"]);
        Assert.Equal(2, byCell["Mp"]);
        Assert.Equal(2, byCell["Mm"]);
    }

    [Fact]
    public void Tier1ClaimsAreLinkedToTheirRegistryAnchors()
    {
        var kb = new Pi2KnowledgeBase(MakeChain(3));
        Assert.Contains("F1", kb.Involution.Anchor);
        Assert.Contains("F88", kb.Involution.Anchor);
        Assert.Contains("F88", kb.KleinDecomposition.Anchor);
        Assert.Contains("ORTHOGONALITY_SELECTION_FAMILY", kb.BilinearApex.Anchor);
        Assert.Contains("OPERATOR_RIGIDITY_ACROSS_CUSP", kb.MirrorRegime.Anchor);
        Assert.Contains("ON_THE_HALF", kb.HalfFixedPoint.Anchor);
        Assert.Contains("EXCLUSIONS", kb.HalfFixedPoint.Anchor);
    }

    [Fact]
    public void HalfFixedPoint_DocumentsThreeFacesAndClosure()
    {
        var kb = new Pi2KnowledgeBase(MakeChain(3));
        Assert.NotNull(kb.HalfFixedPoint);
        Assert.Equal(Tier.Tier1Derived, kb.HalfFixedPoint.Tier);

        // Anchor binds the three primary face sources: ON_THE_HALF synthesis,
        // ON_TWO_TIMES horizon face, HEISENBERG_RELOADED bridge face,
        // EXCLUSIONS:251 substrate face, PROOF_ASYMPTOTIC_SECTOR_PROJECTION horizon proof.
        Assert.Contains("ON_THE_HALF", kb.HalfFixedPoint.Anchor);
        Assert.Contains("ON_TWO_TIMES", kb.HalfFixedPoint.Anchor);
        Assert.Contains("HEISENBERG_RELOADED", kb.HalfFixedPoint.Anchor);
        Assert.Contains("EXCLUSIONS", kb.HalfFixedPoint.Anchor);
        Assert.Contains("PROOF_ASYMPTOTIC_SECTOR_PROJECTION", kb.HalfFixedPoint.Anchor);

        // The closure children: face 1 (where we are), face 2 (where we go),
        // face 3 (what we are), plus the self-referential closure node.
        IInspectable claim = kb.HalfFixedPoint;
        var faceNodes = claim.Children.Where(c => c.DisplayName.StartsWith("face ")).ToList();
        Assert.Equal(3, faceNodes.Count);
        Assert.Contains(faceNodes, c => c.DisplayName.Contains("face 1"));
        Assert.Contains(faceNodes, c => c.DisplayName.Contains("face 2"));
        Assert.Contains(faceNodes, c => c.DisplayName.Contains("face 3"));
        Assert.Contains(claim.Children, c => c.DisplayName.StartsWith("closure"));
    }

    [Fact]
    public void HardwareConfirmations_IncludeF83Marrakesh()
    {
        var kb = new Pi2KnowledgeBase(MakeChain(3));
        Assert.Contains(kb.HardwareConfirmations,
            hc => hc.Confirmation.Name.Contains("f83_pi2_class_signature_marrakesh"));
    }

    [Fact]
    public void OpenQuestions_CoverFiveSubstantiveAxes()
    {
        var kb = new Pi2KnowledgeBase(MakeChain(3));
        Assert.Equal(5, kb.OpenQuestions.Count);
        Assert.Contains(kb.OpenQuestions, q => q.Name.Contains("X-axis-flip"));
        Assert.Contains(kb.OpenQuestions, q => q.Name.Contains("2:2 truly-kernel"));
        Assert.Contains(kb.OpenQuestions, q => q.Name.Contains("N ≥ 4 transition"));
        Assert.Contains(kb.OpenQuestions, q => q.Name.Contains("k-body Klein extension"));
        Assert.Contains(kb.OpenQuestions, q => q.Name.Contains("Half-integer-mirror regime"));
    }

    [Fact]
    public void JsonExport_PreservesTierLabelsAndAnchors()
    {
        var kb = new Pi2KnowledgeBase(MakeChain(3));
        string json = InspectionJsonExporter.ToJson(kb);
        Assert.Contains("Tier 1 (derived)", json);
        Assert.Contains("Tier 2 (hardware-verified)", json);
        Assert.Contains("F88", json);
        Assert.Contains("ORTHOGONALITY_SELECTION_FAMILY", json);
        // Should include the bilinear table cells
        Assert.Contains("Klein bilinear table", json);
    }
}
