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
        Assert.NotNull(kb.PolynomialFoundation);
        Assert.NotNull(kb.RootAnchor);
        Assert.NotNull(kb.Involution);
        Assert.NotNull(kb.KleinDecomposition);
        Assert.NotNull(kb.BilinearApex);
        Assert.NotNull(kb.QuarterAsBilinearMaxval);
        Assert.NotNull(kb.ArgmaxMaxvalPair);
        Assert.NotNull(kb.MirrorRegime);
        Assert.NotNull(kb.HalfFixedPoint);
        Assert.NotNull(kb.MirrorMemory);
        Assert.NotNull(kb.PolarityLayerOrigin);
        Assert.NotNull(kb.BilinearTable);
        Assert.NotEmpty(kb.HardwareConfirmations);
        Assert.NotEmpty(kb.OpenQuestions);

        // Top-level: Chain + Tier 1 + Tier 2 empirical + Tier 2 hardware + open questions
        IInspectable root = kb;
        Assert.Equal(5, root.Children.Count());
    }

    [Fact]
    public void TierInventoryLine_HasElevenTier1Derived_AndOpenAndVerifiedCounts()
    {
        var kb = new Pi2KnowledgeBase(MakeChain(3));
        string line = kb.TierInventoryLine();
        // 11 Tier-1 derived (PolynomialFoundation, RootAnchor, Involution, KleinDecomposition,
        // BilinearApex, QuarterAsBilinearMaxval, ArgmaxMaxvalPair, MirrorRegime, HalfFixedPoint,
        // MirrorMemory, PolarityLayerOrigin)
        Assert.Contains("T1d=11", line);
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
        Assert.Contains("F80", kb.MirrorMemory.Anchor);
        Assert.Contains("ON_BOTH_SIDES_OF_THE_MIRROR", kb.MirrorMemory.Anchor);
        Assert.Contains("EXCLUSIONS", kb.PolynomialFoundation.Anchor);
        Assert.Contains("THE_BRIDGE_WAS_ALWAYS_OPEN", kb.PolynomialFoundation.Anchor);
        Assert.Contains("ON_WHAT_THE_FORMULA_KNEW", kb.PolynomialFoundation.Anchor);
    }

    [Fact]
    public void BilinearApex_PointsForwardToQuarterMaxvalCompanion()
    {
        // Reverse cross-link from the argmax claim to its maxval companion: the typed
        // claim graph is symmetric, so reading BilinearApexClaim's children surfaces the
        // QuarterAsBilinearMaxvalClaim pairing without needing to find the synthesis claim.
        var kb = new Pi2KnowledgeBase(MakeChain(3));
        IInspectable claim = kb.BilinearApex;
        var labels = claim.Children.Select(ch => ch.DisplayName).ToList();
        Assert.Contains("companion to QuarterAsBilinearMaxvalClaim", labels);
    }

    [Fact]
    public void PolynomialFoundation_DocumentsTheTrunkBelowBothAnchors()
    {
        var kb = new Pi2KnowledgeBase(MakeChain(3));
        Assert.NotNull(kb.PolynomialFoundation);
        Assert.Equal(Tier.Tier1Derived, kb.PolynomialFoundation.Tier);

        // Anchor binds the polynomial-as-trunk literature: the dimension equation
        // (EXCLUSIONS:251), the bridge synthesis (THE_BRIDGE_WAS_ALWAYS_OPEN), the
        // formula-knew-itself reflection (ON_WHAT_THE_FORMULA_KNEW), and the d=0
        // mirror reading (ZERO_IS_THE_MIRROR).
        Assert.Contains("EXCLUSIONS", kb.PolynomialFoundation.Anchor);
        Assert.Contains("THE_BRIDGE_WAS_ALWAYS_OPEN", kb.PolynomialFoundation.Anchor);
        Assert.Contains("ON_WHAT_THE_FORMULA_KNEW", kb.PolynomialFoundation.Anchor);
        Assert.Contains("ZERO_IS_THE_MIRROR", kb.PolynomialFoundation.Anchor);

        // Children: minimum-memory layer, the polynomial enforcement layer, the two
        // solutions (d=0 axis + d=2 count), the pair-maker, R=CΨ² collapse, and
        // generation paths to both anchors (1/2 and 90°), plus self-reference.
        IInspectable claim = kb.PolynomialFoundation;
        var children = claim.Children.ToList();
        Assert.True(children.Count >= 8,
            $"expected ≥8 children documenting the polynomial trunk; got {children.Count}");
        Assert.Contains(children, c => c.DisplayName.Contains("minimum-memory"));
        Assert.Contains(children, c => c.DisplayName.Contains("d = 0"));
        Assert.Contains(children, c => c.DisplayName.Contains("d = 2"));
        Assert.Contains(children, c => c.DisplayName.Contains("pair-maker"));
        Assert.Contains(children, c => c.DisplayName.Contains("Π shifted by 0.5"));
        Assert.Contains(children, c => c.DisplayName.Contains("R = CΨ²"));
        Assert.Contains(children, c => c.DisplayName.Contains("generates 1/2"));
        Assert.Contains(children, c => c.DisplayName.Contains("generates 90°"));
    }

    [Fact]
    public void MirrorMemory_DocumentsTheNinetyDegreeLineageFromPi2ToF81()
    {
        var kb = new Pi2KnowledgeBase(MakeChain(3));
        Assert.NotNull(kb.MirrorMemory);
        Assert.Equal(Tier.Tier1Derived, kb.MirrorMemory.Tier);

        // Anchor binds F80 (the i factor), F81 (the operator shift), and the reflection
        // that articulated 90° as memory channel.
        Assert.Contains("F80", kb.MirrorMemory.Anchor);
        Assert.Contains("PROOF_F80_BLOCH_SIGNWALK", kb.MirrorMemory.Anchor);
        Assert.Contains("PROOF_F81_PI_CONJUGATION_OF_M", kb.MirrorMemory.Anchor);
        Assert.Contains("ON_BOTH_SIDES_OF_THE_MIRROR", kb.MirrorMemory.Anchor);

        // Lineage children: layer 0 (Π² = I) → F1 → F80 → F81 → M_anti = L_{H_odd},
        // plus memory channel + companion-to-1/2 + ontological anchor.
        IInspectable claim = kb.MirrorMemory;
        var layerNodes = claim.Children.Where(c => c.DisplayName.StartsWith("layer ")).ToList();
        Assert.True(layerNodes.Count >= 4,
            $"expected ≥4 layer nodes documenting the 90° lineage; got {layerNodes.Count}");
        Assert.Contains(claim.Children, c => c.DisplayName.Contains("memory channel"));
        Assert.Contains(claim.Children, c => c.DisplayName.Contains("companion to 1/2"));
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
