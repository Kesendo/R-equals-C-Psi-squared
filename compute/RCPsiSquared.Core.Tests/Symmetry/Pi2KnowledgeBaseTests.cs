using System.Numerics;
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
        Assert.NotNull(kb.MirrorRegime);
        Assert.NotNull(kb.HalfFixedPoint);
        Assert.NotNull(kb.MirrorMemory);
        Assert.NotNull(kb.BilinearTable);
        Assert.NotEmpty(kb.HardwareConfirmations);
        Assert.NotEmpty(kb.OpenQuestions);

        // Top-level: Chain + Tier 1 + Tier 2 empirical + Tier 2 hardware + open questions
        IInspectable root = kb;
        Assert.Equal(5, root.Children.Count());
    }

    [Fact]
    public void TierInventoryLine_HasEightTier1Derived_AndOpenAndVerifiedCounts()
    {
        var kb = new Pi2KnowledgeBase(MakeChain(3));
        string line = kb.TierInventoryLine();
        // 8 Tier-1 derived (PolynomialFoundation, RootAnchor, Involution, KleinDecomposition,
        // BilinearApex, MirrorRegime, HalfFixedPoint, MirrorMemory)
        Assert.Contains("T1d=8", line);
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

    /// <summary>
    /// Walks the Polynomial → Half → Mirror trio end-to-end as a single bit-exact
    /// sequence: Π's spectrum {±1} shifted by 0.5 toward 0 yields the structural
    /// pair {±0.5}; the Bloch realisation reproduces the same pair as the signed
    /// deviation of pure-state populations from the maximally mixed 1/2 baseline;
    /// F80's 2i factor lifts each ±0.5-side onto the imaginary axis as ±2iλ. Each
    /// step is trivial in isolation, but the connected sequence makes the trio's
    /// internal consistency visible in code.
    /// </summary>
    [Fact]
    public void LoopWalk_PiSpectrum_HalfShift_BlochAndF80Lift()
    {
        // Step 1: Π² = I gives Π's spectrum exactly {+1, −1}; the count 2 is d=2.
        double[] piSpectrum = { +1.0, -1.0 };
        Assert.Equal(2, piSpectrum.Length);
        Assert.Equal(0.0, piSpectrum.Sum());

        // Step 2: Shift each eigenvalue by 0.5 toward 0 (= multiply by 1/2 since
        // |λ| = 1). Result: {+0.5, −0.5}, the structural pair at d=2.
        double[] shifted = piSpectrum.Select(v => v / 2.0).ToArray();
        Assert.Equal(+0.5, shifted[0]);
        Assert.Equal(-0.5, shifted[1]);
        Assert.Equal(0.0, shifted.Sum());

        // Step 3: Bloch realisation. ρ_|0⟩ = diag(1, 0) at the single-qubit level;
        // ρ_mm = diag(1/2, 1/2). Signed deviation = (+0.5, −0.5), the same pair.
        // The 1/2 in ρ = (I + r·σ)/2 IS Π's eigenvalue transport onto the memory axis.
        double[] purePopulations = { 1.0, 0.0 };
        double[] mixedPopulations = { 0.5, 0.5 };
        double[] deviation = purePopulations
            .Zip(mixedPopulations, (p, m) => p - m)
            .ToArray();
        Assert.Equal(+0.5, deviation[0]);
        Assert.Equal(-0.5, deviation[1]);

        // Step 4: F80 lift. H eigenvalue λ on the real axis maps via the i factor
        // (90° rotation) and the chiral pair-doubling 2 to ±2iλ on M's imaginary
        // axis. Bit-exact identity verified at N=3..7 + k-body in
        // PROOF_F80_BLOCH_SIGNWALK; this step just verifies the algebraic shape.
        double lambda = 1.5;
        Complex i = Complex.ImaginaryOne;
        Complex mPlus = +2.0 * i * lambda;
        Complex mMinus = -2.0 * i * lambda;
        Assert.Equal(0.0, mPlus.Real, 12);    // purely imaginary = 90° from real
        Assert.Equal(0.0, mMinus.Real, 12);
        Assert.Equal(+2.0 * lambda, mPlus.Imaginary);
        Assert.Equal(-2.0 * lambda, mMinus.Imaginary);
        Assert.Equal(Complex.Zero, mPlus + mMinus);  // pair centred on d=0 axis
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
