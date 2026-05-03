using RCPsiSquared.Core.F71;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Tests.F71;

public class F71KnowledgeBaseTests
{
    [Theory]
    [InlineData(2)]
    [InlineData(3)]
    [InlineData(5)]
    [InlineData(6)]
    public void MirrorOperator_RSquaredIsIdentity_BitExact(int N)
    {
        var op = new F71MirrorOperator(N);
        Assert.Equal(1 << N, op.R.RowCount);
        Assert.True(op.InvolutionResidual < 1e-12,
            $"R²=I residual at N={N}: {op.InvolutionResidual:E3}");
        Assert.Equal(Tier.Tier1Derived, op.Tier);
    }

    [Theory]
    [InlineData(3, 2, 1, false)]   // N=3: 2 bonds, 1 orbit-pair, no self-paired (2 is even)
    [InlineData(4, 3, 2, true)]    // N=4: 3 bonds, 2 orbits incl. self-paired (3 is odd)
    [InlineData(5, 4, 2, false)]   // N=5: 4 bonds, 2 orbit-pairs, no self-paired
    [InlineData(6, 5, 3, true)]    // N=6: 5 bonds, 3 orbits incl. self-paired
    [InlineData(7, 6, 3, false)]   // N=7: 6 bonds, 3 orbit-pairs, no self-paired
    public void BondOrbitDecomposition_OrbitCountAndCentralStatus(int N, int numBonds, int numOrbits, bool hasSelfPaired)
    {
        var d = new F71BondOrbitDecomposition(N);
        Assert.Equal(numBonds, d.NumBonds);
        Assert.Equal(numOrbits, d.NumOrbits);
        Assert.Equal(hasSelfPaired, d.HasSelfPairedCentralOrbit);
        Assert.Equal(numOrbits, d.Orbits.Count);

        // Verify pair structure: orbit b ↔ N−2−b, with one self-paired iff hasSelfPaired
        int selfPairedCount = d.Orbits.Count(o => o.IsSelfPaired);
        Assert.Equal(hasSelfPaired ? 1 : 0, selfPairedCount);
        foreach (var orbit in d.Orbits.Where(o => !o.IsSelfPaired))
            Assert.Equal(numBonds - 1 - orbit.BondA, orbit.BondB!.Value);
    }

    [Fact]
    public void C1MirrorIdentity_IsTier1DerivedWithKinematicAnchor()
    {
        var c1 = new C1MirrorIdentity();
        Assert.Equal(Tier.Tier1Derived, c1.Tier);
        Assert.Contains("F71", c1.Anchor);
        Assert.Contains("PROOF_C1_MIRROR_SYMMETRY", c1.Anchor);
    }

    [Fact]
    public void F86MirrorGeneralisationLink_PointsToBothF71AndF86Anchors()
    {
        var link = new F86MirrorGeneralisationLink();
        Assert.Equal(Tier.Tier1Derived, link.Tier);
        Assert.Contains("F71", link.Anchor);
        Assert.Contains("PROOF_F86_QPEAK", link.Anchor);
    }

    [Fact]
    public void OpenQuestions_HasFourSubstantiveItemsCoveringNonUniformAndAsymmetric()
    {
        var open = F71OpenQuestions.Standard;
        Assert.Equal(4, open.Count);
        Assert.All(open, q => Assert.Equal(Tier.OpenQuestion, q.Tier));
        // Sanity-check the four named axes are covered.
        Assert.Contains(open, q => q.Name.Contains("non-uniform J_b"));
        Assert.Contains(open, q => q.Name.Contains("non-uniform γ_i"));
        Assert.Contains(open, q => q.Name.Contains("asymmetric initial states"));
        Assert.Contains(open, q => q.Name.Contains("per-F71-orbit substructure"));
    }

    [Fact]
    public void F71KnowledgeBase_ForN5_HasAllStructuralChildren()
    {
        var kb = new F71KnowledgeBase(N: 5);
        Assert.Equal(5, kb.N);
        Assert.NotNull(kb.MirrorOperator);
        Assert.NotNull(kb.BondOrbits);
        Assert.NotNull(kb.C1Identity);
        Assert.NotNull(kb.F86Generalisation);
        Assert.NotNull(kb.OpenQuestions);

        // Top-level tree: N node + Tier 1 derived group + open questions group
        IInspectable root = kb;
        Assert.Equal(3, root.Children.Count());
    }

    [Fact]
    public void F71KnowledgeBase_RejectsTooSmallN()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => new F71KnowledgeBase(N: 1));
        Assert.Throws<ArgumentOutOfRangeException>(() => new F71KnowledgeBase(N: 0));
    }

    [Fact]
    public void F71KnowledgeBase_ClaimsAtTier1Derived_ContainsFourClaims()
    {
        var kb = new F71KnowledgeBase(N: 5);
        var t1 = kb.ClaimsAtTier(Tier.Tier1Derived).ToList();
        // F71MirrorOperator + F71BondOrbitDecomposition + C1MirrorIdentity + F86MirrorGeneralisationLink
        Assert.Equal(4, t1.Count);
    }

    [Fact]
    public void F71KnowledgeBase_TierInventoryLine_HasT1dAndOpenCounts()
    {
        var kb = new F71KnowledgeBase(N: 5);
        string line = kb.TierInventoryLine();
        Assert.Contains("T1d=4", line);
        // 4 OpenQuestion items now appear individually in the tree (no wrapper claim).
        Assert.Contains("open=4", line);
    }

    [Fact]
    public void F71KnowledgeBase_AnchorsReferenced_IsNonEmptyAndDistinct()
    {
        var kb = new F71KnowledgeBase(N: 5);
        var anchors = kb.AnchorsReferenced();
        Assert.NotEmpty(anchors);
        // All four Tier-1 claims have distinct anchors.
        Assert.True(anchors.Count >= 4, $"expected ≥4 distinct anchors, got {anchors.Count}");
    }

    [Fact]
    public void F71KnowledgeBase_JsonExport_PreservesTierLabelsAndAnchors()
    {
        var kb = new F71KnowledgeBase(N: 4);
        string json = InspectionJsonExporter.ToJson(kb);
        Assert.Contains("Tier 1 (derived)", json);
        Assert.Contains("PROOF_C1_MIRROR_SYMMETRY", json);
        // F71MirrorOperator.DisplayName is "R (N=4, dim 16)"; check the dim-16 marker.
        Assert.Contains("R (N=4, dim 16)", json);
    }
}
