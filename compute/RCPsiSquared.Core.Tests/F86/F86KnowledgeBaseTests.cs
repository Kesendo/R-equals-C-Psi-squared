using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.F86;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Resonance;

namespace RCPsiSquared.Core.Tests.F86;

public class F86KnowledgeBaseTests
{
    [Fact]
    public void TPeakLaw_ComputesOneOver4Gamma_AndIsTier1Derived()
    {
        var law = new TPeakLaw(gammaZero: 0.05);
        Assert.Equal(5.0, law.Value, 12);
        Assert.Equal(Tier.Tier1Derived, law.Tier);
    }

    [Fact]
    public void TPeakLaw_RejectsZeroGamma()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => new TPeakLaw(0.0));
    }

    [Theory]
    [InlineData(2.0, 1.0)]
    [InlineData(2.8284, 0.7071)]
    [InlineData(4.0, 0.5)]
    public void QEpLaw_ComputesTwoOverGEff(double gEff, double expectedQEp)
    {
        var law = new QEpLaw(gEff);
        Assert.Equal(expectedQEp, law.Value, 3);
    }

    [Fact]
    public void TwoLevelEpModel_PreEp_AtEp_PostEp_RegimeSwitch()
    {
        // J·g_eff < 2γ₀ → PreEp; J·g_eff = 2γ₀ → AtEp; J·g_eff > 2γ₀ → PostEp
        const double γ0 = 0.05;
        const double gEff = 2.0;
        // EP at J = 2γ₀/g_eff = 0.05.
        var pre = new TwoLevelEpModel(γ0, j: 0.025, gEff);
        var at = new TwoLevelEpModel(γ0, j: 0.05, gEff);
        var post = new TwoLevelEpModel(γ0, j: 0.075, gEff);

        Assert.Equal(EpRegime.PreEp, pre.Regime);
        Assert.Equal(EpRegime.AtEp, at.Regime);
        Assert.Equal(EpRegime.PostEp, post.Regime);

        Assert.True(pre.Discriminant > 0);
        Assert.True(Math.Abs(at.Discriminant) < 1e-10);
        Assert.True(post.Discriminant < 0);
    }

    [Fact]
    public void UniversalShapePrediction_CompareTo_FlagsWithinAndOutside()
    {
        var prediction = new UniversalShapePrediction(
            BondClass.Interior,
            expectedRatio: 0.756,
            tolerance: 0.005,
            witnesses: Array.Empty<UniversalShapeWitness>());

        var within = new PeakResult(QPeak: 1.5, KMax: 0.1, HwhmLeft: 0.756 * 1.5, HwhmRight: null);
        var outside = new PeakResult(QPeak: 1.5, KMax: 0.1, HwhmLeft: 0.74 * 1.5, HwhmRight: null);
        var noHwhm = new PeakResult(QPeak: 1.5, KMax: 0.1, HwhmLeft: null, HwhmRight: null);

        Assert.True(prediction.CompareTo(within).Within);
        Assert.False(prediction.CompareTo(outside).Within);
        Assert.False(prediction.CompareTo(noHwhm).Within);
        Assert.Null(prediction.CompareTo(noHwhm).Actual);
    }

    [Fact]
    public void RetractedClaim_StandardList_ContainsTwoKnownRetractions()
    {
        var retracted = RetractedClaim.Standard;
        Assert.Equal(2, retracted.Count);
        Assert.Contains(retracted, r => r.PreviousFormula.Contains("csc(π/(N+1))"));
        Assert.Contains(retracted, r => r.PreviousFormula.Contains("csc(π/5)"));
        Assert.All(retracted, r => Assert.Equal(Tier.Retracted, r.Tier));
    }

    [Fact]
    public void F86KnowledgeBase_ForC2N5_HasAllStructuralChildren()
    {
        var block = new CoherenceBlock(N: 5, n: 1, gammaZero: 0.05);
        var kb = new F86KnowledgeBase(block);

        Assert.Equal(0.05, kb.Block.GammaZero);
        Assert.Equal(5.0, kb.TPeak.Value, 12);
        Assert.NotNull(kb.QEp);
        Assert.NotNull(kb.Svd);
        Assert.Equal(3, kb.EpTraversal.Count);
        Assert.Equal(BondClass.Interior, kb.InteriorShape.BondClass);
        Assert.Equal(BondClass.Endpoint, kb.EndpointShape.BondClass);
        Assert.Equal(2, kb.Retracted.Count);

        // Top-level tree: Block + Tier1 derived + Tier1 candidate + Tier2 empirical
        // + retracted + open questions + 4-mode note
        IInspectable root = kb;
        Assert.Equal(7, root.Children.Count());
    }

    [Fact]
    public void F86KnowledgeBase_ForC1Block_OmitsInterChannelMachinery()
    {
        // c=1 has only one HD channel (HD=1), so no inter-channel SVD.
        var block = new CoherenceBlock(N: 3, n: 0, gammaZero: 0.05);
        var kb = new F86KnowledgeBase(block);

        Assert.Equal(1, block.C);
        Assert.Null(kb.Svd);
        Assert.Null(kb.QEp);
        Assert.Empty(kb.EpTraversal);
        // t_peak still applies (universal in γ₀ alone).
        Assert.Equal(5.0, kb.TPeak.Value, 12);
    }

    [Fact]
    public void F86KnowledgeBase_CompareToFullScan_FlagsN5InteriorOutsideTolerance()
    {
        // Documented finding 2026-05-02: c=2 N=5 Interior HWHM/Q ≈ 0.7455 is the smallest-N
        // witness, just below the universal 0.756 ± 0.005 tolerance (finite-size effect).
        // Endpoint at N=5 hits 0.770 exactly. Both behaviors are visible via the typed API.
        var block = new CoherenceBlock(N: 5, n: 1, gammaZero: 0.05);
        var kb = new F86KnowledgeBase(block);
        var measuredCurve = new ResonanceScan(block).ComputeKCurve();
        var matches = kb.CompareTo(measuredCurve);

        Assert.Equal(2, matches.Count);
        var interior = matches[0];
        var endpoint = matches[1];
        Assert.Equal(BondClass.Interior, interior.Prediction.BondClass);
        Assert.Equal(BondClass.Endpoint, endpoint.Prediction.BondClass);

        Assert.False(interior.Within,
            $"c=2 N=5 Interior HWHM/Q ({interior.Actual:F4}) should fall OUTSIDE the ±0.005 tolerance " +
            $"around 0.756 — finite-size effect at smallest N.");
        Assert.True(endpoint.Within,
            $"c=2 N=5 Endpoint HWHM/Q ({endpoint.Actual:F4}) should fall WITHIN ±0.005 of 0.770.");
    }

    [Fact]
    public void F86KnowledgeBase_JsonExport_PreservesTierLabelsAndAnchors()
    {
        // Use individual claim JSON exports to avoid triggering all witness computations
        // (which would scan ~18 (c, N) cases up to c=4 N=8, taking minutes).
        var tPeak = new TPeakLaw(0.05);
        var qEp = new QEpLaw(2.83);
        Assert.Contains("Tier 1 (derived)", InspectionJsonExporter.ToJson(tPeak));
        Assert.Contains("docs/ANALYTICAL_FORMULAS.md F86 Statement 1", InspectionJsonExporter.ToJson(qEp));

        var retracted = RetractedClaim.Standard[0];
        var retractedJson = InspectionJsonExporter.ToJson(retracted);
        Assert.Contains("Retracted", retractedJson);
        Assert.Contains("csc(π/(N+1))", retractedJson);

        var openQ = OpenQuestion.Standard[0];
        Assert.Contains("Open question", InspectionJsonExporter.ToJson(openQ));
    }
}
