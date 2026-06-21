using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.F86;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Resonance;

using RCPsiSquared.Core.Knowledge;
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
    public void RetractedClaim_StandardList_ContainsThreeKnownRetractions()
    {
        var retracted = RetractedClaim.Standard;
        Assert.Equal(3, retracted.Count);
        Assert.Contains(retracted, r => r.PreviousFormula.Contains("csc(π/(N+1))"));
        Assert.Contains(retracted, r => r.PreviousFormula.Contains("csc(π/5)"));
        // F86a real-axis-EP mechanism, retracted 2026-06-21 (no eigenvalue coalescence on
        // the real Q axis; the Petermann factor is genuine non-normality, not an EP artifact).
        Assert.Contains(retracted, r => r.Name == "F86a real-axis EP mechanism");
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
        Assert.Equal(3, kb.Retracted.Count);

        Assert.Equal(10, kb.QAnchors.Anchors.Count);
        Assert.NotNull(kb.QAnchors.AnkerAt(1.0));
        Assert.NotNull(kb.QAnchors.AnkerAt(1.5));
        Assert.NotNull(kb.QAnchors.AnkerAt(2.0));
        Assert.NotNull(kb.QAnchors.AnkerAt(2.5));

        Assert.Equal(Tier.Tier1Derived, kb.PolarityPairQPeakDecomposition.Tier);
        Assert.Equal(2.5, PolarityPairQPeakDecompositionClaim.EndpointQPeakSchema);
        Assert.Equal(1.5, PolarityPairQPeakDecompositionClaim.InteriorQPeakSchema);

        // Top-level tree: Block + named Q-anchors + Tier1 derived + Tier1 candidate
        // + Tier2 empirical + retracted + open questions + 4-mode note
        IInspectable root = kb;
        Assert.Equal(8, root.Children.Count());
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
        // OrbitKTable is c=2-only; null for c=1 blocks (parallel to C2UniversalShape gating).
        Assert.Null(kb.OrbitKTable);
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

        var openQ = F86OpenQuestions.Standard[0];
        Assert.Contains("Open question", InspectionJsonExporter.ToJson(openQ));
    }

    [Fact]
    public void F86KnowledgeBase_TierGroups_HoldOnlyClaimsOfMatchingTier()
    {
        // Guards the drift where a Claim is appended to the wrong CollectTierX() bucket and
        // then displays under a tier group whose label contradicts the claim's own Tier.
        // Checks the claims placed directly into each tier group; intentional mixed-tier
        // sub-group tables (EP traversal, per-block Q_peak) are coherent units, not recursed.
        var block = new CoherenceBlock(N: 5, n: 1, gammaZero: 0.05); // c=2 exercises every group
        var kb = new F86KnowledgeBase(block);

        var mismatches = new List<string>();
        foreach (var top in kb.Children)
        {
            if (top is not InspectableNode group) continue;
            Tier[]? allowed = group.DisplayName switch
            {
                "Tier 1 (derived)"   => new[] { Tier.Tier1Derived },
                "Tier 1 (candidate)" => new[] { Tier.Tier1Candidate },
                "Tier 2"             => new[] { Tier.Tier2Verified, Tier.Tier2Empirical },
                _ => null,
            };
            if (allowed is null) continue;
            foreach (var child in group.Children)
                if (child is Claim claim && Array.IndexOf(allowed, claim.Tier) < 0)
                    mismatches.Add($"{claim.GetType().Name} is {claim.Tier} under \"{group.DisplayName}\"");
        }

        Assert.True(mismatches.Count == 0,
            "Claim(s) under a tier group that contradicts their own Tier:\n  - "
            + string.Join("\n  - ", mismatches));
    }
}
