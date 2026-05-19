using RCPsiSquared.Core.Confirmations;
using RCPsiSquared.Core.F1;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Lindblad;

namespace RCPsiSquared.Core.Tests.F1;

public class F1KnowledgeBaseTests
{
    [Fact]
    public void PalindromeIdentity_IsTier1Derived_WithMirrorProofAnchor()
    {
        var f1 = new F1PalindromeIdentity();
        Assert.Equal(Tier.Tier1Derived, f1.Tier);
        Assert.Contains("F1", f1.Anchor);
        Assert.Contains("MIRROR_SYMMETRY_PROOF", f1.Anchor);
    }

    [Theory]
    [InlineData(3, HamiltonianClass.Main, 8.0)]           // (N−1)·4^(N−2) = 2·4 = 8
    [InlineData(4, HamiltonianClass.Main, 48.0)]          // 3·16 = 48
    [InlineData(5, HamiltonianClass.Main, 256.0)]         // 4·64 = 256
    [InlineData(3, HamiltonianClass.SingleBody, 12.0)]    // (2N−3)·4^(N−2) = 3·4 = 12
    [InlineData(4, HamiltonianClass.SingleBody, 80.0)]    // 5·16 = 80
    [InlineData(5, HamiltonianClass.SingleBody, 448.0)]   // 7·64 = 448
    public void Scaling_ChainDefault_MatchesClosedForm(int N, HamiltonianClass cls, double expected)
    {
        var claim = new PalindromeResidualScalingClaim(N, cls);
        Assert.Equal(expected, claim.Factor, 9);
        Assert.Equal(Tier.Tier1Derived, claim.Tier);
    }

    [Fact]
    public void Scaling_GraphAware_RequiresBothOrNeither()
    {
        Assert.Throws<InvalidOperationException>(() =>
            _ = new PalindromeResidualScalingClaim(5, HamiltonianClass.Main, bondCount: 4).Factor);
    }

    [Fact]
    public void Scaling_AdjacentRatio_MatchesClosedForm()
    {
        var main4 = new PalindromeResidualScalingClaim(4, HamiltonianClass.Main);
        // ratio at N=4 → 4·N/(N−1) = 16/3
        Assert.Equal(16.0 / 3.0, main4.AdjacentRatio, 9);

        var single4 = new PalindromeResidualScalingClaim(4, HamiltonianClass.SingleBody);
        // ratio at N=4 → 4·(2N−1)/(2N−3) = 28/5
        Assert.Equal(28.0 / 5.0, single4.AdjacentRatio, 9);
    }

    [Fact]
    public void HardwareConfirmation_WrapsRegistryEntry()
    {
        var conf = ConfirmationsRegistry.Lookup("palindrome_trichotomy");
        Assert.NotNull(conf);
        var claim = new HardwareConfirmationClaim(conf!);
        Assert.Equal(Tier.Tier2Verified, claim.Tier);
        Assert.Contains("palindrome_trichotomy", claim.Name);
        Assert.Contains("ibm_marrakesh", claim.Summary);
    }

    [Fact]
    public void OpenQuestions_IsEmpty()
    {
        // All four items closed 2026-05-18: T1 closed form (Tier 1), depol closed form
        // (Tier 1), non-uniform γ (negative-result closure), general topology
        // (synthesis proof + Tier-2 verification record). See F1OpenQuestions XML doc
        // for the per-item closure references (proof markdown + typed claim). First
        // time the F1 family's open-question collection is empty.
        var open = F1OpenQuestions.Standard;
        Assert.Empty(open);
    }

    [Fact]
    public void F1KnowledgeBase_ForN5_ChainDefault_HasAllStructuralChildren()
    {
        var kb = new F1KnowledgeBase(N: 5);
        Assert.Equal(5, kb.N);
        Assert.Null(kb.BondCount);
        Assert.NotNull(kb.PalindromeIdentity);
        Assert.NotNull(kb.MainScaling);
        Assert.NotNull(kb.SingleBodyScaling);
        Assert.NotNull(kb.T1ResidualClosedForm);
        Assert.NotNull(kb.T1ResidualPi2Decomposition);
        Assert.NotNull(kb.DepolResidualClosedForm);
        Assert.NotNull(kb.F49NonUniformCrossTerm);
        Assert.NotNull(kb.GeneralTopologyVerification);
        Assert.NotNull(kb.KernelDimensionByComponents);
        Assert.NotEmpty(kb.HardwareConfirmations);
        // OpenQuestions is empty as of 2026-05-18 (all four F1 items closed); see
        // F1OpenQuestions XML doc for the per-item closure references.
        Assert.Empty(kb.OpenQuestions);

        // Top-level tree: N node + Tier 1 (derived) group + Tier 2 (verified) group
        // + open questions group. The Tier 1 (candidate) group (added 2026-05-18 with
        // the F4 kernel-dim-by-components bridge) was removed 2026-05-19 when that
        // bridge was promoted to Tier 1 derived (DEGENERACY_PALINDROME Result 2
        // closure of the connected-case upper bound).
        IInspectable root = kb;
        Assert.Equal(4, root.Children.Count());
    }

    [Fact]
    public void F1KnowledgeBase_GeneralTopologyVerification_IsTier2Verified()
    {
        var kb = new F1KnowledgeBase(N: 5);
        Assert.Equal(Tier.Tier2Verified, kb.GeneralTopologyVerification.Tier);
        Assert.Contains("PROOF_F1_GENERAL_TOPOLOGY", kb.GeneralTopologyVerification.Anchor);
        // Spot-check the verified-N set and the topology counts stayed stable across
        // the typed-claim round-trip.
        // N=9 chain was crossed 2026-05-19 via the MklDirect bridge (commit abb2d52);
        // the new frontier at N=10 is memory-pressure rather than LP64 marshalling.
        // See ScaleFrontierBlockedAtN below.
        Assert.Equal(new[] { 5, 6, 7, 8, 9 }, kb.GeneralTopologyVerification.VerifiedNValues);
        Assert.Equal(new[] { 5, 6, 7, 8, 9 }, kb.GeneralTopologyVerification.ScaleUpToN);
        Assert.Equal(10, kb.GeneralTopologyVerification.ScaleFrontierBlockedAtN);
        Assert.False(string.IsNullOrWhiteSpace(kb.GeneralTopologyVerification.ScaleFrontierBlockerReason));
        Assert.True(kb.GeneralTopologyVerification.DisconnectedComponentsVerified);
        Assert.True(kb.GeneralTopologyVerification.WeightedEdgesVerified);
        Assert.True(kb.GeneralTopologyVerification.SingleBodyClassVerified);
        Assert.NotEmpty(kb.GeneralTopologyVerification.SpectrumMetricsDataFiles);
        // chain_N9 landed 2026-05-19 alongside the MklDirect bridge (commit abb2d52);
        // chain_N8 is one of the 4 captured N=8 sweeps.
        Assert.Contains("chain_N8.json", string.Join(";", kb.GeneralTopologyVerification.SpectrumMetricsDataFiles));
        Assert.Contains("chain_N9.json", string.Join(";", kb.GeneralTopologyVerification.SpectrumMetricsDataFiles));
    }

    [Fact]
    public void F1KnowledgeBase_GraphAware_PassesBondCountIntoScaling()
    {
        // Star at N=5: B=4 bonds, D2 = 4² + 4·1² = 20.
        var kb = new F1KnowledgeBase(N: 5, bondCount: 4, degreeSquaredSum: 20);
        // main F = B · 4^(N−2) = 4 · 64 = 256 (same as chain at N=5 because B happens to be same)
        Assert.Equal(256.0, kb.MainScaling.Factor, 9);
        // single-body F = (D2 / 2) · 4^(N−2) = 10 · 64 = 640
        Assert.Equal(640.0, kb.SingleBodyScaling.Factor, 9);
    }

    [Fact]
    public void F1KnowledgeBase_RejectsTooSmallN()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => new F1KnowledgeBase(N: 1));
    }

    [Fact]
    public void F1KnowledgeBase_RejectsHalfGraphArguments()
    {
        Assert.Throws<ArgumentException>(() => new F1KnowledgeBase(N: 5, bondCount: 4));
        Assert.Throws<ArgumentException>(() => new F1KnowledgeBase(N: 5, degreeSquaredSum: 20));
    }

    [Fact]
    public void F1KnowledgeBase_TierInventoryLine_HasT1dT2vAndNoOpen()
    {
        var kb = new F1KnowledgeBase(N: 5);
        string line = kb.TierInventoryLine();
        // 8 Tier-1 derived (was 7 until 2026-05-19): F1 + main + single-body + T1 closed
        // form + T1 Π²-decomposition + depol + F49 non-uniform γ cross-term +
        // F4 kernel-dim-by-components (promoted from Tier 1 candidate to Tier 1 derived
        // 2026-05-19 after DEGENERACY_PALINDROME Result 2 was identified as the closure
        // of the connected-case upper bound).
        // Tier-1 candidate: ZERO after the 2026-05-19 promotion. TierInventoryLine skips
        // empty tiers, so "T1c=" should not appear in the line at all.
        // Tier-2 verified bumped from N hardware confirmations (3) to N+1 with the new
        // F1GeneralTopologyVerifiedClaim general-topology verification record (2026-05-18).
        // Open questions: ZERO after the 2026-05-18 general-topology closure (last F1
        // OpenQuestion resolved via PROOF_F1_GENERAL_TOPOLOGY.md + F1GeneralTopologyVerifiedClaim).
        // TierInventoryLine skips empty tiers, so "open=" should not appear in the line at all.
        Assert.Contains("T1d=8", line);
        Assert.DoesNotContain("T1c=", line);
        Assert.Contains("T2v=", line);
        Assert.DoesNotContain("open=", line);
    }

    [Fact]
    public void F1KnowledgeBase_KernelDimensionByComponents_IsTier1Derived_AndBridgeAnchored()
    {
        // Wire-in spot check for the F4 disconnected-graph bridge: the F1 KB exposes
        // it as a top-level property; its tier is Tier 1 derived (promoted 2026-05-19
        // after DEGENERACY_PALINDROME Result 2 was identified as the closure of the
        // connected-case upper bound); the anchor mentions the proof file, the
        // DEGENERACY_PALINDROME closure source, PROOF_WEIGHT1_DEGENERACY corroborating
        // appendix, and the F1 SLOW_N8 metric data files where the four kernel-dim
        // numbers live.
        var kb = new F1KnowledgeBase(N: 5);
        Assert.Equal(Tier.Tier1Derived, kb.KernelDimensionByComponents.Tier);
        Assert.Contains("PROOF_F4_KERNEL_DIMENSION_BY_COMPONENTS", kb.KernelDimensionByComponents.Anchor);
        Assert.Contains("DEGENERACY_PALINDROME", kb.KernelDimensionByComponents.Anchor);
        Assert.Contains("PROOF_WEIGHT1_DEGENERACY", kb.KernelDimensionByComponents.Anchor);
        Assert.Contains("f1_n8_n9_metrics", kb.KernelDimensionByComponents.Anchor);
        // 4 N=8 anchors bit-exact: chain/ring/star = 9, K_4+disjoint-4-chain = 25.
        Assert.Equal(9, kb.KernelDimensionByComponents.Predict(new[] { 8 }));
        Assert.Equal(25, kb.KernelDimensionByComponents.Predict(new[] { 4, 4 }));
    }

    [Fact]
    public void F1KnowledgeBase_JsonExport_PreservesTierLabelsAndAnchors()
    {
        var kb = new F1KnowledgeBase(N: 4);
        string json = InspectionJsonExporter.ToJson(kb);
        Assert.Contains("Tier 1 (derived)", json);
        Assert.Contains("Tier 2 (hardware-verified)", json);
        Assert.Contains("MIRROR_SYMMETRY_PROOF", json);
        Assert.Contains("ibm_marrakesh", json);
    }
}
