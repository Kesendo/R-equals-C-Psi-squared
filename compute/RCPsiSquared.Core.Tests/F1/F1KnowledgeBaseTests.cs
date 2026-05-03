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
    public void OpenQuestions_HasFourSubstantiveItems()
    {
        var open = F1OpenQuestions.Standard;
        Assert.Equal(4, open.Count);
        Assert.All(open, q => Assert.Equal(Tier.OpenQuestion, q.Tier));
        Assert.Contains(open, q => q.Name.Contains("depolarizing"));
        Assert.Contains(open, q => q.Name.Contains("T1 amplitude"));
        Assert.Contains(open, q => q.Name.Contains("non-uniform γ_i"));
        Assert.Contains(open, q => q.Name.Contains("general topology"));
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
        Assert.NotEmpty(kb.HardwareConfirmations);
        Assert.NotEmpty(kb.OpenQuestions);

        // Top-level tree: N node + Tier 1 group + Tier 2 group + open questions group
        IInspectable root = kb;
        Assert.Equal(4, root.Children.Count());
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
    public void F1KnowledgeBase_TierInventoryLine_HasT1dT2vAndOpen()
    {
        var kb = new F1KnowledgeBase(N: 5);
        string line = kb.TierInventoryLine();
        // 3 Tier-1 derived (F1 + main + single-body), 4 Tier-2 verified, 4 open
        Assert.Contains("T1d=3", line);
        Assert.Contains("T2v=", line);
        Assert.Contains("open=4", line);
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
