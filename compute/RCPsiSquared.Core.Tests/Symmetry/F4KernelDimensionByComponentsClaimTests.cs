using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class F4KernelDimensionByComponentsClaimTests
{
    private static F4KernelDimensionByComponentsClaim BuildClaim() =>
        new F4KernelDimensionByComponentsClaim();

    [Fact]
    public void Tier_IsTier1Derived()
    {
        // Promoted from Tier 1 candidate to Tier 1 derived 2026-05-19: the
        // connected-case upper bound is closed by DEGENERACY_PALINDROME Result 2
        // (magnetization conservation), and the multi-component product follows
        // from standard tensor-sum kernel factorisation. See
        // PROOF_F4_KERNEL_DIMENSION_BY_COMPONENTS § "Upper-bound closure".
        Assert.Equal(Tier.Tier1Derived, BuildClaim().Tier);
    }

    [Fact]
    public void Anchor_NamesProofAndDataFiles()
    {
        var claim = BuildClaim();
        Assert.Contains("PROOF_F4_KERNEL_DIMENSION_BY_COMPONENTS", claim.Anchor);
        Assert.Contains("F4StationaryModeCountPi2Inheritance", claim.Anchor);
        Assert.Contains("F1GeneralTopologyVerifiedClaim", claim.Anchor);
        Assert.Contains("f1_n8_n9_metrics", claim.Anchor);
    }

    [Fact]
    public void Anchor_NamesDegeneracyPalindrome()
    {
        // The 2026-05-19 promotion to Tier 1 derived hinges on
        // experiments/DEGENERACY_PALINDROME.md Result 2 closing the connected-case
        // upper bound. Anchor string must surface it for inspection-time traceability.
        Assert.Contains("DEGENERACY_PALINDROME", BuildClaim().Anchor);
    }

    [Fact]
    public void Anchor_NamesProofWeight1Degeneracy()
    {
        // PROOF_WEIGHT1_DEGENERACY.md Appendix 2026-05-17 provides the per-weight
        // ker breakdown that corroborates the boundary upper-bound across
        // chain/ring/star/K_n at N=3..5.
        Assert.Contains("PROOF_WEIGHT1_DEGENERACY", BuildClaim().Anchor);
    }

    [Theory]
    // The 4 N=8 SLOW_N8 anchors (bit-exact from the f1_n8_n9_metrics/<topology>_N8.json
    // KernelDimension fields, commit 89f725e). Components-as-array.
    // Note: chain N=8, ring N=8, star N=8 all have component-spec [8] → 9; they are
    // three distinct topologies sharing the same prediction. The InlineData rows
    // collapse to two distinct component-spec inputs ([8] and [4, 4]); the four
    // labelled topology rows are spelled out in EmpiricalAnchorsN8_HasFourRows_AllMatchPredict
    // below to keep the topology labels (chain/ring/star) and observed values in the
    // record without the xUnit1025 duplicate-InlineData warning.
    [InlineData(new[] { 8 },        9)]    // chain N=8 = ring N=8 = star N=8 (component spec [8])
    [InlineData(new[] { 4, 4 },     25)]   // K_4 + disjoint 4-chain N=8
    public void Predict_MatchesN8Anchor(int[] componentSizes, int expected)
    {
        Assert.Equal(expected, BuildClaim().Predict(componentSizes));
    }

    [Theory]
    // General Π_c (|c|+1) behaviour beyond the N=8 row:
    [InlineData(new[] { 2 },           3)]      // N=2 connected: 3
    [InlineData(new[] { 5 },           6)]      // N=5 connected: 6 (F4 popcount sectors {0..5})
    [InlineData(new[] { 7 },           8)]      // N=7 connected: 8
    [InlineData(new[] { 1, 1 },        4)]      // two isolated qubits: 2·2 = 4
    [InlineData(new[] { 3, 3 },        16)]     // two disjoint 3-chains at N=6: 4·4 = 16
    [InlineData(new[] { 2, 2, 2 },     27)]     // three disjoint pairs at N=6: 3·3·3 = 27
    [InlineData(new[] { 4, 3, 2, 1 },  120)]    // mixed-size 4 components at N=10: 5·4·3·2 = 120
    public void Predict_GeneralComponentSizes(int[] componentSizes, int expected)
    {
        Assert.Equal(expected, BuildClaim().Predict(componentSizes));
    }

    [Fact]
    public void Predict_Null_Throws()
    {
        Assert.Throws<ArgumentNullException>(() => BuildClaim().Predict(null!));
    }

    [Fact]
    public void Predict_EmptyComponents_Throws()
    {
        Assert.Throws<ArgumentException>(() => BuildClaim().Predict(Array.Empty<int>()));
    }

    [Theory]
    [InlineData(new[] { 0 })]
    [InlineData(new[] { -1 })]
    [InlineData(new[] { 3, 0, 2 })]
    [InlineData(new[] { 3, -2 })]
    public void Predict_NonPositiveSize_Throws(int[] componentSizes)
    {
        Assert.Throws<ArgumentException>(() => BuildClaim().Predict(componentSizes));
    }

    [Fact]
    public void EmpiricalAnchorsN8_HasFourRows_AllMatchPredict()
    {
        var claim = BuildClaim();
        Assert.Equal(4, claim.EmpiricalAnchorsN8.Count);
        foreach (var (topology, sizes, predicted, observed) in claim.EmpiricalAnchorsN8)
        {
            // Predicted column equals the in-line Predict computation.
            Assert.Equal(predicted, claim.Predict(sizes));
            // Observed column equals the predicted column (bit-exact at every N=8 row).
            Assert.Equal(predicted, observed);
            Assert.False(string.IsNullOrWhiteSpace(topology));
        }
    }

    [Fact]
    public void EmpiricalAnchorsN8_TopologyLabels_AreDistinctAndN8()
    {
        var claim = BuildClaim();
        var labels = claim.EmpiricalAnchorsN8.Select(a => a.Topology).ToList();
        Assert.Equal(labels.Count, labels.Distinct().Count());
        Assert.All(labels, label => Assert.Contains("N=8", label));
    }

    [Fact]
    public void AnchorDataFiles_HasFourN8JsonPaths()
    {
        var claim = BuildClaim();
        Assert.Equal(4, claim.AnchorDataFiles.Count);
        Assert.All(claim.AnchorDataFiles, path =>
        {
            Assert.Contains("f1_n8_n9_metrics", path);
            Assert.EndsWith("_N8.json", path);
        });
    }

    [Fact]
    public void ConnectedCase_MatchesDegeneracyPalindromeResult2()
    {
        // DEGENERACY_PALINDROME.md Result 2 (April 2026): d_real(0) = N+1 verified
        // bit-exact at N=2..7 via "rmt" Liouvillian eigenvalue export. This test
        // spot-checks the claim's Predict([N]) prediction against the same
        // closed-form across the wider N=1..10 spec range that this typed claim
        // accepts (no overflow risk, sanity check on the connected specialisation).
        var claim = BuildClaim();
        for (int N = 1; N <= 10; N++)
            Assert.Equal(N + 1, claim.Predict(new[] { N }));
    }

    [Theory]
    // DEGENERACY_PALINDROME.md Result 2 (April 2026): d_real(0) = N+1 verified
    // bit-exact at N=2..7 via "rmt" Liouvillian eigenvalue export. This test
    // cross-checks the claim's Predict([N]) prediction against those independently
    // verified N+1 values (the analytic anchor used to close the connected-case
    // upper bound and lift the claim from Tier1Candidate to Tier1Derived).
    [InlineData(2, 3)]
    [InlineData(3, 4)]
    [InlineData(4, 5)]
    [InlineData(5, 6)]
    [InlineData(6, 7)]
    [InlineData(7, 8)]
    public void ConnectedCase_PredictionMatchesObservedFromDegeneracyPalindrome(
        int N, int observedKernelDimFromDegeneracyPalindrome)
    {
        Assert.Equal(observedKernelDimFromDegeneracyPalindrome, BuildClaim().Predict(new[] { N }));
    }
}
