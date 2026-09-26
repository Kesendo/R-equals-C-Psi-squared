using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using Xunit;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class QuditProductMirrorCapTests
{
    private static (QuditPartialPalindromeCeiling F121, QubitNecessityPi2Inheritance QubitNecessity) MakeParents()
    {
        var qubitNecessity = new QubitNecessityPi2Inheritance(
            new Pi2DyadicLadderClaim(), new Pi2OperatorSpaceMirrorClaim());
        return (new QuditPartialPalindromeCeiling(qubitNecessity), qubitNecessity);
    }

    private static QuditProductMirrorCap MakeClaim()
    {
        var (f121, qubitNecessity) = MakeParents();
        return new QuditProductMirrorCap(f121, qubitNecessity);
    }

    // ------------------------------------------------------------------
    // Claim metadata
    // ------------------------------------------------------------------

    [Fact]
    public void Tier_IsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived, MakeClaim().Tier);
    }

    [Fact]
    public void Anchor_ReferencesProofAndVerifier()
    {
        var claim = MakeClaim();
        Assert.Contains("PROOF_QUDIT_PARTIAL_PALINDROME.md", claim.Anchor);
        Assert.Contains("qudit_product_mirror_cap.py", claim.Anchor);
    }

    [Fact]
    public void TypedParents_AreWired()
    {
        var claim = MakeClaim();
        Assert.NotNull(claim.PartialPalindrome);
        Assert.NotNull(claim.QubitNecessity);
        // The F121 parent carries its own QubitNecessity parent (the typed chain stays intact).
        Assert.NotNull(claim.PartialPalindrome.QubitNecessity);
    }

    [Fact]
    public void Constructor_RejectsNullParents()
    {
        var (f121, qubitNecessity) = MakeParents();
        Assert.Throws<ArgumentNullException>(() => new QuditProductMirrorCap(null!, qubitNecessity));
        Assert.Throws<ArgumentNullException>(() => new QuditProductMirrorCap(f121, null!));
    }

    [Fact]
    public void NotAnIZ2AxisClaim_CubeMapCountsUnchanged()
    {
        // Pinned invariant: the cap is cross-axis structural like AntilinearTriangleClaim
        // and MomentTowerPumpChannelClaim; PolarityCubeMap counts stay untouched.
        Assert.False(typeof(IZ2AxisClaim).IsAssignableFrom(typeof(QuditProductMirrorCap)));
    }

    // ------------------------------------------------------------------
    // Self-check battery (exact integer/permutation arithmetic, built in the ctor)
    // ------------------------------------------------------------------

    [Fact]
    public void Battery_AllCasesPass()
    {
        var claim = MakeClaim();
        Assert.Equal(11, claim.Cases.Count);
        foreach (var c in claim.Cases)
            Assert.True(c.Passes, $"battery case '{c.Name}' failed: expected {c.Expected}, got {c.Actual}");
        Assert.Equal(claim.Cases.Count, claim.PassCount);
    }

    [Fact]
    public void Battery_AlignedSubspace_IsExact()
    {
        var claim = MakeClaim();
        var aligned = claim.Cases.Single(c => c.Name.StartsWith("Π_d exact", StringComparison.Ordinal));
        Assert.True(aligned.Passes, $"aligned-subspace exactness failed: {aligned.Actual}");
    }

    [Fact]
    public void Battery_GroupLaw_Holds()
    {
        var claim = MakeClaim();
        var group = claim.Cases.Single(c => c.Name.StartsWith("mirror group law", StringComparison.Ordinal));
        Assert.True(group.Passes, $"group law failed: {group.Actual}");
    }

    [Fact]
    public void Summary_CarriesTheCapAndTheGroup()
    {
        var claim = MakeClaim();
        Assert.Contains("P(d, N)", claim.Summary);
        Assert.Contains("Z_d ≀ Z₂", claim.Summary);
        Assert.Contains($"{claim.Cases.Count}/{claim.Cases.Count} battery PASS", claim.Summary);
    }

    // ------------------------------------------------------------------
    // The product lemma: P(d, N) = max_m (2d)^(N−2m)·(d³ − d²)^m
    // ------------------------------------------------------------------

    [Fact]
    public void ProductCap_IsTheSwapRankUpToFive_AndAboveItFromSix()
    {
        for (int d = 2; d <= 5; d++)
            for (int n = 1; n <= 6; n++)
                Assert.Equal(QuditProductMirrorCap.ShiftAlignedRank(d, n), QuditProductMirrorCap.ProductCap(d, n));
        Assert.Equal(180L, QuditProductMirrorCap.ProductCap(6, 2));      // d(d² − d) = 6·30
        Assert.Equal(144L, QuditProductMirrorCap.ShiftAlignedRank(6, 2));
        Assert.Equal(12L, QuditProductMirrorCap.ProductCap(6, 1));       // N = 1: only the swap
        Assert.Equal(2160L, QuditProductMirrorCap.ProductCap(6, 3));     // 2d · (d³ − d²)
        Assert.Equal(294L, QuditProductMirrorCap.ProductCap(7, 2));
        Assert.Equal(448L, QuditProductMirrorCap.ProductCap(8, 2));
    }

    [Fact]
    public void ProductCap_IsFullExactlyAtTheQubit_ForEveryD()
    {
        for (int d = 2; d <= 9; d++)
            for (int n = 1; n <= 4; n++)
                Assert.Equal(d == 2, QuditProductMirrorCap.ProductCap(d, n) == QuditPartialPalindromeCeiling.Total(d, n));
    }

    [Fact]
    public void GradeRanks_AreDTwoDAndDSquaredMinusD()
    {
        Assert.Equal(3L, QuditProductMirrorCap.GradeRank(3, 0));
        Assert.Equal(6L, QuditProductMirrorCap.GradeRank(3, 1));
        Assert.Equal(6L, QuditProductMirrorCap.GradeRank(3, 2));
        Assert.Equal(30L, QuditProductMirrorCap.GradeRank(6, 2));
        Assert.Throws<ArgumentOutOfRangeException>(() => QuditProductMirrorCap.GradeRank(3, 3));
    }

    [Fact]
    public void Battery_ProductLemma_IsComputedNotAsserted()
    {
        var claim = MakeClaim();
        var lemma = claim.Cases.Single(c => c.Name.StartsWith("product lemma", StringComparison.Ordinal));
        var optimum = claim.Cases.Single(c => c.Name.StartsWith("product optimum", StringComparison.Ordinal));
        Assert.True(lemma.Passes, lemma.Actual);
        Assert.True(optimum.Passes, optimum.Actual);
    }

    // ------------------------------------------------------------------
    // Translation invariance: exact rank over GF(p)
    // ------------------------------------------------------------------

    [Fact]
    public void TranslationInvariantRank_ReachesTheCeiling_AtThreeThree()
    {
        // 729 × 729 over GF(p); rank_p ≤ rank ≤ ceiling, so 378 is exact.
        int rank = QuditProductMirrorCap.TranslationInvariantRankModP(3, 3);
        Assert.Equal(378, rank);
        Assert.Equal(QuditProductMirrorCap.CombinatorialCeiling(3, 3), rank);
        Assert.True(rank > QuditProductMirrorCap.ProductCap(3, 3));     // 378 > 216
    }

    [Fact]
    public void TranslationInvariantRank_AtTheQubit_IsTheFullSpace()
    {
        Assert.Equal(64, QuditProductMirrorCap.TranslationInvariantRankModP(2, 3));
    }

    [Fact]
    public void ExplicitPermutationBuilder_RejectsDimensionsBeyondIntIndexing()
    {
        Assert.Throws<OverflowException>(() => QuditProductMirrorCap.BuildPiD(3, 20));
    }

    // ------------------------------------------------------------------
    // Static helper contracts
    // ------------------------------------------------------------------

    [Fact]
    public void StaticHelpers_ValidateArguments()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => QuditProductMirrorCap.ProductCap(1, 2));
        Assert.Throws<ArgumentOutOfRangeException>(() => QuditProductMirrorCap.ProductCap(3, 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => QuditProductMirrorCap.CombinatorialCeiling(1, 1));
        Assert.Throws<ArgumentOutOfRangeException>(() => QuditProductMirrorCap.BuildPiD(3, 1, chirality: 2));
        Assert.Throws<ArgumentOutOfRangeException>(() => QuditProductMirrorCap.BuildPiD(1, 1));
    }

    [Fact]
    public void NonProductPart_IsCeilingMinusCap_AndPositiveExactlyWhenNotFull()
    {
        Assert.Equal(18, QuditProductMirrorCap.NonProductPart(3, 2));   // 54 − 36
        Assert.Equal(0, QuditProductMirrorCap.NonProductPart(2, 3));    // d = 2 full
        Assert.Equal(0, QuditProductMirrorCap.NonProductPart(3, 1));    // N = 1 full
        Assert.Equal(64, QuditProductMirrorCap.NonProductPart(4, 2));   // 128 − 64
        Assert.Equal(252, QuditProductMirrorCap.NonProductPart(6, 2));  // 432 − 180, not 432 − 144
        for (int d = 2; d <= 5; d++)
            for (int N = 1; N <= 3; N++)
                Assert.Equal(d > 2 && N >= 2, QuditProductMirrorCap.NonProductPart(d, N) > 0);
    }

    [Fact]
    public void BuildPiD_AtD2N1_IsTheF118Palindromizer()
    {
        // (i, j) ↦ (j, i − 1 mod 2) on pair indices p = 2i + j:
        // (0,0)→(0,1), (0,1)→(1,1), (1,0)→(0,0), (1,1)→(1,0): one 4-cycle, ord = 4 = 2d.
        Assert.Equal(new[] { 1, 3, 0, 2 }, QuditProductMirrorCap.BuildPiD(2, 1));
    }

    [Fact]
    public void BuildPiD_IsAFullRankPermutation_NotTheRestrictedShiftRank()
    {
        var pi = QuditProductMirrorCap.BuildPiD(3, 2);

        Assert.Equal(81, pi.Count);
        Assert.Equal(81, pi.Distinct().Count()); // rank Π_d = d^(2N)
        Assert.Equal(36, QuditProductMirrorCap.ShiftAlignedRank(3, 2)); // rank Π_d P_aligned
    }

    // ------------------------------------------------------------------
    // Direct mathematical spot-check, independent of the claim's battery:
    // the (3, 2) numbers 36/54/81 and the trunk equation.
    // ------------------------------------------------------------------

    [Fact]
    public void SpotCheck_QutritN2_36_54_81_AndTrunkEquation()
    {
        // The three counts at d = 3, N = 2: product cap 36 < ceiling 54 < total 81;
        // non-product gap 18.
        Assert.Equal(36, QuditProductMirrorCap.ProductCap(3, 2));
        Assert.Equal(54, QuditProductMirrorCap.CombinatorialCeiling(3, 2));
        Assert.Equal(81, QuditPartialPalindromeCeiling.Total(3, 2));
        Assert.Equal(18, QuditProductMirrorCap.CombinatorialCeiling(3, 2) - QuditProductMirrorCap.ProductCap(3, 2));

        // Full product mirror ⟺ d² − 2d = 0 ⟺ d = 2: the trunk equation's third appearance.
        for (int d = 2; d <= 5; d++)
            for (int n = 1; n <= 3; n++)
            {
                bool full = QuditProductMirrorCap.ProductCap(d, n) == QuditPartialPalindromeCeiling.Total(d, n);
                Assert.Equal(d * d - 2 * d == 0, full);
                Assert.Equal(d == 2, full);
            }
    }
}
