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
        Assert.Equal(8, claim.Cases.Count);
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
        Assert.Contains("(2d)^N", claim.Summary);
        Assert.Contains("Z_d ≀ Z₂", claim.Summary);
        Assert.Contains($"{claim.Cases.Count}/{claim.Cases.Count} battery PASS", claim.Summary);
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

    // ------------------------------------------------------------------
    // Direct mathematical spot-check, independent of the claim's battery:
    // the (3, 2) numbers 36/54/81 and the trunk equation.
    // ------------------------------------------------------------------

    [Fact]
    public void SpotCheck_QutritN2_36_54_81_AndTrunkEquation()
    {
        // The three counts at d = 3, N = 2: product cap 36 < ceiling 54 < total 81,
        // non-product gap 18.
        Assert.Equal(36, QuditProductMirrorCap.ProductCap(3, 2));
        Assert.Equal(54, QuditProductMirrorCap.CombinatorialCeiling(3, 2));
        Assert.Equal(81, QuditPartialPalindromeCeiling.Total(3, 2));
        Assert.Equal(18, QuditProductMirrorCap.CombinatorialCeiling(3, 2) - QuditProductMirrorCap.ProductCap(3, 2));

        // The trunk equation: the cap is full ⟺ (2d)^N = d^{2N} ⟺ d² − 2d = 0 ⟺ d = 2.
        for (int d = 2; d <= 5; d++)
            for (int n = 1; n <= 3; n++)
            {
                bool full = QuditProductMirrorCap.ProductCap(d, n) == QuditPartialPalindromeCeiling.Total(d, n);
                Assert.Equal(d * d - 2 * d == 0, full);
                Assert.Equal(d == 2, full);
            }
    }
}
