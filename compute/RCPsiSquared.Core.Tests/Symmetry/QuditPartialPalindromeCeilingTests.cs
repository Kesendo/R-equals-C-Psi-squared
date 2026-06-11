using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using Xunit;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class QuditPartialPalindromeCeilingTests
{
    private static QuditPartialPalindromeCeiling MakeClaim() =>
        new QuditPartialPalindromeCeiling(
            new QubitNecessityPi2Inheritance(new Pi2DyadicLadderClaim(), new Pi2OperatorSpaceMirrorClaim()));

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
        Assert.Contains("qutrit_partial_palindrome.py", claim.Anchor);
    }

    [Fact]
    public void TypedParent_IsWired()
    {
        Assert.NotNull(MakeClaim().QubitNecessity);
    }

    [Fact]
    public void Constructor_RejectsNullParent()
    {
        Assert.Throws<ArgumentNullException>(() => new QuditPartialPalindromeCeiling(null!));
    }

    [Fact]
    public void Multiplicity_QutritN2_Is_9_36_36()
    {
        Assert.Equal(9, QuditPartialPalindromeCeiling.Multiplicity(3, 2, 0));
        Assert.Equal(36, QuditPartialPalindromeCeiling.Multiplicity(3, 2, 1));
        Assert.Equal(36, QuditPartialPalindromeCeiling.Multiplicity(3, 2, 2));
    }

    [Fact]
    public void Ceiling_QutritN2_Is_54()
    {
        Assert.Equal(54, QuditPartialPalindromeCeiling.Ceiling(3, 2));
        Assert.Equal(81, QuditPartialPalindromeCeiling.Total(3, 2));
    }

    [Theory]
    [InlineData(2, 1)]
    [InlineData(2, 2)]
    [InlineData(2, 3)]
    public void Ceiling_Qubit_IsFull(int d, int N)
    {
        Assert.Equal(QuditPartialPalindromeCeiling.Total(d, N), QuditPartialPalindromeCeiling.Ceiling(d, N));
    }

    [Theory]
    [InlineData(3, 2)]
    [InlineData(4, 2)]
    [InlineData(3, 3)]
    public void Ceiling_Qudit_IsStrictlyPartial(int d, int N)
    {
        Assert.True(QuditPartialPalindromeCeiling.Ceiling(d, N) < QuditPartialPalindromeCeiling.Total(d, N));
    }

    [Fact]
    public void Battery_AllCasesPass()
    {
        var claim = MakeClaim();
        Assert.Equal(6, claim.Cases.Count);
        foreach (var c in claim.Cases)
            Assert.True(c.Passes, $"battery case '{c.Name}' failed: expected {c.Expected}, got {c.Actual}");
        Assert.Equal(claim.Cases.Count, claim.PassCount);
    }

    [Fact]
    public void Summary_CarriesTheCeilingAndTheQutritNumber()
    {
        var claim = MakeClaim();
        Assert.Contains("(d−1)^k", claim.Summary);
        Assert.Contains("54/81", claim.Summary);
    }
}
