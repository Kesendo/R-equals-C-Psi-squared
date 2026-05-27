using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using Xunit;

namespace RCPsiSquared.Core.Tests.Symmetry;

/// <summary>F112-Y Claim structure tests: BitB axis (shared with F112-Z), Tier1Derived
/// universal N, BitATwinStatus = BitBSpecific (matching F108 Part 3 pattern). F112-Y
/// reuses LindbladBitBPiBalance.BitBParity / IsBitBHomogeneous (no new helpers).</summary>
public class LindbladBitBPiYBalanceTests
{
    private static LindbladBitBPiYBalance Make() =>
        new LindbladBitBPiYBalance(new F108Part3Pi2YEvenAlwaysPalindromic());

    // ============================================================
    // Claim metadata
    // ============================================================

    [Fact]
    public void Z2Axis_Is_BitB() =>
        Assert.Equal(Z2Axis.BitB, Make().Z2Axis);

    [Fact]
    public void Tier_Is_Tier1Derived() =>
        Assert.Equal(Tier.Tier1Derived, Make().Tier);

    [Fact]
    public void BitATwin_Is_Null() =>
        Assert.Null(Make().BitATwin);

    [Fact]
    public void BitATwinStatus_Is_BitBSpecific() =>
        // Y-dephasing is intrinsically a bit_b-axis dephase (Π_Y² grades by bit_b per F38).
        // The BitA-axis sister F112-X (LindbladBitAPiBalance) is on a DIFFERENT axis_d,
        // not the per-letter-mirror twin captured by BitATwinClassification.Filled.
        Assert.Equal(BitATwinClassification.BitBSpecific, Make().BitATwinStatus);

    [Fact]
    public void Theorem_Mentions_BitB_Homogeneity_And_PiY_Polarity()
    {
        string theorem = Make().Theorem;
        Assert.Contains("bit_b-homogeneous", theorem);
        Assert.Contains("Π_Y", theorem);
        Assert.Contains("M_plus_half", theorem);
        Assert.Contains("M_minus_half", theorem);
    }

    [Fact]
    public void DistinctionFromF112Z_Mentions_Klein_V4_D_NonTransport()
    {
        // The distinguishing structural fact: D-conjugation intertwines Π_Z and Π_Y
        // at the operator level but does NOT transport L_Z to a Lindblad-form L_Y.
        string dist = Make().DistinctionFromF112Z;
        Assert.Contains("D-involution", dist);
        Assert.Contains("Welle 13", dist);
        Assert.Contains("Π_Y", dist);
        Assert.Contains("Π_Z", dist);
    }

    [Fact]
    public void Part3_Parent_Is_Wired()
    {
        var claim = Make();
        Assert.NotNull(claim.Part3);
        Assert.IsType<F108Part3Pi2YEvenAlwaysPalindromic>(claim.Part3);
    }
}
