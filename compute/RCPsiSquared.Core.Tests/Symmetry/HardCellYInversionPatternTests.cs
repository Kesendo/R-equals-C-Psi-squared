using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.Symmetry;
using Xunit;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class HardCellYInversionPatternTests
{
    [Fact]
    public void Z2Axis_IsYParity() =>
        Assert.Equal(Z2Axis.YParity, new HardCellYInversionPattern().Z2Axis);

    [Fact]
    public void BitATwin_IsNull() =>
        Assert.Null(new HardCellYInversionPattern().BitATwin);

    [Fact]
    public void BitATwinStatus_IsNotApplicableForThisAxis() =>
        Assert.Equal(BitATwinClassification.NotApplicableForThisAxis,
            ((IZ2AxisClaim)new HardCellYInversionPattern()).BitATwinStatus);

    [Fact]
    public void Tier_IsTier1Candidate() =>
        Assert.Equal(Tier.Tier1Candidate, new HardCellYInversionPattern().Tier);

    // ============================================================
    // Aspect A: diagonal Klein cell membership (closed-form)
    // ============================================================

    [Theory]
    [InlineData(PauliLetter.Z, 0, 1)]
    [InlineData(PauliLetter.X, 1, 0)]
    [InlineData(PauliLetter.Y, 1, 1)]
    public void DiagonalKleinCellForDephase_MatchesF87DissipatorResonance(
        PauliLetter dephase, int expectedBitA, int expectedBitB)
    {
        // F87 dissipator-resonance law: hard appears only when the pair's Klein
        // index matches the dephase letter's own Klein index. Z=(0,1), X=(1,0),
        // Y=(1,1) per the bit_a/bit_b convention of PauliLetter.
        var (bitA, bitB) = HardCellYInversionPattern.DiagonalKleinCellForDephase(dephase);
        Assert.Equal(expectedBitA, bitA);
        Assert.Equal(expectedBitB, bitB);
    }

    [Fact]
    public void DiagonalKleinCellForDephase_RejectsIdentityLetter()
    {
        Assert.Throws<ArgumentException>(() =>
            HardCellYInversionPattern.DiagonalKleinCellForDephase(PauliLetter.I));
    }

    [Fact]
    public void IsDiagonalCell_TrueForMatchingPair()
    {
        Assert.True(HardCellYInversionPattern.IsDiagonalCell((0, 1), PauliLetter.Z));
        Assert.True(HardCellYInversionPattern.IsDiagonalCell((1, 0), PauliLetter.X));
        Assert.True(HardCellYInversionPattern.IsDiagonalCell((1, 1), PauliLetter.Y));
    }

    [Fact]
    public void IsDiagonalCell_FalseForNonMatchingPair()
    {
        // Klein (0,0) is the Mother sector: never diagonal under any dephase.
        Assert.False(HardCellYInversionPattern.IsDiagonalCell((0, 0), PauliLetter.Z));
        Assert.False(HardCellYInversionPattern.IsDiagonalCell((0, 0), PauliLetter.X));
        Assert.False(HardCellYInversionPattern.IsDiagonalCell((0, 0), PauliLetter.Y));
        // Cross-dephase mismatches.
        Assert.False(HardCellYInversionPattern.IsDiagonalCell((0, 1), PauliLetter.X));
        Assert.False(HardCellYInversionPattern.IsDiagonalCell((1, 0), PauliLetter.Y));
        Assert.False(HardCellYInversionPattern.IsDiagonalCell((1, 1), PauliLetter.Z));
    }

    // ============================================================
    // Aspect B: Y-inversion structural reading
    // ============================================================

    [Theory]
    [InlineData(PauliLetter.Z, 0)]
    [InlineData(PauliLetter.X, 0)]
    [InlineData(PauliLetter.Y, 1)]
    public void DominantYParityForDephase_EqualsLetterYParity(
        PauliLetter dephase, int expected)
    {
        // Y-letter has #Y=1 (y_par=1); Z and X have #Y=0 (y_par=0). The "dominant
        // y_par" in the hard cell equals the dephase letter's own y_par; this is
        // the Y-inversion structural reading.
        Assert.Equal(expected, HardCellYInversionPattern.DominantYParityForDephase(dephase));
    }

    [Fact]
    public void DominantYParityForDephase_RejectsIdentityLetter()
    {
        Assert.Throws<ArgumentException>(() =>
            HardCellYInversionPattern.DominantYParityForDephase(PauliLetter.I));
    }

    [Fact]
    public void Theorem_MentionsBothAspectAAndAspectB()
    {
        var claim = new HardCellYInversionPattern();
        Assert.Contains("diagonal Klein", claim.Theorem);
        Assert.Contains("Y-inversion", claim.Theorem);
    }

    [Fact]
    public void F87Corollary_ScopesToDephaseLetterYParity()
    {
        var claim = new HardCellYInversionPattern();
        Assert.Contains("y_par(dephase", claim.F87Corollary);
    }
}
