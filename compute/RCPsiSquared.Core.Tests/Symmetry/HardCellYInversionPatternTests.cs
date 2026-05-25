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
        // Klein (0,0) is the Mother sector — never diagonal under any dephase.
        Assert.False(HardCellYInversionPattern.IsDiagonalCell((0, 0), PauliLetter.Z));
        Assert.False(HardCellYInversionPattern.IsDiagonalCell((0, 0), PauliLetter.X));
        Assert.False(HardCellYInversionPattern.IsDiagonalCell((0, 0), PauliLetter.Y));
        // Cross-dephase mismatches.
        Assert.False(HardCellYInversionPattern.IsDiagonalCell((0, 1), PauliLetter.X));
        Assert.False(HardCellYInversionPattern.IsDiagonalCell((1, 0), PauliLetter.Y));
        Assert.False(HardCellYInversionPattern.IsDiagonalCell((1, 1), PauliLetter.Z));
    }
}
