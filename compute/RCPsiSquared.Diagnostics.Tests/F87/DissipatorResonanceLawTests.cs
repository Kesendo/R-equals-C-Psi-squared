using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Diagnostics.F87;

namespace RCPsiSquared.Diagnostics.Tests.F87;

public class DissipatorResonanceLawTests
{
    [Fact]
    public void Claim_IsTier1Derived()
    {
        var c = new DissipatorResonanceLaw();
        Assert.Equal(Tier.Tier1Derived, c.Tier);
    }

    [Fact]
    public void WitnessTable_HasTwelveEntries_FourCellsTimesThreeLetters()
    {
        var c = new DissipatorResonanceLaw();
        Assert.Equal(12, c.Witnesses.Count);
    }

    [Fact]
    public void MotherCell_IsUniversallyHardFree_AcrossAllThreeLetters()
    {
        var c = new DissipatorResonanceLaw();
        var motherWitnesses = c.Witnesses.Where(w => w.KleinIndex == (0, 0)).ToList();
        Assert.Equal(3, motherWitnesses.Count);
        Assert.All(motherWitnesses, w => Assert.Equal(0, w.HardCount));
        Assert.All(motherWitnesses, w => Assert.Equal(66, w.TotalCount));
    }

    [Theory]
    [InlineData(PauliLetter.Z, 0, 1)]
    [InlineData(PauliLetter.X, 1, 0)]
    [InlineData(PauliLetter.Y, 1, 1)]
    public void DiagonalCells_HaveFiftyHardOfSeventySix_ForEachLetter(
        PauliLetter dephase, int expectedBitA, int expectedBitB)
    {
        var c = new DissipatorResonanceLaw();
        var w = c.Witnesses.Single(x =>
            x.DephaseLetter == dephase && x.KleinIndex == (expectedBitA, expectedBitB));
        Assert.Equal(50, w.HardCount);
        Assert.Equal(76, w.TotalCount);
    }

    [Theory]
    [InlineData(PauliLetter.Z)]
    [InlineData(PauliLetter.X)]
    [InlineData(PauliLetter.Y)]
    public void OffDiagonalCells_AreZeroHard(PauliLetter dephase)
    {
        var c = new DissipatorResonanceLaw();
        var matched = dephase switch
        {
            PauliLetter.Z => (0, 1),
            PauliLetter.X => (1, 0),
            PauliLetter.Y => (1, 1),
            _ => throw new ArgumentOutOfRangeException(nameof(dephase)),
        };
        var offDiagonal = c.Witnesses
            .Where(w => w.DephaseLetter == dephase
                        && w.KleinIndex != (0, 0)
                        && w.KleinIndex != matched)
            .ToList();
        Assert.Equal(2, offDiagonal.Count);
        Assert.All(offDiagonal, w => Assert.Equal(0, w.HardCount));
        Assert.All(offDiagonal, w => Assert.Equal(76, w.TotalCount));
    }

    [Fact]
    public void Anchor_References_TheClosingPythonScript()
    {
        var c = new DissipatorResonanceLaw();
        Assert.Contains("klein_dissipator_resonance.py", c.Anchor);
    }
}
