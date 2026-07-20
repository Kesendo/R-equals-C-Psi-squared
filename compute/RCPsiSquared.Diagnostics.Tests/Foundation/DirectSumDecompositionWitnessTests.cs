using RCPsiSquared.Diagnostics.Foundation;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

public class DirectSumDecompositionWitnessTests
{
    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    public void TheWall_OffParityBlockIsZero_WhileLIsNot(int n)
    {
        var w = new DirectSumDecompositionWitness(n);
        Assert.True(w.OffParityNorm < 1e-10, $"off-parity {w.OffParityNorm:E2}");
        Assert.True(w.FullNorm > 1.0, $"‖L‖ {w.FullNorm:E2} — non-triviality gate");
    }

    [Theory]
    [InlineData(3, 32L)]
    [InlineData(4, 128L)]
    public void SectorDimensions_AreEqualHalves(int n, long expected)
    {
        var w = new DirectSumDecompositionWitness(n);
        Assert.Equal(expected, w.EvenSectorDim);
        Assert.Equal(expected, w.OddSectorDim);
    }

    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    public void PiSectorMap_IsColumnCompleteConsistent(int n)
    {
        Assert.True(new DirectSumDecompositionWitness(n).PiSectorMapConsistent);
    }

    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    public void SectorPalindrome_ResidualMachineZero(int n)
    {
        var w = new DirectSumDecompositionWitness(n);
        Assert.True(w.PalindromeResidualNorm < 1e-9, $"global M {w.PalindromeResidualNorm:E2}");
        Assert.True(w.EvenBlockResidualNorm < 1e-9, $"(even,even) {w.EvenBlockResidualNorm:E2}");
        Assert.True(w.OddBlockResidualNorm < 1e-9, $"(odd,odd) {w.OddBlockResidualNorm:E2}");
    }

    [Fact]
    public void Control_T1BreaksTheMirrorNotTheWall()
    {
        var w = new DirectSumDecompositionWitness(3);
        Assert.True(w.AdOffParityNorm < 1e-10, $"T1 wall {w.AdOffParityNorm:E2} — must survive exactly");
        Assert.True(w.AdPalindromeNorm > 1e-2, $"T1 mirror {w.AdPalindromeNorm:E2} — must break");
    }

    [Fact]
    public void Control_FieldBreaksTheWallNotTheMirror()
    {
        var w = new DirectSumDecompositionWitness(3);
        Assert.True(w.FieldOffParityNorm > 1e-2, $"field wall {w.FieldOffParityNorm:E2} — must break");
        Assert.True(w.FieldPalindromeNorm < 1e-9, $"field mirror {w.FieldPalindromeNorm:E2} — must survive");
    }

    [Fact]
    public void Guard_TooLargeOrTooSmallN_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => new DirectSumDecompositionWitness(1));
        Assert.Throws<ArgumentOutOfRangeException>(() => new DirectSumDecompositionWitness(6));
        Assert.Throws<ArgumentOutOfRangeException>(() => new DirectSumDecompositionWitness(3, gamma: 0.0));
    }
}
