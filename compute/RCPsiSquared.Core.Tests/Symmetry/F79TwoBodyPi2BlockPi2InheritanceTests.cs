using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class F79TwoBodyPi2BlockPi2InheritanceTests
{
    private static F79TwoBodyPi2BlockPi2Inheritance BuildClaim()
    {
        var ladder = new Pi2DyadicLadderClaim();
        var loop = new Pi2I4MemoryLoopClaim();
        var f1 = new F1Pi2Inheritance(new RCPsiSquared.Core.F1.F1PalindromeIdentity(), ladder, loop);
        var klein = new KleinFourCellClaim();
        var mirror = new Pi2OperatorSpaceMirrorClaim();
        return new F79TwoBodyPi2BlockPi2Inheritance(klein, mirror, f1);
    }

    [Fact]
    public void Tier_IsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived, BuildClaim().Tier);
    }

    [Theory]
    // Π²-even: bit_b sum even
    // XX: 0+0 = 0 even, even
    // YY: 1+1 = 2 even
    // ZZ: 1+1 = 2 even
    // YZ: 1+1 = 2 even
    // ZY: 1+1 = 2 even
    [InlineData('X', 'X', 0)]
    [InlineData('Y', 'Y', 0)]
    [InlineData('Z', 'Z', 0)]
    [InlineData('Y', 'Z', 0)]
    [InlineData('Z', 'Y', 0)]
    // Π²-odd: bit_b sum odd
    // XY: 0+1 = 1 odd
    // YX: 1+0 = 1 odd
    // XZ: 0+1 = 1 odd
    // ZX: 1+0 = 1 odd
    [InlineData('X', 'Y', 1)]
    [InlineData('Y', 'X', 1)]
    [InlineData('X', 'Z', 1)]
    [InlineData('Z', 'X', 1)]
    public void Pi2Parity_MatchesBitBSumModTwo(char p, char q, int expected)
    {
        Assert.Equal(expected, BuildClaim().Pi2Parity(p, q));
    }

    [Theory]
    [InlineData('X', 'X', true)]
    [InlineData('Y', 'Z', true)]
    [InlineData('X', 'Y', false)]
    public void IsPi2Even_TrueWhenParityZero(char p, char q, bool expected)
    {
        Assert.Equal(expected, BuildClaim().IsPi2Even(p, q));
    }

    [Theory]
    [InlineData('X', 'Y', true)]
    [InlineData('X', 'Z', true)]
    [InlineData('X', 'X', false)]
    [InlineData('Y', 'Z', false)]
    public void IsPi2Odd_TrueWhenParityOne(char p, char q, bool expected)
    {
        Assert.Equal(expected, BuildClaim().IsPi2Odd(p, q));
    }

    [Theory]
    [InlineData(2, 8.0)]    // 4^2 / 2 = 8
    [InlineData(3, 32.0)]   // 4^3 / 2 = 32
    [InlineData(5, 512.0)]  // 4^5 / 2 = 512
    public void Pi2BlockDimension_EqualsFourToNHalf(int N, double expected)
    {
        Assert.Equal(expected, BuildClaim().Pi2BlockDimension(N), precision: 12);
    }

    [Fact]
    public void EqualBlockDimensionsHold_HoldsAcrossN()
    {
        var f = BuildClaim();
        for (int N = 1; N <= 10; N++)
            Assert.True(f.EqualBlockDimensionsHold(N));
    }

    [Fact]
    public void Pi2OddUniversalityClustersAtN5_HasTwoEntries()
    {
        var clusters = BuildClaim().Pi2OddUniversalityClustersAtN5;
        Assert.Equal(2, clusters.Count);
        Assert.Equal(5.464, clusters[0].SvValue);
        Assert.Equal(512, clusters[0].Multiplicity);
        Assert.Equal(1.464, clusters[1].SvValue);
        Assert.Equal(512, clusters[1].Multiplicity);
    }

    [Fact]
    public void KleinSubcellsForPi2Odd_HasMpAndMm()
    {
        var sub = BuildClaim().KleinSubcellsForPi2Odd;
        Assert.Equal(2, sub.Count);
        Assert.Contains(sub, c => c.Cell.StartsWith("Mp"));
        Assert.Contains(sub, c => c.Cell.StartsWith("Mm"));
    }

    [Fact]
    public void Pi2Parity_InvalidLetter_Throws()
    {
        Assert.Throws<ArgumentException>(() => BuildClaim().Pi2Parity('a', 'X'));
        Assert.Throws<ArgumentException>(() => BuildClaim().Pi2Parity('X', '0'));
    }

    [Fact]
    public void Pi2BlockDimension_NLessThanOne_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().Pi2BlockDimension(N: 0));
    }

    [Fact]
    public void Constructor_NullKlein_Throws()
    {
        var ladder = new Pi2DyadicLadderClaim();
        var loop = new Pi2I4MemoryLoopClaim();
        var f1 = new F1Pi2Inheritance(new RCPsiSquared.Core.F1.F1PalindromeIdentity(), ladder, loop);
        var mirror = new Pi2OperatorSpaceMirrorClaim();
        Assert.Throws<ArgumentNullException>(() =>
            new F79TwoBodyPi2BlockPi2Inheritance(null!, mirror, f1));
    }

    [Fact]
    public void Constructor_NullMirror_Throws()
    {
        var ladder = new Pi2DyadicLadderClaim();
        var loop = new Pi2I4MemoryLoopClaim();
        var f1 = new F1Pi2Inheritance(new RCPsiSquared.Core.F1.F1PalindromeIdentity(), ladder, loop);
        var klein = new KleinFourCellClaim();
        Assert.Throws<ArgumentNullException>(() =>
            new F79TwoBodyPi2BlockPi2Inheritance(klein, null!, f1));
    }

    [Fact]
    public void Constructor_NullF1_Throws()
    {
        var klein = new KleinFourCellClaim();
        var mirror = new Pi2OperatorSpaceMirrorClaim();
        Assert.Throws<ArgumentNullException>(() =>
            new F79TwoBodyPi2BlockPi2Inheritance(klein, mirror, null!));
    }
}
