using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class F63LCommutesPi2Pi2InheritanceTests
{
    private static F63LCommutesPi2Pi2Inheritance BuildClaim()
    {
        var ladder = new Pi2DyadicLadderClaim();
        var mirror = new Pi2OperatorSpaceMirrorClaim();
        var memoryLoop = new Pi2I4MemoryLoopClaim();
        var f38 = new F38Pi2InvolutionPi2Inheritance(ladder, mirror, memoryLoop);
        return new F63LCommutesPi2Pi2Inheritance(f38, ladder);
    }

    [Fact]
    public void Tier_IsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived, BuildClaim().Tier);
    }

    [Fact]
    public void IndependentZ2SymmetryCount_IsTwo()
    {
        Assert.Equal(2, BuildClaim().IndependentZ2SymmetryCount);
    }

    [Fact]
    public void BlockCount_IsExactlyFour()
    {
        Assert.Equal(4.0, BuildClaim().BlockCount, precision: 14);
    }

    [Theory]
    [InlineData(1, 1.0)]    // a_1 = 1 (trivial)
    [InlineData(2, 4.0)]    // a_{-1} = 4
    [InlineData(3, 16.0)]   // a_{-3} = 16
    [InlineData(4, 64.0)]
    [InlineData(5, 256.0)]
    [InlineData(6, 1024.0)]
    public void PerBlockDimension_MatchesLadder(int N, double expected)
    {
        Assert.Equal(expected, BuildClaim().PerBlockDimension(N), precision: 12);
    }

    [Theory]
    [InlineData(1, 1)]
    [InlineData(2, -1)]
    [InlineData(3, -3)]
    [InlineData(7, -11)]
    public void LadderIndexForPerBlock_Is3Minus2N(int N, int expected)
    {
        Assert.Equal(expected, BuildClaim().LadderIndexForPerBlock(N));
    }

    [Theory]
    [InlineData(1, 4.0)]      // 4 · 1
    [InlineData(2, 16.0)]     // 4 · 4 = 4^2
    [InlineData(3, 64.0)]     // 4 · 16 = 4^3
    [InlineData(4, 256.0)]    // 4 · 64 = 4^4
    [InlineData(6, 4096.0)]   // 4 · 1024 = 4^6
    public void FourBlockDimensionsTotal_IsFourPowerN(int N, double expected)
    {
        Assert.Equal(expected, BuildClaim().FourBlockDimensionsTotal(N), precision: 12);
    }

    [Theory]
    // ⌊N/2⌋ + 1 from F63 closed form
    [InlineData(2, 2)]
    [InlineData(3, 2)]
    [InlineData(4, 3)]
    [InlineData(5, 3)]
    [InlineData(6, 4)]
    [InlineData(7, 4)]
    public void Pi2EvenConservedCount_MatchesF63ClosedForm(int N, int expected)
    {
        Assert.Equal(expected, BuildClaim().Pi2EvenConservedCount(N));
    }

    [Theory]
    // ⌈N/2⌉ from F63 closed form
    [InlineData(2, 1)]
    [InlineData(3, 2)]
    [InlineData(4, 2)]
    [InlineData(5, 3)]
    [InlineData(6, 3)]
    [InlineData(7, 4)]
    public void Pi2OddConservedCount_MatchesF63ClosedForm(int N, int expected)
    {
        Assert.Equal(expected, BuildClaim().Pi2OddConservedCount(N));
    }

    [Theory]
    [InlineData(2, 3)]
    [InlineData(3, 4)]
    [InlineData(4, 5)]
    [InlineData(5, 6)]
    [InlineData(7, 8)]
    public void TotalConservedPerSector_IsNPlus1(int N, int expected)
    {
        Assert.Equal(expected, BuildClaim().TotalConservedPerSector(N));
    }

    [Theory]
    [InlineData(2, 8L)]    // 2^3
    [InlineData(3, 32L)]   // 2^5
    [InlineData(4, 128L)]  // 2^7
    [InlineData(5, 512L)]  // 2^9
    public void Pi2PairedBlockDimension_IsTwoPower2NMinus1(int N, long expected)
    {
        Assert.Equal(expected, BuildClaim().Pi2PairedBlockDimension(N));
    }

    [Theory]
    // F63 table: N | sector | cons (e, o) | mirror (e, o)
    // 2 |  8 | (2, 1) | (4, 6)
    // 3 | 32 | (2, 2) | (28, 28)
    // 4 | 128 | (3, 2) | (122, 124)
    // 5 | 512 | (3, 3) | (506, 506)
    [InlineData(2, 4L, 6L)]
    [InlineData(3, 28L, 28L)]
    [InlineData(4, 122L, 124L)]
    [InlineData(5, 506L, 506L)]
    public void MirrorEvenAndOddSector_MatchF63Table(int N, long expectedEven, long expectedOdd)
    {
        var f = BuildClaim();
        Assert.Equal(expectedEven, f.MirrorEvenSector(N));
        Assert.Equal(expectedOdd, f.MirrorOddSector(N));
    }

    [Theory]
    [InlineData(2)]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    public void MatchesF66EndpointMultiplicity_AcrossN(int N)
    {
        // F63 conserved-per-pole count = N + 1 = F66 endpoint multiplicity
        Assert.True(BuildClaim().MatchesF66EndpointMultiplicity(N));
    }

    [Theory]
    [InlineData(1)]
    [InlineData(2)]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    public void PerBlockDimensionAgreesWithF38_AcrossN(int N)
    {
        // F63 per-block dim = 4^(N-1); F38 EigenspaceDimension(N) = 2·4^(N-1)
        // → divide by 2 to recover per-block (Π²-even bit_a sector or Π²-odd bit_a sector)
        Assert.True(BuildClaim().PerBlockDimensionAgreesWithF38(N));
    }

    [Fact]
    public void PerBlockDimension_NLessThanOne_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().PerBlockDimension(0));
    }

    [Fact]
    public void Pi2EvenConservedCount_NLessThanOne_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().Pi2EvenConservedCount(0));
    }

    [Fact]
    public void Constructor_NullF38_Throws()
    {
        Assert.Throws<ArgumentNullException>(() =>
            new F63LCommutesPi2Pi2Inheritance(
                f38: null!,
                ladder: new Pi2DyadicLadderClaim()));
    }

    [Fact]
    public void Constructor_NullLadder_Throws()
    {
        var ladder = new Pi2DyadicLadderClaim();
        var mirror = new Pi2OperatorSpaceMirrorClaim();
        var memoryLoop = new Pi2I4MemoryLoopClaim();
        var f38 = new F38Pi2InvolutionPi2Inheritance(ladder, mirror, memoryLoop);
        Assert.Throws<ArgumentNullException>(() =>
            new F63LCommutesPi2Pi2Inheritance(f38: f38, ladder: null!));
    }
}
