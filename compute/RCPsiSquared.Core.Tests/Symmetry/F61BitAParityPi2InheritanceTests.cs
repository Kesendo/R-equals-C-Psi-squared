using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class F61BitAParityPi2InheritanceTests
{
    private static F61BitAParityPi2Inheritance BuildClaim()
    {
        var ladder = new Pi2DyadicLadderClaim();
        var mirror = new Pi2OperatorSpaceMirrorClaim();
        var memoryLoop = new Pi2I4MemoryLoopClaim();
        var f38 = new F38Pi2InvolutionPi2Inheritance(ladder, mirror, memoryLoop, new HalfAsStructuralFixedPointClaim());
        var f63 = new F63LCommutesPi2Pi2Inheritance(f38, ladder);
        return new F61BitAParityPi2Inheritance(f63, ladder);
    }

    [Fact]
    public void Tier_IsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived, BuildClaim().Tier);
    }

    [Fact]
    public void IndependentZ2SymmetryCount_IsTwoFromF63()
    {
        Assert.Equal(2, BuildClaim().IndependentZ2SymmetryCount);
    }

    [Fact]
    public void BlockCount_IsExactlyFour()
    {
        Assert.Equal(4.0, BuildClaim().BlockCount, precision: 14);
    }

    [Fact]
    public void Z2Axis_IsBitANXY()
    {
        Assert.Equal("bit_a (n_XY)", BuildClaim().Z2Axis);
    }

    [Theory]
    [InlineData(1, 1.0)]
    [InlineData(2, 4.0)]
    [InlineData(3, 16.0)]
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
    [InlineData(1)]
    [InlineData(2)]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    public void PerBlockDimensionAgreesWithF63_AcrossN(int N)
    {
        // F61 and F63 both anchor to a_{3-2N} for per-block dim
        Assert.True(BuildClaim().PerBlockDimensionAgreesWithF63(N));
    }

    [Fact]
    public void BlockCountAgreesWithF63_HoldsExactly()
    {
        Assert.True(BuildClaim().BlockCountAgreesWithF63());
    }

    [Theory]
    [InlineData(0, true)]   // even n_XY → SE-accessible
    [InlineData(1, false)]  // odd n_XY → not SE-accessible
    public void IsSeAccessible_ReturnsTrueForEvenParity(int nXyParity, bool expected)
    {
        Assert.Equal(expected, BuildClaim().IsSeAccessible(nXyParity));
    }

    [Fact]
    public void BreakConditions_IncludesT1AndTransverseFields()
    {
        var conditions = BuildClaim().BreakConditions;
        Assert.Contains(conditions, c => c.Contains("amplitude damping", StringComparison.OrdinalIgnoreCase) || c.Contains("T1", StringComparison.OrdinalIgnoreCase));
        Assert.Contains(conditions, c => c.Contains("h_x"));
        Assert.Contains(conditions, c => c.Contains("h_y"));
    }

    [Fact]
    public void PerBlockDimension_NLessThanOne_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().PerBlockDimension(0));
    }

    [Fact]
    public void Constructor_NullF63_Throws()
    {
        Assert.Throws<ArgumentNullException>(() =>
            new F61BitAParityPi2Inheritance(
                f63: null!,
                ladder: new Pi2DyadicLadderClaim()));
    }

    [Fact]
    public void Constructor_NullLadder_Throws()
    {
        var ladder = new Pi2DyadicLadderClaim();
        var mirror = new Pi2OperatorSpaceMirrorClaim();
        var memoryLoop = new Pi2I4MemoryLoopClaim();
        var f38 = new F38Pi2InvolutionPi2Inheritance(ladder, mirror, memoryLoop, new HalfAsStructuralFixedPointClaim());
        var f63 = new F63LCommutesPi2Pi2Inheritance(f38, ladder);
        Assert.Throws<ArgumentNullException>(() =>
            new F61BitAParityPi2Inheritance(f63: f63, ladder: null!));
    }
}
