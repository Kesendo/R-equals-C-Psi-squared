using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class F38Pi2InvolutionPi2InheritanceTests
{
    private static F38Pi2InvolutionPi2Inheritance BuildClaim() =>
        new F38Pi2InvolutionPi2Inheritance(
            new Pi2DyadicLadderClaim(),
            new Pi2OperatorSpaceMirrorClaim(),
            new Pi2I4MemoryLoopClaim(),
            new HalfAsStructuralFixedPointClaim());

    [Fact]
    public void Tier_IsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived, BuildClaim().Tier);
    }

    [Fact]
    public void CyclicOrder_IsFour()
    {
        Assert.Equal(4, BuildClaim().CyclicOrder);
    }

    [Fact]
    public void Pi2Eigenvalues_AreCanonicalPlusOneMinusOne()
    {
        var eigs = BuildClaim().Pi2Eigenvalues;
        Assert.Equal(new[] { +1, -1 }, eigs);
    }

    [Fact]
    public void Pi2EigenvaluesFromMemoryLoop_AlternatesPlusMinus()
    {
        var eigs = BuildClaim().Pi2EigenvaluesFromMemoryLoop();
        Assert.Equal(new[] { +1, -1, +1, -1 }, eigs);
    }

    [Fact]
    public void MemoryLoopClosesAtFour_HoldsExactly()
    {
        Assert.True(BuildClaim().MemoryLoopClosesAtFour());
    }

    [Theory]
    [InlineData(1, 4L)]
    [InlineData(2, 16L)]
    [InlineData(3, 64L)]
    [InlineData(4, 256L)]
    [InlineData(5, 1024L)]
    [InlineData(6, 4096L)]
    public void FullOperatorSpaceDimension_MatchesFourPowerN(int N, long expected)
    {
        Assert.Equal(expected, BuildClaim().FullOperatorSpaceDimension(N));
    }

    [Theory]
    [InlineData(1, 2L)]
    [InlineData(2, 8L)]
    [InlineData(3, 32L)]
    [InlineData(4, 128L)]
    [InlineData(5, 512L)]
    [InlineData(6, 2048L)]
    public void EigenspaceDimension_IsHalfOfFour_PowerN(int N, long expected)
    {
        Assert.Equal(expected, BuildClaim().EigenspaceDimension(N));
    }

    [Theory]
    [InlineData(1, 2.0)]    // a_0 · a_1 = 2 · 1
    [InlineData(2, 8.0)]    // a_0 · a_{-1} = 2 · 4
    [InlineData(3, 32.0)]   // a_0 · a_{-3} = 2 · 16
    [InlineData(4, 128.0)]  // a_0 · a_{-5} = 2 · 64
    [InlineData(5, 512.0)]  // a_0 · a_{-7} = 2 · 256
    public void EigenspaceDimensionViaLadder_FactorisesAsTwoTimesFourPowerNMinus1(int N, double expected)
    {
        Assert.Equal(expected, BuildClaim().EigenspaceDimensionViaLadder(N), precision: 10);
    }

    [Theory]
    [InlineData(1, 1)]      // 3-2·1 = 1 (a_1 = 1, trivial)
    [InlineData(2, -1)]     // 3-2·2 = -1 (a_{-1} = 4)
    [InlineData(3, -3)]     // 3-2·3 = -3 (a_{-3} = 16)
    [InlineData(4, -5)]
    [InlineData(5, -7)]
    public void LadderIndexForFourPowerNMinus1_Is3Minus2N(int N, int expectedIndex)
    {
        Assert.Equal(expectedIndex, BuildClaim().LadderIndexForFourPowerNMinus1(N));
    }

    [Theory]
    [InlineData(1, 0)]
    [InlineData(2, 1)]
    [InlineData(3, 2)]
    [InlineData(6, 5)]
    public void OperatorSpaceQubitCountFor_IsNMinus1(int N, int expected)
    {
        Assert.Equal(expected, BuildClaim().OperatorSpaceQubitCountFor(N));
    }

    [Theory]
    [InlineData(1)]
    [InlineData(2)]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    public void HalfHalfBalance_IsExactlyOneHalf(int N)
    {
        Assert.Equal(0.5, BuildClaim().HalfHalfBalance(N), precision: 14);
    }

    [Theory]
    [InlineData(1, 4.0)]
    [InlineData(2, 16.0)]
    [InlineData(3, 64.0)]
    [InlineData(4, 256.0)]
    [InlineData(5, 1024.0)]
    [InlineData(6, 4096.0)]
    public void MirrorPinnedFullDimension_AgreesWithFullOperatorSpaceDimension(int N, double expected)
    {
        var f = BuildClaim();
        Assert.Equal(expected, f.MirrorPinnedFullDimension(N), precision: 12);
        Assert.Equal((double)f.FullOperatorSpaceDimension(N), f.MirrorPinnedFullDimension(N), precision: 12);
    }

    [Theory]
    [InlineData(1)]
    [InlineData(2)]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    public void LadderFactorisationHolds_AcrossN(int N)
    {
        Assert.True(BuildClaim().LadderFactorisationHolds(N));
    }

    [Fact]
    public void FullOperatorSpaceDimension_ThrowsForNLessThanOne()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().FullOperatorSpaceDimension(0));
    }

    [Fact]
    public void EigenspaceDimensionViaLadder_ThrowsForNLessThanOne()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().EigenspaceDimensionViaLadder(0));
    }

    [Fact]
    public void Constructor_RejectsNullParents()
    {
        var ladder = new Pi2DyadicLadderClaim();
        var mirror = new Pi2OperatorSpaceMirrorClaim();
        var memoryLoop = new Pi2I4MemoryLoopClaim();
        var half = new HalfAsStructuralFixedPointClaim();
        Assert.Throws<ArgumentNullException>(() =>
            new F38Pi2InvolutionPi2Inheritance(null!, mirror, memoryLoop, half));
        Assert.Throws<ArgumentNullException>(() =>
            new F38Pi2InvolutionPi2Inheritance(ladder, null!, memoryLoop, half));
        Assert.Throws<ArgumentNullException>(() =>
            new F38Pi2InvolutionPi2Inheritance(ladder, mirror, null!, half));
        Assert.Throws<ArgumentNullException>(() =>
            new F38Pi2InvolutionPi2Inheritance(ladder, mirror, memoryLoop, null!));
    }

    [Fact]
    public void TypedParents_HalfIsExposed()
    {
        Assert.NotNull(BuildClaim().Half);
    }
}
