using RCPsiSquared.Core.F1;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class DirectSumDecompositionClaimTests
{
    private static DirectSumDecompositionClaim BuildClaim()
    {
        var ladder = new Pi2DyadicLadderClaim();
        var mirror = new Pi2OperatorSpaceMirrorClaim();
        var memoryLoop = new Pi2I4MemoryLoopClaim();
        var f38 = new F38Pi2InvolutionPi2Inheritance(ladder, mirror, memoryLoop, new HalfAsStructuralFixedPointClaim());
        var f63 = new F63LCommutesPi2Pi2Inheritance(f38, ladder);
        var f61 = new F61BitAParityPi2Inheritance(f63, ladder);
        return new DirectSumDecompositionClaim(new F1PalindromeIdentity(), f61);
    }

    [Fact]
    public void Tier_IsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived, BuildClaim().Tier);
    }

    [Theory]
    [InlineData(2, 8L)]
    [InlineData(3, 32L)]
    [InlineData(4, 128L)]
    [InlineData(5, 512L)]
    [InlineData(6, 2048L)]
    [InlineData(7, 8192L)]
    [InlineData(8, 32768L)]
    public void SectorDimension_IsHalfOf4PowN(int n, long expected)
    {
        var claim = BuildClaim();
        Assert.Equal(expected, claim.SectorDimension(n));
        Assert.Equal(1L << (2 * n), 2 * claim.SectorDimension(n));
    }

    [Theory]
    [InlineData(2, false)]
    [InlineData(3, true)]
    [InlineData(4, false)]
    [InlineData(5, true)]
    public void PiExchangesSectors_IffNOdd(int n, bool expected)
    {
        Assert.Equal(expected, BuildClaim().PiExchangesSectors(n));
    }

    [Fact]
    public void SectorDimension_NLessThanOne_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().SectorDimension(0));
    }

    [Fact]
    public void Anchor_CitesTheProofAndTheProbe()
    {
        var anchor = BuildClaim().Anchor;
        Assert.Contains("DIRECT_SUM_DECOMPOSITION.md", anchor);
        Assert.Contains("direct_sum_scope_probe.py", anchor);
    }

    [Fact]
    public void Parents_AreF1AndF61()
    {
        var claim = BuildClaim();
        Assert.NotNull(claim.F1);
        Assert.NotNull(claim.F61);
    }
}
