using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class F89PathKHbMixedDegreesClaimTests
{
    private static F89PathKHbMixedDegreesClaim BuildClaim() =>
        new F89PathKHbMixedDegreesClaim(new F89TopologyOrbitClosure(new Pi2DyadicLadderClaim()));

    [Fact]
    public void Tier_IsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived, BuildClaim().Tier);
    }

    [Theory]
    [InlineData(1, 2)]      // path-1: N_block=2, dim = 2 · C(2,2) = 2 · 1 = 2
    [InlineData(2, 9)]      // path-2: N_block=3, dim = 3 · C(3,2) = 3 · 3 = 9
    [InlineData(3, 24)]     // path-3: N_block=4, dim = 4 · C(4,2) = 4 · 6 = 24
    [InlineData(4, 50)]
    [InlineData(5, 90)]
    [InlineData(6, 147)]
    public void SeDeFullDimension_MatchesNBlockTimesBinomialNBlock2(int k, int expected)
    {
        Assert.Equal(expected, F89PathKHbMixedDegreesClaim.SeDeFullDimension(k));
    }

    [Theory]
    [InlineData(3, 12)]
    [InlineData(4, 26)]
    [InlineData(5, 45)]
    [InlineData(6, 75)]
    public void S2SymSubBlockDimension_MatchesEmpirical(int k, int expected)
    {
        Assert.Equal(expected, F89PathKHbMixedDegreesClaim.S2SymSubBlockDimension(k));
    }

    [Theory]
    [InlineData(3, 4)]
    [InlineData(4, 8)]
    [InlineData(5, 13)]
    [InlineData(6, 22)]
    public void AtLockedCountInS2Sym_MatchesEmpirical(int k, int expected)
    {
        Assert.Equal(expected, F89PathKHbMixedDegreesClaim.AtLockedCountInS2Sym(k));
    }

    [Theory]
    [InlineData(3, 8)]
    [InlineData(4, 18)]
    [InlineData(5, 32)]
    [InlineData(6, 53)]
    public void HbMixedDegree_EqualsS2SymMinusAtLocked(int k, int expected)
    {
        Assert.Equal(expected, F89PathKHbMixedDegreesClaim.HbMixedSubFactorDegree(k));
    }

    [Fact]
    public void HbMixedDegree_UnsupportedPath_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(
            () => F89PathKHbMixedDegreesClaim.HbMixedSubFactorDegree(2));
    }

    [Fact]
    public void SeDeFullDimension_NegativeKThrows()
    {
        Assert.Throws<ArgumentOutOfRangeException>(
            () => F89PathKHbMixedDegreesClaim.SeDeFullDimension(0));
    }

    [Fact]
    public void S2SymSubBlockDimension_UnsupportedPath_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(
            () => F89PathKHbMixedDegreesClaim.S2SymSubBlockDimension(7));
    }

    [Fact]
    public void AtLockedCountInS2Sym_UnsupportedPath_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(
            () => F89PathKHbMixedDegreesClaim.AtLockedCountInS2Sym(7));
    }
}
