using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class F62WStateBornBelowFoldPi2InheritanceTests
{
    internal static F61BitAParityPi2Inheritance BuildF61(Pi2DyadicLadderClaim ladder)
    {
        var f38 = new F38Pi2InvolutionPi2Inheritance(ladder, new Pi2OperatorSpaceMirrorClaim(), new Pi2I4MemoryLoopClaim(), new HalfAsStructuralFixedPointClaim());
        var f63 = new F63LCommutesPi2Pi2Inheritance(f38, ladder);
        return new F61BitAParityPi2Inheritance(f63, ladder);
    }

    private static F62WStateBornBelowFoldPi2Inheritance BuildClaim()
    {
        var ladder = new Pi2DyadicLadderClaim();
        return new F62WStateBornBelowFoldPi2Inheritance(ladder, new QuarterAsBilinearMaxvalClaim(), BuildF61(ladder));
    }

    [Fact]
    public void Tier_IsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived, BuildClaim().Tier);
    }

    [Fact]
    public void FoldPosition_IsExactlyOneQuarter()
    {
        Assert.Equal(0.25, BuildClaim().FoldPosition, precision: 14);
    }

    [Fact]
    public void NumeratorTwoCoefficient_IsExactlyTwo()
    {
        Assert.Equal(2.0, BuildClaim().NumeratorTwoCoefficient, precision: 14);
    }

    [Theory]
    [InlineData(2, 0.5)]
    [InlineData(3, 1.0 / 3.0)]
    [InlineData(4, 0.25)]
    [InlineData(5, 0.2)]
    [InlineData(10, 0.1)]
    public void PairOffDiagonalElement_Is1OverN(int N, double expected)
    {
        Assert.Equal(expected, BuildClaim().PairOffDiagonalElement(N), precision: 12);
    }

    [Theory]
    // F62 verified table from ANALYTICAL_FORMULAS:
    // N=2: 1/3 = 0.3333 (W_2 = Bell+, above fold)
    // N=3: 10/81 = 0.1235 (below fold)
    // N=5: 26/375 = 0.0693
    // N=10: 68/1500 = 0.04533
    [InlineData(2, 1.0 / 3.0)]
    [InlineData(3, 10.0 / 81.0)]
    [InlineData(5, 26.0 / 375.0)]
    [InlineData(10, 68.0 / 1500.0)]
    public void CPsiAtZeroForWState_MatchesClosedForm(int N, double expected)
    {
        Assert.Equal(expected, BuildClaim().CPsiAtZeroForWState(N), precision: 12);
    }

    [Theory]
    [InlineData(2, false)]   // Bell+ regime, above fold
    [InlineData(3, true)]    // first N below fold
    [InlineData(4, true)]
    [InlineData(5, true)]
    [InlineData(10, true)]
    public void IsBornBelowFold_FollowsClosedForm(int N, bool expected)
    {
        Assert.Equal(expected, BuildClaim().IsBornBelowFold(N));
    }

    [Fact]
    public void SmallestNBelowFold_IsThree()
    {
        Assert.Equal(3, BuildClaim().SmallestNBelowFold);
    }

    [Fact]
    public void BellPlusAboveFold_HoldsAtN2()
    {
        // W_2 = Bell+; CΨ(0) = 1/3 > 1/4
        Assert.True(BuildClaim().BellPlusAboveFold());
    }

    [Fact]
    public void F60_F62_SiblingPair_BothBelowFoldForNAtLeast3()
    {
        // Cross-anchor verification: GHZ_N (F60) and W_N (F62) agree on smallestN = 3.
        var ladder = new Pi2DyadicLadderClaim();
        var f60 = new F60GhzBornBelowFoldPi2Inheritance(
            ladder,
            new PolarityLayerOriginClaim(),
            new QuarterAsBilinearMaxvalClaim(),
            new ArgmaxMaxvalPairClaim());
        var f62 = new F62WStateBornBelowFoldPi2Inheritance(ladder, new QuarterAsBilinearMaxvalClaim(), BuildF61(ladder));

        Assert.Equal(f60.SmallestNBelowFold, f62.SmallestNBelowFold);
        Assert.True(f60.BellPlusAboveFold());
        Assert.True(f62.BellPlusAboveFold());
    }

    [Fact]
    public void CPsiAtZeroForWState_NLessThanTwo_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().CPsiAtZeroForWState(1));
    }

    [Fact]
    public void Constructor_RejectsNullParents()
    {
        var ladder = new Pi2DyadicLadderClaim();
        var quarter = new QuarterAsBilinearMaxvalClaim();
        var f61 = BuildF61(ladder);
        Assert.Throws<ArgumentNullException>(() =>
            new F62WStateBornBelowFoldPi2Inheritance(null!, quarter, f61));
        Assert.Throws<ArgumentNullException>(() =>
            new F62WStateBornBelowFoldPi2Inheritance(ladder, null!, f61));
        Assert.Throws<ArgumentNullException>(() =>
            new F62WStateBornBelowFoldPi2Inheritance(ladder, quarter, null!));
    }

    [Fact]
    public void TypedParents_AreExposed()
    {
        var f = BuildClaim();
        Assert.NotNull(f.Quarter);
        Assert.NotNull(f.F61);
    }
}
