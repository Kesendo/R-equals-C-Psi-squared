using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class F71MirrorSymmetryPi2InheritanceTests
{
    private static F71MirrorSymmetryPi2Inheritance BuildClaim() =>
        new F71MirrorSymmetryPi2Inheritance(new Pi2DyadicLadderClaim());

    [Fact]
    public void Tier_IsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived, BuildClaim().Tier);
    }

    [Theory]
    // Bond pairing b ↔ N−2−b
    // N=3: bonds 0,1 → 1↔0 (pair); 1↔0 (same pair)
    // N=4: bonds 0,1,2 → 0↔2, 1↔1 (self-paired center)
    // N=5: bonds 0..3 → 0↔3, 1↔2 (no center)
    // N=6: bonds 0..4 → 0↔4, 1↔3, 2↔2 (self-paired center)
    [InlineData(3, 0, 1)]
    [InlineData(3, 1, 0)]
    [InlineData(4, 0, 2)]
    [InlineData(4, 1, 1)]
    [InlineData(4, 2, 0)]
    [InlineData(5, 0, 3)]
    [InlineData(5, 1, 2)]
    [InlineData(6, 2, 2)]
    public void MirrorPair_MapsBondToReflection(int N, int b, int expectedPair)
    {
        Assert.Equal(expectedPair, BuildClaim().MirrorPair(N, b));
    }

    [Theory]
    [InlineData(4, 1, true)]    // even N=4, center bond at (4-2)/2 = 1
    [InlineData(6, 2, true)]    // even N=6, center bond at (6-2)/2 = 2
    [InlineData(3, 0, false)]   // odd N=3, no self-paired bond
    [InlineData(5, 1, false)]
    [InlineData(5, 2, false)]
    [InlineData(4, 0, false)]
    public void IsSelfPaired_OnlyForEvenNCenter(int N, int b, bool expected)
    {
        Assert.Equal(expected, BuildClaim().IsSelfPaired(N, b));
    }

    [Theory]
    [InlineData(2, true)]
    [InlineData(3, false)]   // odd N → no center bond
    [InlineData(4, true)]
    [InlineData(5, false)]
    [InlineData(6, true)]
    [InlineData(7, false)]
    [InlineData(8, true)]
    public void HasCenterBond_FollowsNParity(int N, bool expected)
    {
        Assert.Equal(expected, BuildClaim().HasCenterBond(N));
    }

    [Theory]
    [InlineData(2, 0)]   // (2-2)/2 = 0
    [InlineData(4, 1)]   // (4-2)/2 = 1
    [InlineData(6, 2)]
    [InlineData(8, 3)]
    public void CenterBondIndex_ReturnsExpectedValue(int N, int expected)
    {
        Assert.Equal(expected, BuildClaim().CenterBondIndex(N));
    }

    [Fact]
    public void CenterBondIndex_OddN_Throws()
    {
        Assert.Throws<ArgumentException>(() => BuildClaim().CenterBondIndex(3));
        Assert.Throws<ArgumentException>(() => BuildClaim().CenterBondIndex(5));
    }

    [Theory]
    [InlineData(2, 1)]   // 1 bond: bond 0 (self-paired since N-2-0 = 0)
    [InlineData(3, 2)]   // 2 bonds
    [InlineData(4, 3)]
    [InlineData(5, 4)]
    [InlineData(6, 5)]
    public void BondCount_IsNMinus1(int N, int expected)
    {
        Assert.Equal(expected, BuildClaim().BondCount(N));
    }

    [Theory]
    [InlineData(2, 0)]   // 1 bond, that's the center (self-paired); 0 disjoint pairs
    [InlineData(3, 1)]   // 2 bonds, 1 disjoint pair
    [InlineData(4, 1)]   // 3 bonds, 1 disjoint pair + 1 center
    [InlineData(5, 2)]   // 4 bonds, 2 disjoint pairs
    [InlineData(6, 2)]   // 5 bonds, 2 disjoint pairs + 1 center
    [InlineData(7, 3)]
    [InlineData(8, 3)]
    public void PalindromicPairCount_IsFloorBondsOver2(int N, int expected)
    {
        Assert.Equal(expected, BuildClaim().PalindromicPairCount(N));
    }

    [Theory]
    [InlineData(2, 1)]   // 1 component (the center bond)
    [InlineData(3, 1)]   // 1 component (the pair)
    [InlineData(4, 2)]   // 2 components: pair + center
    [InlineData(5, 2)]   // 2 components: 2 pairs
    [InlineData(6, 3)]
    [InlineData(7, 3)]
    [InlineData(8, 4)]
    [InlineData(9, 4)]
    public void IndependentComponentCount_IsFloorNOver2(int N, int expected)
    {
        Assert.Equal(expected, BuildClaim().IndependentComponentCount(N));
    }

    [Theory]
    [InlineData(2, true, 1)]    // 1 = a_1 (self-mirror pivot)
    [InlineData(3, true, 1)]    // 1 = a_1
    [InlineData(4, true, 0)]    // 2 = a_0
    [InlineData(5, true, 0)]    // 2 = a_0
    [InlineData(6, false, null)]   // 3 not on ladder
    [InlineData(7, false, null)]   // 3
    [InlineData(8, true, -1)]   // 4 = a_{-1}
    [InlineData(9, true, -1)]   // 4 = a_{-1}
    [InlineData(10, false, null)]  // 5
    public void LadderIndexForIndependentComponentCount_LandsOnAnchorOrNull(int N, bool onLadder, int? expectedIndex)
    {
        var f = BuildClaim();
        Assert.Equal(onLadder, f.IndependentComponentCountIsLadderAnchor(N));
        Assert.Equal(expectedIndex, f.LadderIndexForIndependentComponentCount(N));
    }

    [Fact]
    public void MirrorPair_NLessThan2_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().MirrorPair(1, 0));
    }

    [Fact]
    public void MirrorPair_BondOutOfRange_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().MirrorPair(5, -1));
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().MirrorPair(5, 4));
    }
}
