using RCPsiSquared.Core.F71;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class F86F71MirrorPi2InheritanceTests
{
    private static F86F71MirrorPi2Inheritance BuildClaim()
    {
        var f71 = new F71MirrorSymmetryPi2Inheritance();
        var f86Link = new F86MirrorGeneralisationLink();
        return new F86F71MirrorPi2Inheritance(f71, f86Link);
    }

    [Fact]
    public void Tier_IsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived, BuildClaim().Tier);
    }

    [Theory]
    [InlineData(5, 0, 3)]    // N=5: b=0 ↔ N-2-b = 3
    [InlineData(5, 1, 2)]    // b=1 ↔ b=2
    [InlineData(6, 2, 2)]    // N=6: b=2 self-paired (N/2-1 = 2)
    [InlineData(7, 1, 4)]    // N=7: b=1 ↔ b=4
    public void MirrorPartnerBond_EqualsNMinus2MinusB(int N, int b, int expected)
    {
        Assert.Equal(expected, BuildClaim().MirrorPartnerBond(N, b));
    }

    [Theory]
    [InlineData(5, 0, false)]   // N=5 odd, no self-paired
    [InlineData(5, 2, false)]   // odd, no
    [InlineData(6, 2, true)]    // N=6 even, b=N/2-1=2 self-paired
    [InlineData(8, 3, true)]    // N=8 even, b=N/2-1=3 self-paired
    public void IsSelfPairedBond_TrueForCenterBondAtEvenN(int N, int b, bool expected)
    {
        Assert.Equal(expected, BuildClaim().IsSelfPairedBond(N, b));
    }

    [Fact]
    public void IsSelfPairedBond_NoneForOddN()
    {
        var f = BuildClaim();
        // N=5: bonds b=0..3 (= N-2). No b satisfies b = N-2-b.
        for (int b = 0; b <= 3; b++)
            Assert.False(f.IsSelfPairedBond(N: 5, b));
    }

    [Theory]
    [InlineData(5, 0)]
    [InlineData(5, 2)]
    [InlineData(7, 3)]
    [InlineData(10, 4)]
    public void MirrorPartnerSumHolds(int N, int b)
    {
        Assert.True(BuildClaim().MirrorPartnerSumHolds(N, b));
    }

    [Fact]
    public void MirrorPartnerBond_OutOfRange_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().MirrorPartnerBond(N: 5, b: -1));
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().MirrorPartnerBond(N: 5, b: 4));   // N-2 = 3 max
    }

    [Fact]
    public void MirrorPartnerBond_NLessThan2_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().MirrorPartnerBond(N: 1, b: 0));
    }

    [Fact]
    public void Constructor_NullF71_Throws()
    {
        var f86Link = new F86MirrorGeneralisationLink();
        Assert.Throws<ArgumentNullException>(() =>
            new F86F71MirrorPi2Inheritance(null!, f86Link));
    }

    [Fact]
    public void Constructor_NullF86Link_Throws()
    {
        var f71 = new F71MirrorSymmetryPi2Inheritance();
        Assert.Throws<ArgumentNullException>(() =>
            new F86F71MirrorPi2Inheritance(f71, null!));
    }
}
