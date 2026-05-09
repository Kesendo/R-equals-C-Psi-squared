using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class F72BlockDiagonalPurityPi2InheritanceTests
{
    private static F72BlockDiagonalPurityPi2Inheritance BuildClaim()
    {
        var ladder = new Pi2DyadicLadderClaim();
        var f70 = new F70DeltaNSelectionRulePi2Inheritance(ladder);
        return new F72BlockDiagonalPurityPi2Inheritance(ladder, f70);
    }

    [Fact]
    public void Tier_IsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived, BuildClaim().Tier);
    }

    [Fact]
    public void BaselineTraceSquared_IsExactlyOneHalf()
    {
        // 1/2 = a_2 = HalfAsStructural; minimum purity for d = 2.
        Assert.Equal(0.5, BuildClaim().BaselineTraceSquared, precision: 14);
    }

    [Fact]
    public void SingleSiteBlockCount_IsExactlyTwo()
    {
        // 2 = a_0 = polynomial root; DD ⊕ CC at k=1.
        Assert.Equal(2.0, BuildClaim().SingleSiteBlockCount, precision: 14);
    }

    [Fact]
    public void BaselineMatchesHalfAsStructural_HoldsExactly()
    {
        Assert.True(BuildClaim().BaselineMatchesHalfAsStructural());
    }

    [Theory]
    [InlineData(1, 2)]    // DD + CC
    [InlineData(2, 3)]    // DD + DC + CC
    [InlineData(3, 4)]
    [InlineData(5, 6)]
    public void SubBlockCountForKLocal_IsKPlusOne(int kLocal, int expected)
    {
        Assert.Equal(expected, BuildClaim().SubBlockCountForKLocal(kLocal));
    }

    [Fact]
    public void SingleSiteBlockCountMatchesF70PlusOne_HoldsExactly()
    {
        // F70.PartialTraceMaxDeltaN(1) = 1; F72 sub-block count at k=1 = 1+1 = 2 = a_0.
        Assert.True(BuildClaim().SingleSiteBlockCountMatchesF70PlusOne());
    }

    [Fact]
    public void SingleSiteBlockNames_AreDdAndCc()
    {
        Assert.Equal(new[] { "DD", "CC" }, BuildClaim().SingleSiteBlockNames);
    }

    [Fact]
    public void PairBlockNames_IncludeDcCross()
    {
        // Pair (k=2) has DC cross specific to pair observables (per F72 generalisation).
        var names = BuildClaim().PairBlockNames;
        Assert.Equal(3, names.Count);
        Assert.Equal("DD", names[0]);
        Assert.Equal("DC", names[1]);
        Assert.Equal("CC", names[2]);
    }

    [Fact]
    public void SubBlockCountForKLocal_KLessThanOne_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().SubBlockCountForKLocal(0));
    }

    [Fact]
    public void Constructor_NullLadder_Throws()
    {
        var ladder = new Pi2DyadicLadderClaim();
        var f70 = new F70DeltaNSelectionRulePi2Inheritance(ladder);
        Assert.Throws<ArgumentNullException>(() =>
            new F72BlockDiagonalPurityPi2Inheritance(null!, f70));
    }

    [Fact]
    public void Constructor_NullF70_Throws()
    {
        var ladder = new Pi2DyadicLadderClaim();
        Assert.Throws<ArgumentNullException>(() =>
            new F72BlockDiagonalPurityPi2Inheritance(ladder, null!));
    }
}
