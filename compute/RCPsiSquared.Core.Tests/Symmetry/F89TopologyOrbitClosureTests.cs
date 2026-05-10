using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class F89TopologyOrbitClosureTests
{
    private static F89TopologyOrbitClosure BuildClaim()
    {
        var ladder = new Pi2DyadicLadderClaim();
        var f70 = new F70DeltaNSelectionRulePi2Inheritance(ladder);
        var f72 = new F72BlockDiagonalPurityPi2Inheritance(ladder, f70);
        var f73 = new F73SpatialSumPurityClosurePi2Inheritance(ladder, f70, f72);
        var f71 = new F71MirrorSymmetryPi2Inheritance();
        return new F89TopologyOrbitClosure(f73, f71);
    }

    [Fact]
    public void Tier_IsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived, BuildClaim().Tier);
    }

    [Fact]
    public void Anchor_PointsToAnalyticalFormulasAndExperimentWriteup()
    {
        var anchor = BuildClaim().Anchor;
        Assert.Contains("docs/ANALYTICAL_FORMULAS.md F89", anchor);
        Assert.Contains("experiments/F89_TOPOLOGY_ORBIT_CLOSURE.md", anchor);
    }

    [Theory]
    [InlineData(2, 0.5)]
    [InlineData(3, 2.0 / 3.0)]
    [InlineData(4, 0.75)]
    [InlineData(7, 6.0 / 7.0)]
    [InlineData(11, 10.0 / 11.0)]
    public void S0ClosedForm_EqualsNMinusOneOverN(int n, double expected)
    {
        Assert.Equal(expected, F89TopologyOrbitClosure.S0ClosedForm(n), precision: 14);
    }

    [Fact]
    public void S0ClosedForm_RejectsNBelowTwo()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => F89TopologyOrbitClosure.S0ClosedForm(1));
        Assert.Throws<ArgumentOutOfRangeException>(() => F89TopologyOrbitClosure.S0ClosedForm(0));
    }

    [Fact]
    public void ChainTopologyClass_EmptyBonds_ReturnsEmptyArray()
    {
        Assert.Equal(Array.Empty<int>(), F89TopologyOrbitClosure.ChainTopologyClass(7, Array.Empty<int>()));
    }

    [Fact]
    public void ChainTopologyClass_SingleBondAtAnyPosition_IsOne()
    {
        for (int b = 0; b <= 5; b++)
            Assert.Equal(new[] { 1 }, F89TopologyOrbitClosure.ChainTopologyClass(7, new[] { b }));
    }

    [Theory]
    [InlineData(new[] { 0, 1 })]
    [InlineData(new[] { 1, 2 })]
    [InlineData(new[] { 2, 3 })]
    [InlineData(new[] { 3, 4 })]
    [InlineData(new[] { 4, 5 })]
    public void ChainTopologyClass_TwoAdjacentBonds_IsTwo(int[] bonds)
    {
        Assert.Equal(new[] { 2 }, F89TopologyOrbitClosure.ChainTopologyClass(7, bonds));
    }

    [Theory]
    [InlineData(new[] { 0, 2 })]
    [InlineData(new[] { 0, 3 })]
    [InlineData(new[] { 0, 5 })]
    [InlineData(new[] { 1, 3 })]
    [InlineData(new[] { 2, 4 })]
    [InlineData(new[] { 3, 5 })]
    public void ChainTopologyClass_TwoDisjointBonds_IsOneOne(int[] bonds)
    {
        Assert.Equal(new[] { 1, 1 }, F89TopologyOrbitClosure.ChainTopologyClass(7, bonds));
    }

    [Theory]
    [InlineData(new[] { 0, 1, 2 })]
    [InlineData(new[] { 1, 2, 3 })]
    [InlineData(new[] { 2, 3, 4 })]
    [InlineData(new[] { 3, 4, 5 })]
    public void ChainTopologyClass_ThreeAdjacentBonds_IsThree(int[] bonds)
    {
        Assert.Equal(new[] { 3 }, F89TopologyOrbitClosure.ChainTopologyClass(7, bonds));
    }

    [Theory]
    [InlineData(new[] { 0, 1, 3 })]
    [InlineData(new[] { 0, 1, 4 })]
    [InlineData(new[] { 0, 1, 5 })]
    [InlineData(new[] { 1, 2, 4 })]
    [InlineData(new[] { 2, 4, 5 })]
    public void ChainTopologyClass_PathTwoPlusIsolated_IsOneTwo(int[] bonds)
    {
        Assert.Equal(new[] { 1, 2 }, F89TopologyOrbitClosure.ChainTopologyClass(7, bonds));
    }

    [Fact]
    public void ChainTopologyClass_AllDisjointTriple_IsOneOneOne()
    {
        Assert.Equal(new[] { 1, 1, 1 }, F89TopologyOrbitClosure.ChainTopologyClass(7, new[] { 0, 2, 4 }));
    }

    [Fact]
    public void ChainTopologyClass_FullSixBondChain_IsSix()
    {
        Assert.Equal(new[] { 6 }, F89TopologyOrbitClosure.ChainTopologyClass(7, new[] { 0, 1, 2, 3, 4, 5 }));
    }

    [Fact]
    public void ChainTopologyClass_RejectsOutOfRangeBond()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            F89TopologyOrbitClosure.ChainTopologyClass(7, new[] { 6 }));
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            F89TopologyOrbitClosure.ChainTopologyClass(7, new[] { -1 }));
    }

    [Fact]
    public void ChainTopologyClass_RejectsDuplicateBonds()
    {
        Assert.Throws<ArgumentException>(() =>
            F89TopologyOrbitClosure.ChainTopologyClass(7, new[] { 0, 0, 1 }));
    }

    [Fact]
    public void ChainTopologyClass_RejectsNullBonds()
    {
        Assert.Throws<ArgumentNullException>(() =>
            F89TopologyOrbitClosure.ChainTopologyClass(7, null!));
    }

    [Fact]
    public void ChainTopologyClass_RejectsNBelowTwo()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            F89TopologyOrbitClosure.ChainTopologyClass(1, new[] { 0 }));
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            F89TopologyOrbitClosure.ChainTopologyClass(0, Array.Empty<int>()));
    }

    [Fact]
    public void AreInSameChainOrbit_TwoAdjacentPairsAtDifferentPositions_AreSameOrbit()
    {
        Assert.True(F89TopologyOrbitClosure.AreInSameChainOrbit(7, new[] { 0, 1 }, new[] { 4, 5 }));
        Assert.True(F89TopologyOrbitClosure.AreInSameChainOrbit(7, new[] { 1, 2 }, new[] { 3, 4 }));
    }

    [Fact]
    public void AreInSameChainOrbit_TwoDisjointAtDifferentGaps_AreSameOrbit()
    {
        Assert.True(F89TopologyOrbitClosure.AreInSameChainOrbit(7, new[] { 0, 2 }, new[] { 0, 5 }));
        Assert.True(F89TopologyOrbitClosure.AreInSameChainOrbit(7, new[] { 1, 3 }, new[] { 2, 5 }));
    }

    [Fact]
    public void AreInSameChainOrbit_AdjacentVsDisjoint_AreDifferentOrbits()
    {
        Assert.False(F89TopologyOrbitClosure.AreInSameChainOrbit(7, new[] { 0, 1 }, new[] { 0, 2 }));
        Assert.False(F89TopologyOrbitClosure.AreInSameChainOrbit(7, new[] { 2, 3 }, new[] { 0, 5 }));
    }

    [Fact]
    public void AreInSameChainOrbit_PathThreeVsAllDisjointTriple_AreDifferentOrbits()
    {
        Assert.False(F89TopologyOrbitClosure.AreInSameChainOrbit(7, new[] { 0, 1, 2 }, new[] { 0, 2, 4 }));
    }

    [Fact]
    public void S0AtN7IsSixSevenths_HoldsExactly()
    {
        Assert.True(BuildClaim().S0AtN7IsSixSevenths());
    }

    [Theory]
    [InlineData(7, 0)]
    [InlineData(7, 1)]
    [InlineData(7, 2)]
    [InlineData(11, 4)]
    public void F71MirrorIsInSameOrbit_HoldsForAllBonds(int n, int b)
    {
        Assert.True(BuildClaim().F71MirrorIsInSameOrbit(n, b));
    }

    [Fact]
    public void F73AnalogConsistent_HoldsExactly()
    {
        Assert.True(BuildClaim().F73AnalogConsistent());
    }

    [Fact]
    public void DisplayName_NamesTopologyOrbitClosure()
    {
        var name = BuildClaim().DisplayName;
        Assert.Contains("F89", name);
        Assert.Contains("topology", name, StringComparison.OrdinalIgnoreCase);
        Assert.Contains("orbit", name, StringComparison.OrdinalIgnoreCase);
    }
}
