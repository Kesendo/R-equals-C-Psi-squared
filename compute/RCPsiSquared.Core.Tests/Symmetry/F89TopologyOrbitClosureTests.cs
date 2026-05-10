using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class F89TopologyOrbitClosureTests
{
    [Fact]
    public void Tier_IsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived, new F89TopologyOrbitClosure().Tier);
    }

    [Fact]
    public void Anchor_PointsToAnalyticalFormulasAndExperimentWriteup()
    {
        var anchor = new F89TopologyOrbitClosure().Anchor;
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
    public void DisplayName_NamesTopologyOrbitClosure()
    {
        var name = new F89TopologyOrbitClosure().DisplayName;
        Assert.Contains("F89", name);
        Assert.Contains("topology", name, StringComparison.OrdinalIgnoreCase);
        Assert.Contains("orbit", name, StringComparison.OrdinalIgnoreCase);
    }

    [Theory]
    [InlineData(7, 1)]
    [InlineData(7, 2)]
    [InlineData(7, 3)]
    [InlineData(11, 5)]
    public void AllIsolatedClosedForm_AtTZero_EqualsSZero(int n, int m)
    {
        double s = F89TopologyOrbitClosure.AllIsolatedClosedForm(n, m, j: 0.075, gammaZero: 0.05, t: 0.0);
        Assert.Equal((double)(n - 1) / n, s, precision: 14);
    }

    [Theory]
    [InlineData(1)]
    [InlineData(2)]
    [InlineData(3)]
    public void AllIsolatedClosedForm_AtInPhaseMoment_IndependentOfM(int m)
    {
        // cos(4Jt) = 1 at t = k·π/(2J); pick k=1 → t = π/(2J).
        const int n = 7;
        const double j = 0.075;
        const double gamma = 0.05;
        double tInPhase = Math.PI / (2.0 * j);

        double s = F89TopologyOrbitClosure.AllIsolatedClosedForm(n, m, j, gamma, tInPhase);
        double sExpected = (double)(n - 1) / n * Math.Exp(-4.0 * gamma * tInPhase);
        Assert.Equal(sExpected, s, precision: 12);
    }

    [Theory]
    [InlineData(7, 1, 0.0, 6.0 / 7.0)]
    [InlineData(7, 1, 21.0, 0.012853, 5)]   // ≈ in-phase moment
    [InlineData(7, 2, 21.0, 0.012853, 5)]   // ≈ in-phase moment
    [InlineData(7, 3, 21.0, 0.012853, 5)]   // ≈ in-phase moment
    public void AllIsolatedClosedForm_MatchesEmpiricalSpotCheck(int n, int m, double t, double expected, int precision = 12)
    {
        double s = F89TopologyOrbitClosure.AllIsolatedClosedForm(n, m, j: 0.075, gammaZero: 0.05, t: t);
        Assert.Equal(expected, s, precision: precision);
    }

    [Fact]
    public void AllIsolatedClosedForm_AsymptoticRateIs4Gamma_AcrossM()
    {
        // For two times t1, t2 both at cos(4Jt)=1, ratio S(t2)/S(t1) = exp(-4γ(t2-t1))
        // independent of m.
        const int n = 7;
        const double j = 0.075;
        const double gamma = 0.05;
        double t1 = Math.PI / (2.0 * j);
        double t2 = 3.0 * Math.PI / (2.0 * j);
        double expectedRatio = Math.Exp(-4.0 * gamma * (t2 - t1));

        for (int m = 1; m <= 3; m++)
        {
            double r = F89TopologyOrbitClosure.AllIsolatedClosedForm(n, m, j, gamma, t2)
                     / F89TopologyOrbitClosure.AllIsolatedClosedForm(n, m, j, gamma, t1);
            Assert.Equal(expectedRatio, r, precision: 12);
        }
    }

    [Fact]
    public void AllIsolatedClosedForm_RejectsBadInputs()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            F89TopologyOrbitClosure.AllIsolatedClosedForm(n: 1, m: 0, j: 1, gammaZero: 0.05, t: 0));
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            F89TopologyOrbitClosure.AllIsolatedClosedForm(n: 7, m: -1, j: 1, gammaZero: 0.05, t: 0));
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            F89TopologyOrbitClosure.AllIsolatedClosedForm(n: 7, m: 4, j: 1, gammaZero: 0.05, t: 0));  // 2m=8 > N=7
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            F89TopologyOrbitClosure.AllIsolatedClosedForm(n: 7, m: 1, j: 1, gammaZero: -0.01, t: 0));
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            F89TopologyOrbitClosure.AllIsolatedClosedForm(n: 7, m: 1, j: 1, gammaZero: 0.05, t: -1));
    }

    [Fact]
    public void AllIsolatedClosedForm_AtMZero_DegeneratesTo_S0_ExpDecay()
    {
        // m=0 means no bonds; pure dephasing on every site (all bare).
        // The closed form should give S(0)·exp(-4γt) with no oscillation.
        const int n = 7;
        const double j = 0.075;
        const double gamma = 0.05;
        double t = 5.0;
        double s = F89TopologyOrbitClosure.AllIsolatedClosedForm(n, m: 0, j, gamma, t);
        Assert.Equal((double)(n - 1) / n * Math.Exp(-4.0 * gamma * t), s, precision: 14);
    }
}
