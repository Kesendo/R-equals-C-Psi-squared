using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>The least Schmidt rank of a thermal element on two sites, the exact half, at the rate point
/// γ↓ = (3, 2), γ↑ = (1, 1), J = 1, ZZ at Δ = 1. The counts are the second route's,
/// simulations/kms_element_least_rank.py (sympy's Gröbner over ℚ(i)): element spaces 36/28, 18 keeping d, 18/10
/// reversing it; 32 charts in each case with 320, 320 and 352 rank-4 equations; the one stage-A survivor t1 = t3 = 1.
/// The witness decides every chart exactly when built, which takes minutes; hence the trait.</summary>
[Trait("Category", "SLOW_LEASTRANK")]
public class ThermalElementLeastRankWitnessTests
{
    private static readonly Lazy<ThermalElementLeastRankWitness> Witness = new(() => new ThermalElementLeastRankWitness());

    [Fact]
    public void TheElementSpaces_HaveTheSecondRoutesDimensions()
    {
        var w = Witness.Value;
        Assert.Equal((36, 28), w.FullElementDimension);
        Assert.Equal(18, w.Cases[0].Dimension);
        Assert.Equal(18, w.Cases[1].Dimension);
        Assert.Equal((18, 10), w.ReversingDimension);
        Assert.Equal(new[] { 1, 2, 1 }, w.SteadyStateSplit);
        Assert.Equal(new[] { 1, 2, 1 }, w.SteadyStateSplitWithZz);
        Assert.Equal((true, true), w.SxyReflects);
        Assert.Equal((true, true), w.ZCorner);
    }

    [Fact]
    public void WithZz_NoElementReversesTheOrder()
    {
        var w = Witness.Value;
        Assert.Equal((false, true), w.OrderOneCoprimeToConjugate);
        Assert.Equal((true, false), w.ReversingTouchesOrderOne);
    }

    [Fact]
    public void RankFour_IsExcludedInEveryCase_ChartByChart()
    {
        var w = Witness.Value;
        Assert.Equal(new[] { 320, 320, 352 }, w.Cases.Select(c => c.EquationCount));
        foreach (var c in w.Cases)
        {
            Assert.True(c.PremisesHold, c.Name);
            Assert.Equal(32, c.ChartCount);
            Assert.Equal(0, c.ChartsSkipped);
            Assert.True(c.RankFourExcluded, $"{c.Name}: {c.ChartsProvenEmpty} of {c.ChartCount} charts proven empty");
            Assert.Equal(new[] { "t1 = 1, t3 = 1" }, c.StageASurvivors);
        }
    }

    [Fact]
    public void TheControl_AProductElementAtOneTemperature_SatisfiesEveryEquation()
    {
        Assert.True(Witness.Value.ControlHolds);
        Assert.Equal(170, Witness.Value.ControlEquationCount);   // the second route's count at (3, 3 | 1, 1)
    }

    [Fact]
    public void AtOneTemperature_ThePatternIsNotForced_AndNothingIsDecided()
    {
        // the steady state is the Gibbs product there: split (0, 1, 0), so rank 4 does not reduce to these equations
        var w = new ThermalElementLeastRankWitness(down: new long[] { 3, 3 }, up: new long[] { 1, 1 });
        Assert.Equal(new[] { 0, 1, 0 }, w.SteadyStateSplit);
        Assert.All(w.Cases, c => Assert.False(c.PremisesHold));
        Assert.All(w.Cases, c => Assert.False(c.RankFourExcluded));
    }
}
