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
    public void WhereRankFourIsRefused_TheOneLitSectorsTogether_NotTheOrderZeroSectorAlone()
    {
        // the second route for the frames without ZZ: simulations/kms_element_why_five_localised.py, part C (modulo one
        // prime, the minors whose entries lie in S; here the equations whose variables lie in S, replayed over ℚ(i));
        // the case with ZZ has this witness alone. d = 0 alone admits the shape, d ≥ 0 alone and d ≤ 0 alone admit it,
        // d ∈ {0, 1, −1} refuses it, the full system refuses it. The equation counts are this witness's own, pinned.
        var w = Witness.Value;
        var equations = new Dictionary<string, (int OrderZero, int Half, int Pair)>
        {
            ["elements keeping d, without ZZ"] = (16, 81, 194),
            ["elements keeping d, with ZZ at Δ = 1"] = (16, 81, 194),
            ["maps U keeping d, i.e. elements U·S_XY reversing d, without ZZ"] = (20, 93, 214),
        };
        foreach (var c in w.Cases)
        {
            Assert.Equal(new[] { "0", "0, 1, 2", "0, −1, −2", "0, 1, −1", "−2, −1, 0, 1, 2" }, c.Localisations.Select(l => l.Sectors));
            var orderZero = c.Localisations[0];
            Assert.Equal(8, orderZero.Coordinates);
            Assert.Equal(equations[c.Name].OrderZero, orderZero.EquationCount);
            Assert.Equal(8, orderZero.ChartCount);
            Assert.Equal(0, orderZero.ChartsSkipped);
            Assert.True(orderZero.PointEvidence, $"{c.Name}: d = 0 alone closed with a unit in every chart");
            var upper = c.Localisations[1];
            Assert.Equal(13, upper.Coordinates);
            Assert.Equal(equations[c.Name].Half, upper.EquationCount);
            Assert.Equal(32, upper.ChartCount);
            Assert.Equal(c.StageASurvivors.Count, upper.ChartsClosedWithoutUnit);   // the survivor chart, pinned against the second route above
            Assert.True(c.StageAReproduced, $"{c.Name}: Localise(0, 1, 2) does not reproduce stage A");
            var lower = c.Localisations[2];
            Assert.Equal(13, lower.Coordinates);
            Assert.Equal(upper.EquationCount, lower.EquationCount);                 // ρ ↦ ρ† carries one relaxation onto the other, sector by sector
            Assert.Equal(32, lower.ChartCount);
            Assert.True(lower.PointEvidence, $"{c.Name}: d ≤ 0 alone closed with a unit in every chart");
            Assert.Equal(1, lower.ChartsClosedWithoutUnit);                        // a regression pin only: the mirror's point can sit in another chart than stage A's, the chart strata depend on the coordinates
            var pair = c.Localisations[3];
            Assert.Equal(16, pair.Coordinates);
            Assert.Equal(equations[c.Name].Pair, pair.EquationCount);
            Assert.Equal(32, pair.ChartCount);
            Assert.True(pair.ProvenEmpty, $"{c.Name}: d in {{0, 1, −1}} proven empty in {pair.ChartsProvenEmpty} of {pair.ChartCount} charts");
            Assert.Empty(pair.PrimesSetAside);
            var all = c.Localisations[4];
            Assert.Equal(c.Dimension, all.Coordinates);
            Assert.True(all.ProvenEmpty);
            Assert.True(c.FifthTermLocalised, c.Name);
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
        Assert.All(w.Cases, c => Assert.False(c.FifthTermLocalised));
        Assert.All(w.Cases, c => Assert.Empty(c.Localisations));
    }
}
