using RCPsiSquared.Core.Calibration;

namespace RCPsiSquared.Core.Tests.Calibration;

/// <summary>Tests the empirical lifecycle classifier, including every strict
/// threshold boundary so equality cannot silently change branch, the Torino
/// calibration record the thresholds were set against, and the Marrakesh 91-day
/// anchor (Q0 PulseStable: below R* on 90 of 91 days, walk 2/90 ≈ 0.022, r mean
/// 0.086).</summary>
public class QubitLifecycleTests
{
    private static Lazy<IReadOnlyDictionary<int, QubitTimeline>> Marrakesh91d =>
        CalibrationFixtures.Marrakesh91d;

    [Fact]
    public void Classify_AlwaysBelow_IsPulseStable()
    {
        var t = CalibrationFixtures.StableTimeline(qid: 0, days: 30, t1Us: 100, t2Us: 30);
        Assert.Equal(LifecycleArchetype.PulseStable, QubitLifecycle.Classify(t));
        Assert.Equal(0.0, QubitLifecycle.RStarBandSwitchRate(t));
        Assert.Equal(1.0, QubitLifecycle.BelowRStarFraction(t), precision: 6);
    }

    [Fact]
    public void Classify_AlwaysAtOrAbove_IsSilentStable()
    {
        var t = CalibrationFixtures.StableTimeline(qid: 0, days: 30, t1Us: 100, t2Us: 80);
        Assert.Equal(LifecycleArchetype.SilentStable, QubitLifecycle.Classify(t));
        Assert.Equal(0.0, QubitLifecycle.RStarBandSwitchRate(t));
        Assert.Equal(0.0, QubitLifecycle.BelowRStarFraction(t));
    }

    [Fact]
    public void Classify_AlternatingDays_IsTwitch()
    {
        var t = CalibrationFixtures.AlternatingTimeline(qid: 0, days: 30);
        Assert.Equal(LifecycleArchetype.Twitch, QubitLifecycle.Classify(t));
        Assert.True(
            QubitLifecycle.RStarBandSwitchRate(t) > QubitLifecycle.HighSwitchRateThreshold);
    }

    [Fact]
    public void Classify_ModerateSwitchRate_IsLifecycle()
    {
        var days = new List<CalibrationDay>();
        for (int i = 0; i < 30; i++)
        {
            int phase = (i / 7) % 2;
            days.Add(new CalibrationDay($"2026-01-{i + 1:D2}", 100, phase == 0 ? 20 : 80));
        }

        var timeline = new QubitTimeline(0, days);
        double rate = QubitLifecycle.RStarBandSwitchRate(timeline);
        Assert.True(rate > QubitLifecycle.ModerateSwitchRateThreshold);
        Assert.True(rate <= QubitLifecycle.HighSwitchRateThreshold);
        Assert.Equal(LifecycleArchetype.Lifecycle, QubitLifecycle.Classify(timeline));
    }

    [Fact]
    public void SwitchRate_EqualityAtPointTwenty_IsLifecycleNotTwitch()
    {
        var timeline = TimelineFromBands(false, false, false, true, true, true);
        Assert.Equal(0.20, QubitLifecycle.RStarBandSwitchRate(timeline), precision: 14);
        Assert.Equal(LifecycleArchetype.Lifecycle, QubitLifecycle.Classify(timeline));
    }

    [Fact]
    public void SwitchRate_EqualityAtPointZeroFive_StaysInStableBranch()
    {
        var bands = Enumerable.Repeat(false, 10).Concat(Enumerable.Repeat(true, 11)).ToArray();
        var timeline = TimelineFromBands(bands);
        Assert.Equal(0.05, QubitLifecycle.RStarBandSwitchRate(timeline), precision: 14);
        Assert.Equal(LifecycleArchetype.ClassicStable, QubitLifecycle.Classify(timeline));
    }

    [Fact]
    public void BelowFraction_EqualityAtPointOne_StaysClassicStable()
    {
        var bands = Enumerable.Repeat(true, 4).Concat(Enumerable.Repeat(false, 36)).ToArray();
        var timeline = TimelineFromBands(bands);
        Assert.Equal(0.10, QubitLifecycle.BelowRStarFraction(timeline), precision: 14);
        Assert.True(QubitLifecycle.RStarBandSwitchRate(timeline) <= 0.05);
        Assert.Equal(LifecycleArchetype.ClassicStable, QubitLifecycle.Classify(timeline));
    }

    [Fact]
    public void BelowFraction_EqualityAtPointSeven_StaysClassicStable()
    {
        var bands = Enumerable.Repeat(true, 28).Concat(Enumerable.Repeat(false, 12)).ToArray();
        var timeline = TimelineFromBands(bands);
        Assert.Equal(0.70, QubitLifecycle.BelowRStarFraction(timeline), precision: 14);
        Assert.True(QubitLifecycle.RStarBandSwitchRate(timeline) <= 0.05);
        Assert.Equal(LifecycleArchetype.ClassicStable, QubitLifecycle.Classify(timeline));
    }

    [Fact]
    public void Classify_FewerThanTwoDays_IsInsufficientData()
    {
        var one = new QubitTimeline(0, new[] { new CalibrationDay("2026-01-01", 100, 50) });
        var none = new QubitTimeline(0, Array.Empty<CalibrationDay>());
        Assert.Equal(LifecycleArchetype.InsufficientData, QubitLifecycle.Classify(one));
        Assert.Equal(LifecycleArchetype.InsufficientData, QubitLifecycle.Classify(none));
    }

    [Fact]
    public void MarrakeshAnchors_ClassifyAsPulseStable()
    {
        var h = Marrakesh91d.Value;
        Assert.Equal(91, h[0].Days.Count);
        Assert.Equal(LifecycleArchetype.PulseStable, QubitLifecycle.Classify(h[0]));
        // Two band switches over 90 day-pairs, one day at or above R*: the counts, exactly.
        Assert.Equal(2.0 / 90, QubitLifecycle.RStarBandSwitchRate(h[0]));
        Assert.Equal(90 / 91.0, QubitLifecycle.BelowRStarFraction(h[0]));
        // Q126/Q127: the only CZ-coupled pair among Marrakesh's stably-below-R* qubits
        // (docs/BOTH_SIDES_VISIBLE.md).
        Assert.Equal(LifecycleArchetype.PulseStable, QubitLifecycle.Classify(h[126]));
        Assert.Equal(LifecycleArchetype.PulseStable, QubitLifecycle.Classify(h[127]));
    }

    [Theory]
    [InlineData(80, 0, LifecycleArchetype.PulseStable)]
    [InlineData(72, 58, LifecycleArchetype.Twitch)]
    [InlineData(98, 53, LifecycleArchetype.Twitch)]
    [InlineData(70, 53, LifecycleArchetype.Twitch)]
    [InlineData(68, 54, LifecycleArchetype.Twitch)]
    [InlineData(105, 15, LifecycleArchetype.Lifecycle)]
    public void TorinoCalibrationRecord_ReproducesTheNamedArchetypes(
        int qubit, int switches, LifecycleArchetype expected)
    {
        // The docs/BOTH_SIDES_VISIBLE.md examples the thresholds were set against:
        // walk 0.000 / 0.322 / 0.294 / 0.083 for Q80 / Q72 / Q98 / Q105, here as switch
        // counts over the 180 day-pairs of the 181-day window.
        var t = CalibrationFixtures.Torino181d.Value[qubit];
        Assert.Equal(181, t.Days.Count);
        Assert.Equal(switches / 180.0, QubitLifecycle.RStarBandSwitchRate(t));
        Assert.Equal(expected, QubitLifecycle.Classify(t));
    }

    private static QubitTimeline TimelineFromBands(params bool[] below)
    {
        var days = below.Select((isBelow, i) => new CalibrationDay(
            $"2026-{i / 28 + 1:D2}-{i % 28 + 1:D2}",
            T1Us: 100,
            T2Us: isBelow ? 20 : 80)).ToArray();
        return new QubitTimeline(0, days);
    }
}
