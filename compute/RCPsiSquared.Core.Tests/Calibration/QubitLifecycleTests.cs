using RCPsiSquared.Core.Calibration;

namespace RCPsiSquared.Core.Tests.Calibration;

/// <summary>Tests the empirical lifecycle classifier, including every strict
/// threshold boundary so equality cannot silently change branch.</summary>
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
    public void MarrakeshAnchors_PreserveMeasuredLifecycleRows()
    {
        var h = Marrakesh91d.Value;
        Assert.Equal(91, h[0].Days.Count);
        Assert.Equal(LifecycleArchetype.PulseStable, QubitLifecycle.Classify(h[0]));
        Assert.True(QubitLifecycle.BelowRStarFraction(h[0]) > 0.95);
        Assert.True(
            QubitLifecycle.RStarBandSwitchRate(h[0]) < QubitLifecycle.ModerateSwitchRateThreshold);
        Assert.Equal(LifecycleArchetype.PulseStable, QubitLifecycle.Classify(h[126]));
        Assert.Equal(LifecycleArchetype.PulseStable, QubitLifecycle.Classify(h[127]));
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
