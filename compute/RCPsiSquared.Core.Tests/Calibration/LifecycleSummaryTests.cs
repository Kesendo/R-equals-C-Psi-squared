using RCPsiSquared.Core.Calibration;

namespace RCPsiSquared.Core.Tests.Calibration;

/// <summary>Tests for <see cref="LifecycleSummary"/>, the path × multi-day-history
/// drift-aware audit. Anchored to the 91-day Marrakesh history.</summary>
public class LifecycleSummaryTests
{
    private static Lazy<IReadOnlyDictionary<int, QubitTimeline>> Marrakesh91d => CalibrationFixtures.Marrakesh91d;

    [Fact]
    public void SoftBreakPath_OnMarrakesh91d_DriftClassified()
    {
        var s = LifecycleSummary.For(Marrakesh91d.Value, new[] { 48, 49, 50 });
        Assert.Equal(3, s.Qubits.Count);
        Assert.False(s.HasMissingHistory);
        // Path [48, 49, 50] has 3 lifecycle qubits per the path-biography review;
        // accept either drift-stable or drift-moderate but not drift-volatile.
        Assert.NotEqual(DriftVerdict.DriftVolatile, s.DriftVerdict);
        Assert.NotEqual(DriftVerdict.InsufficientHistory, s.DriftVerdict);
        Assert.Equal(0, s.TwitchCount);
    }

    [Fact]
    public void FrameworkSnapshotsPath_OnMarrakesh91d_HasQ0PulseStable()
    {
        var s = LifecycleSummary.For(Marrakesh91d.Value, new[] { 0, 1, 2 });
        var q0Stats = s.Qubits[0];
        Assert.Equal(0, q0Stats.Qubit);
        Assert.Equal(LifecycleArchetype.PulseStable, q0Stats.Archetype);
        Assert.True(q0Stats.CrossingFraction > 0.95);
    }

    [Fact]
    public void Q126_Q127_PathOnMarrakesh91d_IsAllStable()
    {
        var s = LifecycleSummary.For(Marrakesh91d.Value, new[] { 126, 127 });
        Assert.True(s.AllStable);
        Assert.Equal(DriftVerdict.DriftStable, s.DriftVerdict);
        Assert.Equal(0, s.TwitchCount);
        Assert.Equal(0, s.LifecycleCount);
        Assert.Equal(2, s.StableCount);
    }

    [Fact]
    public void MissingQubit_FlaggedAsInsufficientData()
    {
        var s = LifecycleSummary.For(Marrakesh91d.Value, new[] { 0, 9999 });
        Assert.Equal(2, s.Qubits.Count);
        Assert.Equal(LifecycleArchetype.InsufficientData, s.Qubits[1].Archetype);
        Assert.True(s.HasMissingHistory);
        Assert.Equal(DriftVerdict.InsufficientHistory, s.DriftVerdict);
    }

    [Fact]
    public void EmptyPath_Throws()
    {
        Assert.Throws<ArgumentException>(() =>
            LifecycleSummary.For(Marrakesh91d.Value, Array.Empty<int>()));
    }

    [Fact]
    public void Headline_IsHumanReadableOneLine()
    {
        var s = LifecycleSummary.For(Marrakesh91d.Value, new[] { 126, 127 });
        string h = s.ToHeadline();
        Assert.Contains("[126, 127]", h);
        Assert.Contains("drift-stable", h);
        Assert.Contains("91 days", h);
    }

    [Fact]
    public void DriftStable_PerSyntheticPath_StableQubits()
    {
        // All-stable synthetic history: 30 days, all r ≈ 0.4 (classical-side, no flips).
        var hist = new Dictionary<int, QubitTimeline>
        {
            [0] = CalibrationFixtures.StableTimeline(0, days: 30, t1Us: 100, t2Us: 80),
            [1] = CalibrationFixtures.StableTimeline(1, days: 30, t1Us: 100, t2Us: 80),
        };
        var s = LifecycleSummary.For(hist, new[] { 0, 1 });
        Assert.True(s.AllStable);
        Assert.Equal(DriftVerdict.DriftStable, s.DriftVerdict);
    }

    [Fact]
    public void DriftVolatile_OneTwitchQubit_FlipsVerdict()
    {
        var hist = new Dictionary<int, QubitTimeline>
        {
            [0] = CalibrationFixtures.StableTimeline(0, days: 30, t1Us: 100, t2Us: 80),  // SilentStable
            [1] = CalibrationFixtures.AlternatingTimeline(1, days: 30),                  // Twitch
        };
        var s = LifecycleSummary.For(hist, new[] { 0, 1 });
        Assert.True(s.AnyTwitch);
        Assert.Equal(DriftVerdict.DriftVolatile, s.DriftVerdict);
        Assert.Equal(1, s.StableCount);
        Assert.Equal(1, s.TwitchCount);
    }

}
