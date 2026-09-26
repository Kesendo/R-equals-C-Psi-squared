using RCPsiSquared.Core.Calibration;

namespace RCPsiSquared.Core.Tests.Calibration;

/// <summary>Tests that history coverage, high-switch flags, and descriptive
/// drift labels remain readings rather than submission advice.</summary>
public class LifecycleSummaryTests
{
    private static Lazy<IReadOnlyDictionary<int, QubitTimeline>> Marrakesh91d =>
        CalibrationFixtures.Marrakesh91d;

    [Fact]
    public void SoftBreakPath_OnMarrakesh91d_HasHistoryAndNoHighSwitchRate()
    {
        var s = LifecycleSummary.For(Marrakesh91d.Value, new[] { 48, 49, 50 });

        Assert.Equal(3, s.Qubits.Count);
        Assert.False(s.HasMissingHistory);
        Assert.NotEqual(DriftVerdict.DriftVolatile, s.DriftVerdict);
        Assert.NotEqual(DriftVerdict.InsufficientHistory, s.DriftVerdict);
        Assert.Equal(0, s.TwitchCount);
        Assert.True(s.HasSufficientHistoryAndNoHighSwitchRate);
    }

    [Fact]
    public void FrameworkSnapshotsPath_PreservesQ0MeasuredHistory()
    {
        var s = LifecycleSummary.For(Marrakesh91d.Value, new[] { 0, 1, 2 });
        var q0 = s.Qubits[0];

        Assert.Equal(0, q0.Qubit);
        Assert.Equal(LifecycleArchetype.PulseStable, q0.Archetype);
        Assert.True(q0.BelowRStarFraction > 0.95);
    }

    [Fact]
    public void Q126Q127_HasStableDescriptiveVerdict()
    {
        var s = LifecycleSummary.For(Marrakesh91d.Value, new[] { 126, 127 });

        Assert.True(s.AllStable);
        Assert.Equal(DriftVerdict.DriftStable, s.DriftVerdict);
        Assert.Equal(0, s.TwitchCount);
        Assert.Equal(0, s.LifecycleCount);
        Assert.Equal(2, s.StableCount);
        Assert.True(s.HasSufficientHistoryAndNoHighSwitchRate);
    }

    [Fact]
    public void MissingQubit_FailsTheHistoryAndSwitchBoolean()
    {
        var s = LifecycleSummary.For(Marrakesh91d.Value, new[] { 0, 9999 });

        Assert.Equal(LifecycleArchetype.InsufficientData, s.Qubits[1].Archetype);
        Assert.True(s.HasMissingHistory);
        Assert.Equal(DriftVerdict.InsufficientHistory, s.DriftVerdict);
        Assert.False(s.HasSufficientHistoryAndNoHighSwitchRate);
    }

    [Fact]
    public void HighSwitchQubit_FailsTheHistoryAndSwitchBoolean()
    {
        var hist = new Dictionary<int, QubitTimeline>
        {
            [0] = CalibrationFixtures.StableTimeline(0, 30, 100, 80),
            [1] = CalibrationFixtures.AlternatingTimeline(1, 30),
        };
        var s = LifecycleSummary.For(hist, new[] { 0, 1 });

        Assert.True(s.AnyTwitch);
        Assert.Equal(DriftVerdict.DriftVolatile, s.DriftVerdict);
        Assert.Equal(1, s.StableCount);
        Assert.Equal(1, s.TwitchCount);
        Assert.False(s.HasSufficientHistoryAndNoHighSwitchRate);
    }

    [Fact]
    public void SyntheticAllStable_IsDriftStable()
    {
        var hist = new Dictionary<int, QubitTimeline>
        {
            [0] = CalibrationFixtures.StableTimeline(0, 30, 100, 80),
            [1] = CalibrationFixtures.StableTimeline(1, 30, 100, 80),
        };
        var s = LifecycleSummary.For(hist, new[] { 0, 1 });

        Assert.True(s.AllStable);
        Assert.Equal(DriftVerdict.DriftStable, s.DriftVerdict);
    }

    [Fact]
    public void SoftBreakPath_IsDriftModerate_WithThreeLifecycleRows()
    {
        // [48, 49, 50] on the 91-day history: 48 Lifecycle (switch 12/90), 49 DriftySilent
        // (below R* on 9 of 91 days, r spread above 0.10), 50 Lifecycle (switch 6/90).
        var s = LifecycleSummary.For(Marrakesh91d.Value, new[] { 48, 49, 50 });

        Assert.Equal(LifecycleArchetype.Lifecycle, s.Qubits[0].Archetype);
        Assert.Equal(LifecycleArchetype.DriftySilent, s.Qubits[1].Archetype);
        Assert.Equal(LifecycleArchetype.Lifecycle, s.Qubits[2].Archetype);
        Assert.Equal(3, s.LifecycleCount);
        Assert.Equal(0, s.StableCount);
        Assert.Equal(DriftVerdict.DriftModerate, s.DriftVerdict);
        Assert.True(s.HasSufficientHistoryAndNoHighSwitchRate);
    }

    [Fact]
    public void InvalidEmptyPath_IsRejected()
    {
        Assert.Throws<ArgumentException>(() =>
            LifecycleSummary.For(Marrakesh91d.Value, Array.Empty<int>()));
    }

    [Fact]
    public void Headline_IsDescriptiveAndOneLine()
    {
        string h = LifecycleSummary.For(Marrakesh91d.Value, new[] { 126, 127 }).ToHeadline();

        Assert.Contains("[126, 127]", h);
        Assert.Contains("drift-stable", h);
        Assert.Contains("91 days", h);
        Assert.DoesNotContain('\n', h);
    }
}
