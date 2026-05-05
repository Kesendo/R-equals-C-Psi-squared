using RCPsiSquared.Core.Calibration;

namespace RCPsiSquared.Core.Tests.Calibration;

/// <summary>Tests for <see cref="LifecycleSummary"/>, the path × multi-day-history
/// drift-aware audit. Anchored to the 91-day Marrakesh history.</summary>
public class LifecycleSummaryTests
{
    private static readonly Lazy<IReadOnlyDictionary<int, QubitTimeline>> Marrakesh91d = new(() =>
        CalibrationHistory.Load(Path.Combine(FindRepoRoot(),
            "data", "ibm_history", "results", "ibm_marrakesh_history.csv")));

    [Fact]
    public void SoftBreakPath_OnMarrakesh91d_DriftClassified()
    {
        var s = LifecycleSummary.For(Marrakesh91d.Value, new[] { 48, 49, 50 });
        Assert.Equal(3, s.Qubits.Count);
        Assert.False(s.HasMissingHistory);
        // Path [48, 49, 50] has 3 lifecycle qubits per the path-biography review;
        // accept either drift-stable or drift-moderate but not drift-volatile.
        Assert.NotEqual("drift-volatile", s.DriftVerdict);
        Assert.NotEqual("insufficient-history", s.DriftVerdict);
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
        Assert.Equal("drift-stable", s.DriftVerdict);
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
        Assert.Equal("insufficient-history", s.DriftVerdict);
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
            [0] = MakeStable(0, days: 30, t1: 100, t2: 80),
            [1] = MakeStable(1, days: 30, t1: 100, t2: 80),
        };
        var s = LifecycleSummary.For(hist, new[] { 0, 1 });
        Assert.True(s.AllStable);
        Assert.Equal("drift-stable", s.DriftVerdict);
    }

    [Fact]
    public void DriftVolatile_OneTwitchQubit_FlipsVerdict()
    {
        var hist = new Dictionary<int, QubitTimeline>
        {
            [0] = MakeStable(0, days: 30, t1: 100, t2: 80),     // SilentStable
            [1] = MakeAlternating(1, days: 30, t1: 100),         // Twitch
        };
        var s = LifecycleSummary.For(hist, new[] { 0, 1 });
        Assert.True(s.AnyTwitch);
        Assert.Equal("drift-volatile", s.DriftVerdict);
        Assert.Equal(1, s.StableCount);
        Assert.Equal(1, s.TwitchCount);
    }

    private static QubitTimeline MakeStable(int qid, int days, double t1, double t2)
    {
        var list = new List<CalibrationDay>(days);
        for (int i = 0; i < days; i++)
            list.Add(new CalibrationDay($"2026-01-{i + 1:D2}", t1, t2));
        return new QubitTimeline(qid, list);
    }

    private static QubitTimeline MakeAlternating(int qid, int days, double t1)
    {
        var list = new List<CalibrationDay>(days);
        for (int i = 0; i < days; i++)
        {
            double t2 = (i % 2 == 0) ? 20.0 : 80.0;
            list.Add(new CalibrationDay($"2026-01-{i + 1:D2}", t1, t2));
        }
        return new QubitTimeline(qid, list);
    }

    private static string FindRepoRoot()
    {
        var dir = new DirectoryInfo(AppContext.BaseDirectory);
        while (dir != null)
        {
            if (File.Exists(Path.Combine(dir.FullName, "MIRROR_THEORY.md"))
             && Directory.Exists(Path.Combine(dir.FullName, "compute")))
                return dir.FullName;
            dir = dir.Parent;
        }
        throw new InvalidOperationException(
            $"could not locate repository root starting from {AppContext.BaseDirectory}");
    }
}
