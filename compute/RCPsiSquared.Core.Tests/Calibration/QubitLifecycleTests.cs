using RCPsiSquared.Core.Calibration;

namespace RCPsiSquared.Core.Tests.Calibration;

/// <summary>Tests for <see cref="QubitLifecycle"/>'s archetype classifier,
/// covering both synthetic timelines (one per archetype) and the Marrakesh
/// 91-day anchor (Q0 should land in PulseStable since the 2026-04-26 path
/// biography review showed walk 0.022 over the window).</summary>
public class QubitLifecycleTests
{
    private static readonly Lazy<IReadOnlyDictionary<int, QubitTimeline>> Marrakesh91d = new(() =>
        CalibrationHistory.Load(Path.Combine(FindRepoRoot(),
            "data", "ibm_history", "results", "ibm_marrakesh_history.csv")));

    [Fact]
    public void Classify_AlwaysQuantumSide_IsPulseStable()
    {
        var t = SyntheticTimeline(qid: 0, days: 30, t1: 100, t2: 30);  // r ≈ 0.15 < R*
        Assert.Equal(LifecycleArchetype.PulseStable, QubitLifecycle.Classify(t));
        Assert.Equal(0.0, QubitLifecycle.WalkRate(t));
        Assert.Equal(1.0, QubitLifecycle.CrossingFraction(t), precision: 6);
    }

    [Fact]
    public void Classify_AlwaysClassical_IsSilentStable()
    {
        var t = SyntheticTimeline(qid: 0, days: 30, t1: 100, t2: 80);  // r = 0.40 > R*
        Assert.Equal(LifecycleArchetype.SilentStable, QubitLifecycle.Classify(t));
        Assert.Equal(0.0, QubitLifecycle.WalkRate(t));
        Assert.Equal(0.0, QubitLifecycle.CrossingFraction(t));
    }

    [Fact]
    public void Classify_AlternatingDays_IsTwitch()
    {
        // 30 days alternating r=0.10 (quantum) and r=0.40 (classical) → walk = 1.0
        var days = new List<CalibrationDay>();
        for (int i = 0; i < 30; i++)
        {
            double t2 = (i % 2 == 0) ? 20.0 : 80.0;
            days.Add(new CalibrationDay($"2026-01-{i + 1:D2}", T1Us: 100, T2Us: t2));
        }
        var t = new QubitTimeline(0, days);
        Assert.Equal(LifecycleArchetype.Twitch, QubitLifecycle.Classify(t));
        Assert.True(QubitLifecycle.WalkRate(t) > QubitLifecycle.TwitchWalkThreshold);
    }

    [Fact]
    public void Classify_LongRunsThenSwitch_IsLifecycle()
    {
        // 15 days quantum-side, 15 days classical-side: walk = 1/29 ≈ 0.034.
        // To land in Lifecycle (walk > 0.05) we need a few more flips.
        var days = new List<CalibrationDay>();
        for (int i = 0; i < 30; i++)
        {
            // 4 flips over 30 days → walk ≈ 4/29 ≈ 0.138 (in lifecycle band)
            int phase = (i / 7) % 2;
            double t2 = phase == 0 ? 20.0 : 80.0;
            days.Add(new CalibrationDay($"2026-01-{i + 1:D2}", T1Us: 100, T2Us: t2));
        }
        var t = new QubitTimeline(0, days);
        double walk = QubitLifecycle.WalkRate(t);
        Assert.True(walk > QubitLifecycle.LifecycleWalkThreshold);
        Assert.True(walk <= QubitLifecycle.TwitchWalkThreshold);
        Assert.Equal(LifecycleArchetype.Lifecycle, QubitLifecycle.Classify(t));
    }

    [Fact]
    public void Classify_OneDay_IsInsufficientData()
    {
        var t = new QubitTimeline(0, new[] { new CalibrationDay("2026-01-01", 100, 50) });
        Assert.Equal(LifecycleArchetype.InsufficientData, QubitLifecycle.Classify(t));
    }

    [Fact]
    public void Classify_EmptyTimeline_IsInsufficientData()
    {
        var t = new QubitTimeline(0, Array.Empty<CalibrationDay>());
        Assert.Equal(LifecycleArchetype.InsufficientData, QubitLifecycle.Classify(t));
    }

    [Fact]
    public void Q0_OnMarrakesh91d_IsPulseStable()
    {
        var t = Marrakesh91d.Value[0];
        Assert.Equal(91, t.Days.Count);
        Assert.Equal(LifecycleArchetype.PulseStable, QubitLifecycle.Classify(t));
        Assert.True(QubitLifecycle.CrossingFraction(t) > 0.95,
            $"Q0 crossing fraction was {QubitLifecycle.CrossingFraction(t):F3}; expected > 0.95");
        Assert.True(QubitLifecycle.WalkRate(t) < QubitLifecycle.LifecycleWalkThreshold,
            $"Q0 walk rate was {QubitLifecycle.WalkRate(t):F3}; expected < {QubitLifecycle.LifecycleWalkThreshold}");
    }

    [Fact]
    public void Q126_Q127_OnMarrakesh91d_AreBothPulseStable()
    {
        // The only CZ-coupled stable-quantum pair on Marrakesh per the
        // 2026-05-05 uniform-quantum search.
        var h = Marrakesh91d.Value;
        Assert.Equal(LifecycleArchetype.PulseStable, QubitLifecycle.Classify(h[126]));
        Assert.Equal(LifecycleArchetype.PulseStable, QubitLifecycle.Classify(h[127]));
    }

    private static QubitTimeline SyntheticTimeline(int qid, int days, double t1, double t2)
    {
        var list = new List<CalibrationDay>(days);
        for (int i = 0; i < days; i++)
            list.Add(new CalibrationDay($"2026-01-{i + 1:D2}", t1, t2));
        return new QubitTimeline(qid, list);
    }

    private static string FindRepoRoot()
    {
        var dir = new DirectoryInfo(AppContext.BaseDirectory);
        while (dir != null)
        {
            if (Directory.Exists(Path.Combine(dir.FullName, "ClaudeTasks"))
             && File.Exists(Path.Combine(dir.FullName, "MIRROR_THEORY.md")))
                return dir.FullName;
            dir = dir.Parent;
        }
        throw new InvalidOperationException(
            $"could not locate repository root starting from {AppContext.BaseDirectory}");
    }
}
