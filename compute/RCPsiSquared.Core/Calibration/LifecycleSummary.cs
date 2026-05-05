namespace RCPsiSquared.Core.Calibration;

/// <summary>Multi-day drift label for a candidate path. Distinct from
/// <see cref="RegimeVerdict"/>: that one classifies the snapshot regime;
/// this one classifies how stable the regime is over the history window.</summary>
public enum DriftVerdict
{
    DriftStable,
    DriftModerate,
    DriftVolatile,
    InsufficientHistory,
}

/// <summary>Kebab-case label conversion for <see cref="DriftVerdict"/>.</summary>
public static class DriftVerdictLabels
{
    public static string Label(this DriftVerdict v) => v switch
    {
        DriftVerdict.DriftStable => "drift-stable",
        DriftVerdict.DriftModerate => "drift-moderate",
        DriftVerdict.DriftVolatile => "drift-volatile",
        DriftVerdict.InsufficientHistory => "insufficient-history",
        _ => throw new ArgumentOutOfRangeException(nameof(v)),
    };
}

/// <summary>One row of a <see cref="LifecycleSummary"/>: a single path qubit
/// with its multi-day boundary statistics and archetype.</summary>
public sealed record QubitLifecycleStats(
    int Qubit,
    int DayCount,
    double RMean,
    double RStdDev,
    double CrossingFraction,
    double WalkRate,
    LifecycleArchetype Archetype);

/// <summary>Drift-aware companion to <see cref="RegimeSummary"/>: instead of
/// classifying a path on a single calibration snapshot, classifies it across
/// a multi-day history window. Answers the recurring pre-submit question
/// "is this path's regime classification *stable* over the experiment-window
/// timescale, or are some qubits drift-volatile?"
///
/// <para>The 2026-04-25 to 2026-04-30 Marrakesh review found Q5 lost 46% of its
/// T2 in 5 days, Q3 gained 30%, Q49 lost 23%. <see cref="RegimeSummary"/> only
/// sees one snapshot at a time; this workflow looks at the trajectory and
/// flags <see cref="LifecycleArchetype.Twitch"/> qubits before they jitter
/// the F88-Lens reading.</para>
///
/// <para>Loaded via <see cref="CalibrationHistory.Load"/> from the daily-
/// calibration CSV produced by <c>data/ibm_history/ibm_history_analysis.py
/// --mode collect</c>. 91-day Marrakesh anchor:
/// <c>data/ibm_history/results/ibm_marrakesh_history.csv</c>.</para></summary>
public sealed record LifecycleSummary(
    IReadOnlyList<int> Path,
    IReadOnlyList<QubitLifecycleStats> Qubits,
    int StableCount,
    int LifecycleCount,
    int TwitchCount,
    int InsufficientDataCount)
{
    /// <summary>True iff every path qubit is in one of the stable archetypes
    /// (PulseStable, SilentStable, ClassicStable). The bar that should be cleared
    /// before spending QPU minutes on a long-time-window experiment.</summary>
    public bool AllStable => StableCount == Qubits.Count;

    /// <summary>True iff at least one path qubit is in
    /// <see cref="LifecycleArchetype.Twitch"/>: the path will read different
    /// physics on different days at the same calibration score.</summary>
    public bool AnyTwitch => TwitchCount > 0;

    /// <summary>True iff the timeline coverage is incomplete (some qubits had
    /// fewer than 2 calibration days). Always interpret a summary with this
    /// flag set as exploratory.</summary>
    public bool HasMissingHistory => InsufficientDataCount > 0;

    /// <summary>Composite drift label, distinct from
    /// <see cref="RegimeSummary.Verdict"/>:
    /// <see cref="Calibration.DriftVerdict.InsufficientHistory"/> if any qubit
    /// lacks &gt; 1 calibration day (always wins);
    /// <see cref="Calibration.DriftVerdict.DriftVolatile"/> if any qubit is
    /// <see cref="LifecycleArchetype.Twitch"/>;
    /// <see cref="Calibration.DriftVerdict.DriftStable"/> if every qubit is
    /// in a stable archetype; otherwise
    /// <see cref="Calibration.DriftVerdict.DriftModerate"/>. Use
    /// <see cref="DriftVerdictLabels.Label"/> for kebab-case strings.</summary>
    public DriftVerdict DriftVerdict =>
        HasMissingHistory ? Calibration.DriftVerdict.InsufficientHistory :
        AnyTwitch ? Calibration.DriftVerdict.DriftVolatile :
        AllStable ? Calibration.DriftVerdict.DriftStable :
        Calibration.DriftVerdict.DriftModerate;

    /// <summary>Build the lifecycle summary for <paramref name="path"/> against
    /// a pre-loaded history map. Qubits missing from the history get
    /// <see cref="LifecycleArchetype.InsufficientData"/>; the call does not throw.</summary>
    public static LifecycleSummary For(
        IReadOnlyDictionary<int, QubitTimeline> history,
        IReadOnlyList<int> path)
    {
        if (path.Count == 0) throw new ArgumentException("path must not be empty", nameof(path));

        var qubits = new List<QubitLifecycleStats>(path.Count);
        int stable = 0, life = 0, twitch = 0, missing = 0;

        foreach (int qid in path)
        {
            if (!history.TryGetValue(qid, out var timeline) || timeline.Days.Count < 2)
            {
                qubits.Add(new QubitLifecycleStats(
                    Qubit: qid,
                    DayCount: timeline?.Days.Count ?? 0,
                    RMean: 0,
                    RStdDev: 0,
                    CrossingFraction: 0,
                    WalkRate: 0,
                    Archetype: LifecycleArchetype.InsufficientData));
                missing++;
                continue;
            }

            var arch = QubitLifecycle.Classify(timeline);
            qubits.Add(new QubitLifecycleStats(
                Qubit: qid,
                DayCount: timeline.Days.Count,
                RMean: QubitLifecycle.RMean(timeline),
                RStdDev: QubitLifecycle.RStdDev(timeline),
                CrossingFraction: QubitLifecycle.CrossingFraction(timeline),
                WalkRate: QubitLifecycle.WalkRate(timeline),
                Archetype: arch));

            switch (arch)
            {
                case LifecycleArchetype.PulseStable:
                case LifecycleArchetype.SilentStable:
                case LifecycleArchetype.ClassicStable:
                    stable++; break;
                case LifecycleArchetype.Lifecycle:
                case LifecycleArchetype.DriftySilent:
                    life++; break;
                case LifecycleArchetype.Twitch:
                    twitch++; break;
            }
        }

        return new LifecycleSummary(path, qubits, stable, life, twitch, missing);
    }

    /// <summary>One-line human-readable summary, e.g.
    /// <c>"path [48, 49, 50] | drift-stable (3 stable, 0 lifecycle, 0 twitch) over 91 days"</c>.</summary>
    public string ToHeadline()
    {
        int days = Qubits.Count == 0 ? 0 : Qubits.Max(q => q.DayCount);
        string pathStr = "[" + string.Join(", ", Path) + "]";
        return $"path {pathStr} | {DriftVerdict.Label()} ({StableCount} stable, "
             + $"{LifecycleCount} lifecycle, {TwitchCount} twitch) over {days} days";
    }
}
