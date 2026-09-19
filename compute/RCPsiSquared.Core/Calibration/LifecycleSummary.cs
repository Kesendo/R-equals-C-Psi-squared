namespace RCPsiSquared.Core.Calibration;

/// <summary>Descriptive multi-day drift label for a path.</summary>
public enum DriftVerdict
{
    DriftStable,
    DriftModerate,
    DriftVolatile,
    InsufficientHistory,
}

/// <summary>Stable text labels for <see cref="DriftVerdict"/>.</summary>
public static class DriftVerdictLabels
{
    public static string Label(this DriftVerdict value) => value switch
    {
        DriftVerdict.DriftStable => "drift-stable",
        DriftVerdict.DriftModerate => "drift-moderate",
        DriftVerdict.DriftVolatile => "drift-volatile",
        DriftVerdict.InsufficientHistory => "insufficient-history",
        _ => throw new ArgumentOutOfRangeException(nameof(value)),
    };
}

/// <summary>One qubit's descriptive multi-day statistics.</summary>
public sealed record QubitLifecycleStats(
    int Qubit,
    int DayCount,
    double RMean,
    double RStdDev,
    double BelowRStarFraction,
    double RStarBandSwitchRate,
    LifecycleArchetype Archetype);

/// <summary>History coverage, switch behavior, and drift labels for a path.
/// The object reports what the supplied time series contains; it does not make
/// a hardware submission recommendation.</summary>
public sealed record LifecycleSummary(
    IReadOnlyList<int> Path,
    IReadOnlyList<QubitLifecycleStats> Qubits,
    int StableCount,
    int LifecycleCount,
    int TwitchCount,
    int InsufficientDataCount)
{
    public bool AllStable => StableCount == Qubits.Count;
    public bool AnyTwitch => TwitchCount > 0;
    public bool HasMissingHistory => InsufficientDataCount > 0;

    /// <summary>Composite descriptive label. Missing history takes precedence,
    /// followed by high-switch behavior, full stability, then moderate drift.</summary>
    public DriftVerdict DriftVerdict =>
        HasMissingHistory ? Calibration.DriftVerdict.InsufficientHistory :
        AnyTwitch ? Calibration.DriftVerdict.DriftVolatile :
        AllStable ? Calibration.DriftVerdict.DriftStable :
        Calibration.DriftVerdict.DriftModerate;

    /// <summary>True exactly when every path row has at least two days and no
    /// row exceeds the pinned high-switch-rate threshold.</summary>
    public bool HasSufficientHistoryAndNoHighSwitchRate =>
        !HasMissingHistory && !AnyTwitch;

    /// <summary>Build a descriptive summary. Missing rows are retained as
    /// <see cref="LifecycleArchetype.InsufficientData"/>.</summary>
    public static LifecycleSummary For(
        IReadOnlyDictionary<int, QubitTimeline> history,
        IReadOnlyList<int> path)
    {
        if (path.Count == 0)
            throw new ArgumentException("path must not be empty", nameof(path));

        var qubits = new List<QubitLifecycleStats>(path.Count);
        int stable = 0, lifecycle = 0, twitch = 0, missing = 0;

        foreach (int qid in path)
        {
            if (!history.TryGetValue(qid, out var timeline) || timeline.Days.Count < 2)
            {
                qubits.Add(new QubitLifecycleStats(
                    qid,
                    timeline?.Days.Count ?? 0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    LifecycleArchetype.InsufficientData));
                missing++;
                continue;
            }

            LifecycleArchetype archetype = QubitLifecycle.Classify(timeline);
            qubits.Add(new QubitLifecycleStats(
                qid,
                timeline.Days.Count,
                QubitLifecycle.RMean(timeline),
                QubitLifecycle.RStdDev(timeline),
                QubitLifecycle.BelowRStarFraction(timeline),
                QubitLifecycle.RStarBandSwitchRate(timeline),
                archetype));

            switch (archetype)
            {
                case LifecycleArchetype.PulseStable:
                case LifecycleArchetype.SilentStable:
                case LifecycleArchetype.ClassicStable:
                    stable++;
                    break;
                case LifecycleArchetype.Lifecycle:
                case LifecycleArchetype.DriftySilent:
                    lifecycle++;
                    break;
                case LifecycleArchetype.Twitch:
                    twitch++;
                    break;
            }
        }

        return new LifecycleSummary(path, qubits, stable, lifecycle, twitch, missing);
    }

    /// <summary>One-line descriptive summary for logs.</summary>
    public string ToHeadline()
    {
        int days = Qubits.Count == 0 ? 0 : Qubits.Max(q => q.DayCount);
        string path = "[" + string.Join(", ", Path) + "]";
        return $"path {path} | {DriftVerdict.Label()} ({StableCount} stable, "
             + $"{LifecycleCount} lifecycle, {TwitchCount} twitch) over {days} days";
    }
}
