namespace RCPsiSquared.Core.Calibration;

/// <summary>Empirical archetype for a multi-day trajectory of the calibration
/// ratio r around R*. These names summarize a time series; they do not assign
/// a universal physical phase.</summary>
public enum LifecycleArchetype
{
    /// <summary>Low switch rate and more than 70% of days below R*.</summary>
    PulseStable,
    /// <summary>Low switch rate, fewer than 10% below R*, and small r spread.</summary>
    SilentStable,
    /// <summary>Low switch rate and a middle-band fraction below R*.</summary>
    ClassicStable,
    /// <summary>Switch rate above 0.05 and at most 0.20.</summary>
    Lifecycle,
    /// <summary>Switch rate strictly above 0.20.</summary>
    Twitch,
    /// <summary>Low switch rate, fewer than 10% below R*, and large r spread.</summary>
    DriftySilent,
    /// <summary>Fewer than two calibration days.</summary>
    InsufficientData,
}

/// <summary>Pure descriptive statistics and an empirical classifier for a
/// <see cref="QubitTimeline"/>.</summary>
public static class QubitLifecycle
{
    public const double MostlyBelowRStarFractionThreshold = 0.7;
    public const double MaximumBelowRStarFractionForMostlyAbove = 0.1;
    public const double ModerateSwitchRateThreshold = 0.05;
    public const double HighSwitchRateThreshold = 0.20;
    public const double DriftySilentStdDev = 0.10;

    /// <summary>Mean r over the available days.</summary>
    public static double RMean(QubitTimeline timeline) =>
        timeline.Days.Count == 0 ? 0.0 : timeline.Days.Average(d => d.RParam);

    /// <summary>Population standard deviation of r over the available days.</summary>
    public static double RStdDev(QubitTimeline timeline)
    {
        if (timeline.Days.Count < 2) return 0.0;
        double mean = RMean(timeline);
        double sumSq = 0.0;
        foreach (var day in timeline.Days)
        {
            double delta = day.RParam - mean;
            sumSq += delta * delta;
        }
        return Math.Sqrt(sumSq / timeline.Days.Count);
    }

    /// <summary>Fraction of days for which r &lt; R*.</summary>
    public static double BelowRStarFraction(QubitTimeline timeline) =>
        timeline.Days.Count == 0
            ? 0.0
            : timeline.Days.Count(d => d.RStarBand == Regime.BelowRStar)
              / (double)timeline.Days.Count;

    /// <summary>Fraction of consecutive day-pairs whose binary R* band differs.</summary>
    public static double RStarBandSwitchRate(QubitTimeline timeline)
    {
        var days = timeline.Days;
        if (days.Count < 2) return 0.0;

        int switches = 0;
        Regime previous = days[0].RStarBand;
        for (int i = 1; i < days.Count; i++)
        {
            Regime current = days[i].RStarBand;
            if (current != previous) switches++;
            previous = current;
        }
        return switches / (double)(days.Count - 1);
    }

    /// <summary>Apply the pinned empirical thresholds. The comparisons are
    /// intentionally strict: 0.20 is Lifecycle, 0.05 remains on the stable
    /// branch, and fractions exactly 0.10 or 0.70 are ClassicStable.</summary>
    public static LifecycleArchetype Classify(QubitTimeline timeline)
    {
        if (timeline.Days.Count < 2) return LifecycleArchetype.InsufficientData;

        double switchRate = RStarBandSwitchRate(timeline);
        double belowFraction = BelowRStarFraction(timeline);

        if (switchRate > HighSwitchRateThreshold) return LifecycleArchetype.Twitch;
        if (switchRate > ModerateSwitchRateThreshold) return LifecycleArchetype.Lifecycle;
        if (belowFraction > MostlyBelowRStarFractionThreshold)
            return LifecycleArchetype.PulseStable;
        if (belowFraction < MaximumBelowRStarFractionForMostlyAbove)
            return RStdDev(timeline) < DriftySilentStdDev
                ? LifecycleArchetype.SilentStable
                : LifecycleArchetype.DriftySilent;
        return LifecycleArchetype.ClassicStable;
    }
}
