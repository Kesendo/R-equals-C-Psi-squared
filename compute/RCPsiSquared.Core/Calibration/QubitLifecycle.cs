namespace RCPsiSquared.Core.Calibration;

/// <summary>Coarse archetype label for a qubit's multi-day boundary trajectory.
/// Distinct from the snapshot-only <see cref="Regime"/> classification
/// (<see cref="QubitRegime"/>): an archetype captures whether a qubit *stays* on
/// one side of R*, drifts slowly, or twitches rapidly.</summary>
public enum LifecycleArchetype
{
    /// <summary>Walk &lt; <see cref="QubitLifecycle.LifecycleWalkThreshold"/>,
    /// crossing &gt; <see cref="QubitLifecycle.PulseStableCrossing"/>: stably
    /// quantum-side over the whole window. The Q80 archetype.</summary>
    PulseStable,
    /// <summary>Walk &lt; <see cref="QubitLifecycle.LifecycleWalkThreshold"/>,
    /// crossing &lt; <see cref="QubitLifecycle.SilentStableCrossing"/>, low
    /// std: stably classical-side, never approaches the boundary.</summary>
    SilentStable,
    /// <summary>Walk &lt; <see cref="QubitLifecycle.LifecycleWalkThreshold"/>,
    /// crossing in the middle band: stably classical but near the boundary.</summary>
    ClassicStable,
    /// <summary>Walk in [<see cref="QubitLifecycle.LifecycleWalkThreshold"/>,
    /// <see cref="QubitLifecycle.TwitchWalkThreshold"/>): drifts slowly across
    /// the boundary with long runs on each side. The Q105 archetype.</summary>
    Lifecycle,
    /// <summary>Walk &gt; <see cref="QubitLifecycle.TwitchWalkThreshold"/>:
    /// rapidly flipping daily. The Q72/Q98 archetype.</summary>
    Twitch,
    /// <summary>Walk ≤ <see cref="QubitLifecycle.LifecycleWalkThreshold"/>,
    /// crossing &lt; <see cref="QubitLifecycle.SilentStableCrossing"/>, and
    /// std ≥ <see cref="QubitLifecycle.DriftySilentStdDev"/>: mostly above
    /// the boundary, with large excursions when it does deviate.</summary>
    DriftySilent,
    /// <summary>Fewer than 2 calibration days available, walk rate undefined.</summary>
    InsufficientData,
}

/// <summary>Multi-day boundary-trajectory classifier. Operates on a
/// <see cref="QubitTimeline"/> and returns the qubit's
/// <see cref="LifecycleArchetype"/> plus its component statistics.
///
/// <para>Calibrated against the BOTH_SIDES_VISIBLE.md named examples on the
/// 180-day Torino history: Q80 (consistent crosser, walk 0.000) →
/// <see cref="LifecycleArchetype.PulseStable"/>; Q72 (rhythmic, walk 0.322)
/// and Q98 (lifecycle, walk 0.294) → <see cref="LifecycleArchetype.Twitch"/>;
/// Q105 (long active then silent, walk 0.083) →
/// <see cref="LifecycleArchetype.Lifecycle"/>. The day-level twitch reading
/// of Q98 is consistent with the doc's week-level lifecycle reading; the two
/// scales describe the same trajectory.</para></summary>
public static class QubitLifecycle
{
    /// <summary>Walk rate above which a qubit is classified as
    /// <see cref="LifecycleArchetype.Twitch"/>. 0.20 reproduces Q72/Q98/Q70/Q68
    /// from BOTH_SIDES_VISIBLE.md as Twitch and Q80 as PulseStable.</summary>
    public const double TwitchWalkThreshold = 0.20;

    /// <summary>Walk rate above which a qubit is at least
    /// <see cref="LifecycleArchetype.Lifecycle"/> (slow drift across the
    /// boundary). 0.05 reproduces Q105 from BOTH_SIDES_VISIBLE.md as Lifecycle.</summary>
    public const double LifecycleWalkThreshold = 0.05;

    /// <summary>Crossing fraction above which a stable (low-walk) qubit is
    /// <see cref="LifecycleArchetype.PulseStable"/>.</summary>
    public const double PulseStableCrossing = 0.7;

    /// <summary>Crossing fraction below which a stable (low-walk) qubit is
    /// either <see cref="LifecycleArchetype.SilentStable"/> or
    /// <see cref="LifecycleArchetype.DriftySilent"/>.</summary>
    public const double SilentStableCrossing = 0.1;

    /// <summary>Std-dev threshold separating <see cref="LifecycleArchetype.SilentStable"/>
    /// from <see cref="LifecycleArchetype.DriftySilent"/> in the low-walk
    /// low-crossing regime.</summary>
    public const double DriftySilentStdDev = 0.10;

    /// <summary>Mean of r = T2/(2·T1) over the timeline.</summary>
    public static double RMean(QubitTimeline timeline) =>
        timeline.Days.Count == 0 ? 0.0 : timeline.Days.Average(d => d.RParam);

    /// <summary>Population std-dev of r over the timeline.</summary>
    public static double RStdDev(QubitTimeline timeline)
    {
        if (timeline.Days.Count < 2) return 0.0;
        double mean = RMean(timeline);
        double sumSq = 0;
        foreach (var d in timeline.Days) { double diff = d.RParam - mean; sumSq += diff * diff; }
        return Math.Sqrt(sumSq / timeline.Days.Count);
    }

    /// <summary>Fraction of days the qubit was in the
    /// <see cref="Regime.QuantumSide"/> regime (r &lt; R*).</summary>
    public static double CrossingFraction(QubitTimeline timeline) =>
        timeline.Days.Count == 0
            ? 0.0
            : timeline.Days.Count(d => d.Regime == Regime.QuantumSide)
              / (double)timeline.Days.Count;

    /// <summary>Walk rate: fraction of consecutive day-pairs where the qubit
    /// flipped across R* (i.e. the <see cref="Regime"/> classification of
    /// <see cref="CalibrationDay.Regime"/> changed). Primary empirical indicator
    /// of boundary volatility.</summary>
    public static double WalkRate(QubitTimeline timeline)
    {
        var days = timeline.Days;
        if (days.Count < 2) return 0.0;
        int flips = 0;
        Regime prev = days[0].Regime;
        for (int i = 1; i < days.Count; i++)
        {
            Regime curr = days[i].Regime;
            if (curr != prev) flips++;
            prev = curr;
        }
        return (double)flips / (days.Count - 1);
    }

    /// <summary>Classify the qubit's lifecycle archetype over the timeline window.</summary>
    public static LifecycleArchetype Classify(QubitTimeline timeline)
    {
        if (timeline.Days.Count < 2) return LifecycleArchetype.InsufficientData;
        double walk = WalkRate(timeline);
        double crossing = CrossingFraction(timeline);

        if (walk > TwitchWalkThreshold) return LifecycleArchetype.Twitch;
        if (walk > LifecycleWalkThreshold) return LifecycleArchetype.Lifecycle;
        if (crossing > PulseStableCrossing) return LifecycleArchetype.PulseStable;
        if (crossing < SilentStableCrossing)
            return RStdDev(timeline) < DriftySilentStdDev
                ? LifecycleArchetype.SilentStable
                : LifecycleArchetype.DriftySilent;
        return LifecycleArchetype.ClassicStable;
    }
}
