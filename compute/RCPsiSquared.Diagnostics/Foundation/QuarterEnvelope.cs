using System;
using System.Collections.Generic;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>A finite CΨ-envelope reading: local maxima (parabolic-apex heights), predecessor rises,
/// and the last downward threshold crossing after which the sampled suffix stays below. This reader
/// describes one supplied window; it does not turn a finite zero-rise sample into a theorem. See
/// docs/superpowers/specs/2026-06-12-symphony-envelope-witness-design.md and
/// docs/proofs/PROOF_MONOTONICITY_CPSI.md.</summary>
public readonly record struct EnvelopeReading(
    IReadOnlyList<(double Time, double ApexValue)> Maxima,
    int RiseCount,
    bool IsNonIncreasing,
    double? FirstRiseTime,
    double MaxRiseMagnitude,
    double? LastStayBelowCrossingTime);

/// <summary>Pure envelope analysis of a CΨ(t) curve. Peak heights use a local three-point quadratic
/// fit on the supplied time coordinates. It is exact for a quadratic; for a generic smooth curve its
/// local error is O(h³), so refinement remains part of the reading. A flat sampled top uses the last
/// plateau sample as the local-maximum index; the fitted apex is therefore an estimator convention,
/// not an intrinsic feature. Rises compare each fitted maximum with its predecessor.</summary>
public static class QuarterEnvelope
{
    /// <summary>Analyse a CΨ(t) curve into its envelope reading.</summary>
    /// <param name="cpsi">The CΨ samples; same length as <paramref name="tGrid"/>.</param>
    /// <param name="tGrid">The strictly-ascending time grid the samples sit on.</param>
    /// <param name="threshold">The threshold (default 0.25) for the last-stay-below crossing.</param>
    /// <param name="riseTol">A predecessor-rise counts only if it exceeds this (default 1e-9).</param>
    public static EnvelopeReading Of(double[] cpsi, double[] tGrid,
                                     double threshold = 0.25, double riseTol = 1e-9)
    {
        ArgumentNullException.ThrowIfNull(cpsi);
        ArgumentNullException.ThrowIfNull(tGrid);
        if (cpsi.Length != tGrid.Length)
            throw new ArgumentException("CΨ samples and time grid must have equal length.");
        if (!double.IsFinite(threshold))
            throw new ArgumentOutOfRangeException(nameof(threshold), "threshold must be finite.");
        if (!double.IsFinite(riseTol) || riseTol < 0.0)
            throw new ArgumentOutOfRangeException(nameof(riseTol), "rise tolerance must be finite and non-negative.");

        var crossings = ThresholdCrossings.Extract(cpsi, tGrid, threshold);
        int n = cpsi.Length;

        // 1. local maxima indices: left endpoint if the curve starts descending; interior peaks
        //    (rise-or-flat in, strict fall out); never the right endpoint.
        var idx = new List<int>();
        if (n >= 2 && cpsi[0] > cpsi[1]) idx.Add(0);
        for (int i = 1; i < n - 1; i++)
            if (cpsi[i] >= cpsi[i - 1] && cpsi[i] > cpsi[i + 1]) idx.Add(i);

        // 2. fitted apex for each maximum (endpoints keep their raw time and value).
        var maxima = new List<(double Time, double ApexValue)>(idx.Count);
        foreach (int i in idx) maxima.Add(ParabolicApex(cpsi, tGrid, i));

        // 3. Raw predecessor order and maximum are independent of the reporting bar. Only the
        //    count and first reported time use riseTol.
        int riseCount = 0;
        double? firstRise = null;
        double maxRise = 0.0;
        bool isNonIncreasing = true;
        for (int k = 1; k < maxima.Count; k++)
        {
            double delta = maxima[k].ApexValue - maxima[k - 1].ApexValue;
            if (delta > 0.0)
            {
                isNonIncreasing = false;
                if (delta > maxRise) maxRise = delta;
            }
            if (delta > riseTol)
            {
                riseCount++;
                firstRise ??= maxima[k].Time;
            }
        }

        // 4. LAST downward crossing whose post-crossing sampled suffix is strictly below. The shared
        //    extractor collapses equality plateaux; its time is the last equality sample before exit.
        double? lastStayBelow = null;
        foreach (var crossing in crossings)
        {
            if (crossing.Direction >= 0) continue;
            bool staysBelow = true;
            for (int k = crossing.AfterIndex; k < n; k++)
                if (cpsi[k] >= threshold) { staysBelow = false; break; }
            if (staysBelow) lastStayBelow = crossing.Time;
        }

        return new EnvelopeReading(maxima, riseCount, isNonIncreasing, firstRise, maxRise, lastStayBelow);
    }

    /// <summary>The vertex of the quadratic through the actual three points (i−1, i, i+1).
    /// Endpoints, convex/linear fits, non-finite arithmetic, a vertex outside the bracket, or a
    /// fitted value below the selected sample fall back to that sample. No absolute curvature cutoff
    /// is used because curvature has units.</summary>
    private static (double Time, double ApexValue) ParabolicApex(double[] v, double[] t, int i)
    {
        if (i == 0 || i == v.Length - 1) return (t[i], v[i]);
        double x0 = t[i - 1], x1 = t[i], x2 = t[i + 1];
        double y0 = v[i - 1], y1 = v[i], y2 = v[i + 1];
        double s01 = (y1 - y0) / (x1 - x0);
        double s12 = (y2 - y1) / (x2 - x1);
        double a = (s12 - s01) / (x2 - x0);
        double h0 = x0 - x1;
        double b = s01 - a * h0;
        double u = -b / (2.0 * a);
        double apexTime = x1 + u;
        double apexValue = y1 + 0.5 * b * u;
        if (double.IsFinite(a) && double.IsFinite(b) && double.IsFinite(u)
            && double.IsFinite(apexTime) && double.IsFinite(apexValue)
            && a < 0.0 && x0 < apexTime && apexTime < x2 && apexValue >= y1)
            return (apexTime, apexValue);
        return (x1, y1);
    }
}

/// <summary>One threshold event with its time, direction, and first strict sample after the event.
/// Equality is a set/interval: an interior equality run is one event only when its strict neighbors
/// lie on opposite sides. Its scalar convention is the last equality sample before the exit edge.</summary>
internal readonly record struct ThresholdCrossing(double Time, int Direction, int AfterIndex);

/// <summary>The single structural threshold-event extractor shared by Symphony projections and the
/// envelope suffix reader.</summary>
internal static class ThresholdCrossings
{
    internal static IReadOnlyList<ThresholdCrossing> Extract(
        double[] values, double[] times, double threshold)
    {
        ArgumentNullException.ThrowIfNull(values);
        ArgumentNullException.ThrowIfNull(times);
        if (values.Length != times.Length)
            throw new ArgumentException("sample and time arrays must have equal length.");
        if (!double.IsFinite(threshold))
            throw new ArgumentOutOfRangeException(nameof(threshold), "threshold must be finite.");
        for (int i = 0; i < values.Length; i++)
        {
            if (!double.IsFinite(values[i]) || !double.IsFinite(times[i]))
                throw new ArgumentException("samples and times must be finite.");
            if (i > 0 && times[i] <= times[i - 1])
                throw new ArgumentException("time grid must be strictly increasing.", nameof(times));
        }

        var events = new List<ThresholdCrossing>();
        int index = 0;
        while (index < values.Length)
        {
            double current = values[index] - threshold;
            if (current == 0.0)
            {
                int firstEqual = index;
                int lastEqual = index;
                while (lastEqual + 1 < values.Length && values[lastEqual + 1] == threshold)
                    lastEqual++;
                if (firstEqual > 0 && lastEqual + 1 < values.Length)
                {
                    double left = values[firstEqual - 1] - threshold;
                    double right = values[lastEqual + 1] - threshold;
                    if ((left > 0.0 && right < 0.0) || (left < 0.0 && right > 0.0))
                        events.Add(new ThresholdCrossing(
                            times[lastEqual], left > 0.0 ? -1 : +1, lastEqual + 1));
                }
                index = lastEqual + 1;
                continue;
            }

            if (index + 1 < values.Length)
            {
                double next = values[index + 1] - threshold;
                if ((current > 0.0 && next < 0.0) || (current < 0.0 && next > 0.0))
                {
                    double absCurrent = Math.Abs(current);
                    double absNext = Math.Abs(next);
                    double scale = Math.Max(absCurrent, absNext);
                    double fraction = (absCurrent / scale) /
                        ((absCurrent / scale) + (absNext / scale));
                    double time = times[index] + fraction * (times[index + 1] - times[index]);
                    events.Add(new ThresholdCrossing(
                        time, current > 0.0 ? -1 : +1, index + 1));
                }
            }
            index++;
        }
        return events;
    }
}
