using System;
using System.Collections.Generic;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The CΨ envelope reading: the local maxima (parabolic-apex heights), whether the maxima
/// form a non-increasing sequence (predecessor semantics — the Envelope Theorem's claim), the first
/// and largest predecessor-rise, and the absorbing ¼ envelope fold. See
/// docs/superpowers/specs/2026-06-12-symphony-envelope-witness-design.md and
/// docs/proofs/PROOF_MONOTONICITY_CPSI.md.</summary>
public readonly record struct EnvelopeReading(
    IReadOnlyList<(double Time, double ApexValue)> Maxima,
    int RiseCount,
    bool IsNonIncreasing,
    double? FirstRiseTime,
    double MaxRiseMagnitude,
    double? EnvelopeFoldTime);

/// <summary>Pure envelope analysis of a CΨ(t) curve. Two pieces of rigor are baked in (the Fable
/// round-2 review): peak heights are compared at their PARABOLIC APEX (3-point fit, O(Δt⁴) error,
/// so ~1% real beating rises are not buried under ~24% grid clipping), and rises are PREDECESSOR
/// rises (each peak vs the one before — the theorem's "non-increasing sequence", NOT a running max).</summary>
public static class QuarterEnvelope
{
    /// <summary>Analyse a CΨ(t) curve into its envelope reading.</summary>
    /// <param name="cpsi">The CΨ samples; same length as <paramref name="tGrid"/>.</param>
    /// <param name="tGrid">The strictly-ascending time grid the samples sit on.</param>
    /// <param name="threshold">The ¼ boundary (default 0.25) for the envelope fold.</param>
    /// <param name="riseTol">A predecessor-rise counts only if it exceeds this (default 1e-9).</param>
    public static EnvelopeReading Of(double[] cpsi, double[] tGrid,
                                     double threshold = 0.25, double riseTol = 1e-9)
    {
        int n = cpsi.Length;

        // 1. local maxima indices: left endpoint if the curve starts descending; interior peaks
        //    (rise-or-flat in, strict fall out); never the right endpoint.
        var idx = new List<int>();
        if (n >= 2 && cpsi[0] > cpsi[1]) idx.Add(0);
        for (int i = 1; i < n - 1; i++)
            if (cpsi[i] >= cpsi[i - 1] && cpsi[i] > cpsi[i + 1]) idx.Add(i);

        // 2. parabolic-apex height for each maximum (endpoints keep their raw value).
        var maxima = new List<(double Time, double ApexValue)>(idx.Count);
        foreach (int i in idx) maxima.Add((tGrid[i], ParabolicApex(cpsi, i)));

        // 3. predecessor rises above riseTol.
        int riseCount = 0;
        double? firstRise = null;
        double maxRise = 0.0;
        for (int k = 1; k < maxima.Count; k++)
        {
            double d = maxima[k].ApexValue - maxima[k - 1].ApexValue;
            if (d > riseTol)
            {
                riseCount++;
                firstRise ??= maxima[k].Time;
                if (d > maxRise) maxRise = d;
            }
        }

        // 4. envelope fold: the LAST downward ¼-crossing after which the curve stays below ¼.
        double? fold = null;
        for (int s = 1; s < n; s++)
        {
            if (cpsi[s - 1] > threshold && cpsi[s] < threshold)
            {
                bool staysBelow = true;
                for (int k = s + 1; k < n; k++)
                    if (cpsi[k] >= threshold) { staysBelow = false; break; }
                if (staysBelow)
                {
                    // At most one downward crossing can satisfy stays-below (once below for good, there is
                    // no later down-crossing); the overwrite simply keeps that last absorbing one.
                    double frac = (cpsi[s - 1] - threshold) / (cpsi[s - 1] - cpsi[s]);
                    fold = tGrid[s - 1] + frac * (tGrid[s] - tGrid[s - 1]);
                }
            }
        }

        return new EnvelopeReading(maxima, riseCount, riseCount == 0, firstRise, maxRise, fold);
    }

    /// <summary>The apex value of the parabola through (i−1, i, i+1). Endpoints keep their raw value;
    /// a (near-)degenerate parabola falls back to the sample.</summary>
    private static double ParabolicApex(double[] v, int i)
    {
        if (i == 0 || i == v.Length - 1) return v[i];
        double y0 = v[i - 1], y1 = v[i], y2 = v[i + 1];
        double den = y0 - 2.0 * y1 + y2;
        if (Math.Abs(den) < 1e-15) return y1;
        double delta = 0.5 * (y0 - y2) / den;
        return y1 - 0.25 * (y0 - y2) * delta;
    }
}
