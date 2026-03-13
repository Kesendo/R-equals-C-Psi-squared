using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Complex;

namespace RCPsiSquared.Compute;

/// <summary>
/// Mirror symmetry analysis of decay rate spectra.
/// </summary>
public static class MirrorAnalysis
{
    public record MirrorResult(
        double Score,
        int Below,
        int Above,
        int AtCenter,
        int Matched,
        double Center
    );

    /// <summary>
    /// Check mirror symmetry of rates around a center value.
    /// Returns score from 0 (no symmetry) to 1 (perfect).
    /// </summary>
    public static MirrorResult CheckSymmetry(List<double> rates, double center, double tol = 0.005)
    {
        var below = rates.Where(r => r < center - tol).ToList();
        var above = rates.Where(r => r > center + tol).ToList();
        var atCenter = rates.Count(r => Math.Abs(r - center) < tol);

        int matched = 0;
        foreach (var r in below)
        {
            double mirror = 2 * center - r;
            double closest = above.Count > 0 ? above.MinBy(a => Math.Abs(a - mirror)) : 999;
            if (Math.Abs(closest - mirror) < tol)
                matched++;
        }

        double score = below.Count > 0 ? (double)matched / below.Count : 1.0;
        return new MirrorResult(score, below.Count, above.Count, atCenter, matched, center);
    }

    /// <summary>
    /// Full spectrum analysis: boundary check, mirror check, statistics.
    /// </summary>
    public static void Analyze(List<double> rates, int nQubits, double gamma, TextWriter output)
    {
        if (rates.Count == 0)
        {
            output.WriteLine("    No oscillatory rates found.");
            return;
        }

        double min = rates.Min();
        double max = rates.Max();
        double center = nQubits * gamma;
        var mirror = CheckSymmetry(rates, center);

        output.WriteLine($"    Rates: {rates.Count}");
        output.WriteLine($"    Min: {min:F6} ({min / gamma:F4}g)");
        output.WriteLine($"    Max: {max:F6} ({max / gamma:F4}g)");
        output.WriteLine($"    BW:  {(max - min) / gamma:F4}g");
        output.WriteLine($"    Predicted max 2(N-1)g: {2 * (nQubits - 1) * gamma:F6} ({2 * (nQubits - 1):F1}g)");
        output.WriteLine($"    Min matches 2g: {Math.Abs(min / gamma - 2) < 0.01}");
        output.WriteLine($"    Max matches 2(N-1)g: {Math.Abs(max / gamma - 2 * (nQubits - 1)) < 0.01}");
        output.WriteLine($"    Mirror symmetry: {mirror.Score:P1} ({mirror.Matched}/{mirror.Below})");
    }
}
