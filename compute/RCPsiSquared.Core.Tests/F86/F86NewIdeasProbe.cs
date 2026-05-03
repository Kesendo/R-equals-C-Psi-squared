using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.Decomposition;
using RCPsiSquared.Core.F86;
using RCPsiSquared.Core.Resonance;
using Xunit.Abstractions;

namespace RCPsiSquared.Core.Tests.F86;

/// <summary>Diagnostic probes for new F86 structural ideas. Each is a [Skip = "diagnostic"]
/// fact that can be unblocked to print exploration data into the test output. Promote to
/// real tests + typed claims once the pattern is verified.
/// </summary>
public class F86NewIdeasProbe(ITestOutputHelper output)
{
    [Fact(Skip = "diagnostic — full-grid scan; remove Skip to reproduce 2026-05-03 exploration data")]
    public void Probe_QPeakOverQEpRatio_AcrossWitnesses()
    {
        // Idea (a): Q_peak / Q_EP per (c, N, BondClass). Q_peak is from the K-curve
        // observable; Q_EP = 2/σ_0 from the inter-channel SVD. The ratio should be class-
        // specific (Interior ≠ Endpoint) and N-stable if it's a real structural fingerprint.
        output.WriteLine("c | N | bond     | Q_peak  | σ_0     | Q_EP    | Q_peak/Q_EP");
        output.WriteLine("--+---+----------+---------+---------+---------+-------------");

        foreach (var (c, N) in new[] { (2, 5), (2, 6), (2, 7), (2, 8),
                                       (3, 5), (3, 6), (3, 7), (3, 8) })
        {
            int n = c - 1;
            var block = new CoherenceBlock(N, n, gammaZero: 0.05);
            var svd = InterChannelSvd.Build(block, hd1: 1, hd2: 3);
            double qEp = 2.0 / svd.Sigma0;

            // Use a wide grid (40 points to 6.0) so Endpoint Q_peak (~2.5) is in range.
            var qGrid = ResonanceScan.LinearQGrid(0.20, 6.00, 40);
            var curve = new ResonanceScan(block).ComputeKCurve(qGrid);

            var interior = curve.Peak(BondClass.Interior);
            var endpoint = curve.Peak(BondClass.Endpoint);

            output.WriteLine($"{c} | {N} | Interior | {interior.QPeak,7:F4} | {svd.Sigma0,7:F4} | {qEp,7:F4} | {interior.QPeak / qEp,11:F4}");
            output.WriteLine($"{c} | {N} | Endpoint | {endpoint.QPeak,7:F4} | {svd.Sigma0,7:F4} | {qEp,7:F4} | {endpoint.QPeak / qEp,11:F4}");
        }
    }

    [Fact(Skip = "diagnostic — full-grid scan; remove Skip to reproduce 2026-05-03 exploration data")]
    public void Probe_Sigma0AcrossChromaticities()
    {
        // Idea (b): σ_0 trajectory for c=2, c=3, c=4. Look for pattern in the asymptote.
        // σ_0(c=2) → 2√2 ≈ 2.828 (already known). σ_0(c=3) and σ_0(c=4) asymptotes are open.
        output.WriteLine("c | N | n | σ_0     | σ_0/√(2(c−1))");
        output.WriteLine("--+---+---+---------+--------------");

        foreach (var (c, N) in new[] { (2, 5), (2, 6), (2, 7), (2, 8),
                                       (3, 5), (3, 6), (3, 7), (3, 8),
                                       (4, 7), (4, 8) })
        {
            int n = c - 1;
            var block = new CoherenceBlock(N, n, gammaZero: 0.05);
            var svd = InterChannelSvd.Build(block, hd1: 1, hd2: 3);
            double normalised = svd.Sigma0 / Math.Sqrt(2.0 * (c - 1));
            output.WriteLine($"{c} | {N} | {n} | {svd.Sigma0,7:F4} | {normalised,12:F4}");
        }
    }

    [Fact]
    public void Probe_InterChannelSvd_AllSingularValues_ForC2()
    {
        // The high-Q peak at b=1 might come from coupling to a non-top SVD mode.
        // For each c=2 N case, print all singular values + their corresponding
        // Q_EP = 2/σ_k to see if any σ_k matches the high-Q observation.
        output.WriteLine("c=2 InterChannelSvd singular values + corresponding Q_EP = 2/σ_k:");
        output.WriteLine("(N | σ_0 | σ_1 | σ_2 | ... | Q_EP_0 | Q_EP_1 | Q_EP_2 | ...)");
        foreach (int N in new[] { 5, 6, 7, 8 })
        {
            var block = new CoherenceBlock(N, 1, gammaZero: 0.05);
            var svd = RCPsiSquared.Core.Decomposition.InterChannelSvd.Build(block, hd1: 1, hd2: 3);
            var sigmas = svd.SingularValues.Where(s => s > 1e-10).ToArray();
            var sigmaStr = string.Join(", ", sigmas.Select(s => s.ToString("F4")));
            var qEpStr = string.Join(", ", sigmas.Select(s => (2.0 / s).ToString("F3")));
            output.WriteLine($"N={N}: σ = [{sigmaStr}]  →  Q_EP = [{qEpStr}]");
        }
    }

    [Fact]
    public void Probe_HighQFeature_c2N7_b1_ExtendedGrid()
    {
        // Investigate the c=2 N=7 high-Q observation: bonds b=1 and b=5 (orbit 1) showed
        // Q_peak >6 in the [0.2, 6.0] grid. Extend to [0.2, 20] to see if there's a real
        // peak at higher Q or if the K-curve is monotonic / has another structure.
        var block = new CoherenceBlock(N: 7, n: 1, gammaZero: 0.05);
        var qGrid = ResonanceScan.LinearQGrid(0.20, 20.0, 100);
        var curve = new ResonanceScan(block).ComputeKCurve(qGrid);

        output.WriteLine($"c=2 N=7, grid [0.2, 20.0] × 100 pts:");
        output.WriteLine("bond | Q_peak | K_max  | first 5 Q points where K decreases");

        for (int b = 0; b < 6; b++)
        {
            var peak = curve.PeakAtBond(b);
            var classCurve = new double[qGrid.Length];
            for (int q = 0; q < qGrid.Length; q++) classCurve[q] = curve.KByBond[b, q];

            // Find local maxima (where K[q-1] < K[q] > K[q+1])
            var localMaxima = new List<(double Q, double K)>();
            for (int q = 1; q < qGrid.Length - 1; q++)
            {
                if (classCurve[q - 1] < classCurve[q] && classCurve[q] > classCurve[q + 1])
                    localMaxima.Add((qGrid[q], classCurve[q]));
            }

            string maxList = string.Join(", ",
                localMaxima.Take(3).Select(m => $"({m.Q:F2}, K={m.K:F4})"));
            output.WriteLine($"  b{b} | {peak.QPeak:F3} | {peak.KMax:F4} | local maxima: {maxList}");
        }
    }

    [Fact(Skip = "diagnostic — full-grid scan; remove Skip + skip c=4 N=8 to avoid OOM. Data captured in PerF71OrbitObservation.")]
    public void Probe_PerBondQPeak_F71MirrorInvariance()
    {
        // F71 spatial-mirror: bond b and N−2−b are partners. Group bonds by F71 orbit
        // and print per-orbit Q_peak values to detect substructure (Endpoint vs Interior
        // vs central self-paired vs flanking).
        output.WriteLine("Per-F71-orbit Q_peak (γ₀ = 0.05, grid [0.2, 6.0] × 60 pts):");
        output.WriteLine("orbit indexed from outside (0 = endpoint pair) inward; * = self-paired");
        output.WriteLine("");

        var cases = new[]
        {
            (2, 5), (2, 6), (2, 7), (2, 8),
            (3, 5), (3, 6), (3, 7), (3, 8),
            (4, 7), (4, 8),
        };

        foreach (var (c, N) in cases)
        {
            int n = c - 1;
            var block = new CoherenceBlock(N, n, gammaZero: 0.05);
            var qGrid = ResonanceScan.LinearQGrid(0.20, 6.00, 60);
            var curve = new ResonanceScan(block).ComputeKCurve(qGrid);

            int numBonds = N - 1;
            int numOrbits = (numBonds + 1) / 2;
            var line = new System.Text.StringBuilder();
            line.Append($"c={c} N={N} ({numBonds} bonds, {numOrbits} orbits): ");
            for (int orbit = 0; orbit < numOrbits; orbit++)
            {
                int bA = orbit;
                int bB = numBonds - 1 - orbit;
                var peakA = curve.PeakAtBond(bA);
                if (bA == bB)
                {
                    line.Append($"orbit{orbit}*={peakA.QPeak:F3} ");
                }
                else
                {
                    var peakB = curve.PeakAtBond(bB);
                    line.Append($"orbit{orbit}={peakA.QPeak:F3}/{peakB.QPeak:F3} ");
                }
            }
            output.WriteLine(line.ToString());
        }
    }
}
