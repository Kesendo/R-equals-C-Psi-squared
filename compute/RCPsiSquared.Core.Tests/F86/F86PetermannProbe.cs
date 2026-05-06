using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.F86;
using RCPsiSquared.Core.Resonance;
using Xunit.Abstractions;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Tests.F86;

/// <summary>F86 Petermann-factor analysis at deep post-EP. The two-regime K-curve structure
/// at c=2 N≥7 (orbit-1 saturation plateau at Q ≈ 7-8) suggested a near-EP region; Petermann
/// factor K_n = 1/|⟨l_n|r_n⟩|² (unit-norm left/right eigvecs) diverges at exact EP.
///
/// <para>Connection to <see href="hypotheses/FRAGILE_BRIDGE.md">FRAGILE_BRIDGE</see>: there
/// the Hopf bifurcation in the complex γ plane signals an EP via Petermann K = 403 above
/// γ_crit. Same diagnostic, different parameter axis (Q vs γ).</para>
/// </summary>
public class F86PetermannProbe(ITestOutputHelper output)
{
    [Fact]
    public void Probe_PetermannArgMax_VsN_AtC2()
    {
        output.WriteLine("c=2: Q_EP_max from Petermann factor max K_n, vs N");
        output.WriteLine("(N | Q_EP_Interior | K_Interior | Q_EP_Endpoint | K_Endpoint)");

        foreach (int N in new[] { 5, 6, 7, 8 })
        {
            var block = new CoherenceBlock(N, 1, gammaZero: 0.05);
            var decomp = block.Decomposition;

            double bestQInt = 0, bestKInt = 0, bestQEnd = 0, bestKEnd = 0;
            var qGrid = ResonanceScan.LinearQGrid(0.5, 4.0, 70);
            foreach (double q in qGrid)
            {
                double j = q * block.GammaZero;
                var L = decomp.D + (Complex)j * decomp.MhTotal;
                double maxK = ComputeMaxPetermannK(L);

                // Interior region (Q < 2), Endpoint region (Q ≥ 2).
                if (q < 2.0)
                {
                    if (maxK > bestKInt) { bestKInt = maxK; bestQInt = q; }
                }
                else
                {
                    if (maxK > bestKEnd) { bestKEnd = maxK; bestQEnd = q; }
                }
            }

            output.WriteLine($"N={N}: Q_EP_Int={bestQInt:F3} (K={bestKInt:F0})  |  Q_EP_End={bestQEnd:F3} (K={bestKEnd:F0})");
        }
    }

    /// <summary>Finer-grid Petermann sweep across N=5..8 at c=2, with explicit Interior/Endpoint
    /// peak tracking. Counterpart to <see cref="Probe_PetermannArgMax_VsN_AtC2"/> with twice the
    /// Q-resolution near the empirical Q_peak ≈ 1.5-2.5 region. Used to decide whether K grows
    /// monotonically with N (real-axis hit of FRAGILE_BRIDGE's complex-γ EP) or stays bounded
    /// (off-axis siblings). The σ_0 R-even/R-odd degeneracy at even N (A3 finding) predicts a
    /// parity asymmetry in the K-spike pattern across N.
    /// </summary>
    [Fact(Skip = "Expensive ~30s probe; run via --filter when investigating local↔global EP link or refreshing LocalGlobalEpLink witnesses")]
    public void Probe_PetermannFineGrid_C2_VsN()
    {
        output.WriteLine($"c=2 Petermann fine-grid sweep, N=5..8");
        output.WriteLine($"Q-grid: {LocalGlobalEpLink.SweepQPoints} points uniform on " +
                         $"[{LocalGlobalEpLink.SweepQMin:F2}, {LocalGlobalEpLink.SweepQMax:F2}]");
        output.WriteLine("(N | dim | maxK_global | argmax_Q | maxK_Int | argmaxQ_Int | maxK_End | argmaxQ_End | parity)");
        output.WriteLine("--------------------------------------------------------------------------------------------");

        var qGrid = ResonanceScan.LinearQGrid(
            LocalGlobalEpLink.SweepQMin,
            LocalGlobalEpLink.SweepQMax,
            LocalGlobalEpLink.SweepQPoints);

        // Single-pass EVD: compute per-Q slice arrays once, then derive both the
        // header summary (max-K + Interior/Endpoint peaks) and the printed slices
        // from the same cached arrays. Halves the EVD work compared to the earlier
        // two-pass version.
        var slices = new Dictionary<int, double[]>();
        var summaries = new List<(int N, int Dim, double MaxKGlobal, double QAtMaxGlobal,
                                  double MaxKInt, double QAtMaxInt, double MaxKEnd, double QAtMaxEnd)>();

        foreach (int N in new[] { 5, 6, 7, 8 })
        {
            var block = new CoherenceBlock(N, 1, gammaZero: LocalGlobalEpLink.SweepGammaZero);
            var decomp = block.Decomposition;
            int dim = decomp.D.RowCount;

            double maxKGlobal = 0, qAtMaxGlobal = 0;
            double maxKInt = 0, qAtMaxInt = 0;
            double maxKEnd = 0, qAtMaxEnd = 0;
            var arr = new double[qGrid.Length];

            for (int iQ = 0; iQ < qGrid.Length; iQ++)
            {
                double q = qGrid[iQ];
                double j = q * block.GammaZero;
                var L = decomp.D + (Complex)j * decomp.MhTotal;
                double maxK = ComputeMaxPetermannK(L);
                arr[iQ] = maxK;

                if (maxK > maxKGlobal) { maxKGlobal = maxK; qAtMaxGlobal = q; }
                if (q < 2.0)
                {
                    if (maxK > maxKInt) { maxKInt = maxK; qAtMaxInt = q; }
                }
                else
                {
                    if (maxK > maxKEnd) { maxKEnd = maxK; qAtMaxEnd = q; }
                }
            }

            slices[N] = arr;
            summaries.Add((N, dim, maxKGlobal, qAtMaxGlobal, maxKInt, qAtMaxInt, maxKEnd, qAtMaxEnd));
        }

        // Pass 2: format the cached results — no recomputation.
        foreach (var s in summaries)
        {
            string parity = (s.N % 2 == 0) ? "even" : "odd";
            output.WriteLine(
                $"N={s.N} | dim={s.Dim,5} | maxK={s.MaxKGlobal,9:F1} at Q={s.QAtMaxGlobal:F3} " +
                $"| Int: K={s.MaxKInt,9:F1} at Q={s.QAtMaxInt:F3} " +
                $"| End: K={s.MaxKEnd,9:F1} at Q={s.QAtMaxEnd:F3} | {parity}");
        }

        output.WriteLine("");
        output.WriteLine("Per-Q maxK slices (every 4th grid point):");
        output.WriteLine("Q       | N=5         | N=6         | N=7         | N=8");
        output.WriteLine("--------+-------------+-------------+-------------+-------------");

        for (int iQ = 0; iQ < qGrid.Length; iQ += 4)
        {
            output.WriteLine(
                $"{qGrid[iQ],6:F3}  | {slices[5][iQ],11:F1} | {slices[6][iQ],11:F1} " +
                $"| {slices[7][iQ],11:F1} | {slices[8][iQ],11:F1}");
        }
    }

    [Fact]
    public void Probe_PetermannFactor_C2N7_VsQ()
    {
        var block = new CoherenceBlock(N: 7, n: 1, gammaZero: 0.05);
        var decomp = block.Decomposition;

        // Pre-EP, around canonical Q_EP, plateau region, deep post-EP.
        var qValues = new[]
        {
            0.3, 0.5, 0.71, 1.0, 1.5, 2.0, 3.0, 5.0,
            6.0, 6.5, 7.0, 7.2, 7.5, 8.0, 8.5, 9.0,
            10.0, 12.0, 15.0, 20.0,
        };

        output.WriteLine($"c=2 N=7, block-L dim {decomp.D.RowCount}");
        output.WriteLine("Q       | max K_n | mode index | top-3 K_n values");
        output.WriteLine("--------+---------+------------+--------------------");

        foreach (double q in qValues)
        {
            double j = q * block.GammaZero;
            var L = decomp.D + (Complex)j * decomp.MhTotal;
            var (maxK, iMax, allK) = ComputeAllPetermannK(L);

            var top3 = allK
                .Select((k, idx) => (K: k, Index: idx))
                .OrderByDescending(x => x.K)
                .Take(3)
                .ToArray();
            string top3Str = string.Join(", ", top3.Select(x => $"K[{x.Index}]={x.K:F1}"));

            output.WriteLine($"{q,6:F2}  | {maxK,7:F1} | {iMax,10} | {top3Str}");
        }
    }

    /// <summary>Petermann factor K_n = ‖r_n‖²·‖l_n‖² (in the EVD biorthogonal basis where
    /// ⟨l_n|r_n⟩ = 1) for every eigenmode of <paramref name="L"/>. Returns the per-mode
    /// array plus the (max, argmax) for convenience.</summary>
    private static (double Max, int IMax, double[] All) ComputeAllPetermannK(ComplexMatrix L)
    {
        var evd = L.Evd();
        var R = evd.EigenVectors;
        var Rinv = R.Inverse();
        int dim = R.RowCount;

        var k_n = new double[dim];
        for (int n = 0; n < dim; n++)
        {
            double rNorm2 = 0, lNorm2 = 0;
            for (int i = 0; i < dim; i++)
            {
                double rMag = R[i, n].Magnitude;
                rNorm2 += rMag * rMag;
                double lMag = Rinv[n, i].Magnitude;
                lNorm2 += lMag * lMag;
            }
            k_n[n] = rNorm2 * lNorm2;
        }

        int iMax = 0;
        for (int n = 1; n < dim; n++) if (k_n[n] > k_n[iMax]) iMax = n;
        return (k_n[iMax], iMax, k_n);
    }

    private static double ComputeMaxPetermannK(ComplexMatrix L) => ComputeAllPetermannK(L).Max;
}
