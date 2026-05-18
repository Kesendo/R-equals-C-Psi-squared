using System.Numerics;
using System.Text.Json;
using RCPsiSquared.Core.BlockSpectrum;

namespace RCPsiSquared.Core.F1;

/// <summary>F1 spectrum-statistics utility for the N=8 / N=9 F1 palindromic-pairing
/// dogfood tests. Captures five metric groups in a single immutable record:
/// wall-time profile (Group 1), pairing-precision histogram (Group 2),
/// spectrum-structure invariants (Group 3), block-decomposition cost picture
/// (Group 4), and the Hamiltonian / dissipator reproducibility anchor (Group 5).
///
/// <para>The record serialises to pretty-printed JSON via <see cref="ToJson"/>, suitable
/// for committing under <c>simulations/results/f1_n8_n9_metrics/</c>. Both the F1 N=8 and
/// the F1 N=9 chain test classes call <see cref="Compute"/> on the per-block spectrum
/// produced by <see cref="LiouvillianBlockSpectrum.ComputeSpectrumPerBlock"/>, log the
/// metrics through <c>ITestOutputHelper</c>, and persist the JSON.</para>
///
/// <para>Anchors:
/// <see cref="F1GeneralTopologyVerifiedClaim"/>,
/// <c>compute/RCPsiSquared.Core.Tests/F1/F1GeneralTopologyN8BlockSpectrumTests.cs</c>,
/// <c>compute/RCPsiSquared.Core.Tests/F1/F1GeneralTopologyN9BlockSpectrumChainTests.cs</c>,
/// <c>docs/proofs/PROOF_F1_GENERAL_TOPOLOGY.md</c> (verification table populated from
/// the JSON outputs).</para></summary>
public static class F1SpectrumStatistics
{
    /// <summary>Tolerance (in absolute units of Re/Im) that classifies a coordinate as
    /// machine zero. Used for kernel-dimension, pure-imaginary, and real-eigenvalue counts.
    /// 1e-9 is loose enough to absorb MKL Evd accumulation across the 81 (N=8) or 100 (N=9)
    /// joint-popcount block diagonalisations while tight enough to separate the physical
    /// dissipation gap from machine zero in all systems exercised here.</summary>
    public const double MachineZeroTolerance = 1e-9;

    /// <summary>Bin width for the "distinct binned eigenvalues" degeneracy proxy. Each
    /// eigenvalue is rounded to the nearest <c>BinWidth</c> grid point in both real and
    /// imaginary parts; the count of distinct (binnedRe, binnedIm) pairs is the proxy.
    /// 1e-9 matches <see cref="MachineZeroTolerance"/>.</summary>
    public const double BinWidth = 1e-9;

    /// <summary>Outlier definition for pairing-distance histogram: pairs whose distance
    /// exceeds <c>OutlierMultiplier · median</c> are flagged. 100 is permissive enough
    /// that uniform-precision spectra report zero outliers; tight enough to surface the
    /// few worst-case Evd accumulations on the largest sectors.</summary>
    public const double OutlierMultiplier = 100.0;

    /// <summary>Bond record for JSON serialisation. Mirrors <see cref="ChainSystems.Bond"/>
    /// but is a record-type-local payload so the JSON file is fully self-describing
    /// without pulling in core Pauli types.</summary>
    public sealed record BondRecord(int Site1, int Site2, double J);

    /// <summary>Full TopologyMetrics payload covering Groups 1-5. Immutable record;
    /// one instance per (N, topology, H, γ) system. Serialised by
    /// <see cref="ToJson"/>.</summary>
    public sealed record TopologyMetrics(
        // ----- Group 1: wall-time profile (dogfood evidence) -----
        int N,
        string TopologyName,
        double TotalWallSeconds,
        double ComputeSpectrumWallSeconds,
        double EffectiveSpeedupOverDense,
        // ----- Group 2: palindromic-pairing precision -----
        double MaxPairingDistance,
        double MeanPairingDistance,
        double MedianPairingDistance,
        double P99PairingDistance,
        int OutlierPairCount,
        double MinPairingDistance,
        // ----- Group 3: spectrum structure -----
        int SpectrumSize,
        double MinReal,
        double MaxReal,
        double MinImag,
        double MaxImag,
        double DissipationGap,
        int KernelDimension,
        int PureImaginaryCount,
        int RealEigenvalueCount,
        int DistinctBinnedEigenvalueCount,
        // ----- Group 4: block structure -----
        int SectorCount,
        int PrimarySectorCount,
        int MaxBlockSize,
        int MaxBlockSectorPCol,
        int MaxBlockSectorPRow,
        IReadOnlyList<int> Top3BlockSizes,
        long TotalBlockCubicCost,
        // ----- Group 5: Hamiltonian + dissipator setup -----
        double JValue,
        double GammaValue,
        double SigmaShift,
        string HamiltonianClass,
        IReadOnlyList<BondRecord> Bonds);

    /// <summary>Compute every metric group from a per-block spectrum + the wall-time
    /// pair already captured by the test scaffold. The method does NOT itself invoke
    /// the eigensolver; it only post-processes the eigenvalue array and reads block
    /// structure from <see cref="JointPopcountSectorBuilder.Build"/>.</summary>
    /// <param name="N">Qubit count.</param>
    /// <param name="topologyName">Human-readable topology label written verbatim to JSON.</param>
    /// <param name="hamiltonianClass">Hamiltonian class label (e.g. "Heisenberg XXX (XX+YY+ZZ)").</param>
    /// <param name="jValue">Hamiltonian coupling J (J=1 in the standard sweep).</param>
    /// <param name="gammaValue">Per-site Z-dephasing rate γ (γ=0.5 in the standard sweep).</param>
    /// <param name="bonds">Bond list; serialised to JSON for reproducibility.</param>
    /// <param name="spectrum">Per-block Liouvillian spectrum, length 4^N.</param>
    /// <param name="totalWallSeconds">Total test wall time including scaffolding and assertions.</param>
    /// <param name="computeSpectrumWallSeconds">Wall time of the
    /// <see cref="LiouvillianBlockSpectrum.ComputeSpectrumPerBlock"/> call only.</param>
    public static TopologyMetrics Compute(
        int N,
        string topologyName,
        string hamiltonianClass,
        double jValue,
        double gammaValue,
        IReadOnlyList<(int i, int j)> bonds,
        Complex[] spectrum,
        double totalWallSeconds,
        double computeSpectrumWallSeconds)
    {
        if (N < 1) throw new ArgumentOutOfRangeException(nameof(N), N, "N must be >= 1");
        if (spectrum is null) throw new ArgumentNullException(nameof(spectrum));
        int expected = 1 << (2 * N);
        if (spectrum.Length != expected)
            throw new ArgumentException(
                $"spectrum length {spectrum.Length} != 4^N = {expected} at N={N}",
                nameof(spectrum));

        double sigma = N * gammaValue;
        double sigmaShift = -2.0 * sigma;

        // ----- Group 2: pairing precision via greedy nearest-neighbour matching -----
        // Identical algorithm to MultisetAssert.NearestNeighbourEqual; here we keep all
        // distances rather than fail-fast so the histogram can be summarised.
        var shifted = new Complex[expected];
        for (int i = 0; i < expected; i++) shifted[i] = -2.0 * sigma - spectrum[i];

        var distances = NearestNeighbourDistances(spectrum, shifted);
        var sortedDistances = (double[])distances.Clone();
        Array.Sort(sortedDistances);

        double maxDist = sortedDistances[^1];
        double minDist = sortedDistances[0];
        double meanDist = sortedDistances.Average();
        double medianDist = Percentile(sortedDistances, 0.50);
        double p99Dist = Percentile(sortedDistances, 0.99);
        double outlierThreshold = medianDist * OutlierMultiplier;
        int outlierCount = 0;
        for (int i = 0; i < sortedDistances.Length; i++)
            if (sortedDistances[i] > outlierThreshold) outlierCount++;

        // ----- Group 3: spectrum-structure invariants -----
        double minRe = double.PositiveInfinity, maxRe = double.NegativeInfinity;
        double minIm = double.PositiveInfinity, maxIm = double.NegativeInfinity;
        int kernelDim = 0;
        int pureImagCount = 0;
        int realCount = 0;
        double dissipationGap = double.PositiveInfinity;
        for (int i = 0; i < expected; i++)
        {
            var l = spectrum[i];
            double re = l.Real, im = l.Imaginary;
            if (re < minRe) minRe = re;
            if (re > maxRe) maxRe = re;
            if (im < minIm) minIm = im;
            if (im > maxIm) maxIm = im;

            double absRe = Math.Abs(re);
            double absIm = Math.Abs(im);
            bool reIsZero = absRe < MachineZeroTolerance;
            bool imIsZero = absIm < MachineZeroTolerance;

            // Kernel: |λ| < tol (both Re and Im machine-zero).
            if (reIsZero && imIsZero) kernelDim++;
            // Pure-imaginary undamped oscillation: Re ~ 0 AND Im not ~ 0.
            else if (reIsZero) pureImagCount++;
            // Real eigenvalue: Im ~ 0 (per the plan's definition; kernel is the subset
            // with Re also ~ 0).
            if (imIsZero) realCount++;

            // Dissipation gap: smallest non-trivial |Re|. The kernel is allowed (|Re|<tol
            // excluded); we want the slowest non-stationary decay rate.
            if (absRe > MachineZeroTolerance && absRe < dissipationGap)
                dissipationGap = absRe;
        }
        if (double.IsPositiveInfinity(dissipationGap)) dissipationGap = 0.0;

        // Distinct binned eigenvalue count (degeneracy proxy at 1e-9 grid resolution).
        var binned = new HashSet<(long binRe, long binIm)>(capacity: expected / 4);
        double invBin = 1.0 / BinWidth;
        for (int i = 0; i < expected; i++)
        {
            var l = spectrum[i];
            long bRe = (long)Math.Round(l.Real * invBin);
            long bIm = (long)Math.Round(l.Imaginary * invBin);
            binned.Add((bRe, bIm));
        }
        int distinctBinned = binned.Count;

        // ----- Group 4: block structure (LiouvillianBlockSpectrum-specific) -----
        var decomp = JointPopcountSectorBuilder.Build(N);
        int sectorCount = decomp.SectorRanges.Count;

        // Primary-sector count: X⊗N pairing collapses paired sectors onto a single eig
        // call. Count distinct "primaries" (lex-smaller of each pair plus self-paired).
        int primarySectorCount = CountPrimarySectorsByXNPairing(N, decomp.SectorRanges);

        // Max-block size + sector label + Top-3 block sizes by descending size.
        var sectorsByDesc = decomp.SectorRanges
            .OrderByDescending(s => s.Size)
            .ToArray();
        var maxSector = sectorsByDesc[0];
        int maxBlock = maxSector.Size;
        var top3 = sectorsByDesc.Take(3).Select(s => s.Size).ToArray();

        // Total per-block cubic cost (sum of block-size³); cross-N scaling figure.
        long totalCubic = 0L;
        foreach (var s in decomp.SectorRanges)
            totalCubic += (long)s.Size * s.Size * s.Size;

        // Group 1: effective speedup over naive dense (4^N)³.
        // Cast to double early — at N=9, (4^9)³ = 1.8e16 overflows long.
        double denseCubic = Math.Pow(expected, 3.0);
        double effectiveSpeedup = totalCubic > 0 ? denseCubic / totalCubic : 0.0;

        // ----- Group 5: Hamiltonian + dissipator setup payload -----
        var bondRecords = bonds.Select(b => new BondRecord(b.i, b.j, jValue)).ToArray();

        return new TopologyMetrics(
            N: N,
            TopologyName: topologyName,
            TotalWallSeconds: totalWallSeconds,
            ComputeSpectrumWallSeconds: computeSpectrumWallSeconds,
            EffectiveSpeedupOverDense: effectiveSpeedup,
            MaxPairingDistance: maxDist,
            MeanPairingDistance: meanDist,
            MedianPairingDistance: medianDist,
            P99PairingDistance: p99Dist,
            OutlierPairCount: outlierCount,
            MinPairingDistance: minDist,
            SpectrumSize: expected,
            MinReal: minRe,
            MaxReal: maxRe,
            MinImag: minIm,
            MaxImag: maxIm,
            DissipationGap: dissipationGap,
            KernelDimension: kernelDim,
            PureImaginaryCount: pureImagCount,
            RealEigenvalueCount: realCount,
            DistinctBinnedEigenvalueCount: distinctBinned,
            SectorCount: sectorCount,
            PrimarySectorCount: primarySectorCount,
            MaxBlockSize: maxBlock,
            MaxBlockSectorPCol: maxSector.PCol,
            MaxBlockSectorPRow: maxSector.PRow,
            Top3BlockSizes: top3,
            TotalBlockCubicCost: totalCubic,
            JValue: jValue,
            GammaValue: gammaValue,
            SigmaShift: sigmaShift,
            HamiltonianClass: hamiltonianClass,
            Bonds: bondRecords);
    }

    /// <summary>Serialise metrics to pretty-printed JSON. Determinist filename pattern
    /// is chosen by the caller (e.g. <c>&lt;topology&gt;_N&lt;N&gt;.json</c>).</summary>
    public static string ToJson(TopologyMetrics metrics)
    {
        var opts = new JsonSerializerOptions { WriteIndented = true };
        return JsonSerializer.Serialize(metrics, opts);
    }

    /// <summary>Append-or-replace JSON write into the given directory; creates the
    /// directory if missing. Used by both the N=8 and N=9 test classes.</summary>
    public static void WriteJson(TopologyMetrics metrics, string directoryPath, string fileName)
    {
        Directory.CreateDirectory(directoryPath);
        var fullPath = Path.Combine(directoryPath, fileName);
        File.WriteAllText(fullPath, ToJson(metrics), System.Text.Encoding.UTF8);
    }

    /// <summary>Greedy nearest-neighbour matching: for each <c>actual[i]</c>, find the
    /// closest still-unmatched <c>expected[j]</c> by Euclidean (Magnitude) distance and
    /// record that distance. Returns the array of N matched distances. Same algorithm as
    /// <c>MultisetAssert.NearestNeighbourEqual</c>; here we expose the full distance
    /// vector instead of fail-fast on first violation.</summary>
    private static double[] NearestNeighbourDistances(Complex[] actual, Complex[] expected)
    {
        if (actual.Length != expected.Length)
            throw new ArgumentException("actual / expected length mismatch");
        int n = actual.Length;
        var taken = new bool[n];
        var distances = new double[n];
        for (int i = 0; i < n; i++)
        {
            var x = actual[i];
            int bestIdx = -1;
            double bestDist = double.MaxValue;
            for (int j = 0; j < n; j++)
            {
                if (taken[j]) continue;
                double d = (x - expected[j]).Magnitude;
                if (d < bestDist) { bestDist = d; bestIdx = j; }
            }
            if (bestIdx < 0) throw new InvalidOperationException("no candidate");
            taken[bestIdx] = true;
            distances[i] = bestDist;
        }
        return distances;
    }

    /// <summary>Linear-interpolated percentile from a sorted ascending array. Returns
    /// <c>sorted[(int)((n-1)·q)]</c> with linear fractional interpolation to the next
    /// point. q must be in [0, 1].</summary>
    private static double Percentile(double[] sortedAsc, double q)
    {
        if (sortedAsc.Length == 0) return 0.0;
        if (q <= 0.0) return sortedAsc[0];
        if (q >= 1.0) return sortedAsc[^1];
        double pos = q * (sortedAsc.Length - 1);
        int lo = (int)Math.Floor(pos);
        int hi = (int)Math.Ceiling(pos);
        if (lo == hi) return sortedAsc[lo];
        double frac = pos - lo;
        return sortedAsc[lo] * (1.0 - frac) + sortedAsc[hi] * frac;
    }

    /// <summary>Count primary sectors under the X⊗N pairing (p_c, p_r) ↔ (N-p_c, N-p_r).
    /// Self-paired sectors (p_c == N-p_c AND p_r == N-p_r, only possible at even N with
    /// p_c = p_r = N/2) count themselves; non-self-paired sectors count exactly once each
    /// when the lex-smaller representative is encountered first.</summary>
    private static int CountPrimarySectorsByXNPairing(
        int N, IReadOnlyList<JointPopcountSectorBuilder.SectorRange> sectors)
    {
        var primaries = new HashSet<(int pc, int pr)>(capacity: sectors.Count);
        foreach (var s in sectors)
        {
            int pcMirror = N - s.PCol;
            int prMirror = N - s.PRow;
            // Lex-smaller of {(pc, pr), (N-pc, N-pr)} is the primary representative.
            var primary = (s.PCol, s.PRow).CompareTo((pcMirror, prMirror)) <= 0
                ? (s.PCol, s.PRow)
                : (pcMirror, prMirror);
            primaries.Add(primary);
        }
        return primaries.Count;
    }
}
