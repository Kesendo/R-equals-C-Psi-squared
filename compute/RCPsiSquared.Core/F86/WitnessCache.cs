using System.Collections.Concurrent;
using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.Decomposition;
using RCPsiSquared.Core.F86.Item1Derivation;
using RCPsiSquared.Core.Resonance;

namespace RCPsiSquared.Core.F86;

/// <summary>Per-(c, N, γ₀) cache of <see cref="KCurve"/> and <see cref="InterChannelSvd"/>
/// objects so multiple witness types (HWHM ratio, Q_peak position, shape-function values,
/// σ_0 scaling) at the same (c, N) share a single computation instead of re-running it.
///
/// <para>Construction is on-demand via <see cref="GetOrCompute"/> and <see cref="GetOrComputeSvd"/>:
/// the first witness asking for (c=3, N=7) triggers the underlying scan/SVD; subsequent
/// witnesses at the same point reuse the cached result. Thread-safe via
/// <see cref="ConcurrentDictionary{TKey, TValue}"/>.</para>
/// </summary>
public sealed class WitnessCache
{
    public IReadOnlyList<double>? QGrid { get; }

    private readonly ConcurrentDictionary<(int C, int N, double GammaZero), KCurve> _kCurves = new();
    private readonly ConcurrentDictionary<(int C, int N, double GammaZero), InterChannelSvd> _svds = new();
    private readonly ConcurrentDictionary<(int N, int LowerPopcount, double GammaZero), C2HwhmRatio> _c2HwhmRatios = new();

    public WitnessCache(IReadOnlyList<double>? qGrid = null)
    {
        QGrid = qGrid;
    }

    /// <summary>Get or compute the K-curve for the (c, N, γ₀) point. The (c, N) → n mapping
    /// uses the smallest n giving chromaticity c (i.e. n = c − 1, valid when N ≥ 2c − 1).</summary>
    public KCurve GetOrCompute(int c, int N, double gammaZero) =>
        _kCurves.GetOrAdd((c, N, gammaZero), key =>
        {
            var block = BuildBlock(key.C, key.N, key.GammaZero);
            var scan = new ResonanceScan(block);
            return QGrid is null ? scan.ComputeKCurve() : scan.ComputeKCurve(QGrid);
        });

    /// <summary>Get or compute the HD=1 ↔ HD=3 inter-channel SVD for (c, N, γ₀). Used by
    /// <see cref="SigmaZeroChromaticityScaling"/> witnesses and any other code that needs
    /// σ_0 = top singular value.</summary>
    public InterChannelSvd GetOrComputeSvd(int c, int N, double gammaZero) =>
        _svds.GetOrAdd((c, N, gammaZero), key =>
            InterChannelSvd.Build(BuildBlock(key.C, key.N, key.GammaZero), hd1: 1, hd2: 3));

    /// <summary>Get or compute <see cref="C2HwhmRatio"/> for a c=2 block, keyed on
    /// (N, n, γ₀). Shared between <see cref="C2UniversalShapeDerivation"/> and
    /// <see cref="PerF71OrbitKTable"/> so an F86KnowledgeBase inspection that pulls both
    /// pays the full-block Q-scan once instead of twice.</summary>
    public C2HwhmRatio GetOrComputeC2HwhmRatio(CoherenceBlock block)
    {
        if (block.C != 2)
            throw new ArgumentException(
                $"GetOrComputeC2HwhmRatio applies only to c=2 blocks; got c={block.C} (N={block.N}, n={block.LowerPopcount}).",
                nameof(block));
        return _c2HwhmRatios.GetOrAdd(
            (block.N, block.LowerPopcount, block.GammaZero),
            _ => C2HwhmRatio.Build(block));
    }

    private static CoherenceBlock BuildBlock(int c, int N, double gammaZero)
    {
        if (N < 2 * c - 1)
            throw new ArgumentException($"N={N} too small for chromaticity c={c} (need N ≥ {2 * c - 1})");
        int n = c - 1;
        return new CoherenceBlock(N, n, gammaZero);
    }

    /// <summary>Process-wide default cache used when witnesses are constructed without an
    /// explicit cache. Lazy: never instantiated until something pulls from it.</summary>
    public static WitnessCache Default { get; } = new();
}
