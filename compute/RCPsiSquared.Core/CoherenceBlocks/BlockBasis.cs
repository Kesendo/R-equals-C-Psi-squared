using System.Numerics;

namespace RCPsiSquared.Core.CoherenceBlocks;

/// <summary>Computational basis for the (n, n+1) coherence block of an N-qubit chain.
///
/// Big-endian: site 0 = MSB. The flat index of a basis pair |p⟩⟨q| with popcount(p)=n,
/// popcount(q)=n+1 is <c>IndexP(p) · Mq + IndexQ(q)</c>.
///
/// <para>States are <c>long</c>: at N ≥ 32 a popcount-1 state can set bit 31 or above,
/// beyond <c>int</c> range. Indices (IndexP/IndexQ/FlatIndex results, Mp/Mq/MTotal) stay
/// <c>int</c> — there are only Mp·Mq of them.</para>
/// </summary>
public sealed class BlockBasis
{
    public int N { get; }
    public int LowerPopcount { get; }
    public IReadOnlyList<long> StatesP { get; }
    public IReadOnlyList<long> StatesQ { get; }

    private readonly Dictionary<long, int> _pToIdx;
    private readonly Dictionary<long, int> _qToIdx;

    public int Mp => StatesP.Count;
    public int Mq => StatesQ.Count;
    public int MTotal => Mp * Mq;

    public BlockBasis(int N, int n)
    {
        if (N < 1) throw new ArgumentOutOfRangeException(nameof(N), $"N must be ≥ 1; got {N}.");
        if (n < 0 || n >= N) throw new ArgumentOutOfRangeException(nameof(n),
            $"n must be in [0, N-1]; got n={n}, N={N}.");
        this.N = N;
        LowerPopcount = n;
        StatesP = PopcountStates(N, n);
        StatesQ = PopcountStates(N, n + 1);
        _pToIdx = new Dictionary<long, int>(StatesP.Count);
        for (int i = 0; i < StatesP.Count; i++) _pToIdx[StatesP[i]] = i;
        _qToIdx = new Dictionary<long, int>(StatesQ.Count);
        for (int i = 0; i < StatesQ.Count; i++) _qToIdx[StatesQ[i]] = i;
    }

    public int IndexP(long p) => _pToIdx[p];
    public int IndexQ(long q) => _qToIdx[q];
    public int FlatIndex(long p, long q) => _pToIdx[p] * Mq + _qToIdx[q];

    /// <summary>Computational-basis states with given popcount, sorted by integer value.
    /// The <c>1L &lt;&lt; N</c> upper bound keeps the scan correct through N = 62.</summary>
    public static IReadOnlyList<long> PopcountStates(int N, int popcount)
    {
        var result = new List<long>();
        long upper = 1L << N;
        for (long x = 0; x < upper; x++)
            if (BitOperations.PopCount((ulong)x) == popcount)
                result.Add(x);
        return result;
    }
}
