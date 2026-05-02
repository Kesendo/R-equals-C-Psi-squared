using System.Numerics;

namespace RCPsiSquared.Core.CoherenceBlocks;

/// <summary>Computational basis for the (n, n+1) coherence block of an N-qubit chain.
///
/// Big-endian: site 0 = MSB. The flat index of a basis pair |p⟩⟨q| with popcount(p)=n,
/// popcount(q)=n+1 is <c>IndexP(p) · Mq + IndexQ(q)</c>.
/// </summary>
public sealed class BlockBasis
{
    public int N { get; }
    public int LowerPopcount { get; }
    public IReadOnlyList<int> StatesP { get; }
    public IReadOnlyList<int> StatesQ { get; }

    private readonly Dictionary<int, int> _pToIdx;
    private readonly Dictionary<int, int> _qToIdx;

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
        _pToIdx = new Dictionary<int, int>(StatesP.Count);
        for (int i = 0; i < StatesP.Count; i++) _pToIdx[StatesP[i]] = i;
        _qToIdx = new Dictionary<int, int>(StatesQ.Count);
        for (int i = 0; i < StatesQ.Count; i++) _qToIdx[StatesQ[i]] = i;
    }

    public int IndexP(int p) => _pToIdx[p];
    public int IndexQ(int q) => _qToIdx[q];
    public int FlatIndex(int p, int q) => _pToIdx[p] * Mq + _qToIdx[q];

    /// <summary>Computational-basis states with given popcount, sorted by integer value.</summary>
    public static IReadOnlyList<int> PopcountStates(int N, int popcount)
    {
        var result = new List<int>();
        for (int x = 0; x < (1 << N); x++)
            if (BitOperations.PopCount((uint)x) == popcount)
                result.Add(x);
        return result;
    }
}
