namespace RCPsiSquared.Core.CoherenceBlocks;

/// <summary>F74: number of distinct pure-dephasing rates in the (n, n+1) coherence block.</summary>
public static class Chromaticity
{
    /// <summary>c(n, N) = min(n, N − 1 − n) + 1. Pure rates are 2γ₀·HD with HD ∈ {1, 3, …, 2c−1}.</summary>
    public static int Compute(int N, int n)
    {
        if (N < 1) throw new ArgumentOutOfRangeException(nameof(N), $"N must be ≥ 1; got {N}.");
        if (n < 0 || n >= N) throw new ArgumentOutOfRangeException(nameof(n),
            $"n must be in [0, N-1]; got n={n}, N={N}.");
        return Math.Min(n, N - 1 - n) + 1;
    }

    /// <summary>HD values realised in the (n, n+1) block: {1, 3, …, 2c−1}.</summary>
    public static IReadOnlyList<int> HammingDistances(int N, int n)
    {
        int c = Compute(N, n);
        return Enumerable.Range(0, c).Select(k => 2 * k + 1).ToArray();
    }
}
