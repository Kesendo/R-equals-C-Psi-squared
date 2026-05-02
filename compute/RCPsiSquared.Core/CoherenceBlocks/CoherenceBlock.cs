namespace RCPsiSquared.Core.CoherenceBlocks;

/// <summary>The (n, n+1) popcount coherence block of an N-qubit chain under Z-dephasing γ₀.
///
/// <para>The block carries N, n, γ₀, the basis (popcount-state lists + index maps), and a
/// lazy-built <see cref="BlockLDecomposition"/> for L = D + Σ_b J_b · M_H_per_bond[b]. The
/// chromaticity is c(n, N) = min(n, N − 1 − n) + 1 (F74); pure rates are 2γ₀·HD with
/// HD ∈ {1, 3, …, 2c − 1}.</para>
/// </summary>
public sealed class CoherenceBlock
{
    public int N { get; }
    public int LowerPopcount { get; }
    public double GammaZero { get; }
    public BlockBasis Basis { get; }
    public int C => Chromaticity.Compute(N, LowerPopcount);
    public int NumBonds => N - 1;

    public CoherenceBlock(int N, int n, double gammaZero)
    {
        if (gammaZero < 0) throw new ArgumentOutOfRangeException(nameof(gammaZero),
            $"γ₀ must be ≥ 0; got {gammaZero}.");
        this.N = N;
        LowerPopcount = n;
        GammaZero = gammaZero;
        Basis = new BlockBasis(N, n);
    }

    private BlockLDecomposition? _decomp;
    public BlockLDecomposition Decomposition => _decomp ??= BlockLDecomposition.Build(this);
}
