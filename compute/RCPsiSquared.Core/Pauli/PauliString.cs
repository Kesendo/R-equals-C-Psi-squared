using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Pauli;

/// <summary>Multi-qubit Pauli operators: tensor products of single-letter Pauli matrices.
///
/// Convention: site 0 = leftmost Kronecker factor (most-significant bit in the
/// 2^N basis index). This matches the rest of <c>Core</c> (BlockBasis, Lindblad).
/// </summary>
public static class PauliString
{
    /// <summary>σ_{letters[0]} ⊗ σ_{letters[1]} ⊗ … ⊗ σ_{letters[N−1]} as a 2^N × 2^N matrix.</summary>
    public static ComplexMatrix Build(IReadOnlyList<PauliLetter> letters)
    {
        if (letters.Count == 0)
            throw new ArgumentException("at least one letter required");
        var result = PauliMatrix.Of(letters[0]);
        for (int i = 1; i < letters.Count; i++)
            result = result.KroneckerProduct(PauliMatrix.Of(letters[i]));
        return result;
    }

    /// <summary>Single-Pauli operator on <paramref name="site"/>, identity elsewhere.
    /// 2^N × 2^N matrix.</summary>
    public static ComplexMatrix SiteOp(int N, int site, PauliLetter letter)
    {
        if (site < 0 || site >= N)
            throw new ArgumentOutOfRangeException(nameof(site), $"site must be in [0, N-1]; got {site}, N={N}");
        var letters = new PauliLetter[N];
        for (int i = 0; i < N; i++) letters[i] = PauliLetter.I;
        letters[site] = letter;
        return Build(letters);
    }

    /// <summary>Per-site (X_i, Y_i, Z_i) operators. Useful for Bloch-component readouts.
    /// Cached by N — matrices are immutable structurally.</summary>
    public static IReadOnlyList<(ComplexMatrix X, ComplexMatrix Y, ComplexMatrix Z)> SitePaulis(int N) =>
        _sitePauliCache.GetOrAdd(N, BuildSitePaulis);

    private static readonly System.Collections.Concurrent.ConcurrentDictionary<int,
        IReadOnlyList<(ComplexMatrix X, ComplexMatrix Y, ComplexMatrix Z)>> _sitePauliCache = new();

    private static IReadOnlyList<(ComplexMatrix X, ComplexMatrix Y, ComplexMatrix Z)> BuildSitePaulis(int N)
    {
        var result = new (ComplexMatrix, ComplexMatrix, ComplexMatrix)[N];
        for (int i = 0; i < N; i++)
            result[i] = (SiteOp(N, i, PauliLetter.X), SiteOp(N, i, PauliLetter.Y), SiteOp(N, i, PauliLetter.Z));
        return result;
    }
}
