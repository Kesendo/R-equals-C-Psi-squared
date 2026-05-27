using System.Numerics;
using RCPsiSquared.Core.Pauli;

namespace RCPsiSquared.Diagnostics.Polarity;

/// <summary>Sparse Pauli-basis representation of L_σ = -i[σ, ·] for a single
/// Pauli string σ at chain length N. Per the structural Pauli algebra:
/// L_σ[τ] = -2i · ε(σ,τ) · (σ·τ) where ε(σ,τ) ∈ {0,1} is the anticommutation
/// indicator. So each column τ of the 4^N × 4^N matrix has at most ONE nonzero
/// entry, at row index = index(σ·τ as Pauli string) with value -2i · (phase
/// of σ·τ). Total nonzero count for non-identity σ: 2·4^(N-1).
///
/// <para>Memory: ~12 KB per L_σ at N=5 (vs 16 MB dense), ~50 KB at N=6 (vs 256 MB),
/// ~200 KB at N=7 (vs 4 GB). Unlocks N=6 enumeration as empirical anchor for the
/// universal-N F112 closure (Welle 11).</para>
///
/// <para>Internal storage: parallel arrays (ColIndices[k], RowIndices[k], Values[k])
/// hold the k-th nonzero. ColIndices is the Pauli-basis index of the input τ;
/// RowIndices is the Pauli-basis index of the output σ·τ; Values is the complex
/// coefficient. The arrays are sorted by ColIndices ascending so that
/// FrobeniusInnerSparse can do a column-index merge in O(min(nnz_A, nnz_B)).</para></summary>
public sealed record SparseLSigma(
    int N,
    int Dim,             // 4^N (precomputed for convenience)
    int[] ColIndices,    // length = nnz, sorted ascending
    int[] RowIndices,    // length = nnz, parallel to ColIndices
    Complex[] Values)    // length = nnz, parallel
{
    public int Nnz => ColIndices.Length;
}
