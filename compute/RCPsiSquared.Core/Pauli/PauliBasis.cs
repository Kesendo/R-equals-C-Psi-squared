using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;
using ComplexVector = MathNet.Numerics.LinearAlgebra.Vector<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Pauli;

/// <summary>4^N Pauli-basis decomposition of density matrices.
///
/// The Pauli strings σ_α (α ∈ [0, 4^N)) form an orthogonal basis of d×d matrices
/// (d = 2^N) under the Hilbert-Schmidt inner product: Tr(σ_α σ_β) = d · δ_αβ.
/// Any d×d operator A decomposes as
///   A = Σ_α vec[α] · σ_α,    vec[α] = Tr(σ_α · A) / d.
///
/// For a 1-qubit density matrix this reduces to the standard Bloch-vector form
/// ρ = (I + ⟨X⟩ X + ⟨Y⟩ Y + ⟨Z⟩ Z) / 2 with vec[I] = 1/2, vec[X] = ⟨X⟩/2, etc.
///
/// <para>Also exposes <see cref="VecToPauliBasisTransform"/>, the d² × 4^N matrix M with
/// columns vec_F(σ_α) (column-major flattening). Cached by N. M† M = 2^N · I. Used by
/// Symmetry.PalindromeResidual, Observables.PiProtectedObservables, and the F81 / DZero
/// Π-decomposition diagnostics — the canonical Pauli↔vec round-trip in matrix form.</para>
/// </summary>
public static class PauliBasis
{
    private static readonly System.Collections.Concurrent.ConcurrentDictionary<int, ComplexMatrix> _transformCache = new();

    /// <summary>The d² × 4^N matrix M with M[:, α] = vec_F(σ_α) column-stacked. Cached by N.
    /// Callers should treat the returned matrix as read-only.
    ///
    /// <para><b>Convention note (Welle 10d audit, 2026-05-27)</b>: T columns are vec_F
    /// (column-major) per the BuildVecToPauliBasisTransform loop. The standard codebase
    /// pattern combines T with a vec_R-style commutator superoperator L_vec = -i(H ⊗ I −
    /// I ⊗ H^T) (e.g. Liouvillian.cs, LindbladianBuilder.cs, BondPerturbation.cs,
    /// PauliDephasingDissipator.cs, PalindromeResidual.cs, PiDecomposition.cs,
    /// F112NonHermitianBasisEnumeration.cs). Reading T columns as vec_R implicitly maps
    /// σ_k ↦ σ_k^T, so the resulting Pauli-basis matrix is D · L_natural · D where
    /// D = diag((−1)^n_Y(k)) is a real diagonal unitary involution. All current callers
    /// compute D-conjugation-invariant quantities (Frobenius norms, inner products,
    /// eigenvalues, zero patterns), so the convention twist is invisible — but a new
    /// caller that reads individual matrix-entry SIGNS or compares against a
    /// third-party natural L would need to apply the (−1)^(n_Y(row) + n_Y(col)) sign
    /// correction (see <c>F112NonHermitianBasisEnumeration.BuildSparseLSigma</c> for the
    /// applied form). A bonus structural identity surfaced: D · Π_Z · D = Π_Y, so the
    /// twist coincides with the Z↔Y dephasing-letter swap on the operator space.</para></summary>
    public static ComplexMatrix VecToPauliBasisTransform(int N) =>
        _transformCache.GetOrAdd(N, BuildVecToPauliBasisTransform);

    private static ComplexMatrix BuildVecToPauliBasisTransform(int N)
    {
        int d = 1 << N;
        long d2 = 1L << (2 * N);
        var M = Matrix<Complex>.Build.Dense(d * d, (int)d2);
        for (long k = 0; k < d2; k++)
        {
            var letters = PauliIndex.FromFlat(k, N);
            var sigma = PauliString.Build(letters);
            for (int j = 0; j < d; j++)
                for (int i = 0; i < d; i++)
                    M[j * d + i, (int)k] = sigma[i, j];
        }
        return M;
    }

    /// <summary>Decompose ρ into Pauli-basis coefficients: vec[α] = Tr(σ_α · ρ) / 2^N.</summary>
    public static ComplexVector ToPauliVector(ComplexMatrix rho, int N)
    {
        long d2 = 1L << (2 * N);
        double inv = 1.0 / (1 << N);
        var vec = MathNet.Numerics.LinearAlgebra.Vector<Complex>.Build.Dense((int)d2);
        for (long k = 0; k < d2; k++)
        {
            var letters = PauliIndex.FromFlat(k, N);
            var sigma = PauliString.Build(letters);
            vec[(int)k] = (sigma * rho).Trace() * inv;
        }
        return vec;
    }

    /// <summary>Reconstruct A = Σ_α vec[α] · σ_α from its Pauli-basis coefficients.</summary>
    public static ComplexMatrix FromPauliVector(ComplexVector pauliVec, int N)
    {
        long d2 = 1L << (2 * N);
        if (pauliVec.Count != d2)
            throw new ArgumentException($"expected length 4^N = {d2}; got {pauliVec.Count}");
        int d = 1 << N;
        var rho = Matrix<Complex>.Build.Dense(d, d);
        for (long k = 0; k < d2; k++)
        {
            if (pauliVec[(int)k] == Complex.Zero) continue;
            var letters = PauliIndex.FromFlat(k, N);
            rho = rho + pauliVec[(int)k] * PauliString.Build(letters);
        }
        return rho;
    }
}
