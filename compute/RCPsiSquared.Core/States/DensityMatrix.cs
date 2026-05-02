using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;
using ComplexVector = MathNet.Numerics.LinearAlgebra.Vector<System.Numerics.Complex>;

namespace RCPsiSquared.Core.States;

/// <summary>Coercion between pure-state vectors and density matrices.
///
/// For a state vector |ψ⟩ of length 2^N, the corresponding density matrix is
/// ρ = |ψ⟩⟨ψ|. For an existing 2^N × 2^N matrix, this is a no-op identity.
/// </summary>
public static class DensityMatrix
{
    /// <summary>Outer product |ψ⟩⟨ψ|.</summary>
    public static ComplexMatrix FromStateVector(ComplexVector psi)
    {
        int d = psi.Count;
        var rho = Matrix<Complex>.Build.Dense(d, d);
        for (int i = 0; i < d; i++)
            for (int j = 0; j < d; j++)
                rho[i, j] = psi[i] * Complex.Conjugate(psi[j]);
        return rho;
    }

    /// <summary>Infer N from a 2^N × 2^N matrix or 2^N vector.</summary>
    public static int InferN(int dimension)
    {
        int n = (int)Math.Round(Math.Log2(dimension));
        if ((1 << n) != dimension)
            throw new ArgumentException($"dimension {dimension} is not a power of 2");
        return n;
    }
}
