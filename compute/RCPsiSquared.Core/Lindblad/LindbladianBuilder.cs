using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Lindblad;

/// <summary>Build the d²×d² Liouvillian for a system with Hamiltonian H and Lindblad
/// dissipator operators c_k.
///
/// L = −i [H, ·] + Σ_k ( c_k(·)c_k† − ½ {c_k†c_k, ·} )
///
/// The vec form follows the framework convention
///   L_matrix = −i (H ⊗ I − I ⊗ H^T)  +  Σ_k ( c_k ⊗ c_k* − ½ (c_k†c_k ⊗ I) − ½ (I ⊗ (c_k†c_k)^T) )
/// matching <c>simulations/framework/lindblad.py</c>'s <c>lindbladian_general</c>. The same
/// vec convention is used by <c>palindrome_residual</c> and the Pauli-basis transform.
/// </summary>
public static class LindbladianBuilder
{
    public static ComplexMatrix Build(ComplexMatrix H, IReadOnlyList<ComplexMatrix> cOps)
    {
        AssertHermitian(H);
        int d = H.RowCount;
        var I = Matrix<Complex>.Build.DenseIdentity(d);

        var hKronI = H.KroneckerProduct(I);
        var iKronHt = I.KroneckerProduct(H.Transpose());
        var L = -Complex.ImaginaryOne * (hKronI - iKronHt);

        foreach (var c in cOps)
        {
            var cConj = c.Conjugate();
            var cDag = c.ConjugateTranspose();
            var cDagC = cDag * c;
            L = L
                + c.KroneckerProduct(cConj)
                - 0.5 * cDagC.KroneckerProduct(I)
                - 0.5 * I.KroneckerProduct(cDagC.Transpose());
        }
        return L;
    }

    private static void AssertHermitian(ComplexMatrix H)
    {
        var diff = H - H.ConjugateTranspose();
        if (diff.FrobeniusNorm() > 1e-10)
            throw new ArgumentException("Hamiltonian H must be Hermitian.");
    }
}
