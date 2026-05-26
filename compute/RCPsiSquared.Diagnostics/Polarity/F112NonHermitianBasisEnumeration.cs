using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.Pauli;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Diagnostics.Polarity;

/// <summary>F112 non-Hermitian extension verifier: enumerates the Pauli-string basis
/// at length N and checks the open identity Im⟨L_{σ_α,-i}, L_{σ_β,-i}⟩ = 0 for every
/// pair. If all pairs vanish bit-exact (within tolerance), F112 non-Hermitian extension
/// is Tier1Derived at this N by bilinearity + basis spanning.
///
/// <para>C# port of <c>simulations/_f112_open_identity_basis_enum.py</c>. Existing Python
/// anchor: N=2, 3, 4 bit-exact (35,112 distinct ordered pairs total). This primitive
/// extends to N=5 (524,800 pairs, ~16 GB working memory).</para>
///
/// <para>Per pair cost: one Frobenius inner product over a 4^N × 4^N complex matrix.
/// Per-α pre-cache cost: 4^N L matrices each of size (4^N)² complex = 16·4^(2N) bytes.
/// At N=5: 1024 matrices × 16 MB each = 16 GB total cache.</para></summary>
public static class F112NonHermitianBasisEnumeration
{
    /// <summary>Build L_H = -i[H, ·] in the Pauli basis. Mirrors the Python
    /// <c>build_L_H_pauli</c>: L_vec = -i (H ⊗ I − I ⊗ H^T) is the d² × d² vec-basis
    /// commutator superoperator; rotated into the 4^N-dim Pauli basis via
    /// T^† L_vec T / 2^N where T = <see cref="PauliBasis.VecToPauliBasisTransform"/>.
    /// Returns a 4^N × 4^N complex matrix.</summary>
    public static ComplexMatrix BuildLHInPauliBasis(ComplexMatrix H, int N)
    {
        if (H is null) throw new ArgumentNullException(nameof(H));
        int d = 1 << N;
        if (H.RowCount != d || H.ColumnCount != d)
            throw new ArgumentException($"H must be {d}x{d} (2^N x 2^N); got {H.RowCount}x{H.ColumnCount}", nameof(H));

        var Id = Matrix<Complex>.Build.DenseIdentity(d);
        // L_vec = -i (H ⊗ I − I ⊗ H^T)
        var lVec = -Complex.ImaginaryOne * (H.KroneckerProduct(Id) - Id.KroneckerProduct(H.Transpose()));
        var T = PauliBasis.VecToPauliBasisTransform(N);
        return T.ConjugateTranspose() * lVec * T / (1 << N); // 2^N
    }

    /// <summary>Project a 4^N × 4^N operator-space matrix M onto the eigenspace of the
    /// Π-conjugation map X ↦ Π X Π⁻¹ at the requested eigenvalue λ ∈ {+1, -1, +i, -i}.
    /// Π is order-4 on operator space (Π⁴ = I) so the standard idempotent projector is
    /// P_λ(M) = (1/4) Σ_{k=0..3} λ^{-k} Π^k M Π^{-k}.</summary>
    public static ComplexMatrix ProjectOntoPiEigenspace(ComplexMatrix M, ComplexMatrix pi, Complex targetEigenvalue)
    {
        if (M is null) throw new ArgumentNullException(nameof(M));
        if (pi is null) throw new ArgumentNullException(nameof(pi));

        var piInv = pi.ConjugateTranspose();
        var result = Matrix<Complex>.Build.Dense(M.RowCount, M.ColumnCount);
        var curPi = Matrix<Complex>.Build.DenseIdentity(pi.RowCount);
        var curPiInv = Matrix<Complex>.Build.DenseIdentity(pi.RowCount);
        Complex lambdaToK = Complex.One;
        for (int k = 0; k < 4; k++)
        {
            var coef = (1.0 / lambdaToK) / 4.0;
            result = result + coef * (curPi * M * curPiInv);
            curPi = curPi * pi;
            curPiInv = curPiInv * piInv;
            lambdaToK = lambdaToK * targetEigenvalue;
        }
        return result;
    }

    /// <summary>Frobenius inner product ⟨A, B⟩ = Σ A[i,j]* · B[i,j] = Tr(A† B).</summary>
    public static Complex FrobeniusInner(ComplexMatrix A, ComplexMatrix B)
    {
        if (A is null) throw new ArgumentNullException(nameof(A));
        if (B is null) throw new ArgumentNullException(nameof(B));
        if (A.RowCount != B.RowCount || A.ColumnCount != B.ColumnCount)
            throw new ArgumentException($"shape mismatch: A is {A.RowCount}x{A.ColumnCount}, B is {B.RowCount}x{B.ColumnCount}");

        Complex sum = Complex.Zero;
        for (int i = 0; i < A.RowCount; i++)
            for (int j = 0; j < A.ColumnCount; j++)
                sum += Complex.Conjugate(A[i, j]) * B[i, j];
        return sum;
    }
}
