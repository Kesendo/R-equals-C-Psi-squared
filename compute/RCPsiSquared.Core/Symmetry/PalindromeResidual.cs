using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.Pauli;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>The F1-palindrome residual M = Π · L_Pauli · Π⁻¹ + L_Pauli + 2σ·I.
///
/// For "truly" Hamiltonians and pure dephasing, the residual is exactly zero
/// (machine precision): Π·L·Π⁻¹ + L = −2σ·I, the F1 palindrome identity. For
/// non-truly Hamiltonians or T1 dissipators, ‖M‖²_F has the closed-form scaling
/// given by <see cref="Lindblad.PalindromeResidualScaling"/>.
///
/// L is supplied in the framework's vec form (the same as
/// <see cref="Lindblad.LindbladianBuilder"/>'s output). Internally the residual is
/// computed in the 4^N Pauli-string basis via L_pauli = M_vec_to_pauli^† · L · M_vec_to_pauli / 2^N.
/// The 1/2^N normalisation is correct because <c>M_vec_to_pauli</c> satisfies
/// M† · M = 2^N · I, so M⁻¹ = M† / 2^N and the similarity transform is L · M⁻¹.
/// </summary>
public static class PalindromeResidual
{
    public static ComplexMatrix Build(ComplexMatrix LVec, int N, double sigmaGamma,
        PauliLetter dephaseLetter = PauliLetter.Z)
    {
        var transform = PauliBasis.VecToPauliBasisTransform(N);
        var lPauli = (transform.ConjugateTranspose() * LVec * transform) / Math.Pow(2, N);
        var pi = PiOperator.BuildFull(N, dephaseLetter);
        var piInv = pi.ConjugateTranspose(); // Π is unitary signed permutation
        long d2 = 1L << (2 * N);
        var identity = Matrix<Complex>.Build.SparseIdentity((int)d2);
        return pi * lPauli * piInv + lPauli + (Complex)(2.0 * sigmaGamma) * identity;
    }
}
