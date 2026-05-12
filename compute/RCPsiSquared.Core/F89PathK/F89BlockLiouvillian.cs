using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.Pauli;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.F89PathK;

/// <summary>F89 block Liouvillian super-operator on the column-major vec convention
/// <c>vec(M)[b·d + a] = M[a, b]</c>.
///
/// <para><c>L_block = -i[H_B, ·] + Σ_l γ (Z_l ρ Z_l - ρ)</c> on n_block-qubit Hilbert space,
/// with <c>H_B = J · Σ (X_b X_{b+1} + Y_b Y_{b+1})</c>. Returned as a (d² × d²) dense matrix
/// where d = 2^n_block.</para>
///
/// <para>Column-major vec form: <c>L = -i (I_d ⊗ H - H^T ⊗ I_d) + Σ_l γ (Z_l^T ⊗ Z_l - I ⊗ I)</c>.
/// This matches Python <c>simulations/_f89_pathk_lib.py</c> bit-for-bit. To compare against
/// <see cref="Lindblad.LindbladianBuilder"/> (row-major vec), build dρ/dt directly and check
/// both forms reproduce it; see <see cref="F89PathKConvention"/>.</para>
/// </summary>
public static class F89BlockLiouvillian
{
    public static ComplexMatrix BuildBlockL(double J, double gamma, int nBlock)
    {
        if (nBlock < 2)
            throw new ArgumentOutOfRangeException(nameof(nBlock), nBlock, "nBlock must be ≥ 2.");
        if (gamma < 0)
            throw new ArgumentOutOfRangeException(nameof(gamma), gamma, "gamma must be ≥ 0.");

        int d = 1 << nBlock;
        var H = F89BlockHamiltonian.BuildBlockH(J, nBlock);
        var Id = Matrix<Complex>.Build.DenseIdentity(d);

        var idKronH = Id.KroneckerProduct(H);
        var hTKronId = H.Transpose().KroneckerProduct(Id);
        var L = -Complex.ImaginaryOne * (idKronH - hTKronId);

        var idKronId = Id.KroneckerProduct(Id);
        for (int l = 0; l < nBlock; l++)
        {
            var Zl = PauliString.SiteOp(nBlock, l, PauliLetter.Z);
            var zlTKronZl = Zl.Transpose().KroneckerProduct(Zl);
            L = L + (Complex)gamma * (zlTKronZl - idKronId);
        }
        return L;
    }
}
