using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.Pauli;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.F89PathK;

/// <summary>F89 XY block Hamiltonian on n_block sites: <c>H_B = J · Σ_b (X_b X_{b+1} + Y_b Y_{b+1})</c>.
///
/// <para>F89 convention has no 1/2 prefactor (contrast F86 <see cref="Pauli.PauliHamiltonian.XYChain"/>).
/// See <see cref="F89PathKConvention"/> for the J/2 vs J distinction.</para>
/// </summary>
public static class F89BlockHamiltonian
{
    public static ComplexMatrix BuildBlockH(double J, int nBlock)
    {
        if (nBlock < 2)
            throw new ArgumentOutOfRangeException(nameof(nBlock), nBlock, "nBlock must be ≥ 2.");
        int d = 1 << nBlock;
        var H = Matrix<Complex>.Build.Dense(d, d);
        Complex jc = J;
        for (int b = 0; b < nBlock - 1; b++)
        {
            var xb = PauliString.SiteOp(nBlock, b, PauliLetter.X);
            var xbp = PauliString.SiteOp(nBlock, b + 1, PauliLetter.X);
            var yb = PauliString.SiteOp(nBlock, b, PauliLetter.Y);
            var ybp = PauliString.SiteOp(nBlock, b + 1, PauliLetter.Y);
            H = H + jc * (xb * xbp + yb * ybp);
        }
        return H;
    }
}
