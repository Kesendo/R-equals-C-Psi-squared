using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.Pauli;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.F89PathK;

/// <summary>Builds the F89 block Liouvillian super-operator
/// <c>L = -i (I_d ⊗ H - H^T ⊗ I_d) + Σ_l γ (Z_l^T ⊗ Z_l - I ⊗ I)</c> on the column-major
/// vec convention. See <see cref="F89PathKConvention"/>.</summary>
public static class F89BlockLiouvillian
{
    public static ComplexMatrix BuildBlockL(double J, double gamma, int nBlock)
    {
        if (nBlock < 2)
            throw new ArgumentOutOfRangeException(nameof(nBlock), nBlock, "nBlock must be ≥ 2.");
        if (nBlock > 6)
            throw new ArgumentOutOfRangeException(nameof(nBlock), nBlock,
                "BuildBlockL exceeds .NET 2 GB single-array limit at n_block ≥ 7; " +
                "use a sparse / matrix-free path for larger blocks (see RCPsiSquared.Compute Liouvillian.BuildDirectNative or RCPsiSquared.Propagate MatrixFreePropagator).");
        if (gamma < 0)
            throw new ArgumentOutOfRangeException(nameof(gamma), gamma, "gamma must be ≥ 0.");

        int d = 1 << nBlock;
        var H = F89BlockHamiltonian.BuildBlockH(J, nBlock);
        var Id = Matrix<Complex>.Build.DenseIdentity(d);

        var idKronH = Id.KroneckerProduct(H);
        var hTKronId = H.Transpose().KroneckerProduct(Id);
        var L = -Complex.ImaginaryOne * (idKronH - hTKronId);

        var siteOps = PauliString.SitePaulis(nBlock);
        var zSum = Matrix<Complex>.Build.Dense(d * d, d * d);
        for (int l = 0; l < nBlock; l++)
        {
            var Zl = siteOps[l].Z;
            zSum = zSum + Zl.Transpose().KroneckerProduct(Zl);
        }
        L = L + (Complex)gamma * zSum;

        // Subtract n_block · γ from the diagonal in place (avoids allocating Id ⊗ Id).
        double diagShift = nBlock * gamma;
        int dd = d * d;
        for (int i = 0; i < dd; i++) L[i, i] -= diagShift;

        return L;
    }
}
