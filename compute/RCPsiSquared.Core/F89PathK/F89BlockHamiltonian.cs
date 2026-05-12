using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.Pauli;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.F89PathK;

/// <summary>Builds the F89-convention XY chain Hamiltonian
/// <c>H = J·Σ_b (X_b X_{b+1} + Y_b Y_{b+1})</c> on n_block sites.
/// See <see cref="F89PathKConvention"/> for the J convention.</summary>
public static class F89BlockHamiltonian
{
    public static ComplexMatrix BuildBlockH(double J, int nBlock)
    {
        if (nBlock < 2)
            throw new ArgumentOutOfRangeException(nameof(nBlock), nBlock, "nBlock must be ≥ 2.");
        return PauliHamiltonian.XYChain(nBlock, 2.0 * J).ToMatrix();
    }
}
