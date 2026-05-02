using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.Pauli;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>K-sublattice (chiral) operator: K_full = ⊗_{odd i} Z_i.
///
/// K anti-commutes with H_xy = (X_a X_{a+1} + Y_a Y_{a+1}) / 2: K H K = −H. Z-dephasing
/// commutes with K trivially. When K anti-commutes with H, K is a super-operator symmetry
/// of L. K_full has K² = I (involutory).
///
/// In the Altland-Zirnbauer scheme this is the chiral / sublattice symmetry of class BDI
/// — see project_chiral_partnership memory and PT_SYMMETRY_ANALYSIS for the broader
/// classification context.
/// </summary>
public static class ChiralK
{
    public enum Classification
    {
        /// <summary>K H K = +H (K-symmetric).</summary>
        KEven,
        /// <summary>K H K = −H (K-anti-symmetric: e.g. XY/Heisenberg hopping).</summary>
        KOdd,
        /// <summary>Neither: K-mixed.</summary>
        KMixed,
    }

    public static ComplexMatrix BuildFull(int N)
    {
        var I2 = Matrix<Complex>.Build.DenseIdentity(2);
        var z = PauliMatrix.Of(PauliLetter.Z);
        ComplexMatrix? result = null;
        for (int i = 0; i < N; i++)
        {
            var factor = (i % 2 == 1) ? z : I2;
            result = result is null ? factor : result.KroneckerProduct(factor);
        }
        return result!;
    }

    public static Classification ClassifyHamiltonian(ComplexMatrix H, int N, ComplexMatrix? kFull = null)
    {
        kFull ??= BuildFull(N);
        var khk = kFull * H * kFull;
        if ((khk - H).FrobeniusNorm() < 1e-10) return Classification.KEven;
        if ((khk + H).FrobeniusNorm() < 1e-10) return Classification.KOdd;
        return Classification.KMixed;
    }

    /// <summary>Per-Pauli-string K-eigenvalue (+1, −1, or 0 for mixed).</summary>
    public static int ClassifyPauliString(IReadOnlyList<PauliLetter> letters, int N, ComplexMatrix? kFull = null)
    {
        kFull ??= BuildFull(N);
        var sigma = PauliString.Build(letters);
        var ksk = kFull * sigma * kFull;
        if ((ksk - sigma).FrobeniusNorm() < 1e-10) return +1;
        if ((ksk + sigma).FrobeniusNorm() < 1e-10) return -1;
        return 0;
    }
}
