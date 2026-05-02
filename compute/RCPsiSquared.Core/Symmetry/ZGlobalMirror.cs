using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.Pauli;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;
using ComplexVector = MathNet.Numerics.LinearAlgebra.Vector<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>Global Z-string operator Z⊗N = ⊗_l Z_l.
///
/// Z⊗N satisfies Z⊗N · σ_α · Z⊗N = (−1)^n_XY · σ_α for any Pauli string σ_α (n_XY counts
/// X+Y letters = bit_a-parity in our convention). For an X-basis Néel state |+−+−…⟩,
/// Z⊗N gives |−+−+…⟩ — the AFM-mirror partner.
///
/// Operators with even n_XY per term (XX, YY, ZZ, II, σ⁻σ⁺ pairs, Z-detuning δ_l Z_l) commute
/// with Z⊗N. Single-site X or Y (transverse fields) anti-commute and break Z⊗N.
/// </summary>
public static class ZGlobalMirror
{
    public static ComplexMatrix Build(int N)
    {
        var z = PauliMatrix.Of(PauliLetter.Z);
        ComplexMatrix result = z;
        for (int i = 1; i < N; i++) result = result.KroneckerProduct(z);
        return result;
    }

    public static ComplexVector Apply(ComplexVector psi, int N) => Build(N) * psi;

    public static ComplexMatrix Conjugate(ComplexMatrix rho, int N)
    {
        var z = Build(N);
        return z * rho * z;
    }
}
