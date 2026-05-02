using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.Pauli;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Lindblad;

/// <summary>L = −i[H, ·] + Σ_l γ_l (P_l ρ P_l − ρ) for P ∈ {X, Y, Z} dephasing.
///
/// Z-dephasing is the F-framework canonical case: satisfies the **F1** palindrome identity
/// Π·L·Π⁻¹ + L + 2Σγ·I = 0 bit-exactly when H is "truly" (XY/Heisenberg). See
/// docs/ANALYTICAL_FORMULAS.md F1 and <see cref="Symmetry.PalindromeResidual"/>.
///
/// Dissipator-resonance law (**F77** + dephasing-axis dependence, verified at N=4 k=3 over
/// 294 Z₂³-homogeneous pairs, 2026-05-01): F77-hardness lives in the Klein cell that matches
/// the dephase letter's Klein index — Z=(0,1), X=(1,0), Y=(1,1) in the bit_a/bit_b convention
/// of <see cref="Pauli.PauliLetter"/>.
/// </summary>
public static class PauliDephasingDissipator
{
    public static ComplexMatrix Build(ComplexMatrix H, IReadOnlyList<double> gammaPerSite, PauliLetter dephaseLetter = PauliLetter.Z)
    {
        if (dephaseLetter == PauliLetter.I)
            throw new ArgumentException("dephase_letter must be one of X, Y, Z.");
        int d = H.RowCount;
        int N = (int)Math.Round(Math.Log2(d));
        if ((1 << N) != d) throw new ArgumentException($"H dim {d} not a power of 2");
        if (gammaPerSite.Count != N)
            throw new ArgumentException($"gamma list has length {gammaPerSite.Count}, expected N={N}");

        var I = Matrix<Complex>.Build.DenseIdentity(d);
        var L = -Complex.ImaginaryOne * (H.KroneckerProduct(I) - I.KroneckerProduct(H.Transpose()));
        for (int l = 0; l < N; l++)
        {
            double gamma = gammaPerSite[l];
            if (gamma == 0) continue;
            var Pl = PauliString.SiteOp(N, l, dephaseLetter);
            L = L + (Complex)gamma * (Pl.KroneckerProduct(Pl.Conjugate()) - I.KroneckerProduct(I));
        }
        return L;
    }

    /// <summary>Pure Z-dephasing convenience: L = −i[H, ·] + Σ_l γ_l (Z_l ρ Z_l − ρ).</summary>
    public static ComplexMatrix BuildZ(ComplexMatrix H, IReadOnlyList<double> gammaPerSite) =>
        Build(H, gammaPerSite, PauliLetter.Z);
}
