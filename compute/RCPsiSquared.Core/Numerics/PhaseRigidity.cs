using System.Collections.Generic;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using ComplexVector = MathNet.Numerics.LinearAlgebra.Vector<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Numerics;

/// <summary>Phase rigidity of each eigenmode of a (generally non-normal) complex matrix L. For right
/// eigenvectors |R_i⟩ (L|R_i⟩ = λ_i|R_i⟩) and left eigenvectors ⟨L_i| (⟨L_i|L = λ_i⟨L_i|),
/// r_i = |⟨L_i|R_i⟩| / (‖L_i‖·‖R_i‖) with the Hermitian inner product. r_i = 1 for a normal/isolated
/// mode; r_i → 0 at an exceptional point (the coalescing left/right eigenvectors become orthogonal,
/// the Petermann factor 1/r_i² diverges). Computed in Petermann form r_i = 1/(‖R⁻¹_row_i‖·‖R_col_i‖)
/// from a single Evd(L) (R = eigenvector matrix): near an EP the matching R⁻¹ row diverges so r → 0,
/// while a merely-degenerate non-defective mode keeps a bounded R⁻¹ row so r stays > 0. This is
/// basis-robust at a degeneracy, where eigenvalue-matched left/right overlaps mis-pair and fake r → 0.</summary>
public static class PhaseRigidity
{
    /// <summary>One eigenmode: its eigenvalue, its phase rigidity, and its right eigenvector.</summary>
    public readonly record struct Mode(Complex Lambda, double Rigidity, ComplexVector Right);

    /// <summary>Per-eigenvalue phase rigidity of L, returned in L's Evd eigenvalue order.</summary>
    public static IReadOnlyList<Mode> Compute(Matrix<Complex> L)
    {
        var evd = L.Evd();
        var lam = evd.EigenValues;
        var r = evd.EigenVectors;     // columns = right eigenvectors
        var rInv = r.Inverse();       // rows = left eigenvectors (biorthogonal duals)
        int n = lam.Count;

        var modes = new List<Mode>(n);
        for (int i = 0; i < n; i++)
        {
            var right = r.Column(i);
            // Petermann form: ‖left‖·‖right‖ with left = R⁻¹ row i. Diverges at an EP ⟹ rigidity → 0.
            double rigidity = 1.0 / (rInv.Row(i).L2Norm() * right.L2Norm());
            modes.Add(new Mode(lam[i], rigidity, right));
        }
        return modes;
    }
}
