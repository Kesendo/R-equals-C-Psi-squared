using System.Collections.Generic;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using ComplexVector = MathNet.Numerics.LinearAlgebra.Vector<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Numerics;

/// <summary>Phase rigidity of each eigenmode of a (generally non-normal) complex matrix L. For right
/// eigenvectors |R_i⟩ (L|R_i⟩ = λ_i|R_i⟩) and left eigenvectors ⟨L_i| (⟨L_i|L = λ_i⟨L_i|),
/// r_i = |⟨L_i|R_i⟩| / (‖L_i‖·‖R_i‖) with the Hermitian inner product. r_i = 1 for a mode of a normal matrix.
/// The rigidity tends to zero at an exceptional point (the coalescing left/right eigenvectors become orthogonal,
/// the Petermann factor 1/r_i² diverges). Computed in Petermann form r_i = 1/(‖R⁻¹_row_i‖·‖R_col_i‖)
/// from a single Evd(L) (R = eigenvector matrix): near an EP the matching R⁻¹ row diverges so r → 0.
/// For a simple isolated eigenvalue, 0 < r_i ≤ 1 is the usual per-mode Petermann rigidity.
/// This construction avoids the separate-eigensolver left/right matching failure that can fake r → 0.
/// At an exact repeated semisimple eigenvalue, however, the individual values depend on the eigenbasis chosen inside its invariant subspace; use a spectral-projector or Jordan diagnostic there.</summary>
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
        var rigidities = RigiditiesFromRightEigenvectors(r);
        int n = lam.Count;

        var modes = new List<Mode>(n);
        for (int i = 0; i < n; i++)
        {
            var right = r.Column(i);
            modes.Add(new Mode(lam[i], rigidities[i], right));
        }
        return modes;
    }

    /// <summary>Petermann rigidities for one chosen right-eigenvector matrix. At an exact repeated
    /// eigenvalue these are properties of that chosen eigenbasis, not invariants of each mode.</summary>
    internal static double[] RigiditiesFromRightEigenvectors(Matrix<Complex> rightEigenvectors)
    {
        var inverse = rightEigenvectors.Inverse(); // rows = biorthogonal left duals
        int n = rightEigenvectors.ColumnCount;
        var rigidities = new double[n];
        for (int i = 0; i < n; i++)
        {
            rigidities[i] = 1.0 /
                (inverse.Row(i).L2Norm() * rightEigenvectors.Column(i).L2Norm());
        }
        return rigidities;
    }
}
