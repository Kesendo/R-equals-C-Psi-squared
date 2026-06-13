using System.Collections.Generic;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using ComplexVector = MathNet.Numerics.LinearAlgebra.Vector<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Numerics;

/// <summary>Phase rigidity of each eigenmode of a (generally non-normal) complex matrix L. For right
/// eigenvectors |R_i⟩ (L|R_i⟩ = λ_i|R_i⟩) and left eigenvectors ⟨L_i| (⟨L_i|L = λ_i⟨L_i|),
/// r_i = |⟨L_i|R_i⟩| / (‖L_i‖·‖R_i‖) with the Hermitian inner product. r_i = 1 for a normal/isolated
/// mode; r_i → 0 at an exceptional point (the coalescing left/right eigenvectors become orthogonal,
/// the Petermann factor 1/r_i² diverges). Computed from Evd(L) and Evd(Lᴴ) matched by eigenvalue,
/// which is numerically stabler near an EP than inverting the (there-singular) eigenvector matrix.</summary>
public static class PhaseRigidity
{
    /// <summary>One eigenmode: its eigenvalue, its phase rigidity, and its right eigenvector.</summary>
    public readonly record struct Mode(Complex Lambda, double Rigidity, ComplexVector Right);

    /// <summary>Per-eigenvalue phase rigidity of L, returned in L's Evd eigenvalue order.</summary>
    public static IReadOnlyList<Mode> Compute(Matrix<Complex> L)
    {
        var evdR = L.Evd();
        var evdL = L.ConjugateTranspose().Evd();
        var lamR = evdR.EigenValues;
        var lamL = evdL.EigenValues;
        var vR = evdR.EigenVectors;
        var vL = evdL.EigenVectors;
        int n = lamR.Count;

        var used = new bool[n];
        var modes = new List<Mode>(n);
        for (int i = 0; i < n; i++)
        {
            // L's left eigenvector for λ_i is Lᴴ's right eigenvector for conj(λ_i): match by that.
            Complex target = Complex.Conjugate(lamR[i]);
            int best = -1;
            double bestDist = double.PositiveInfinity;
            for (int j = 0; j < n; j++)
            {
                if (used[j]) continue;
                double dist = (lamL[j] - target).Magnitude;
                if (dist < bestDist) { bestDist = dist; best = j; }
            }
            used[best] = true;

            var right = vR.Column(i);
            var left = vL.Column(best);
            // Hermitian inner product ⟨left|right⟩, both normalized.
            Complex ip = left.ConjugateDotProduct(right);
            double r = ip.Magnitude / (left.L2Norm() * right.L2Norm());
            modes.Add(new Mode(lamR[i], r, right));
        }
        return modes;
    }
}
