using System.Numerics;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Diagnostics.Ptf;

/// <summary>The first-order mixing of the slow modes under a bond perturbation: the matrix
/// ⟨W_s | V_L | M_{s'}⟩ on the slowest modes of an unperturbed Liouvillian L_A. This is the
/// operator-algebra face of "the in-between is more than a rotation".
///
/// <para>On the crossover axis the in-between is a pure rotation Ad_{R_z} (a unitary on operator
/// space) and every eigenvalue is frozen. On the J-defect axis it is this mixing: the diagonal entry
/// ⟨W_s|V_L|M_s⟩ is the first-order eigenvalue shift δλ_s, and the off-diagonal entries drive the
/// eigenvector mixing δM_s = Σ_{s'≠s} [⟨W_{s'}|V_L|M_s⟩ / (λ_s − λ_{s'})] M_{s'}, the PTF α_i.</para>
///
/// <para>The shift profile is the telling reading. The kernel (the steady states, Re λ ≈ 0) is
/// protected to machine zero by U(1) excitation conservation: V_L (a hopping bilinear) preserves the
/// excitation sector, so its diagonal on the excitation-projector modes vanishes exactly. The slow
/// coherences (Re λ ≈ −2γ) are NOT protected: their diagonal is the bond-energy expectation
/// ⟨ψ|H_bond|ψ⟩ and they shift at first order. So the J-defect protects the kernel but moves the
/// rest of the spectrum (which nonetheless stays palindromic), and the off-diagonal mixing is alive
/// throughout, unlike the crossover's rigid rotation.</para>
///
/// <para>Slow modes are the <c>slowCount</c> eigenvalues of smallest |Re λ| (slowest decay). Left
/// covectors W are the rows of R⁻¹ (R the right-eigenvector matrix); this assumes L_A is
/// diagonalizable, which holds for dephasing plus a generic Hamiltonian.</para></summary>
public static class SlowModeMixing
{
    /// <summary>The mixing reading: the matrix, the per-mode first-order shift magnitudes
    /// |⟨W_s|V_L|M_s⟩| (≈ 0 on the protected kernel, alive on the moving coherences), the off-diagonal
    /// Frobenius mass (the eigenvector mixing), and the slow eigenvalues.</summary>
    public sealed record Reading(
        ComplexMatrix Mixing,
        IReadOnlyList<double> DiagonalShiftMagnitudes,
        double OffDiagonalMass,
        IReadOnlyList<Complex> SlowEigenvalues);

    /// <summary>Diagonalize <paramref name="lA"/>, take the <paramref name="slowCount"/> slowest modes
    /// (smallest |Re λ|), and form ⟨W_s | V_L | M_{s'}⟩ with V_L = <paramref name="vL"/>.</summary>
    public static Reading Compute(ComplexMatrix lA, ComplexMatrix vL, int slowCount)
    {
        var evd = lA.Evd();
        var evals = evd.EigenValues;
        var r = evd.EigenVectors;          // columns = right eigenvectors M_s
        int total = evals.Count;
        int k = Math.Min(slowCount, total);

        // Slowest modes: smallest |Re λ| (Re λ ≤ 0, so this is the least-negative / closest to 0).
        var idx = Enumerable.Range(0, total).OrderBy(i => Math.Abs(evals[i].Real)).Take(k).ToArray();

        var w = r.Inverse();               // rows = left covectors W_s (dual basis)
        var rSlow = ComplexMatrix.Build.Dense(r.RowCount, k, (row, j) => r[row, idx[j]]);
        var wSlow = ComplexMatrix.Build.Dense(k, w.ColumnCount, (i, col) => w[idx[i], col]);
        var mixing = PerturbationMatrixElements.Compute(rSlow, wSlow, vL); // k × k

        var diagonalShifts = new double[k];
        double off2 = 0.0;
        for (int i = 0; i < k; i++)
            for (int j = 0; j < k; j++)
            {
                double mag = mixing[i, j].Magnitude;
                if (i == j) diagonalShifts[i] = mag;
                else off2 += mag * mag;
            }

        var slowEvals = idx.Select(i => evals[i]).ToList();
        return new Reading(mixing, diagonalShifts, Math.Sqrt(off2), slowEvals);
    }
}
