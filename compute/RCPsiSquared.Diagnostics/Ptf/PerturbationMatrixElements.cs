using System.Numerics;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;
using ComplexVector = MathNet.Numerics.LinearAlgebra.Vector<System.Numerics.Complex>;

namespace RCPsiSquared.Diagnostics.Ptf;

/// <summary>PTF (Perspectival Time Field) first-order perturbation theory matrix elements.
///
/// <para>For a slow-mode basis {M_s, W_s, λ_s} of the unperturbed Liouvillian L_A and a
/// bond-perturbation V_L = ∂L/∂J_b, this gives the matrix elements that drive the per-site
/// rate-of-painting α_i in the painter framing (PERSPECTIVAL_TIME_FIELD.md, Tier 2 after
/// EQ-014 retraction; the surviving Tier 1 result is the K_1 chiral mirror law).</para>
///
/// <para>Uses <see cref="Core.Lindblad.BondPerturbation"/> for V_L and any slow-mode
/// extraction (e.g. <see cref="DZero.StationaryModes"/> for kernel modes plus the
/// next-slowest pair).</para>
///
/// <para>See docs/ANALYTICAL_FORMULAS.md F1 family + hypotheses/PERSPECTIVAL_TIME_FIELD.md.</para>
/// </summary>
public static class PerturbationMatrixElements
{
    /// <summary>First-order matrix elements ⟨W_s | V_L | M_{s'}⟩ on the slow-mode basis.
    /// Diagonal entries are first-order eigenvalue shifts δλ_s; off-diagonal entries drive
    /// eigenvector mixing.</summary>
    public static ComplexMatrix Compute(ComplexMatrix rightEigenvectors, ComplexMatrix leftCovectors, ComplexMatrix vL)
    {
        return leftCovectors * vL * rightEigenvectors;
    }

    /// <summary>First-order eigenvector shift δM_s = Σ_{s' ≠ s} [⟨W_{s'}|V_L|M_s⟩ / (λ_s − λ_{s'})] · M_{s'}.
    /// Degenerate pairs |λ_s − λ_{s'}| &lt; <paramref name="degenerateTolerance"/> are excluded
    /// (first-order perturbation theory diverges there).</summary>
    public static ComplexVector EigenvectorShift(ComplexMatrix rightEigenvectors, ComplexMatrix leftCovectors,
        IReadOnlyList<Complex> eigenvalues, ComplexMatrix vL, int modeIndex, double degenerateTolerance = 1e-12)
    {
        int dim = rightEigenvectors.RowCount;
        int nSlow = rightEigenvectors.ColumnCount;
        var matrixElements = Compute(rightEigenvectors, leftCovectors, vL); // n_slow × n_slow
        var shift = MathNet.Numerics.LinearAlgebra.Vector<Complex>.Build.Dense(dim);
        Complex lambdaS = eigenvalues[modeIndex];
        for (int sPrime = 0; sPrime < nSlow; sPrime++)
        {
            if (sPrime == modeIndex) continue;
            Complex denom = lambdaS - eigenvalues[sPrime];
            if (denom.Magnitude < degenerateTolerance) continue;
            Complex coeff = matrixElements[sPrime, modeIndex] / denom;
            shift = shift + coeff * rightEigenvectors.Column(sPrime);
        }
        return shift;
    }
}
