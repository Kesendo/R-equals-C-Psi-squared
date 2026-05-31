using System.Numerics;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;
using ComplexVector = MathNet.Numerics.LinearAlgebra.Vector<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Propagation;

/// <summary>The loop, as a reusable primitive: eigendecompose a Liouvillian L once, then
/// propagate a vectorized density matrix vec(ρ) over a τ-grid and read observables off it.
///
/// <para>exp(L·τ)·v = R·diag(exp(λ_i·τ))·R⁻¹·v, with (λ, R) = Evd(L). Each observable is a
/// covector w with ⟨O⟩(τ) = Σ_k w[k]·vec(ρ(τ))[k] (the unconjugated bilinear sum). The vec
/// convention (row-major / C-order vs column-major / F-order) lives entirely in the caller:
/// this primitive is convention-agnostic, it only needs vec(ρ₀) and the covectors to agree
/// with the L that was handed in.</para>
///
/// <para>This is the eig-once propagation idiom otherwise inlined in
/// <c>BlockCpsiTrajectory.BuildCore</c> and <c>StationaryModes</c>; extracted so the post-EP
/// flow does not add a third copy. Diagnostic-scale dense Evd (small N).</para></summary>
public static class SpectralPropagator
{
    /// <summary>Returns <c>result[observable][τ] = Re(Σ_k covector[k]·vec(ρ(τ))[k])</c>.</summary>
    public static double[][] Evolve(
        ComplexMatrix L,
        ComplexVector rho0Vec,
        IReadOnlyList<ComplexVector> observableCovectors,
        IReadOnlyList<double> tauGrid)
    {
        if (L is null) throw new ArgumentNullException(nameof(L));
        if (rho0Vec is null) throw new ArgumentNullException(nameof(rho0Vec));
        if (rho0Vec.Count != L.RowCount)
            throw new ArgumentException($"rho0Vec dim {rho0Vec.Count} != L dim {L.RowCount}", nameof(rho0Vec));

        var evd = L.Evd();
        var R = evd.EigenVectors;
        var lambda = evd.EigenValues;
        var c0 = R.Solve(rho0Vec);            // R·c0 = vec(ρ₀)  ⇒  mode coefficients
        int d2 = L.RowCount;

        var result = new double[observableCovectors.Count][];
        for (int o = 0; o < result.Length; o++) result[o] = new double[tauGrid.Count];

        var expc = ComplexVector.Build.Dense(d2);
        for (int ti = 0; ti < tauGrid.Count; ti++)
        {
            double tau = tauGrid[ti];
            for (int i = 0; i < d2; i++) expc[i] = Complex.Exp(lambda[i] * tau) * c0[i];
            var vt = R * expc;                // vec(ρ(τ))
            for (int o = 0; o < observableCovectors.Count; o++)
            {
                var w = observableCovectors[o];
                Complex acc = Complex.Zero;
                for (int k = 0; k < d2; k++) acc += w[k] * vt[k];
                result[o][ti] = acc.Real;
            }
        }
        return result;
    }
}
