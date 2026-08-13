namespace RCPsiSquared.Core.Lindblad;

/// <summary>F26: closed form CΨ(t) for the Bell+ state under arbitrary 3-axis Pauli noise.
///
/// Bell+ = (|00⟩ + |11⟩)/√2 evolves under c_x = √γ_x · X_l, c_y = √γ_y · Y_l, c_z = √γ_z · Z_l
/// dissipators on each of the 2 qubits. Then
///
///   CΨ(t) = max(u, v) · (1 + u² + v² + w²) / 12
///   u = exp(−α·t),  α = 4·(γ_y + γ_z)
///   v = exp(−β·t),  β = 4·(γ_x + γ_z)
///   w = exp(−δ·t),  δ = 4·(γ_x + γ_y)
///
/// The l₁-coherence prefactor is max(u, v), NOT u: the proof's u·(…) form holds under the
/// WLOG re-sort α ≤ β (PROOF_MONOTONICITY_CPSI.md Part 2). For pure Y noise (α = 4γ, β = 0)
/// the physical rates violate that WLOG, and dropping the re-sort silently yields the pure-Z
/// form and the wrong K_Y = 0.0374 (the reverted 2026-04-29 error; K_Y = K_X = ln(2)/8).
///
/// Tier-1 proven (PROOF_MONOTONICITY_CPSI.md). Replaces a Lindblad-master-equation solver
/// for multi-axis Pauli noise on Bell+ states.
/// </summary>
public static class CpsiBellPlus
{
    public static double At(double gammaX, double gammaY, double gammaZ, double t)
    {
        double alpha = 4 * (gammaY + gammaZ);
        double beta = 4 * (gammaX + gammaZ);
        double delta = 4 * (gammaX + gammaY);
        double u = Math.Exp(-alpha * t);
        double v = Math.Exp(-beta * t);
        double w = Math.Exp(-delta * t);
        double l1 = Math.Max(u, v);
        return l1 * (1 + u * u + v * v + w * w) / 12.0;
    }

    /// <summary>K = γ·t at the cusp where CΨ first crosses 1/4. Channel-specific, γ-invariant.</summary>
    public static class CuspK
    {
        /// <summary>Pure Z: u³ + u = 3/2 at the crossing with u = e^{−4K}. The cubic has negative
        /// discriminant, so its real root is reachable by real radicals and K_Z = −¼·ln(u*) is
        /// exact too; 0.03735013… is a four-place reading of it, not a missing closed form. Kept
        /// as a decimal const because Symphony.KFoldN2 chains off it at compile time.</summary>
        public const double PureZ = 0.0374;

        /// <summary>Pure X: the l₁-coherence NORM max(|c_x|, |c_y|) is pinned (the XX correlation is
        /// conserved; the |00⟩⟨11| element itself is not pinned, it runs ½ → ¼), so CΨ = (1+v²)/6,
        /// which is 1/4 at v² = e^{−8K} = 1/2, giving K = ln(2)/8 = 0.0866433… A rational multiple
        /// of ln 2, so it is held as the expression and never as a decimal.</summary>
        public static readonly double PureX = Math.Log(2.0) / 8.0;

        /// <summary>Pure Y = PureX, and by a symmetry rather than a coincidence: the local Clifford
        /// S ⊗ S sends X → Y, carrying (Bell+ under pure Y) to (Bell− under pure X), and CΨ sees
        /// only |c_x|, |c_y|, |c_z|.</summary>
        public static readonly double PureY = Math.Log(2.0) / 8.0;

        /// <summary>Depolarizing (γ/3 on each axis): u³ + u/3 = 1 with u = e^{−8K/3}, again a cubic
        /// with negative discriminant, so K_depol = −⅜·ln(u*) is exact as well (0.04395476…).
        /// Kept a decimal const for the same reason as PureZ: const consumers chain at compile
        /// time, so converting this to a computed expression would break them.</summary>
        public const double Depolarizing = 0.0440;
    }
}
