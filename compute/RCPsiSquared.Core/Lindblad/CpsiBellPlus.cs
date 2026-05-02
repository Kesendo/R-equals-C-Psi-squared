namespace RCPsiSquared.Core.Lindblad;

/// <summary>F26: closed form CΨ(t) for the Bell+ state under arbitrary 3-axis Pauli noise.
///
/// Bell+ = (|00⟩ + |11⟩)/√2 evolves under c_x = √γ_x · X_l, c_y = √γ_y · Y_l, c_z = √γ_z · Z_l
/// dissipators on each of the 2 qubits. Then
///
///   CΨ(t) = u · (1 + u² + v² + w²) / 12
///   u = exp(−α·t),  α = 4·(γ_y + γ_z)
///   v = exp(−β·t),  β = 4·(γ_x + γ_z)
///   w = exp(−δ·t),  δ = 4·(γ_x + γ_y)
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
        return u * (1 + u * u + v * v + w * w) / 12.0;
    }

    /// <summary>K = γ·t at the cusp where CΨ first crosses 1/4. Channel-specific, γ-invariant.</summary>
    public static class CuspK
    {
        public const double PureZ = 0.0374;
        public const double PureX = 0.0867;
        public const double PureY = 0.0374; // = PureZ by F26 functional symmetry
        public const double Depolarizing = 0.0440; // γ/3 on each axis
    }
}
