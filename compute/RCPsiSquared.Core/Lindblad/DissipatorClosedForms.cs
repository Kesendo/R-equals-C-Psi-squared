using System.Numerics;

namespace RCPsiSquared.Core.Lindblad;

/// <summary>Closed-form constants for the F1-palindrome residual contribution of single-class
/// and cross-class dissipators given their Pauli decomposition.
///
/// For c = α·X + β·Y + δ·Z + ε·I (the ε·I component is Π-trivial), the per-class contribution
/// to ‖M(c)‖² at H = 0 is
///
///   ‖M(c)‖²_F = 4^(N−1) · [ c1 · Σ γ_l² + c2 · (Σ γ_l)² ]
///
/// with closed forms verified empirically on 25+ test classes (pure Paulis, σ⁻, σ⁺, real and
/// complex superpositions, scaled operators).
/// </summary>
public static class DissipatorClosedForms
{
    /// <summary>(c1, c2) for a single-class dissipator decomposed as c = α·X + β·Y + δ·Z (+ ε·I).
    /// Closed forms reflect the Π-bit_b-parity structure (X is bit_b-even, Y/Z are bit_b-odd).</summary>
    public static (double C1, double C2) C1C2FromPauli(Complex alpha, Complex beta, Complex delta)
    {
        double a2 = alpha.Magnitude * alpha.Magnitude;
        double b2 = beta.Magnitude * beta.Magnitude;
        double d2 = delta.Magnitude * delta.Magnitude;
        double normSq = a2 + b2 + d2;
        double imAB = (alpha * Complex.Conjugate(beta)).Imaginary;
        double imAD = (alpha * Complex.Conjugate(delta)).Imaginary;
        double c1 = 16 * a2 * normSq + 32 * b2 * d2 + 16 * imAB * imAB + 16 * imAD * imAD;
        double c2 = 16 * normSq * normSq;
        return (c1, c2);
    }

    /// <summary>Cross-term d2 for two single-class dissipators c_k = α_k·X + β_k·Y + δ_k·Z + ε_k·I.
    /// Universal closed form (verified across 196 class-pair combinations, all topologies):
    ///   d2 = 32 · ‖c1_traceless‖² · ‖c2_traceless‖²</summary>
    public static double D2FromPauli(Complex alpha1, Complex beta1, Complex delta1,
        Complex alpha2, Complex beta2, Complex delta2)
    {
        double n1 = alpha1.Magnitude * alpha1.Magnitude + beta1.Magnitude * beta1.Magnitude + delta1.Magnitude * delta1.Magnitude;
        double n2 = alpha2.Magnitude * alpha2.Magnitude + beta2.Magnitude * beta2.Magnitude + delta2.Magnitude * delta2.Magnitude;
        return 32.0 * n1 * n2;
    }
}
