using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.F1;

/// <summary>Π²-orthogonal Pythagorean decomposition of the F1 T1-residual Frobenius norm
/// (Tier 1 derived; pure-T1 case, no Hamiltonian, Z-dephasing absorbed by the +2Σγ·I shift).
///
/// <para>The total residual norm of <see cref="F1T1ResidualClosedForm"/> splits orthogonally
/// into Π²-antisymmetric and Π²-symmetric components,</para>
///
/// <code>
///     ‖M_anti(T1)‖²_F  =  4^(N−1) · Σ_l γ²_T1_l                       (F82/F84 amplitude-damping content)
///     ‖M_sym(T1)‖²_F   =  4^(N−1) · [ 2·Σ_l γ²_T1_l  +  4·(Σ_l γ_T1_l)² ]   (Π²-even complement)
/// </code>
///
/// <para>By Π²-orthogonality the two pieces sum bit-exact to the parent total
/// <see cref="F1T1ResidualClosedForm.Predict"/>:</para>
///
/// <code>
///     ‖M(T1)‖²_F = ‖M_anti(T1)‖²_F + ‖M_sym(T1)‖²_F
///                = 4^(N−1) · [ 3·Σγ² + 4·(Σγ)² ]
/// </code>
///
/// <para>Anti side. The Π²-antisymmetric part is exactly the F82/F84 amplitude-damping
/// signature <c>D_{T1, odd}</c>: a single per-site (Z, I) Pauli-basis entry of weight γ_T1_l
/// at each site. For pure T1 (no Hamiltonian) the F81 identity Π·M·Π⁻¹ = M − 2·D_{T1, odd}
/// collapses M_anti directly onto D_{T1, odd}, so ‖M_anti‖²_F = ‖D_{T1, odd}‖²_F =
/// √(Σγ²)·2^(N−1) squared = 4^(N−1)·Σγ² (matches F82 closed form). Both anti coefficients
/// are local: (local, cross) = (1, 0).</para>
///
/// <para>Sym side. The Π²-symmetric complement carries every Π²-even entry of the per-site
/// kernel M_l = Π·D_T1·Π⁻¹ + D_T1: the diagonal (X, X), (Y, Y), (Z, Z) entries that produce
/// the cooperative |tr(M_l)|² · 4^(N−2) cross-site sum, plus the residual local content not
/// absorbed by the (Z, I) anti channel. The (local, cross) split is (2, 4).</para>
///
/// <para>Anchor: <c>docs/proofs/PROOF_F1_T1_RESIDUAL_CLOSED_FORM.md</c> Step 7 (Pythagorean
/// split derived from F82's M_anti = D_{T1, odd} identity); parent total at
/// <c>compute/RCPsiSquared.Core/F1/F1T1ResidualClosedForm.cs</c>; registered alongside the
/// parent on <c>compute/RCPsiSquared.Core/F1/F1KnowledgeBase.cs</c>. F82 anchor:
/// <c>docs/proofs/PROOF_F82_T1_DISSIPATOR_CORRECTION.md</c>; F84 thermal generalisation:
/// <c>docs/proofs/PROOF_F84_AMPLITUDE_DAMPING.md</c>.</para>
/// </summary>
public sealed class F1T1ResidualPi2Decomposition : Claim
{
    /// <summary>Coefficient of Σ γ²_T1 in ‖M_anti(T1)‖²_F: 1 (the entire anti side is local).</summary>
    public const double AntisymmetricLocalCoefficient = 1.0;

    /// <summary>Coefficient of (Σ γ_T1)² in ‖M_anti(T1)‖²_F: 0 (no cross-site cooperative anti content).</summary>
    public const double AntisymmetricCrossCoefficient = 0.0;

    /// <summary>Coefficient of Σ γ²_T1 in ‖M_sym(T1)‖²_F: 2 = (parent total 3) − (anti local 1).</summary>
    public const double SymmetricLocalCoefficient = 2.0;

    /// <summary>Coefficient of (Σ γ_T1)² in ‖M_sym(T1)‖²_F: 4 (whole cooperative piece is Π²-even).</summary>
    public const double SymmetricCrossCoefficient = 4.0;

    public F1T1ResidualPi2Decomposition()
        : base("F1 T1-residual Π²-decomposition: ‖M_anti‖² = 4^(N−1)·Σγ², ‖M_sym‖² = 4^(N−1)·[2·Σγ² + 4·(Σγ)²]",
               Tier.Tier1Derived,
               "docs/proofs/PROOF_F1_T1_RESIDUAL_CLOSED_FORM.md Step 7 + " +
               "compute/RCPsiSquared.Core/F1/F1T1ResidualClosedForm.cs (parent total) + " +
               "compute/RCPsiSquared.Core/F1/F1KnowledgeBase.cs")
    { }

    /// <summary>Predicted ‖M_anti(T1)‖²_F (Π²-antisymmetric / F82-F84 amplitude-damping side).</summary>
    public static double PredictAntisymmetric(int N, IReadOnlyList<double> gammaT1) =>
        PredictWithCoefficients(N, gammaT1, AntisymmetricLocalCoefficient, AntisymmetricCrossCoefficient);

    /// <summary>Convenience overload for uniform γ_T1: ‖M_anti(T1)‖²_F = 4^(N−1) · N · γ²_T1.</summary>
    public static double PredictAntisymmetricUniform(int N, double gammaT1) =>
        PredictUniformWithCoefficients(N, gammaT1, AntisymmetricLocalCoefficient, AntisymmetricCrossCoefficient);

    /// <summary>Predicted ‖M_sym(T1)‖²_F (Π²-symmetric complement; cooperative cross-site piece lives here).</summary>
    public static double PredictSymmetric(int N, IReadOnlyList<double> gammaT1) =>
        PredictWithCoefficients(N, gammaT1, SymmetricLocalCoefficient, SymmetricCrossCoefficient);

    /// <summary>Convenience overload for uniform γ_T1: ‖M_sym(T1)‖²_F = 4^(N−1) · γ² · (2N + 4N²).</summary>
    public static double PredictSymmetricUniform(int N, double gammaT1) =>
        PredictUniformWithCoefficients(N, gammaT1, SymmetricLocalCoefficient, SymmetricCrossCoefficient);

    /// <summary>Shared kernel for the explicit-γ Predict variants: 4^(N−1) · (localCoeff·Σγ² + crossCoeff·(Σγ)²).
    /// Validates N ≥ 2 and gammaT1.Count == N before computing.</summary>
    private static double PredictWithCoefficients(
        int N, IReadOnlyList<double> gammaT1, double localCoeff, double crossCoeff)
    {
        if (N < 2) throw new ArgumentOutOfRangeException(nameof(N), $"N must be ≥ 2; got {N}");
        if (gammaT1.Count != N)
            throw new ArgumentException(
                $"gammaT1 length ({gammaT1.Count}) must equal N ({N}).", nameof(gammaT1));
        double sumSq = 0.0;
        double sum = 0.0;
        for (int l = 0; l < N; l++)
        {
            sumSq += gammaT1[l] * gammaT1[l];
            sum += gammaT1[l];
        }
        return Math.Pow(4, N - 1) * (localCoeff * sumSq + crossCoeff * sum * sum);
    }

    /// <summary>Shared kernel for the uniform-γ Predict variants: 4^(N−1) · γ² · (localCoeff·N + crossCoeff·N²).
    /// Validates N ≥ 2.</summary>
    private static double PredictUniformWithCoefficients(
        int N, double gammaT1, double localCoeff, double crossCoeff)
    {
        if (N < 2) throw new ArgumentOutOfRangeException(nameof(N), $"N must be ≥ 2; got {N}");
        // Σγ² = N·γ², (Σγ)² = N²·γ² → 4^(N−1) · γ² · (localCoeff·N + crossCoeff·N²).
        return Math.Pow(4, N - 1) * gammaT1 * gammaT1 * (localCoeff * N + crossCoeff * N * N);
    }

    public override string DisplayName =>
        "F1 T1-residual Π²-decomposition: ‖M_anti‖² = 4^(N−1)·Σγ²; ‖M_sym‖² = 4^(N−1)·[2·Σγ²+4·(Σγ)²]";

    public override string Summary =>
        "Π²-orthogonal Pythagorean split of the F1 T1-residual; anti side = F82/F84 " +
        "amplitude-damping content (local-only, coeffs (1, 0)); sym side carries the cooperative " +
        "cross-site piece (coeffs (2, 4)); the two sum bit-exact to F1T1ResidualClosedForm";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("statement (anti)",
                summary: "‖M_anti(T1)‖²_F = 4^(N−1) · Σ_l γ²_T1_l  (F82/F84 amplitude-damping content)");
            yield return new InspectableNode("statement (sym)",
                summary: "‖M_sym(T1)‖²_F  = 4^(N−1) · [ 2·Σ_l γ²_T1_l + 4·(Σ_l γ_T1_l)² ]  (Π²-even complement)");
            yield return InspectableNode.RealScalar("anti local coefficient (Σγ²)", AntisymmetricLocalCoefficient);
            yield return InspectableNode.RealScalar("anti cross-site coefficient ((Σγ)²)", AntisymmetricCrossCoefficient);
            yield return InspectableNode.RealScalar("sym local coefficient (Σγ²)", SymmetricLocalCoefficient);
            yield return InspectableNode.RealScalar("sym cross-site coefficient ((Σγ)²)", SymmetricCrossCoefficient);
            yield return new InspectableNode("Pythagorean closure",
                summary: "anti + sym = (1+2)·Σγ² + (0+4)·(Σγ)² = 3·Σγ² + 4·(Σγ)² (parent F1T1ResidualClosedForm total)");
            yield return new InspectableNode("anti side identification",
                summary: "for pure T1, M_anti = D_{T1, odd} via F81's Π·M·Π⁻¹ = M − 2·D_{T1, odd}; " +
                         "‖D_{T1, odd}‖_F = √(Σγ²)·2^(N−1) (F82 closed form)");
            yield return new InspectableNode("F82/F84 cross-reference",
                summary: "anti side maps to F82 σ⁻ (vacuum) closed form; F84 generalises to thermal Δγ = γ_↓ − γ_↑");
            yield return new InspectableNode("derivation",
                summary: "Step 7 of PROOF_F1_T1_RESIDUAL_CLOSED_FORM.md: M_anti collapses onto D_{T1, odd}; " +
                         "the Π²-orthogonal complement carries the remainder");
            // Example values are computed (not hardcoded) so they stay in sync with the formulas above.
            const int exampleN = 3;
            const double exampleGamma = 0.1;
            double exampleAnti = PredictAntisymmetricUniform(exampleN, exampleGamma);
            double exampleSym = PredictSymmetricUniform(exampleN, exampleGamma);
            double exampleTotal = F1T1ResidualClosedForm.PredictUniform(exampleN, exampleGamma);
            yield return new InspectableNode("verification",
                summary: $"section 6 of simulations/_f1_t1_residual_verify.py: N={exampleN} uniform γ_T1={exampleGamma} → " +
                         $"‖M‖²={exampleTotal:F2}, ‖M_anti‖²={exampleAnti:F2}, ‖M_sym‖²={exampleSym:F2} (machine precision)");
        }
    }
}
