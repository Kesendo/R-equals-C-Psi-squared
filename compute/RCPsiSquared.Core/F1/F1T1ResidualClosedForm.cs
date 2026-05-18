using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.F1;

/// <summary>F1 palindrome-residual closed form for the T1 amplitude-damping block
/// (Tier 1 derived; bit-exact at N = 2..5).
///
/// <code>
///     ‖M(T1)‖²_F  =  4^(N−1) · [ 3·Σ_l γ²_T1_l  +  4·(Σ_l γ_T1_l)² ]
/// </code>
///
/// <para>In the framework's orthonormal Pauli-string basis. The (3, 4) pair derives
/// from the per-site kernel M_l = Π·D_T1·Π⁻¹ + D_T1 with ‖M_l‖²_F = 7 and
/// |tr(M_l)|² = 16: multi-site assembly via tr(M_l† M_l′) = |tr(M_l)|² · 4^(N−2)
/// for l ≠ l′ gives (7 − 4)·Σγ² + 4·(Σγ)².</para>
///
/// <para>Structural properties:
/// <list type="bullet">
///   <item>Hamiltonian-independent: T1 block is Frobenius-orthogonal to the H block.</item>
///   <item>γ_Z-independent: Z-dephasing is absorbed by the +2Σγ·I shift.</item>
///   <item>Per-site (3·Σγ²) + cross-site (4·(Σγ)²) split; cooperative piece dominates
///         as N grows (ratio (4/3)·N for uniform γ_T1).</item>
/// </list></para>
///
/// <para>Anchor: <c>docs/proofs/PROOF_F1_T1_RESIDUAL_CLOSED_FORM.md</c> (Steps 1-5).
/// Verification: <c>simulations/_f1_t1_residual_verify.py</c> (six sections,
/// bit-exact match to machine precision N = 2..5 across uniform and non-uniform γ_T1).</para>
/// </summary>
public sealed class F1T1ResidualClosedForm : Claim
{
    /// <summary>Per-site Frobenius norm squared of M_l = Π·D_T1·Π⁻¹ + D_T1 at γ = 1.</summary>
    public const double PerSiteFrobeniusSquared = 7.0;

    /// <summary>|tr(M_l)|² at γ = 1; drives the cross-site (cooperative) coefficient.</summary>
    public const double PerSiteTraceSquared = 16.0;

    /// <summary>Coefficient of Σ γ²_T1 in the closed form: c1 = ‖M_l‖² − |tr(M_l)|²/4 = 7 − 4 = 3.</summary>
    public const double LocalCoefficient = 3.0;

    /// <summary>Coefficient of (Σ γ_T1)² in the closed form: c2 = |tr(M_l)|²/4 = 16/4 = 4.</summary>
    public const double CrossSiteCoefficient = 4.0;

    public F1T1ResidualClosedForm()
        : base("F1 T1-residual closed form: ‖M(T1)‖² = 4^(N−1) · [3·Σγ² + 4·(Σγ)²]",
               Tier.Tier1Derived,
               "docs/proofs/PROOF_F1_T1_RESIDUAL_CLOSED_FORM.md")
    { }

    /// <summary>Predicted ‖M(T1)‖²_F at the given chain length and per-site γ_T1 rates.</summary>
    public static double Predict(int N, IReadOnlyList<double> gammaT1)
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
        return Math.Pow(4, N - 1) * (LocalCoefficient * sumSq + CrossSiteCoefficient * sum * sum);
    }

    /// <summary>Convenience overload: uniform γ_T1 across all N sites.</summary>
    public static double PredictUniform(int N, double gammaT1)
    {
        if (N < 2) throw new ArgumentOutOfRangeException(nameof(N), $"N must be ≥ 2; got {N}");
        // Σγ² = N·γ², (Σγ)² = N²·γ² → 4^(N−1) · γ² · (3N + 4N²).
        return Math.Pow(4, N - 1) * gammaT1 * gammaT1 * (LocalCoefficient * N + CrossSiteCoefficient * N * N);
    }

    public override string DisplayName =>
        "F1 T1-residual: ‖M(T1)‖² = 4^(N−1) · [3·Σγ²_T1 + 4·(Σγ_T1)²]";

    public override string Summary =>
        "T1-block closed form in framework Pauli basis; H-independent, γ_Z-independent, " +
        "Frobenius-orthogonal to H and Z blocks; (3, 4) from per-site ‖M_l‖² = 7, |tr(M_l)|² = 16";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("statement",
                summary: "‖M(T1)‖²_F = 4^(N−1) · [3·Σ_l γ²_T1_l + 4·(Σ_l γ_T1_l)²]");
            yield return InspectableNode.RealScalar("per-site ‖M_l‖²_F (γ=1)", PerSiteFrobeniusSquared);
            yield return InspectableNode.RealScalar("per-site |tr(M_l)|² (γ=1)", PerSiteTraceSquared);
            yield return InspectableNode.RealScalar("local coefficient (Σγ²)", LocalCoefficient);
            yield return InspectableNode.RealScalar("cross-site coefficient ((Σγ)²)", CrossSiteCoefficient);
            yield return new InspectableNode("derivation",
                summary: "per-site M_l = Π·D_T1·Π⁻¹ + D_T1, tensor-assembled via " +
                         "tr(M_l† M_l′) = |tr(M_l)|² · 4^(N−2) for l ≠ l′");
            yield return new InspectableNode("orthogonality",
                summary: "T1 block Frobenius-orthogonal to H and Z blocks; cross-terms vanish identically");
            yield return new InspectableNode("verification",
                summary: "bit-exact at N = 2..5, uniform and non-uniform γ_T1 (simulations/_f1_t1_residual_verify.py)");
        }
    }
}
