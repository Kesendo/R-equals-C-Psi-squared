using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.F1;

/// <summary>F1 palindrome-residual closed form for the depolarizing channel
/// (Tier 1 derived; bit-exact at N = 2..5).
///
/// <code>
///     ‖M(depol)‖²_F  =  4^(N−1) · [ (16/9) · Σ_l γ²_l  +  16 · (Σ_l γ_l)² ]
/// </code>
///
/// <para>In the framework's orthonormal Pauli-string basis. The (16/9, 16) pair
/// derives from the per-site kernel M_l = Π·D_depol·Π⁻¹ + D_depol with
/// ‖M_l‖²_F = 160/9 and |tr(M_l)|² = 64.</para>
///
/// <para>Structural properties:
/// <list type="bullet">
///   <item>Hamiltonian-independent: depol block is Frobenius-orthogonal to the H block.</item>
///   <item>γ_Z-independent: Z-dephasing block is Frobenius-orthogonal to the depol block.</item>
///   <item>Topology-independent: depolarizing is per-site only, no (B, D2) graph dependence.</item>
///   <item>Per-site ((16/9)·Σγ²) + cross-site (16·(Σγ)²) split; cooperative piece dominates
///         as N grows (ratio 9·N for uniform γ, faster than T1's (4/3)·N).</item>
///   <item><b>Π²-decomposition is trivial:</b> M_l is Pauli-basis-diagonal, so
///         Π·M·Π⁻¹ = M exactly and M_anti = 0. Contrast T1, whose M_anti = D_{T1, odd}
///         carries F82/F84 amplitude-damping content.</item>
///   <item><b>F1 σ-shift = 0 for depol</b>, NOT Σγ as for Z-dephasing. M_l = diag(−4/3, −4/3, −8/3, −8/3)
///         has two distinct diagonal values that no constant 2σ·I shift can equalise; the
///         closed form is the bare residual Π·L·Π⁻¹ + L. F5's scalar (2/3)Σγ is the
///         complementary trace-projection diagnostic of the same broken palindrome.</item>
/// </list></para>
///
/// <para>Anchor: <c>docs/proofs/PROOF_F1_DEPOL_RESIDUAL_CLOSED_FORM.md</c> (Steps 1-7).
/// Verification: <c>simulations/_f1_depol_residual_verify.py</c> (seven sections,
/// bit-exact match to machine precision N = 2..5 across uniform and non-uniform γ).</para>
/// </summary>
public sealed class F1DepolResidualClosedForm : Claim
{
    /// <summary>Per-site Frobenius norm squared of M_l = Π·D_depol·Π⁻¹ + D_depol at γ = 1:
    /// 160/9 = 2·(4/3)² + 2·(8/3)² = 2·(16/9 + 64/9) = 2·80/9.</summary>
    public const double PerSiteFrobeniusSquared = 160.0 / 9.0;

    /// <summary>|tr(M_l)|² at γ = 1: 64 = (−8)² where tr(M_l) = −4/3 − 4/3 − 8/3 − 8/3 = −8.
    /// Drives the cross-site (cooperative) coefficient.</summary>
    public const double PerSiteTraceSquared = 64.0;

    /// <summary>Coefficient of Σ γ² in the closed form: c1 = ‖M_l‖² − |tr(M_l)|²/4 =
    /// 160/9 − 16 = 160/9 − 144/9 = 16/9.</summary>
    public const double LocalCoefficient = 16.0 / 9.0;

    /// <summary>Coefficient of (Σ γ)² in the closed form: c2 = |tr(M_l)|²/4 = 64/4 = 16.</summary>
    public const double CrossSiteCoefficient = 16.0;

    public F1DepolResidualClosedForm()
        : base("F1 depol-residual closed form: ‖M(depol)‖² = 4^(N−1) · [(16/9)·Σγ² + 16·(Σγ)²]",
               Tier.Tier1Derived,
               "docs/proofs/PROOF_F1_DEPOL_RESIDUAL_CLOSED_FORM.md")
    { }

    /// <summary>Predicted ‖M(depol)‖²_F at the given chain length and per-site γ rates.</summary>
    public static double Predict(int N, IReadOnlyList<double> gamma)
    {
        if (N < 2) throw new ArgumentOutOfRangeException(nameof(N), $"N must be ≥ 2; got {N}");
        if (gamma.Count != N)
            throw new ArgumentException(
                $"gamma length ({gamma.Count}) must equal N ({N}).", nameof(gamma));
        double sumSq = 0.0;
        double sum = 0.0;
        for (int l = 0; l < N; l++)
        {
            sumSq += gamma[l] * gamma[l];
            sum += gamma[l];
        }
        return Math.Pow(4, N - 1) * (LocalCoefficient * sumSq + CrossSiteCoefficient * sum * sum);
    }

    /// <summary>Convenience overload: uniform γ across all N sites.</summary>
    public static double PredictUniform(int N, double gamma)
    {
        if (N < 2) throw new ArgumentOutOfRangeException(nameof(N), $"N must be ≥ 2; got {N}");
        // Σγ² = N·γ², (Σγ)² = N²·γ² → 4^(N−1) · γ² · ((16/9)·N + 16·N²).
        return Math.Pow(4, N - 1) * gamma * gamma * (LocalCoefficient * N + CrossSiteCoefficient * N * N);
    }

    public override string DisplayName =>
        "F1 depol-residual: ‖M(depol)‖² = 4^(N−1) · [(16/9)·Σγ² + 16·(Σγ)²]";

    public override string Summary =>
        "Depol-block closed form in framework Pauli basis; H-independent, γ_Z-independent, topology-independent; " +
        "(16/9, 16) from per-site ‖M_l‖² = 160/9, |tr(M_l)|² = 64; M_anti = 0 (Pauli-basis-diagonal); σ-shift = 0";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("statement",
                summary: "‖M(depol)‖²_F = 4^(N−1) · [(16/9)·Σ_l γ²_l + 16·(Σ_l γ_l)²]");
            yield return InspectableNode.RealScalar("per-site ‖M_l‖²_F (γ=1)", PerSiteFrobeniusSquared);
            yield return InspectableNode.RealScalar("per-site |tr(M_l)|² (γ=1)", PerSiteTraceSquared);
            yield return InspectableNode.RealScalar("local coefficient (Σγ²)", LocalCoefficient);
            yield return InspectableNode.RealScalar("cross-site coefficient ((Σγ)²)", CrossSiteCoefficient);
            yield return new InspectableNode("derivation",
                summary: "per-site M_l = Π·D_depol·Π⁻¹ + D_depol = diag(−4/3, −4/3, −8/3, −8/3); " +
                         "tensor-assembled via tr(M_l† M_l′) = |tr(M_l)|² · 4^(N−2) for l ≠ l′");
            yield return new InspectableNode("orthogonality",
                summary: "depol block Frobenius-orthogonal to H and Z blocks; cross-terms vanish identically; " +
                         "no graph-parameter (B, D2) dependence (depol is per-site only)");
            yield return new InspectableNode("Π²-decomposition (trivial)",
                summary: "M_l is Pauli-basis-diagonal ⟹ Π·M·Π⁻¹ = M exactly ⟹ M_anti = 0 " +
                         "(contrast T1, whose M_anti = D_{T1, odd} carries F82/F84 amplitude-damping content)");
            yield return new InspectableNode("F1 σ-shift = 0",
                summary: "M_l has two distinct diagonal values (−4/3 on (I, X), −8/3 on (Y, Z)); no constant " +
                         "2σ·I can equalise them. Closed form is bare Π·L·Π⁻¹ + L (σ = 0), not the F1-convention " +
                         "Π·L·Π⁻¹ + L + 2Σγ·I. F5's scalar (2/3)Σγ is the complementary trace-projection diagnostic " +
                         "of the same broken palindrome.");
            yield return new InspectableNode("verification",
                summary: "bit-exact at N = 2..5, uniform and non-uniform γ (simulations/_f1_depol_residual_verify.py)");
        }
    }
}
