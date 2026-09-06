using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.F1;

/// <summary>Centered F1 palindrome-residual closed form for the depolarizing channel
/// (Tier 1 derived; verified to machine precision at N = 2..5).
///
/// <code>
///     ‖M_F1(depol)‖²_F  =  4^(N−1) · (16/9) · Σ_l γ²_l
/// </code>
///
/// <para>In the framework's orthonormal Pauli-string basis, with the F1 shift
/// σ = Σ_l γ_l. The bare per-site kernel has ‖M_l‖²_F = 160/9 and trace −8;
/// adding 2γ_l I makes it traceless with norm squared 16/9, so all cross-site
/// terms vanish.</para>
///
/// <para>Structural properties:
/// <list type="bullet">
///   <item>Hamiltonian-independent: depol block is Frobenius-orthogonal to the H block.</item>
///   <item>γ_Z-independent: Z-dephasing block is Frobenius-orthogonal to the depol block.</item>
///   <item>Topology-independent: depolarizing is per-site only, no (B, D2) graph dependence.</item>
///   <item>Purely per-site: (16/9)·Σγ²; F1 centering removes the bare residual's
///         cooperative 16·(Σγ)² term exactly.</item>
///   <item><b>Π²-decomposition is trivial:</b> M_l is Pauli-basis-diagonal, so
///         Π·M·Π⁻¹ = M exactly and M_anti = 0. Contrast T1, whose M_anti = D_{T1, odd}
///         carries F82/F84 amplitude-damping content.</item>
///   <item><b>F1 shift σ = Σγ.</b> It removes the diagonal mean but not the ±2γ/3
///         split. F5's (2/3)Σγ is the extreme pair-sum shortfall, equivalently
///         the spectral norm of the centered diagonal residual.</item>
/// </list></para>
///
/// <para>Anchor: <c>docs/proofs/PROOF_F1_DEPOL_RESIDUAL_CLOSED_FORM.md</c> (Steps 1-7).
/// Verification: <c>simulations/f1_depol_residual_verify.py</c> (seven sections,
/// match to machine precision N = 2..5 across uniform and non-uniform γ).</para>
/// </summary>
public sealed class F1DepolResidualClosedForm : Claim
{
    /// <summary>Per-site Frobenius norm squared of M_l = Π·D_depol·Π⁻¹ + D_depol at γ = 1:
    /// 160/9 = 2·(4/3)² + 2·(8/3)² = 2·(16/9 + 64/9) = 2·80/9.</summary>
    public const double PerSiteFrobeniusSquared = 160.0 / 9.0;

    /// <summary>|tr(M_l)|² at γ = 1: 64 = (−8)² where tr(M_l) = −4/3 − 4/3 − 8/3 − 8/3 = −8.
    /// Drives the bare residual's cross-site term, which F1 centering removes.</summary>
    public const double PerSiteTraceSquared = 64.0;

    /// <summary>Coefficient of Σ γ² in the closed form: c1 = ‖M_l‖² − |tr(M_l)|²/4 =
    /// 160/9 − 16 = 160/9 − 144/9 = 16/9.</summary>
    public const double LocalCoefficient = 16.0 / 9.0;

    /// <summary>Coefficient of (Σ γ)² after F1 centering. It is zero because each
    /// centered local kernel is traceless.</summary>
    public const double CrossSiteCoefficient = 0.0;

    public F1DepolResidualClosedForm()
        : base("F1 depol-residual closed form: ‖M_F1(depol)‖² = 4^(N−1) · (16/9)·Σγ²",
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
        for (int l = 0; l < N; l++)
        {
            sumSq += gamma[l] * gamma[l];
        }
        return Math.Pow(4, N - 1) * LocalCoefficient * sumSq;
    }

    /// <summary>Convenience overload: uniform γ across all N sites.</summary>
    public static double PredictUniform(int N, double gamma)
    {
        if (N < 2) throw new ArgumentOutOfRangeException(nameof(N), $"N must be ≥ 2; got {N}");
        // Σγ² = N·γ² → 4^(N−1) · γ² · (16/9)·N.
        return Math.Pow(4, N - 1) * gamma * gamma * LocalCoefficient * N;
    }

    public override string DisplayName =>
        "F1 depol-residual: ‖M_F1(depol)‖² = 4^(N−1) · (16/9)·Σγ²";

    public override string Summary =>
        "Centered depol-block closed form in framework Pauli basis; H-independent, topology-independent; " +
        "local coefficient 16/9, cross-site coefficient 0; M_anti = 0; F1 shift σ = Σγ";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("statement",
                summary: "‖M_F1(depol)‖²_F = 4^(N−1) · (16/9)·Σ_l γ²_l, with σ = Σ_l γ_l");
            yield return InspectableNode.RealScalar("per-site ‖M_l‖²_F (γ=1)", PerSiteFrobeniusSquared);
            yield return InspectableNode.RealScalar("per-site |tr(M_l)|² (γ=1)", PerSiteTraceSquared);
            yield return InspectableNode.RealScalar("local coefficient (Σγ²)", LocalCoefficient);
            yield return InspectableNode.RealScalar("centered cross-site coefficient ((Σγ)²)", CrossSiteCoefficient);
            yield return new InspectableNode("derivation",
                summary: "per-site M_l = Π·D_depol·Π⁻¹ + D_depol = diag(−4/3, −4/3, −8/3, −8/3); " +
                         "tensor-assembled via tr(M_l† M_l′) = |tr(M_l)|² · 4^(N−2) for l ≠ l′");
            yield return new InspectableNode("orthogonality",
                summary: "depol block Frobenius-orthogonal to H and Z blocks; cross-terms vanish identically; " +
                         "no graph-parameter (B, D2) dependence (depol is per-site only)");
            yield return new InspectableNode("Π²-decomposition (trivial)",
                summary: "M_l is Pauli-basis-diagonal ⟹ Π·M·Π⁻¹ = M exactly ⟹ M_anti = 0 " +
                         "(contrast T1, whose M_anti = D_{T1, odd} carries F82/F84 amplitude-damping content)");
            yield return new InspectableNode("F1-centering shift σ = Σγ",
                summary: "Adding 2Σγ·I removes the bare residual's mean and cross-site term, leaving local " +
                         "entries ±2γ_l/3. F5's (2/3)Σγ is the extreme pair-sum shortfall / spectral norm, " +
                         "not a trace projection of the bare residual.");
            yield return new InspectableNode("verification",
                summary: "verified to machine precision at N = 2..5, uniform and non-uniform γ (simulations/f1_depol_residual_verify.py)");
        }
    }
}
