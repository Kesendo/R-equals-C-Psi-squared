using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.F86.Item1Derivation;

/// <summary>F86 c=2 bare-doubled-PTF closed form for the K_b observable (Tier-1-Derived
/// 2026-05-16). Provides the explicit analytical expression that underlies the universal
/// constants <see cref="C2HwhmRatio.BareDoubledPtfXPeak"/> (= 2.196910329331) and
/// <see cref="C2HwhmRatio.BareDoubledPtfHwhmRatio"/> (= 0.671535517861, rounded as
/// 0.671535 in <see cref="C2HwhmRatio"/>).
///
/// <para>Setting (γ₀ = 1 dimensionless, t_peak = 1/4 = 1/(4γ₀)):</para>
/// <code>
///   L_2(x) = [[ -2, +2ix ],
///             [ +2ix, -6  ]]               x = Q/Q_EP = Q·g_eff/2
///   ρ_0    = (1, 0)        S = I (probe-block identity)
///   K_b(x) = 2 Re⟨ρ(t_peak)| dρ/dJ⟩
/// </code>
///
/// <para>Eigenvalues split at the EP x=1 into post-EP (complex conjugate pair) and
/// pre-EP (two real). Eigenvector geometry: |v_±|² = 2 universal. Duhamel kernel at
/// t = 1/(4γ₀) reduces via direct sympy algebra (see
/// <c>simulations/_f86_doubled_ptf_bare_floor_derivation.py</c>) to two regime-specific
/// closed forms:</para>
///
/// <code>
///   Post-EP (x > 1, ξ = √(x²-1)):
///     K_b(x) = e⁻² · x · [(ξ² + 2)·cos(ξ) − 2] / ξ⁴
///
///   Pre-EP (x &lt; 1, μ = √(1-x²)):
///     K_b(x) = e⁻² · x · [(2 − μ²)·cosh(μ) − 2] / μ⁴
///
///   EP limit (x = 1):
///     K_b(1) = −5·e⁻² / 12
/// </code>
///
/// <para>The post-EP and pre-EP forms are analytic continuations of each other under
/// ξ ↔ iμ (cos → cosh, ξ² → −μ²). Both have the same removable singularity at the EP
/// where (ξ² + 2)·cos(ξ) − 2 ~ −5ξ⁴/12 + 11ξ⁶/360 + O(ξ⁸) — the numerator vanishes at
/// the same order as the denominator. Stable evaluation near the EP uses the Taylor
/// expansion to avoid cancellation.</para>
///
/// <para><b>x_peak and HWHM_left as implicit transcendental equations:</b></para>
///
/// <code>
///   x_peak:  ξ ≥ 1 satisfying dK_b/dξ = 0 (post-EP regime).
///            Numerical solution: ξ_peak = 1.956122438683, x_peak = √(ξ²+1) = 2.196910329331.
///   x_half:  |K_b(x_half)| = |K_b(x_peak)| / 2, with x_half &lt; x_peak.
///            Numerical solution: x_half = 0.721607013629 (PRE-EP regime).
///   ratio:   HWHM_left/x_peak = (x_peak − x_half) / x_peak = 0.671535517861.
/// </code>
///
/// <para><b>The HWHM_left line crosses the EP</b>: x_half &lt; 1 &lt; x_peak. This is
/// the structural reason naïve post-EP-only analyses miss the 0.671535 value — both
/// the post-EP closed form (for K_b at x_peak) AND the pre-EP closed form (for K_b at
/// x_half) are required. Implemented in <see cref="EvaluateKb"/> which dispatches by
/// regime.</para>
///
/// <para><b>No clean algebraic closed form found for the ratio 0.671535517861</b>;
/// candidates tested against include (π−1)/π, 1 − 1/(π + 1/π), e⁻¹ + (1 − e⁻¹)/2,
/// 2/√(2π/e), √(1 − e⁻²), erf(1); closest miss is (π−1)/π = 0.6817 (off by 10⁻²).
/// The transcendental value follows from the implicit equations above and is bit-exact
/// reproducible via the explicit K_b formulas; the rounded <c>0.671535</c> in
/// <see cref="C2HwhmRatio.BareDoubledPtfHwhmRatio"/> is the 6-place truncation of this
/// constant.</para>
///
/// <para>Anchors: <see cref="C2HwhmRatio.BareDoubledPtfXPeak"/>,
/// <see cref="C2HwhmRatio.BareDoubledPtfHwhmRatio"/>,
/// <c>docs/superpowers/syntheses/2026-05-06-direction-b-attempt.md</c>,
/// <c>simulations/_f86_doubled_ptf_bare_floor_derivation.py</c> (sympy derivation +
/// brute Duhamel verification: closed-form vs eigendecomposition residual ≤ 7·10⁻¹⁷ at
/// 7 sample x values spanning both regimes).</para>
///
/// <para><b>Out of scope</b>: the per-bond-class lift from 0.671535 floor to empirical
/// values 0.7506 (Interior) and 0.7728 (Endpoint). That lift is +0.08/+0.10 and is
/// L4-blocked at 4-mode reduction level (see <c>docs/proofs/PROOF_F86B_OBSTRUCTION.md</c>);
/// requires full block-L analytical structure. Tracked separately in
/// <see cref="F86HwhmClosedFormClaim"/> via fitted (α_subclass, β_subclass) parameters.</para></summary>
public sealed class C2BareDoubledPtfClosedForm : Claim
{
    /// <summary>x_peak to 12 digits of precision: ξ_peak² + 1 with ξ_peak = 1.956122438683
    /// satisfying dK_b/dξ = 0. The 6-place rounding stored in
    /// <see cref="C2HwhmRatio.BareDoubledPtfXPeak"/> is the truncation of this value.</summary>
    public const double XPeakPrecise = 2.196910329331;

    /// <summary>HWHM_left/x_peak to 12 digits of precision. The 6-place rounding stored
    /// in <see cref="C2HwhmRatio.BareDoubledPtfHwhmRatio"/> is the truncation of this
    /// value. No clean algebraic closed form among standard π/e/√n candidates matches.</summary>
    public const double HwhmLeftOverXPeakPrecise = 0.671535517861;

    /// <summary>EP limit value: K_b(x=1) = −5·e⁻²/12. The removable-singularity limit of
    /// both regime closed forms at the EP.</summary>
    public static readonly double KbAtEp = -5.0 * Math.Exp(-2.0) / 12.0;

    /// <summary>Evaluate K_b(x) via the appropriate regime closed form. Stable across
    /// the EP via Taylor expansion when |ξ| or |μ| &lt; 0.05.</summary>
    public static double EvaluateKb(double x)
    {
        if (x <= 0) throw new ArgumentOutOfRangeException(nameof(x), x, "x must be > 0.");
        const double eMinus2 = 0.1353352832366127; // = e^(-2)
        const double small = 0.05;

        if (x > 1.0)
        {
            double xi = Math.Sqrt(x * x - 1.0);
            double bracket;
            if (xi < small)
            {
                double xi2 = xi * xi;
                double xi4 = xi2 * xi2;
                double xi6 = xi4 * xi2;
                bracket = -5.0 * xi4 / 12.0 + 11.0 * xi6 / 360.0;
            }
            else
            {
                bracket = (xi * xi + 2.0) * Math.Cos(xi) - 2.0;
            }
            double xi_4 = xi * xi * xi * xi;
            return eMinus2 * x * bracket / xi_4;
        }

        if (x < 1.0)
        {
            double mu = Math.Sqrt(1.0 - x * x);
            double bracket;
            if (mu < small)
            {
                double mu2 = mu * mu;
                double mu4 = mu2 * mu2;
                double mu6 = mu4 * mu2;
                bracket = -5.0 * mu4 / 12.0 - 11.0 * mu6 / 360.0;
            }
            else
            {
                bracket = (2.0 - mu * mu) * Math.Cosh(mu) - 2.0;
            }
            double mu_4 = mu * mu * mu * mu;
            return eMinus2 * x * bracket / mu_4;
        }

        return KbAtEp; // x == 1
    }

    public C2BareDoubledPtfClosedForm()
        : base("F86 c=2 bare-doubled-PTF closed form: K_b(x) explicit analytical expression in two regimes (post-EP cos, pre-EP cosh, EP limit −5·e⁻²/12); underlies the universal constants x_peak = 2.196910329331 (= C2HwhmRatio.BareDoubledPtfXPeak) and HWHM_left/x_peak = 0.671535517861 (= C2HwhmRatio.BareDoubledPtfHwhmRatio rounded to 6 places); x_half = 0.721607 sits in pre-EP regime so both closed forms are required for the HWHM evaluation",
               Tier.Tier1Derived,
               "compute/RCPsiSquared.Core/F86/Item1Derivation/C2HwhmRatio.cs (universal constants), docs/superpowers/syntheses/2026-05-06-direction-b-attempt.md (model setup), simulations/_f86_doubled_ptf_bare_floor_derivation.py (sympy derivation + brute Duhamel verification, residual ≤ 7e-17 at 7 sample x); transcendental constants reproducible via the explicit closed forms; per-bond-class lift to 0.7506/0.7728 remains L4-blocked (PROOF_F86B_OBSTRUCTION.md) and out of scope")
    { }

    public override string DisplayName => "F86 c=2 bare-doubled-PTF K_b closed form (post-EP cos + pre-EP cosh)";

    public override string Summary =>
        $"K_b(x) = e⁻²·x·[(ξ²+2)·cos(ξ) − 2]/ξ⁴ post-EP (ξ=√(x²−1)); = e⁻²·x·[(2−μ²)·cosh(μ) − 2]/μ⁴ pre-EP (μ=√(1−x²)); EP: K_b(1)=−5e⁻²/12. x_peak={XPeakPrecise}, HWHM/x_peak={HwhmLeftOverXPeakPrecise} ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("Closed form post-EP",
                summary: "K_b(x) = e⁻² · x · [(ξ² + 2)·cos(ξ) − 2] / ξ⁴ with ξ = √(x²−1), valid for x > 1");
            yield return new InspectableNode("Closed form pre-EP",
                summary: "K_b(x) = e⁻² · x · [(2 − μ²)·cosh(μ) − 2] / μ⁴ with μ = √(1−x²), valid for x ∈ (0, 1)");
            yield return new InspectableNode("EP limit",
                summary: "K_b(1) = −5·e⁻² / 12 ≈ −0.05636967 (removable singularity, common limit of both regimes)");
            yield return new InspectableNode("x_peak",
                summary: $"{XPeakPrecise:F12} (ξ_peak = 1.956122438683 from dK/dξ = 0 transcendental)");
            yield return new InspectableNode("x_half (HWHM crossing)",
                summary: "0.721607013629 (pre-EP! HWHM_left line crosses the EP)");
            yield return new InspectableNode("HWHM_left / x_peak",
                summary: $"{HwhmLeftOverXPeakPrecise:F12} (transcendental, no clean algebraic match)");
            yield return new InspectableNode("Lift out of scope",
                summary: "Per-bond-class +0.08/+0.10 to 0.7506/0.7728 is L4-blocked at 4-mode (PROOF_F86B_OBSTRUCTION.md), tracked separately in F86HwhmClosedFormClaim via fitted (α, β) parameters");
            yield return new InspectableNode("Verification",
                summary: "Bit-exact match (residual ≤ 7e-17) against brute Duhamel evaluator at x ∈ {0.5, 0.722, 1.0, 1.5, 2.197, 3.0, 5.0} spanning both regimes; see C2BareDoubledPtfClosedFormTests");
        }
    }
}
