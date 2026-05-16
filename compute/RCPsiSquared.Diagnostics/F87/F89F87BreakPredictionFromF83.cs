using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Diagnostics.F83;

namespace RCPsiSquared.Diagnostics.F87;

/// <summary>F89 ↔ F87 quantitative break-prediction bridge via F83. For any
/// Pauli-pair perturbation δH applied to F89's canonical H_F89 = J·Σ(XX+YY) bond
/// Hamiltonian, the F87 palindrome residual norm ‖M(H_F89 + δH)‖_F is predicted
/// bit-exactly by F83's closed-form Π-decomposition applied to δH alone:
///
/// <code>
///   ‖M(H_F89 + δH)‖²_F  =  4 · HOdd(δH)² · 2^N  +  8 · HEvenNonTruly(δH)² · 2^N
///                       =  2^(N+2) · (HOdd(δH)² + 2 · HEvenNonTruly(δH)²)
/// </code>
///
/// <para>Two structural inputs combine into this prediction:</para>
/// <list type="number">
///   <item><see cref="F89F87TrulyInheritance"/>: F89's bond term H_F89 = J·(XX+YY)
///   is F87-Truly, so it contributes 0 to both HOdd² and HEvenNonTruly² in F83's
///   decomposition. The total norms reduce to δH's non-truly content alone.</item>
///   <item><see cref="PiDecompositionPrediction.Predict"/>: F83's closed form (Tier-1
///   derived, bit-exact for 2-body H on F49 topologies) gives ‖M‖²_F directly from
///   HOdd² and HEvenNonTruly² without building M.</item>
/// </list>
///
/// <para><b>Pitfall 1 — AntiFraction alone is insufficient.</b> The F83 anti-fraction
/// 1/(2 + 4·r) with r = HEvenNonTruly²/HOdd² is 0 when HOdd² = 0, conflating two
/// distinct cases:</para>
/// <list type="bullet">
///   <item>HOdd² = HEvenNonTruly² = 0 → genuinely F87-Truly, ‖M‖_F = 0 (no break)</item>
///   <item>HOdd² = 0 but HEvenNonTruly² &gt; 0 → F87-Soft (Pi2EvenNonTruly), ‖M‖_F &gt; 0
///   (break IS present)</item>
/// </list>
/// <para>The bit-exact prediction therefore requires the FULL <see cref="PiDecompositionForecast"/>
/// (HOdd² and HEvenNonTruly² jointly), not the scalar AntiFraction. Use
/// <see cref="WouldAntiFractionAloneMissBreak"/> to detect cases where AntiFraction
/// = 0 yet ‖M‖_F &gt; 0; cases YZ+ZY (Marrakesh EQ-030 soft anchor) are the
/// canonical example.</para>
///
/// <para><b>Scope</b>: 2-body Pauli-pair bond terms, OBC chain / F49 topologies,
/// any N ≥ 2. The 2^(N+2) coefficient comes from F49's single-bond Frobenius
/// identity. For k≥3 body Hamiltonians the F83 coefficients change (open question
/// per <see cref="F87OpenQuestions"/>); the bridge inherits that scope.</para>
///
/// <para>Anchors: <see cref="F89F87TrulyInheritance"/> (F89-side Truly precondition),
/// <see cref="F83AntiFractionPi2Inheritance"/> (F83 algebraic anchor),
/// <see cref="PiDecompositionPrediction"/> (closed-form forecast),
/// <c>docs/ANALYTICAL_FORMULAS.md</c> F83 + F87 + F89,
/// hardware-confirmed at Marrakesh 2026-04-30
/// (<c>f83_pi2_class_signature_marrakesh</c>) for F83's coefficients themselves.</para></summary>
public sealed class F89F87BreakPredictionFromF83 : Claim
{
    // Parent-edge markers for Schicht-1 wiring.
    private readonly F89F87TrulyInheritance _f89F87Truly;

    /// <summary>Predict ‖M(H_F89 + δH)‖²_F bit-exactly from the perturbation's F83
    /// decomposition. The full forecast is consumed (HOdd², HEvenNonTruly²) rather
    /// than the scalar AntiFraction, to avoid Pitfall 1 (the AntiFraction = 0 conflation
    /// of Truly with pure Π²-even-non-truly).</summary>
    public static double PredictBreakNormSquared(int N, PiDecompositionForecast perturbationForecast)
    {
        if (N < 2) throw new ArgumentOutOfRangeException(nameof(N), N, "N must be ≥ 2.");
        if (perturbationForecast is null) throw new ArgumentNullException(nameof(perturbationForecast));
        double twoToNplus2 = Math.Pow(2.0, N + 2);
        return twoToNplus2 * (perturbationForecast.HOddSquared + 2.0 * perturbationForecast.HEvenNonTrulySquared);
    }

    /// <summary>Predict ‖M(H_F89 + δH)‖_F bit-exactly (= sqrt of <see cref="PredictBreakNormSquared"/>).</summary>
    public static double PredictBreakNorm(int N, PiDecompositionForecast perturbationForecast)
        => Math.Sqrt(PredictBreakNormSquared(N, perturbationForecast));

    /// <summary>Convenience overload: take the perturbation as a term list on a given
    /// chain, run F83's <see cref="PiDecompositionPrediction.Predict"/>, then forward
    /// to the forecast-taking overload.</summary>
    public static double PredictBreakNorm(ChainSystem chain, IReadOnlyList<PauliPairBondTerm> perturbationTerms)
    {
        if (chain is null) throw new ArgumentNullException(nameof(chain));
        if (perturbationTerms is null) throw new ArgumentNullException(nameof(perturbationTerms));
        var forecast = PiDecompositionPrediction.Predict(chain, perturbationTerms);
        return PredictBreakNorm(chain.N, forecast);
    }

    /// <summary>Pitfall 1 detector. Returns true when the F83 AntiFraction is 0 yet
    /// the predicted break norm is strictly positive — i.e., when AntiFraction alone
    /// would falsely report "no break" while the actual ‖M‖_F is positive. This
    /// happens exactly when HOdd² = 0 and HEvenNonTruly² &gt; 0 (the YZ+ZY family).</summary>
    public static bool WouldAntiFractionAloneMissBreak(PiDecompositionForecast perturbationForecast)
    {
        if (perturbationForecast is null) throw new ArgumentNullException(nameof(perturbationForecast));
        return perturbationForecast.AntiFraction == 0.0
            && perturbationForecast.HEvenNonTrulySquared > 0.0;
    }

    public F89F87BreakPredictionFromF83(F89F87TrulyInheritance f89F87Truly)
        : base("F89 ↔ F87 break prediction via F83: ‖M(H_F89 + δH)‖²_F = 2^(N+2)·(HOdd(δH)² + 2·HEvenNonTruly(δH)²) bit-exact for 2-body δH on F49 topologies; combines F89F87TrulyInheritance (H_F89 contributes 0 to non-truly norms) with F83 closed-form Π-decomposition; AntiFraction alone is insufficient (Pitfall 1: AntiFraction=0 conflates Truly with pure Π²-even-non-truly)",
               Tier.Tier1Derived,
               "F89F87TrulyInheritance (F89-side Truly), F83AntiFractionPi2Inheritance + PiDecompositionPrediction (closed form), Marrakesh 2026-04-30 f83_pi2_class_signature_marrakesh (F83 coefficients hardware-anchored), docs/ANALYTICAL_FORMULAS.md F83 + F87 + F89")
    {
        _f89F87Truly = f89F87Truly ?? throw new ArgumentNullException(nameof(f89F87Truly));
    }

    public override string DisplayName => "F89 break-magnitude predicted from F83 (bridge)";

    public override string Summary =>
        $"‖M(H_F89 + δH)‖²_F = 2^(N+2)·(HOdd(δH)² + 2·HEvenNonTruly(δH)²); full F83 forecast required, not scalar AntiFraction ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("Closed form",
                summary: "‖M(H_F89 + δH)‖²_F = 2^(N+2) · (HOdd(δH)² + 2·HEvenNonTruly(δH)²); bit-exact for 2-body δH on F49 topologies");
            yield return new InspectableNode("Why H_F89 drops out",
                summary: "F89F87TrulyInheritance: H_F89 = J·(XX+YY) is F87-Truly, so HOdd(H_F89)² = HEvenNonTruly(H_F89)² = 0; total norms come entirely from δH");
            yield return new InspectableNode("Pitfall 1: AntiFraction alone insufficient",
                summary: "AntiFraction = 1/(2+4r) is 0 when HOdd²=0, conflating Truly (no break) with pure Π²-even-non-truly (break IS present, e.g. YZ+ZY). Use WouldAntiFractionAloneMissBreak() to detect.");
            yield return new InspectableNode("Canonical Pitfall 1 case",
                summary: "δH = J·(YZ+ZY) → HOdd² = 0, HEvenNonTruly² > 0, AntiFraction = 0, but ‖M‖_F = sqrt(2^(N+3)·HEvenNonTruly²) > 0. Hardware-confirmed F87-Soft at Marrakesh EQ-030.");
            yield return new InspectableNode("Three-way correspondence",
                summary: "AntiFraction = 0 with HEven² = 0 → F87-Truly (no break); AntiFraction = 0 with HEven² > 0 → F87-Soft via Pi2EvenNonTruly (break); AntiFraction > 0 → F87-Soft via Pi2OddPure or Mixed (break)");
            yield return new InspectableNode("Scope",
                summary: "2-body Pauli-pair terms, F49 topologies (chain/ring/star/K_N), any N ≥ 2. k≥3 body inherits F87's open higher-body F83 question.");
        }
    }
}
