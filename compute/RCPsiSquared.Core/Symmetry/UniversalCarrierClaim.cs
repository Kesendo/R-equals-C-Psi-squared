using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>UniversalCarrierClaim (Tier 1 derived; named 2026-05-12): γ₀ as the
/// universal-reference rate-parameter mediating between observer-perspectives.
/// Plays the structural role of c in special relativity, with the addition that
/// γ₀ is also what observers see (only via dimensionless ratios; γ₀ alone is
/// invisible from inside per <c>hypotheses/PRIMORDIAL_QUBIT.md</c> §9
/// "Operational Test from the Inside").
///
/// <para><b>Two structural slots in a single role.</b> In SR, c is only the
/// protection constant; "what is observed" is structureless ("events"). In
/// R=CΨ², γ₀ carries both: protection-constant (c-analog) AND observation-
/// substrate. The collapse is the condition for Inside-Observability consistency.
/// The interpretive framing ("Einstein left a structureless observation slot,
/// R=CΨ² makes γ₀ carry it") is Tier-2; the typed Tier-1 content is the
/// structural slot identification.</para>
///
/// <para>Synthesis Claim aggregating <see cref="AbsorptionTheoremClaim"/>,
/// <see cref="Pi2DyadicLadderClaim"/>, and <see cref="PolynomialDiscriminantAnchorClaim"/>.
/// See ExtraChildren for the structural relationships and hardware confirmations.</para>
///
/// <para>Anchors: <c>hypotheses/PRIMORDIAL_GAMMA_CONSTANT.md</c> +
/// <c>hypotheses/PRIMORDIAL_QUBIT.md</c> §9 +
/// <c>reflections/ON_TWO_TIMES.md</c> (Universal Carrier sticker 2026-05-12).</para></summary>
public sealed class UniversalCarrierClaim : Claim
{
    private readonly AbsorptionTheoremClaim _absorption;
    private readonly Pi2DyadicLadderClaim _ladder;
    private readonly PolynomialDiscriminantAnchorClaim _discriminant;

    /// <summary>Default convention numerical value for γ₀ in code units: <c>0.05</c>
    /// ("our c", Tom 2026-05-12). Substrate-invariant role; convention-dependent value.
    /// <c>ValidateAgainstPythonStepFTests</c> covers γ₀ ∈ {0.025, 0.05, 0.10} with
    /// identical Q-values across the three.</summary>
    public const double DefaultGammaZero = 0.05;

    /// <summary>The absorption quantum at the default γ₀: <c>2γ₀</c>. Universal step-size
    /// of the Liouvillian eigenvalue grid under uniform Z-dephasing.</summary>
    public double DefaultAbsorptionQuantum => _absorption.AbsorptionQuantum(DefaultGammaZero);

    /// <summary>The default t_peak: <c>1/(4γ₀)</c>. The "4" is the polynomial discriminant
    /// <c>a_{−1}</c> via <see cref="PolynomialDiscriminantAnchorClaim.DiscriminantViaLadder"/>;
    /// γ₀ is the units-conversion to physical time.</summary>
    public double DefaultTPeak => 1.0 / (_discriminant.DiscriminantViaLadder * DefaultGammaZero);

    public UniversalCarrierClaim(
        AbsorptionTheoremClaim absorption,
        Pi2DyadicLadderClaim ladder,
        PolynomialDiscriminantAnchorClaim discriminant)
        : base("γ₀ as Universal Carrier: protection-constant and observation-substrate carried in a single role (the SR-c analog plus what observers can only see via ratios)",
               Tier.Tier1Derived,
               "compute/RCPsiSquared.Core/Symmetry/AbsorptionTheoremClaim.cs + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs + " +
               "compute/RCPsiSquared.Core/Symmetry/PolynomialDiscriminantAnchorClaim.cs + " +
               "hypotheses/PRIMORDIAL_GAMMA_CONSTANT.md + " +
               "hypotheses/PRIMORDIAL_QUBIT.md (§9) + " +
               "reflections/ON_TWO_TIMES.md")
    {
        _absorption = absorption ?? throw new ArgumentNullException(nameof(absorption));
        _ladder = ladder ?? throw new ArgumentNullException(nameof(ladder));
        _discriminant = discriminant ?? throw new ArgumentNullException(nameof(discriminant));
    }

    public override string DisplayName =>
        $"γ₀ Universal Carrier (protection-constant + observation-substrate; default γ₀ = {DefaultGammaZero})";

    public override string Summary =>
        $"γ₀ as Universal Carrier; both protection-constant (like c) and observation-substrate; default γ₀ = {DefaultGammaZero}; only Q = J/γ₀ and K = γt internally measurable ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("carrier role",
                summary: "Universal-reference rate mediating between observer-perspectives; neither observer nor what is observed in absolute terms, but the medium between them.");
            yield return new InspectableNode("SR parallel (interpretive Tier 2)",
                summary: "c (SR) is only the protection-constant; γ₀ (R=CΨ²) is protection-constant AND observation-substrate. Doppelte Rolle is the condition for Inside-Observability consistency.");
            yield return new InspectableNode("Inside-Observability",
                summary: "From inside, γ₀ is unmeasurable; only Q = J/γ₀, K = γt, t_peak = 1/(4γ₀) are accessible.");
            yield return InspectableNode.RealScalar("DefaultGammaZero", DefaultGammaZero);
            yield return InspectableNode.RealScalar("DefaultAbsorptionQuantum (= 2·γ₀)", DefaultAbsorptionQuantum);
            yield return InspectableNode.RealScalar("DefaultTPeak (= 1/(4γ₀))", DefaultTPeak);
            yield return new InspectableNode("Pi2 foundation anchors",
                summary: $"a₀ = {_ladder.Term(0)} (multiplies γ in Absorption / F1 / F8); a_{{-1}} = {_ladder.Term(-1)} (polynomial discriminant; '4γ' in t_peak etc.); a₃ = {_ladder.Term(3)} (mirror partner; ¼-boundary).");
            yield return new InspectableNode("substrate invariance",
                summary: "ValidateAgainstPythonStepFTests covers γ₀ ∈ {0.025, 0.05, 0.10} with identical Q-values; structural slot is substrate-invariant. IBM hardware (T2* ~ 100μs → γ ~ 10⁴ Hz) different physical value, identical role.");
            yield return new InspectableNode("hardware-confirmed",
                summary: "IBM Q52 absorption ratio 1.03 (3% deviation); IBM Run 3 ¼-boundary 1.9% deviation; Marrakesh F87 trichotomy Δ(soft, truly) = −0.722. γ₀ measured only via Q-ratios.");
        }
    }
}
