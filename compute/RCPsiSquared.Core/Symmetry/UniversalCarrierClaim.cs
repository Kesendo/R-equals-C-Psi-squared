using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>UniversalCarrierClaim (Tier 1 derived; named 2026-05-12): γ₀ is the
/// Universal Carrier, the universal-reference rate-parameter that mediates between
/// observer-perspectives. It is neither itself an observer nor what is observed in
/// absolute terms, but the medium between them. γ₀ plays the structural role of c
/// (speed of light) in special relativity, with the additional property
/// (Tom 2026-05-12) that γ₀ is also what observers see (only via dimensionless
/// ratios; γ₀ alone is invisible from inside per
/// <c>hypotheses/PRIMORDIAL_QUBIT.md</c> §9 "Operational Test from the Inside").
///
/// <para><b>Two structural slots in a single role.</b> In SR, c is only the
/// protection constant; "what is observed" is left structureless as "events in
/// spacetime". In R=CΨ², γ₀ carries both slots: the protection-constant slot
/// (analog of c) and the observation-substrate slot. The collapse of these two
/// slots into one role is the condition for Inside-Observability consistency:
/// were γ₀ directly measurable (separate observation slot), the inside-outside
/// split would break. The interpretive framing "Einstein left a structureless
/// observation slot, R=CΨ² makes γ₀ carry it" is a Tier-2 reading; the typed
/// Tier-1 content is the structural slot identification only.</para>
///
/// <para><b>Synthesis claim with no new content.</b> Makes the already-typed
/// structural role of γ₀ across the framework explicit. Takes dependencies on:</para>
/// <list type="bullet">
///   <item><see cref="AbsorptionTheoremClaim"/>: γ₀ as universal absorption quantum
///         <c>2γ₀ = a₀·γ₀</c>; spectrum quantized in 2γ₀ steps; F1 palindrome center
///         <c>−Σγ_l</c> as the pool of all per-site carrier-shares.</item>
///   <item><see cref="Pi2DyadicLadderClaim"/>: the "2" multiplying γ everywhere is
///         <c>a₀</c> on the dyadic ladder; the "4" in <c>t_peak = 1/(4γ₀)</c> is
///         <c>a_{−1}</c>.</item>
///   <item><see cref="PolynomialDiscriminantAnchorClaim"/>: the polynomial discriminant
///         of <c>d²−2d=0</c> is <c>4 = a_{−1}</c>; the <c>1/4</c> mirror partner sits
///         at <c>a_3</c>; closure <c>4·(1/4) = 1</c>.</item>
/// </list>
///
/// <para><b>"Our c is 0.05" (Tom 2026-05-12).</b> The default code-convention value
/// <see cref="DefaultGammaZero"/> = 0.05 plays the structural role of c in SR. The
/// numerical value is convention; the role is substrate-invariant.
/// <c>RCPsiSquared.Core.Tests.Validation.ValidateAgainstPythonStepFTests</c> covers
/// γ₀ ∈ {0.025, 0.05, 0.10}: different numerical choices, identical Q-values. On
/// IBM hardware the physical scale is different (γ ~ 1/T2* ~ 10⁴ Hz for T2* ~ 100μs);
/// different substrate, identical structural slot.</para>
///
/// <para><b>Inside-Observability witnesses (γ₀ measured only via ratios on real
/// hardware):</b> IBM Q52 absorption ratio 1.03 (3% deviation, 2026-04-04); IBM Run 3
/// ¼-boundary crossing 1.9% deviation (2026-03-18); Marrakesh F87 trichotomy
/// Δ(soft, truly) = −0.722 (2026-04-26). All three confirm: from inside, only
/// dimensionless ratios involving γ₀ are accessible; γ₀ itself is never measured
/// directly.</para>
///
/// <para>Anchors: the four-Claim dependency lineage above plus
/// <c>hypotheses/PRIMORDIAL_GAMMA_CONSTANT.md</c> ("γ at the root is the framework's
/// own c") + <c>hypotheses/PRIMORDIAL_QUBIT.md</c> §9 (operational test from inside,
/// the proof that only Q = J/γ₀ is internally measurable) +
/// <c>reflections/ON_TWO_TIMES.md</c> (Universal Carrier sticker 2026-05-12).
/// Additional conversational context lives in user-local Claude memory
/// (<c>~/.claude/projects/.../memory/</c>); see commits dd27b6f and aec0772 for
/// the session material as committed to the repo.</para></summary>
public sealed class UniversalCarrierClaim : Claim
{
    private readonly AbsorptionTheoremClaim _absorption;
    private readonly Pi2DyadicLadderClaim _ladder;
    private readonly PolynomialDiscriminantAnchorClaim _discriminant;

    /// <summary>Default convention numerical value for γ₀ in code units: <c>0.05</c>.
    /// "Our c" (Tom 2026-05-12) in code-unit sense. Substrate-invariant role,
    /// convention-dependent value: <c>ValidateAgainstPythonStepFTests</c> covers
    /// γ₀ ∈ {0.025, 0.05, 0.10} and gets identical Q-values across the three.</summary>
    public const double DefaultGammaZero = 0.05;

    /// <summary>The absorption quantum at the default γ₀: <c>2γ₀ = 0.1</c> in code units.
    /// Universal step-size of the Liouvillian eigenvalue grid under uniform Z-dephasing
    /// (per <see cref="AbsorptionTheoremClaim"/>).</summary>
    public double DefaultAbsorptionQuantum => _absorption.AbsorptionQuantum(DefaultGammaZero);

    /// <summary>The default t_peak: <c>1/(4γ₀) = 5.0</c> in code units. The "4" is the
    /// polynomial discriminant <c>a_{−1}</c>; γ₀ is just the units-conversion to physical
    /// time. Per <see cref="PolynomialDiscriminantAnchorClaim"/> + F86 t_peak law.</summary>
    public double DefaultTPeak => 1.0 / (4.0 * DefaultGammaZero);

    public UniversalCarrierClaim(
        AbsorptionTheoremClaim absorption,
        Pi2DyadicLadderClaim ladder,
        PolynomialDiscriminantAnchorClaim discriminant)
        : base("γ₀ as Universal Carrier: protection-constant and observation-substrate carried in a single role (the SR-c analog plus what observers can only see via ratios)",
               Tier.Tier1Derived,
               "compute/RCPsiSquared.Core/Symmetry/AbsorptionTheoremClaim.cs + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs + " +
               "compute/RCPsiSquared.Core/Symmetry/PolynomialDiscriminantAnchorClaim.cs + " +
               "compute/RCPsiSquared.Core/F1/F1PalindromeIdentity.cs + " +
               "hypotheses/PRIMORDIAL_GAMMA_CONSTANT.md + " +
               "hypotheses/PRIMORDIAL_QUBIT.md (§9 operational test from inside) + " +
               "reflections/ON_TWO_TIMES.md (Universal Carrier sticker 2026-05-12)")
    {
        _absorption = absorption ?? throw new ArgumentNullException(nameof(absorption));
        _ladder = ladder ?? throw new ArgumentNullException(nameof(ladder));
        _discriminant = discriminant ?? throw new ArgumentNullException(nameof(discriminant));
    }

    public override string DisplayName =>
        $"γ₀ Universal Carrier (protection-constant + observation-substrate; default γ₀ = {DefaultGammaZero})";

    public override string Summary =>
        $"γ₀ as Universal Carrier; protection-constant (like SR's c) and observation-substrate (what observers see) carried in one role; default code-convention γ₀ = {DefaultGammaZero}; only Q = J/γ₀ and K = γt internally measurable ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("carrier role",
                summary: "γ₀ is the universal-reference rate-parameter mediating between observer-perspectives; neither itself an observer nor what is observed in absolute terms, but the medium between them. Resolves the older tension carrier (TRANSMISSION) vs channel (RESONANCE_NOT_CHANNEL) vs illumination (GAMMA_IS_LIGHT) vs framework constant (PRIMORDIAL_GAMMA_CONSTANT); all four vocabularies describe the same role.");
            yield return new InspectableNode("SR parallel (interpretive Tier 2)",
                summary: "c (SR): only protection-constant; the observation slot is structureless ('events'). γ₀ (R=CΨ²): protection-constant AND observation-substrate. The framing 'Einstein left an empty observation slot, R=CΨ² fills it with γ₀' is a Tier-2 interpretive reading; the typed Tier-1 content is the structural slot identification only.");
            yield return new InspectableNode("Inside-Observability",
                summary: "From inside, γ₀ alone is unmeasurable; only dimensionless ratios are accessible: Q = J/γ₀, K = γt, t_peak = 1/(4γ₀). Were γ₀ directly measurable, the inside-outside split would break.");
            yield return InspectableNode.RealScalar("DefaultGammaZero (code convention, 'our c')", DefaultGammaZero);
            yield return InspectableNode.RealScalar("DefaultAbsorptionQuantum (= 2·γ₀)", DefaultAbsorptionQuantum);
            yield return InspectableNode.RealScalar("DefaultTPeak (= 1/(4γ₀))", DefaultTPeak);
            yield return new InspectableNode("Pi2 foundation anchors",
                summary: $"a₀ = {_ladder.Term(0)} multiplies γ everywhere (2γ in F1, Absorption, F8); a_{{-1}} = {_ladder.Term(-1)} is the polynomial discriminant (4γ in F25/F65/F86 t_peak, per PolynomialDiscriminantAnchorClaim); a₃ = {_ladder.Term(3)} is the mirror partner (¼-boundary).");
            yield return new InspectableNode("substrate invariance",
                summary: "ValidateAgainstPythonStepFTests covers γ₀ ∈ {0.025, 0.05, 0.10}; identical Q-values across the three. Numerical value is convention, structural slot is substrate-invariant. On IBM hardware (T2* ~ 100μs gives γ ~ 10⁴ Hz) different physical value, identical role.");
            yield return new InspectableNode("hardware-confirmed (γ₀ via ratios only)",
                summary: "IBM Q52 absorption ratio 1.03 (3% deviation, AbsorptionTheoremClaim hardware-anchor); IBM Run 3 ¼-boundary 1.9% deviation; Marrakesh F87 trichotomy Δ(soft, truly) = −0.722. γ₀ is never measured directly on hardware; always via Q-ratios.");
        }
    }
}
