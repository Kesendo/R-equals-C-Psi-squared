using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>UniversalCarrierClaim (Tier 1 derived; named 2026-05-12): γ₀ is the
/// Universal Carrier — the universal-reference rate-parameter that mediates
/// between observer-perspectives, neither itself an observer nor what is observed
/// in absolute terms, but the medium between them. It plays the structural role of
/// the speed of light c in special relativity, with the addition (Tom 2026-05-12)
/// that γ₀ is also what observers see (only via dimensionless ratios; γ₀ alone
/// is invisible from inside per <c>hypotheses/PRIMORDIAL_QUBIT.md</c> §9
/// Inside-Observability).
///
/// <para><b>Schutz-Konstante UND Beobachtetes in einer Rolle:</b> Einstein lieferte
/// "Zeit ist relativ zum Beobachter" plus c als Schutz-Konstante, ließ aber implizit
/// was eigentlich beobachtet wird (in SR-Praxis: "Events in der Raumzeit", strukturlos).
/// R=CΨ² füllt die Lücke: γ₀ trägt beide Rollen. Würde γ₀ direkt messbar sein
/// (separater Beobachtungs-Slot), wäre Inside-Observability verletzt; dass γ₀ zugleich
/// Konstante und Beobachtetes ist, ist die Bedingung dafür dass die Inside-Outside-
/// Trennung (per <c>docs/proofs/INCOMPLETENESS_PROOF.md</c>) algebraisch konsistent
/// bleibt.</para>
///
/// <para><b>Synthesis claim — introduces no new content:</b> makes the already-typed
/// structural role of γ₀ across the framework explicit. Takes dependencies on:</para>
/// <list type="bullet">
///   <item><see cref="AbsorptionTheoremClaim"/>: γ₀ as universal absorption quantum
///         <c>2γ₀ = a₀·γ₀</c>; spectrum quantized in 2γ₀ steps; F1 palindrome center
///         <c>−Σγ_l</c> (pool of all per-site carrier-shares).</item>
///   <item><see cref="Pi2DyadicLadderClaim"/>: the "2" multiplying γ everywhere is
///         <c>a₀</c> on the dyadic ladder; the "4" in <c>t_peak = 1/(4γ₀)</c> is
///         <c>a_{−1}</c>.</item>
///   <item><see cref="PolynomialDiscriminantAnchorClaim"/>: the polynomial discriminant
///         of <c>d²−2d=0</c> is <c>4 = a_{−1}</c>; the <c>1/4</c> mirror partner sits
///         at <c>a_3</c>; closure <c>4·(1/4) = 1</c>.</item>
/// </list>
///
/// <para><b>"Unsere c ist 0.05" (Tom 2026-05-12):</b> the default code-convention
/// value <see cref="DefaultGammaZero"/> = 0.05 plays the structural role of c in SR.
/// The numerical value is convention; the role is substrate-invariant.
/// <c>F86KnowledgeBase</c> tests γ₀ ∈ {0.025, 0.05, 0.10} for γ₀-invariance — different
/// numerical choices, identical Q-values. On IBM hardware the physical scale is
/// different (γ ~ 1/T2* ~ 10⁴ Hz for T2* ~ 100μs); different substrate, identical
/// structural slot.</para>
///
/// <para><b>Inside-Observability witnesses (γ₀ measured only via ratios on real
/// hardware):</b> IBM Q52 absorption ratio 1.03 (3% deviation, Apr 4 2026), IBM
/// Run 3 ¼-boundary crossing 1.9% deviation (Mar 18 2026), Marrakesh F87 trichotomy
/// Δ(soft − truly) = −0.722 (Apr 26 2026). All three confirm: from inside, only
/// dimensionless ratios involving γ₀ are accessible; γ₀ itself never measured directly.</para>
///
/// <para>Anchors: the four-Claim dependency lineage above plus
/// <c>hypotheses/PRIMORDIAL_GAMMA_CONSTANT.md</c> ("γ at the root is the framework's
/// own c") + <c>hypotheses/PRIMORDIAL_QUBIT.md</c> §9 (Inside-Observability proof) +
/// <c>reflections/ON_TWO_TIMES.md</c> (Universal Carrier sticker 2026-05-12) +
/// <c>memory/project_einstein_gap_filled_by_gamma0.md</c> (Einstein-Lücke
/// vervollständigt) + <c>memory/project_time_is_gamma0_observer.md</c> v3
/// (perspektiv-relative Zeit) + <c>memory/project_ptf_is_enacted_not_defined.md</c>
/// (PTF wird getan, nicht definiert).</para></summary>
public sealed class UniversalCarrierClaim : Claim
{
    private readonly AbsorptionTheoremClaim _absorption;
    private readonly Pi2DyadicLadderClaim _ladder;
    private readonly PolynomialDiscriminantAnchorClaim _discriminant;

    /// <summary>Default convention numerical value for γ₀ in code units: <c>0.05</c>.
    /// "Unsere c" (Tom 2026-05-12) in code-unit sense. Substrate-invariant role,
    /// convention-dependent value: <c>F86KnowledgeBase</c> tests γ₀ ∈ {0.025, 0.05, 0.10}
    /// for γ₀-invariance and gets identical Q-values across the three.</summary>
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
        : base("γ₀ als Universal Carrier: Schutz-Konstante UND Beobachtetes in einer Rolle (Einstein-Lücke vervollständigt)",
               Tier.Tier1Derived,
               "compute/RCPsiSquared.Core/Symmetry/AbsorptionTheoremClaim.cs + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs + " +
               "compute/RCPsiSquared.Core/Symmetry/PolynomialDiscriminantAnchorClaim.cs + " +
               "compute/RCPsiSquared.Core/F1/F1PalindromeIdentity.cs + " +
               "hypotheses/PRIMORDIAL_GAMMA_CONSTANT.md + " +
               "hypotheses/PRIMORDIAL_QUBIT.md (§9 Inside-Observability) + " +
               "reflections/ON_TWO_TIMES.md (Universal Carrier sticker 2026-05-12)")
    {
        _absorption = absorption ?? throw new ArgumentNullException(nameof(absorption));
        _ladder = ladder ?? throw new ArgumentNullException(nameof(ladder));
        _discriminant = discriminant ?? throw new ArgumentNullException(nameof(discriminant));
    }

    public override string DisplayName =>
        $"γ₀ Universal Carrier (Schutz-Konstante UND Beobachtetes; default γ₀ = {DefaultGammaZero})";

    public override string Summary =>
        $"γ₀ als Universal Carrier; Schutz-Konstante (wie SR's c) UND was beobachtet wird (Einstein-Lücke vervollständigt); default code-convention γ₀ = {DefaultGammaZero}; nur Q = J/γ₀ und K = γt von innen messbar ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("carrier role",
                summary: "γ₀ ist das universal-reference rate-parameter das zwischen Beobachter-Perspektiven vermittelt; weder selbst Beobachter noch Beobachtetes-im-absoluten-Sinn, sondern das Medium dazwischen. Resolves the older tension carrier (TRANSMISSION) vs channel (RESONANCE_NOT_CHANNEL) vs illumination (GAMMA_IS_LIGHT) vs framework constant (PRIMORDIAL_GAMMA_CONSTANT) — alle vier Vokabulare beschreiben dieselbe Rolle.");
            yield return new InspectableNode("Einstein parallel (SR ↔ R=CΨ²)",
                summary: "c (SR): nur Schutz-Konstante; was beobachtet wird = strukturlose 'Events'. γ₀ (R=CΨ²): Schutz-Konstante UND was beobachtet wird (Einstein ließ den Beobachtungs-Slot strukturlos; R=CΨ² füllt ihn). Doppelrolle ist Bedingung für Inside-Observability-Konsistenz.");
            yield return new InspectableNode("Inside-Observability",
                summary: "von innen ist γ₀ alleine unmessbar; nur dimensionslose Verhältnisse zugänglich: Q = J/γ₀, K = γt, t_peak = 1/(4γ₀). Würde γ₀ direkt messbar sein, wäre Inside-Outside-Trennung verletzt.");
            yield return InspectableNode.RealScalar("DefaultGammaZero (code convention, 'unsere c')", DefaultGammaZero);
            yield return InspectableNode.RealScalar("DefaultAbsorptionQuantum (= 2·γ₀)", DefaultAbsorptionQuantum);
            yield return InspectableNode.RealScalar("DefaultTPeak (= 1/(4γ₀))", DefaultTPeak);
            yield return new InspectableNode("Pi2-Foundation anchors",
                summary: $"a₀ = {_ladder.Term(0)} multipliziert γ überall (2γ in F1, F-Absorption, F8); a_{{-1}} = {_ladder.Term(-1)} ist die Polynom-Diskriminante (4γ in F25/F65/F86 t_peak per PolynomialDiscriminantAnchorClaim); a₃ = {_ladder.Term(3)} ist Mirror-Partner (¼-boundary).");
            yield return new InspectableNode("substrate invariance",
                summary: "F86KnowledgeBase testet γ₀ ∈ {0.025, 0.05, 0.10} → identische Q-Werte; numerischer Wert ist Konvention, struktureller Slot ist substrate-invariant. Auf IBM Hardware (T2* ~ 100μs → γ ~ 10⁴ Hz) anderer physikalischer Wert, identische Rolle.");
            yield return new InspectableNode("hardware-confirmed (γ₀ via ratios only)",
                summary: "IBM Q52 absorption ratio 1.03 (3% deviation, AbsorptionTheoremClaim hardware-anchor) + IBM Run 3 ¼-boundary 1.9% deviation + Marrakesh F87 trichotomy Δ(soft−truly) = −0.722. γ₀ wird nie direkt gemessen; immer über Q-Verhältnisse.");
            yield return new InspectableNode("memory anchors",
                summary: "memory/project_einstein_gap_filled_by_gamma0.md (2026-05-12 Tom) + memory/project_time_is_gamma0_observer.md v3 (perspektiv-relativ; γ₀, wir, vermutlich Bath haben jeweils eigene Zeit) + memory/project_ptf_is_enacted_not_defined.md (PTF wird getan, nicht definiert; Iteration ist der Wechsel)");
        }
    }
}
