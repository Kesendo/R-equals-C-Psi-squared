using RCPsiSquared.Core.F86.Item1Derivation;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>The cusp-approach family (our "family of approach shapes") wired into the typed graph. The
/// partial-entanglement initial state |ψ(α)⟩ = cosα|00⟩ + sinα|11⟩ under Z-dephasing has the coherence
/// CΨ(α,t) = w₀·e^(−4γt) + w₁·e^(−12γt), w₀ = s(1−s²/2)/3, w₁ = s³/6, s = sin2α (Tier-1, verified
/// bit-exact against the Lindblad evolution; implemented in
/// <c>RCPsiSquared.Diagnostics.Foundation.OddHarmonicApproach</c> / <c>ApproachFamilyField</c>, the
/// <c>--axis approach</c> eyepiece). This Claim places that family in the existing typed graph instead of
/// leaving it isolated, via four edges.
///
/// <para>(1) <see cref="UniversalCarrierClaim"/>: every member shares the slowest mode, the carrier rate
/// 4γ₀, and collapses onto it at late time. The "4" is the polynomial discriminant carried inside the
/// Universal Carrier; the carrier the family shares IS γ₀ in its universal-carrier role.</para>
///
/// <para>(2) <see cref="C2BareDoubledPtfClosedForm"/>: a c=2 doubled-PTF kinship. The family is a
/// two-mode structure (the 4γ carrier + a 12γ harmonic, the 3:1 odd-harmonic ratio), the state-space
/// DECAY observable; C2's K_b is the parameter-space SUSCEPTIBILITY observable carrying the same 3:1
/// ratio and the same carrier role. Siblings on the shared structure, NOT a hidden identity: the family's
/// 12γ harmonic is a purity×coherence cross term, while the block / K_b harmonic is the HD=3 sector mode,
/// and the two live at different (intensity vs amplitude) levels. The 3:1 ratio and the carrier are what
/// they share; a viewpoint, not a thing (as with the cusp/EP F95 siblinghood).</para>
///
/// <para>(3) <see cref="TwoReadingsClaim"/>: the closed form (algebra) and the Lindblad trajectory
/// (dynamics) are two readings of the one approach. "The slowing is ours" is this pair: the carrier
/// reading is steady, the apparent slowing lives in the observable.</para>
///
/// <para>(4) <see cref="F25CPsiBellPlusPi2Inheritance"/>: the Bell+ member (s = 1, weights 1/6, 1/6)
/// reproduces F25's two-exponential exactly; the family is F25's one-parameter generalization.</para>
///
/// <para>Tier-1 derived for the closed form and the four edges; the "decay-face sibling of the
/// susceptibility-face" framing is a viewpoint (the shared two-mode carrier structure), not an asserted
/// same-object identity. The cusp's other typed home is <see cref="TransitionBridgeF95SiblingClaim"/>
/// (the F95 angle); these carrier / PTF edges are orthogonal to that one.</para>
///
/// <para>Anchors: <c>compute/RCPsiSquared.Diagnostics/Foundation/OddHarmonicApproach.cs</c> +
/// <c>docs/NAVIGATING_THE_DIMENSIONS.md</c> ("The family of approach shapes") +
/// <c>simulations/approach_family.py</c>.</para></summary>
public sealed class ApproachFamilyCarrierClaim : Claim
{
    /// <summary>Edge 1: the shared carrier 4γ₀ (γ₀ in its universal-carrier role; the "4" the discriminant).</summary>
    public UniversalCarrierClaim Carrier { get; }

    /// <summary>Edge 2: the c=2 doubled-PTF kinship (shared 3:1 ratio + carrier; decay vs susceptibility; kinship not identity).</summary>
    public C2BareDoubledPtfClosedForm C2Ptf { get; }

    /// <summary>Edge 3: algebra (closed form) vs dynamics (Lindblad), two readings of the one approach.</summary>
    public TwoReadingsClaim TwoReadings { get; }

    /// <summary>Edge 4: the Bell+ member (s=1) reproduces F25; the family generalizes F25.</summary>
    public F25CPsiBellPlusPi2Inheritance F25 { get; }

    /// <summary>The carrier rate coefficient: the slowest mode decays at 4γ (the polynomial discriminant × γ).</summary>
    public const double CarrierRateCoefficient = 4.0;

    /// <summary>The harmonic rate coefficient: the fast mode decays at 12γ = 3 × the carrier (the 3:1 ratio).</summary>
    public const double HarmonicRateCoefficient = 12.0;

    /// <summary>The entanglement threshold s = 3/4: the approach crosses ¼ iff s &gt; 3/4 (CΨ(0)=s/3 &gt; ¼).</summary>
    public const double CrossingThresholdS = 0.75;

    public ApproachFamilyCarrierClaim(
        UniversalCarrierClaim carrier,
        C2BareDoubledPtfClosedForm c2Ptf,
        TwoReadingsClaim twoReadings,
        F25CPsiBellPlusPi2Inheritance f25)
        : base("The cusp-approach family CΨ(α,t)=w₀e^(−4γt)+w₁e^(−12γt) (|ψ(α)⟩=cosα|00⟩+sinα|11⟩) shares the " +
               "universal carrier 4γ₀ (every member collapses onto it), is a c=2 two-mode decay sibling of the F86 " +
               "K_b susceptibility (C2 bare-doubled-PTF) on the shared 4γ/12γ 3:1 skeleton (kinship, not identity), " +
               "its algebra and Lindblad dynamics are a two-readings pair, and the Bell+ member is F25",
               Tier.Tier1Derived,
               "compute/RCPsiSquared.Diagnostics/Foundation/OddHarmonicApproach.cs + " +
               "compute/RCPsiSquared.Core/Symmetry/UniversalCarrierClaim.cs + " +
               "compute/RCPsiSquared.Core/F86/Item1Derivation/C2BareDoubledPtfClosedForm.cs + " +
               "compute/RCPsiSquared.Core/Symmetry/TwoReadingsClaim.cs + " +
               "compute/RCPsiSquared.Core/Symmetry/F25CPsiBellPlusPi2Inheritance.cs + " +
               "docs/NAVIGATING_THE_DIMENSIONS.md + simulations/approach_family.py")
    {
        Carrier = carrier ?? throw new ArgumentNullException(nameof(carrier));
        C2Ptf = c2Ptf ?? throw new ArgumentNullException(nameof(c2Ptf));
        TwoReadings = twoReadings ?? throw new ArgumentNullException(nameof(twoReadings));
        F25 = f25 ?? throw new ArgumentNullException(nameof(f25));
    }

    /// <summary>Builds the claim with fresh parent chains (for standalone use; the registry wires via b.Get).
    /// None of the four parents expose a Build()/Shared factory, so the parent chains are constructed
    /// directly (TransitionBridgeF95SiblingClaim.Build() style). A single Pi2DyadicLadderClaim and a single
    /// PolynomialFoundationClaim are shared across the sub-chains that need them (the ladder feeds Absorption,
    /// the discriminant, the Universal Carrier, and F25; the polynomial feeds the discriminant and TwoReadings).</summary>
    public static ApproachFamilyCarrierClaim Build()
    {
        // Shared Pi2-Foundation roots (parameterless ctors).
        var ladder = new Pi2DyadicLadderClaim();
        var polynomial = new PolynomialFoundationClaim();
        var qubit = new QubitDimensionalAnchorClaim();
        var quarter = new QuarterAsBilinearMaxvalClaim();

        // Edge 1: UniversalCarrierClaim(AbsorptionTheoremClaim, Pi2DyadicLadderClaim, PolynomialDiscriminantAnchorClaim).
        var absorption = new AbsorptionTheoremClaim(ladder);
        var discriminant = new PolynomialDiscriminantAnchorClaim(polynomial, qubit, ladder);
        var carrier = new UniversalCarrierClaim(absorption, ladder, discriminant);

        // Edge 2: C2BareDoubledPtfClosedForm (parameterless ctor).
        var c2Ptf = new C2BareDoubledPtfClosedForm();

        // Edge 3: TwoReadingsClaim(PolynomialFoundationClaim).
        var twoReadings = new TwoReadingsClaim(polynomial);

        // Edge 4: F25CPsiBellPlusPi2Inheritance(Pi2DyadicLadderClaim, QuarterAsBilinearMaxvalClaim).
        var f25 = new F25CPsiBellPlusPi2Inheritance(ladder, quarter);

        return new ApproachFamilyCarrierClaim(carrier, c2Ptf, twoReadings, f25);
    }

    /// <summary>Shared singleton; the claim is a structural synthesis, block-independent.</summary>
    public static ApproachFamilyCarrierClaim Shared { get; } = Build();

    public override string DisplayName =>
        "Approach family carrier wiring (cusp-approach as a carried c=2 decay reading; Bell+ = F25)";

    public override string Summary =>
        "the cusp-approach family CΨ(α,t)=w₀e^(−4γt)+w₁e^(−12γt) shares the universal carrier 4γ₀, is a c=2 " +
        "two-mode decay sibling of the F86 K_b susceptibility (kinship not identity), its algebra/dynamics are a " +
        $"two-readings pair, and the Bell+ member (s=1) is F25; crosses ¼ iff s>3/4 ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("the shared carrier (edge 1)",
                summary: "every member collapses onto the slowest mode 4γ₀; the '4' is the polynomial discriminant carried inside the Universal Carrier; the carrier the family shares is γ₀ in its universal-carrier role.");
            yield return new InspectableNode("the c=2 doubled-PTF kinship (edge 2, a viewpoint not an identity)",
                summary: "the family (decay observable) and C2's K_b (susceptibility observable) share the 4γ/12γ 3:1 odd-harmonic ratio and the carrier role; the 12γ harmonic is a purity×coherence cross term here, an HD=3 sector mode there, at different (intensity vs amplitude) levels. Siblings on the shared structure, not the same object.");
            yield return new InspectableNode("the two readings (edge 3)",
                summary: "the closed form (algebra) and the Lindblad trajectory (dynamics) are two readings of the one approach; 'the slowing is ours' is this pair, the carrier reading steady, the apparent slowing in the observable.");
            yield return new InspectableNode("the Bell+ member = F25 (edge 4)",
                summary: "s=1 gives weights (1/6, 1/6), so CΨ(1,t)=(1/6)e^(−4γt)+(1/6)e^(−12γt)=f(1+f²)/6, exactly F25; the family is F25's one-parameter (entanglement) generalization.");
            yield return InspectableNode.RealScalar("carrier rate coefficient (×γ)", CarrierRateCoefficient);
            yield return InspectableNode.RealScalar("harmonic rate coefficient (×γ, = 3×carrier)", HarmonicRateCoefficient);
            yield return InspectableNode.RealScalar("entanglement crossing threshold s = sin2α", CrossingThresholdS);
            yield return Carrier;
            yield return C2Ptf;
            yield return TwoReadings;
            yield return F25;
        }
    }
}
