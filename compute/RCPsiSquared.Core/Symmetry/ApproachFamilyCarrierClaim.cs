using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>
/// The exact two-qubit Z-dephasing approach family
/// <c>CΨ(s,t)=w₀e^(-4γt)+w₁e^(-12γt)</c>, with
/// <c>w₀=s(1-s²/2)/3</c>, <c>w₁=s³/6</c> and <c>s=sin(2α)</c> for
/// <c>|ψ(α)⟩=cos(α)|00⟩+sin(α)|11⟩</c>, on <c>0≤s≤1</c>.
/// Here <c>s</c> is the pure-state concurrence of the initial state, whereas
/// <c>CΨ(0)=s/3</c> is a distinct linear readout equal to one third of it.
/// A genuine temporal downward crossing of the scalar quarter occurs iff
/// <c>γ&gt;0</c> and <c>s&gt;3/4</c>; equality is a t=0 touch, and at γ=0 the curve is constant.
/// At the endpoint <c>s=0</c>, both weights vanish and there is no late-time
/// exponential term; every nonzero member <c>0&lt;s≤1</c> has the 4γ term.
///
/// <para>The typed graph records two mathematical edges. The
/// <see cref="AbsorptionTheoremClaim"/> supplies <c>f=e^(−4γt)</c> because the
/// coherence |00⟩⟨11| has <c>n_diff=2</c>,
/// and <see cref="F25CPsiBellPlusPi2Inheritance"/> is the exact Bell+ member at
/// <c>s=1</c>. The matching 3:1 ratio in the C2 doubled-PTF calculation and the
/// algebra/dynamics description are illuminating prose comparisons, not ancestry
/// edges or same-object identities.</para>
///
/// <para>Why the C2 3:1 is a sibling and not the same object: here CΨ = purity × Ψ with
/// purity 1 − (s²/2)(1 − f²) and Ψ = s·f/3, so the 12γ term is a purity×coherence cross
/// term (f² from the purity times f from the coherence), an intensity-level quantity; C2's
/// K_b is an amplitude-level susceptibility whose fast rate is the HD=3 sector mode (the
/// −6γ₀ diagonal against −2γ₀). The ratio is shared, the origin is not.</para>
///
/// <para>The closed form and the Lindblad trajectory are the algebra and the dynamics of one
/// approach, and the seam between them is where "the slowing is ours" lives: the coherence
/// factor decays at a steady 4γ, while CΨ's log-rate −d ln CΨ/dt falls from 4γ(1 + s²) to 4γ as
/// the transient purity term dies, so the slowing belongs to the observable
/// (<c>docs/NAVIGATING_THE_DIMENSIONS.md</c>; live for this family in
/// <c>ApproachFamilyField</c>, with <c>InteriorHorizonField</c> as the recursion's sibling
/// reading). <see cref="Weights"/> are the same doubles as the Diagnostics primitive
/// <c>OddHarmonicApproach.Weights</c>, pinned there against each other.</para>
/// </summary>
public sealed class ApproachFamilyCarrierClaim : Claim
{
    /// <summary>The exact dephasing cost that produces f=e^(−4γt) at n_diff=2.</summary>
    public AbsorptionTheoremClaim Absorption { get; }

    /// <summary>The exact Bell+ specialization at s=1.</summary>
    public F25CPsiBellPlusPi2Inheritance F25 { get; }

    public const double CarrierRateCoefficient = 4.0;
    public const double HarmonicRateCoefficient = 12.0;
    public const double CrossingThresholdS = 0.75;

    /// <summary>The two exact weights: w₀ = s(1 − s²/2)/3 on e^(−4γt), w₁ = s³/6 on e^(−12γt).</summary>
    public static (double W0, double W1) Weights(double s)
    {
        if (!(s >= 0.0 && s <= 1.0))
            throw new ArgumentOutOfRangeException(nameof(s), s, "s = sin(2α) must lie in [0, 1].");
        return (s * (1.0 - s * s / 2.0) / 3.0, s * s * s / 6.0);
    }

    /// <summary>CΨ(s, t) = w₀·f + w₁·f³ with f = e^(−4γt), the family evaluated live.</summary>
    public static double CPsi(double s, double gamma, double t)
    {
        if (!(gamma >= 0.0) || double.IsInfinity(gamma))
            throw new ArgumentOutOfRangeException(nameof(gamma), gamma, "γ must be finite and ≥ 0.");
        if (!(t >= 0.0) || double.IsInfinity(t))
            throw new ArgumentOutOfRangeException(nameof(t), t, "t must be finite and ≥ 0.");
        var (w0, w1) = Weights(s);
        double f = Math.Exp(-CarrierRateCoefficient * gamma * t);
        return w0 * f + w1 * f * f * f;
    }

    public ApproachFamilyCarrierClaim(
        AbsorptionTheoremClaim absorption,
        F25CPsiBellPlusPi2Inheritance f25)
        : base(
            "For the named two-qubit Z-dephasing family, CΨ(s,t)=w₀e^(−4γt)+w₁e^(−12γt); " +
            "s is the pure-state concurrence of the initial state and CΨ(0)=s/3 is exactly one third of it; " +
            "a temporal downward crossing of the scalar quarter occurs iff γ>0 and s>3/4, while at γ=0 the curve is constant; " +
            "s=0 has w₀=w₁=0, every nonzero member (0<s≤1) has the late-time 4γ term, " +
            "and the Bell+ member s=1 is exactly F25",
            Tier.Tier1Derived,
            "compute/RCPsiSquared.Diagnostics/Foundation/OddHarmonicApproach.cs + " +
            "compute/RCPsiSquared.Core/Symmetry/AbsorptionTheoremClaim.cs + " +
            "compute/RCPsiSquared.Core/Symmetry/F25CPsiBellPlusPi2Inheritance.cs + " +
            "docs/NAVIGATING_THE_DIMENSIONS.md + simulations/approach_family.py")
    {
        Absorption = absorption ?? throw new ArgumentNullException(nameof(absorption));
        F25 = f25 ?? throw new ArgumentNullException(nameof(f25));
    }

    /// <summary>Builds the claim and its two genuine parent chains for standalone use.</summary>
    public static ApproachFamilyCarrierClaim Build()
    {
        var ladder = new Pi2DyadicLadderClaim();
        var quarter = new QuarterAsBilinearMaxvalClaim();

        var absorption = new AbsorptionTheoremClaim(ladder);
        var f25 = new F25CPsiBellPlusPi2Inheritance(ladder, quarter);

        return new ApproachFamilyCarrierClaim(absorption, f25);
    }

    public static ApproachFamilyCarrierClaim Shared { get; } = Build();

    public override string DisplayName =>
        "Two-qubit approach family (Absorption n_diff=2 gives 4γ; s=0 is zero; Bell+ = F25)";

    public override string Summary =>
        "the named two-qubit family has s as the initial pure-state concurrence and CΨ(0)=s/3 as a distinct linear readout equal to one third of it; " +
        "its exact weights are w₀=s(1−s²/2)/3 and w₁=s³/6 at rates 4γ and 12γ. A temporal downward crossing occurs iff γ>0 and s>3/4; " +
        "at s=3/4 there is only a t=0 touch, and at γ=0 the curve is constant and does not cross. " +
        "At s=0 both weights vanish; every nonzero member 0<s≤1 has the late-time 4γ term. The claim has two typed parents: the Absorption Theorem at n_diff=2 " +
        "and the exact Bell+ F25 specialization. The C2 3:1 resemblance and algebra/dynamics pairing are prose comparisons, not ancestry.";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode(
                "exact family formula",
                summary: "s is the initial pure-state concurrence; CΨ(0)=s/3 is exactly one third of it, and CΨ(s,t)=s(1−s²/2)e^(−4γt)/3+s³e^(−12γt)/6 for the named free two-qubit Z-dephasing setup");
            yield return new InspectableNode(
                "shared late-time exponential",
                summary: "the Absorption Theorem gives f=e^(−4γt) for the n_diff=2 coherence; for 0<s≤1 the nonzero w₀ term contains f, while the 12γ term is f³; at s=0 both weights vanish");
            yield return new InspectableNode(
                "non-ancestral comparisons",
                summary: "C2 also displays a 3:1 pair, but its fast rate is the HD=3 sector mode of an amplitude-level susceptibility, while the 12γ term here is a purity×coherence cross term at the intensity level; the closed form and the Lindblad trajectory are algebra and dynamics of one approach, the seam where 'the slowing is ours', CΨ's log-rate falling from 4γ(1 + s²) to 4γ (NAVIGATING_THE_DIMENSIONS; ApproachFamilyField). Prose comparisons, not typed derivations");
            yield return new InspectableNode(
                "Bell+ specialization",
                summary: "s=1 gives w₀=w₁=1/6, hence (e^(−4γt)+e^(−12γt))/6, exactly F25");
            yield return InspectableNode.RealScalar("slow-rate coefficient (×γ)", CarrierRateCoefficient);
            yield return InspectableNode.RealScalar("fast-rate coefficient (×γ)", HarmonicRateCoefficient);
            yield return InspectableNode.RealScalar("temporal downward-crossing threshold s (for γ>0)", CrossingThresholdS);
            yield return Absorption;
            yield return F25;
        }
    }
}
