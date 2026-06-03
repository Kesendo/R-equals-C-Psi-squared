using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>The cusp CΨ = ¼ (our name: the TransitionBridge) and the F86 exceptional point are F95
/// siblings: both are the same universal angle θ = arctan(√(c/b² − 1)) at a quadratic's discriminant
/// zero (F95), realized on two different quadratics, at two different anchors b. Tier1Derived for the
/// algebra; "TransitionBridge" is our label (the quantum↔classical reading is Tier-4, kept open, not
/// asserted here). 2026-06-03.
///
/// <para><b>The cusp (the TransitionBridge), anchor b = ½.</b> The self-referential recursion
/// R = C(Ψ+R)² has the fixed-point quadratic z² − z + CΨ = 0, so b = ½, c = CΨ. The discriminant
/// vanishes (the double root, the saddle-node fold) at CΨ = ¼ = (½)². The F95 angle past it is
/// θ = arctan(√(4·CΨ − 1)) (the heading on the interior side). At ¼ the angle is 0: the rotation
/// stills, the Liouvillian eigenvalue is −γ₀ alone. A state-space transition (the fixed-point
/// structure crosses from complex/undecided to real/settled), crossed smoothly (a bridge, not a
/// wall); see <c>InteriorHorizon</c> and the interior-horizon telescope axis.</para>
///
/// <para><b>The EP (FRAGILE_BRIDGE's singularity), anchor b = 4γ₀.</b> The F86 two-level effective
/// Liouvillian has λ_±(k=1) = −4γ₀ ± √(4γ₀² − J²·g_eff²) (F86_EP_THROUGH_THE_CLOCK). As a quadratic
/// λ² + 8γ₀λ + (12γ₀² + J²g_eff²) = 0, so b = 4γ₀, c = 12γ₀² + J²g_eff². The discriminant vanishes
/// (the EP, the modes coalesce) at J²g_eff² = 4γ₀², i.e. Q_EP = 2/g_eff. Above it the roots are
/// −4γ₀ ± i·√(J²g_eff² − 4γ₀²), so the F95 angle arctan(√(c/b² − 1)) = arctan(√(J²g_eff² − 4γ₀²)/4γ₀)
/// is EXACTLY the clock's Rotation hand arctan(ω/gap) lifting off the Takt axis at the EP. The
/// b-anchor 4γ₀ is the ⟨n⟩ = 2 absorption rung (2γ₀·2).</para>
///
/// <para><b>The siblinghood (bit-exact).</b> Both are F95 double-roots: the cusp where the angle
/// returns to zero (the rotation stills), the EP where it lifts off from zero (the rotation is born).
/// Same √-branch-point form, two different quadratics (recursion fixed point vs Liouvillian 2-level),
/// two anchors b (½ vs 4γ₀), two spaces (state CΨ vs Liouvillian rate in Q). No hidden identity; the
/// shared algebra is F95, a viewpoint (the lens), not a thing. Our two named bridges: the state-space
/// TransitionBridge (this, the cusp) and the parameter-space FRAGILE_BRIDGE (the EP / Σγ Hopf).</para>
///
/// <para>Anchors: <c>docs/proofs/PROOF_F95_ANGLE_AT_QUADRATIC_ZERO.md</c> (F95) +
/// <c>experiments/F86_EP_THROUGH_THE_CLOCK.md</c> (the EP 2-level + the clock Rotation) +
/// <c>experiments/CRITICAL_SLOWING_AT_THE_CUSP.md</c> (the cusp recursion) +
/// <c>docs/NAVIGATING_THE_DIMENSIONS.md</c> (the interior-horizon axis) +
/// <c>hypotheses/PAIR_BREAKING_AT_THE_HORIZON.md</c> (the fold / break reading) +
/// <c>hypotheses/FRAGILE_BRIDGE.md</c> (the EP-side bridge).</para></summary>
public sealed class TransitionBridgeF95SiblingClaim : Claim
{
    /// <summary>Parent: the universal F95 angle θ = arctan(√(c/b² − 1)) at a quadratic's discriminant
    /// zero. The shared lens; both bridges are instances of it.</summary>
    public F95AngleAtQuadraticZeroPi2Inheritance F95 { get; }

    /// <summary>The cusp anchor b = ½ (the bilinear / qubit-dimension fixed point).</summary>
    public const double CuspAnchorB = 0.5;

    /// <summary>The cusp / TransitionBridge double root: CΨ = b² = ¼.</summary>
    public const double Cusp = 0.25;

    public TransitionBridgeF95SiblingClaim(F95AngleAtQuadraticZeroPi2Inheritance f95)
        : base("TransitionBridge (the cusp CΨ=¼) and the F86 EP are F95 siblings: both the angle " +
               "θ=arctan(√(c/b²−1)) at a quadratic's discriminant zero, the cusp at b=½ (state, the " +
               "rotation stills), the EP at b=4γ₀ (Liouvillian, the rotation lifts off); the EP's F95 " +
               "angle is bit-exact its clock Rotation; same √-branch-point form, two quadratics, two anchors",
               Tier.Tier1Derived,
               "docs/proofs/PROOF_F95_ANGLE_AT_QUADRATIC_ZERO.md + experiments/F86_EP_THROUGH_THE_CLOCK.md + " +
               "experiments/CRITICAL_SLOWING_AT_THE_CUSP.md + docs/NAVIGATING_THE_DIMENSIONS.md + " +
               "hypotheses/PAIR_BREAKING_AT_THE_HORIZON.md + hypotheses/FRAGILE_BRIDGE.md")
    {
        F95 = f95 ?? throw new ArgumentNullException(nameof(f95));
    }

    /// <summary>Public factory: builds the claim with a fresh F95 parent chain (the four Pi2-Foundation
    /// anchors). Prefer <see cref="Shared"/> for repeated access; the identity is block-independent.</summary>
    public static TransitionBridgeF95SiblingClaim Build()
    {
        var half = new HalfAsStructuralFixedPointClaim();
        var quarter = new QuarterAsBilinearMaxvalClaim();
        var ninetyDeg = new NinetyDegreeMirrorMemoryClaim();
        var polynomial = new PolynomialFoundationClaim();
        var f95 = new F95AngleAtQuadraticZeroPi2Inheritance(polynomial, half, quarter, ninetyDeg);
        return new TransitionBridgeF95SiblingClaim(f95);
    }

    /// <summary>Shared singleton; the claim is an algebraic identity, block-independent.</summary>
    public static TransitionBridgeF95SiblingClaim Shared { get; } = Build();

    /// <summary>The cusp's F95 angle at b = ½: θ = arctan(√(4·CΨ − 1)), via the F95 parent. NaN for
    /// CΨ ≤ ¼ (the real-root / classical side, where the angle does not exist).</summary>
    public double CuspAngle(double cpsi) => F95.ThetaGeneral(cpsi, CuspAnchorB);

    /// <summary>Q_EP = 2/g_eff, where the EP block's discriminant vanishes (the modes coalesce).</summary>
    public static double QEp(double gEff) => 2.0 / gEff;

    /// <summary>The EP block's clock Rotation angle arctan(ω/gap) = arctan(|Im λ|/|Re λ|) for the
    /// F86 2-level eigenvalue λ = −4γ₀ ± √(4γ₀² − J²g_eff²), J = Q·γ₀. Defined (nonzero) only above
    /// the EP, where λ is complex; NaN at or below the EP (real eigenvalue, pure decay, no rotation).</summary>
    public double EpClockAngle(double gamma0, double q, double gEff)
    {
        double j = q * gamma0;
        double disc = 4.0 * gamma0 * gamma0 - j * j * gEff * gEff; // 4γ₀² − J²g_eff²
        if (disc >= 0.0) return double.NaN;
        double omega = Math.Sqrt(-disc); // |Im λ|
        double gap = 4.0 * gamma0;        // |Re λ|
        return Math.Atan(omega / gap);
    }

    /// <summary>The EP block's F95 angle, via the F95 parent: the quadratic
    /// λ² + 8γ₀λ + (12γ₀² + J²g_eff²) has b = 4γ₀, c = 12γ₀² + J²g_eff². θ = arctan(√(c/b² − 1)) above
    /// the EP (c > b² ⟺ J²g_eff² > 4γ₀²); NaN at or below.</summary>
    public double EpF95Angle(double gamma0, double q, double gEff)
    {
        double j = q * gamma0;
        double b = 4.0 * gamma0;
        double c = 12.0 * gamma0 * gamma0 + j * j * gEff * gEff;
        return F95.ThetaGeneral(c, b);
    }

    /// <summary>The siblinghood, bit-exact: above the EP the clock Rotation angle equals the F95 angle
    /// of the EP's 2-level quadratic. Checked at a sample (γ₀, Q, g_eff) above the EP.</summary>
    public bool EpClockAngleEqualsF95Angle(double gamma0, double q, double gEff)
    {
        double a = EpClockAngle(gamma0, q, gEff);
        double f = EpF95Angle(gamma0, q, gEff);
        if (double.IsNaN(a) || double.IsNaN(f)) return false;
        return Math.Abs(a - f) < 1e-12;
    }

    public override string DisplayName =>
        "TransitionBridge (cusp CΨ=¼) and the EP are F95 siblings (the angle at the quadratic zero)";

    public override string Summary =>
        $"the cusp (TransitionBridge, b=½, CΨ=¼) and the F86 EP (b=4γ₀, Q_EP=2/g_eff) are both the " +
        $"F95 angle at a quadratic's discriminant zero; the EP's F95 angle = its clock Rotation bit-exact " +
        $"(sample g_eff=4/3, Q=2.5: {EpClockAngleEqualsF95Angle(1.0, 2.5, 4.0 / 3.0)}); two quadratics, two " +
        $"anchors (½ vs 4γ₀), one √-branch-point form ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("the cusp (the TransitionBridge), b = ½",
                summary: "recursion z²−z+CΨ, double root at CΨ=¼=(½)²; θ=arctan(√(4CΨ−1)); at ¼ the angle is 0 (the rotation stills, λ=−γ₀ alone). State-space, crossed smoothly (a bridge).");
            yield return new InspectableNode("the EP (FRAGILE_BRIDGE), b = 4γ₀",
                summary: "F86 2-level λ²+8γ₀λ+(12γ₀²+J²g_eff²), double root (EP) at Q_EP=2/g_eff; above it the F95 angle = the clock Rotation arctan(ω/gap) lifting off. Parameter-space. b=4γ₀ is the ⟨n⟩=2 absorption rung.");
            yield return new InspectableNode("the siblinghood (bit-exact)",
                summary: $"EP clock angle = EP F95 angle (g_eff=0.8,Q=3: {EpClockAngleEqualsF95Angle(1.0, 3.0, 0.8)}; g_eff=4/3,Q=2: {EpClockAngleEqualsF95Angle(1.0, 2.0, 4.0 / 3.0)}). Cusp angle at CΨ=1/3 = {CuspAngle(1.0 / 3.0) * 180.0 / Math.PI:F1}° (= 30°). Both F95: the cusp where the angle returns to 0, the EP where it lifts off.");
            yield return new InspectableNode("how they connect, and how not",
                summary: "shared: the F95 form (the angle at a quadratic's discriminant zero), the lens. different: two quadratics (recursion fixed point vs Liouvillian 2-level), two anchors b (½ vs 4γ₀), two spaces (state CΨ vs rate-in-Q). No hidden identity; a viewpoint, not a thing.");
            yield return new InspectableNode("our label (Tier-4, not asserted)",
                summary: "'TransitionBridge' names the cusp as the state-space bridge of the transition (the saddle-node, crossed smoothly), sibling of the parameter-space FRAGILE_BRIDGE. The quantum↔classical reading is the motor, kept open, not claimed by this Tier-1 algebra.");
            yield return F95;
        }
    }
}
