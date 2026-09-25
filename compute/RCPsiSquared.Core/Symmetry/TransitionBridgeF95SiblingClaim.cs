using System.Numerics;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>Two distinct positive-b quadratics can be evaluated through F95's coordinate
/// θ = arctan(√(c/b² − 1)): the dimensionless recurrence at c = ¼, and the specified F86
/// two-level decay polynomial at its exceptional point. They have different variables, anchors,
/// and physical status. Tier1Derived covers only the two algebraic substitutions and their
/// numerical angle equality; the name TransitionBridge is the page's image, not a physical
/// classification.
///
/// <para><b>The cusp (the TransitionBridge), anchor b = ½.</b> The self-referential recursion
/// R = C(Ψ+R)² has the fixed-point quadratic z² − z + CΨ = 0, so b = ½, c = CΨ. The discriminant
/// vanishes (the double root, the saddle-node fold) at CΨ = ¼ = (½)². The F95 angle past it is
/// θ = arctan(√(4·CΨ − 1)) on the complex-root side. At ¼ the two recurrence roots meet and the
/// angle is zero. This statement is about the dimensionless recurrence variable <c>z_rec</c>, not
/// about a Liouvillian eigenvalue or a universal physical transition.</para>
///
/// <para><b>The genuine F86 toy 2×2 EP, anchor b = 4γ₀.</b> The F86 two-level effective
/// Liouvillian has λ_±(k=1) = −4γ₀ ± √(4γ₀² − J²·g_eff²) (F86_EP_THROUGH_THE_CLOCK). In the positive
/// decay variable z=−λ, its quadratic is z² − 8γ₀z + (12γ₀² + J²g_eff²) = 0, so b = 4γ₀ &gt; 0 and
/// c = 12γ₀² + J²g_eff². The anchor is an absorption rung: at J = 0 the two modes decay at 2γ₀ and
/// 6γ₀, the rungs ⟨n_XY⟩ = 1 and 3 of α = 2γ₀⟨n_XY⟩, and b = 4γ₀ is their midpoint, the rung
/// ⟨n_XY⟩ = 2. The discriminant vanishes
/// (the EP, the modes coalesce) at J²g_eff² = 4γ₀², i.e. Q_EP = 2/g_eff. Above it the roots are
/// 4γ₀ ± i·√(J²g_eff² − 4γ₀²) in <c>z_decay=-lambda</c>, so the F95 angle
/// arctan(√(c/b² − 1)) = arctan(√(J²g_eff² − 4γ₀²)/4γ₀) is the same real-arithmetic
/// expression as the toy clock's Rotation hand arctan(ω/gap).</para>
///
/// <para><b>The comparison.</b> These are two distinct positive-b quadratic applications:
/// the recurrence angle returns to zero at its discriminant boundary, whereas the F86 angle lifts
/// off from zero at its exceptional point. The variables are <c>z_rec</c> and
/// <c>z_decay=-lambda</c>; the anchors are ½ and 4γ₀. F95 supplies only the coordinate used to
/// compare them, not a shared physical object. The separate Σγ=0 FRAGILE_BRIDGE system (two
/// qubits per chain) has, at generic couplings, an exceptional point of its own, a defective EP2 on the real γ axis at its
/// threshold γ_crit, where, in the first popcount block of its Liouvillian to go unstable, two
/// eigenvalues on the imaginary axis meet at zero decay (at five exact couplings its threshold is
/// zero instead, with no EP). It is not an instance of this quadratic, whose double root sits at
/// the positive decay z_decay = b = 4γ₀.</para>
///
/// <para>Anchors: <c>docs/proofs/PROOF_F95_ANGLE_AT_QUADRATIC_ZERO.md</c> (F95) +
/// <c>experiments/F86_EP_THROUGH_THE_CLOCK.md</c> (the EP 2-level + the clock Rotation) +
/// <c>experiments/CRITICAL_SLOWING_AT_THE_CUSP.md</c> (the cusp recursion) +
/// <c>docs/NAVIGATING_THE_DIMENSIONS.md</c> (the interior-horizon axis) +
/// <c>hypotheses/PAIR_BREAKING_AT_THE_HORIZON.md</c> (the fold / break reading) +
/// <c>hypotheses/FRAGILE_BRIDGE.md</c> (the separate gain-loss system, whose threshold at two qubits per chain is generically its own EP2 on the real γ axis).</para></summary>
public sealed class TransitionBridgeF95SiblingClaim : Claim
{
    /// <summary>Parent: the positive-b F95 angle θ = arctan(√(c/b² − 1)) at a quadratic's discriminant
    /// zero. It is a reusable coordinate; it does not identify its applications.</summary>
    public F95AngleAtQuadraticZeroPi2Inheritance F95 { get; }

    /// <summary>The cusp anchor b = ½ (the bilinear / qubit-dimension fixed point).</summary>
    public const double CuspAnchorB = 0.5;

    /// <summary>The cusp / TransitionBridge double root: CΨ = b² = ¼.</summary>
    public const double Cusp = 0.25;

    public TransitionBridgeF95SiblingClaim(F95AngleAtQuadraticZeroPi2Inheritance f95)
        : base("Two distinct positive-b quadratic applications of F95: the dimensionless recurrence " +
               "z_rec²−z_rec+CΨ at b=½, and the F86 decay polynomial in z_decay=−λ at b=4γ₀; " +
               "the two F86 angle evaluations agree within numerical tolerance, without identifying the objects",
               Tier.Tier1Derived,
               "docs/proofs/PROOF_F95_ANGLE_AT_QUADRATIC_ZERO.md + experiments/F86_EP_THROUGH_THE_CLOCK.md + " +
               "experiments/CRITICAL_SLOWING_AT_THE_CUSP.md + docs/NAVIGATING_THE_DIMENSIONS.md + " +
               "hypotheses/PAIR_BREAKING_AT_THE_HORIZON.md + hypotheses/FRAGILE_BRIDGE.md")
    {
        F95 = f95 ?? throw new ArgumentNullException(nameof(f95));
    }

    /// <summary>Public factory: builds the claim with a fresh, parentless F95 claim.</summary>
    public static TransitionBridgeF95SiblingClaim Build()
    {
        var f95 = new F95AngleAtQuadraticZeroPi2Inheritance();
        return new TransitionBridgeF95SiblingClaim(f95);
    }

    /// <summary>Shared singleton; the claim is an algebraic identity, block-independent.</summary>
    public static TransitionBridgeF95SiblingClaim Shared { get; } = Build();

    /// <summary>The cusp's F95 angle at b = ½: θ = arctan(√(4·CΨ − 1)), via the F95 parent. NaN for
    /// CΨ &lt; ¼; exactly zero at the discriminant-zero boundary CΨ = ¼.</summary>
    public double CuspAngle(double cpsi) => F95.ThetaGeneral(cpsi, CuspAnchorB);

    /// <summary>Q_EP = 2/g_eff, where the EP block's discriminant vanishes (the modes coalesce).</summary>
    public static double QEp(double gEff)
    {
        RequirePositiveFiniteGEff(gEff);
        return 2.0 / gEff;
    }

    /// <summary>Positive quadratic anchor b = 4γ₀ in z_decay = −λ, the absorption rung ⟨n_XY⟩ = 2.</summary>
    public double EpAnchorB(double gamma0)
    {
        RequirePositiveGammaZero(gamma0);
        return 4.0 * gamma0;
    }

    /// <summary>F86 toy two-level Liouvillian roots λ = −4γ₀ ± √(4γ₀² − J²g_eff²).</summary>
    public (Complex Plus, Complex Minus) EpLiouvillianRoots(
        double gamma0,
        double q,
        double gEff)
    {
        RequirePositiveGammaZero(gamma0);
        RequireValidQ(q);
        RequirePositiveFiniteGEff(gEff);
        double j = q * gamma0;
        double discriminantQuarter = 4.0 * gamma0 * gamma0 - j * j * gEff * gEff;
        Complex branch = Complex.Sqrt(new Complex(discriminantQuarter, 0.0));
        return (
            new Complex(-4.0 * gamma0, 0.0) + branch,
            new Complex(-4.0 * gamma0, 0.0) - branch);
    }

    /// <summary>Positive-decay roots constructed branch by branch as z_decay = −λ. Their sum is 8γ₀
    /// and their midpoint is b = 4γ₀.</summary>
    public (Complex FromLambdaPlus, Complex FromLambdaMinus) EpDecayRoots(
        double gamma0,
        double q,
        double gEff)
    {
        var lambda = EpLiouvillianRoots(gamma0, q, gEff);
        return (-lambda.Plus, -lambda.Minus);
    }

    /// <summary>The EP block's clock Rotation angle arctan(ω/gap) = arctan(|Im λ|/|Re λ|) for the
    /// F86 2-level eigenvalue λ = −4γ₀ ± √(4γ₀² − J²g_eff²), J = Q·γ₀. Nonzero above
    /// the EP, zero at the EP, and NaN below the EP.</summary>
    public double EpClockAngle(double gamma0, double q, double gEff)
    {
        RequirePositiveGammaZero(gamma0);
        RequireValidQ(q);
        RequirePositiveFiniteGEff(gEff);
        if (q == QEp(gEff)) return 0.0;
        double j = q * gamma0;
        double disc = 4.0 * gamma0 * gamma0 - j * j * gEff * gEff; // 4γ₀² − J²g_eff²
        if (disc > 0.0) return double.NaN;
        var z = EpDecayRoots(gamma0, q, gEff).FromLambdaPlus;
        return Math.Atan(Math.Abs(z.Imaginary) / z.Real);
    }

    /// <summary>The EP block's F95 angle, via the F95 parent: the quadratic
    /// z² − 8γ₀z + (12γ₀² + J²g_eff²), for z=−λ, has b = 4γ₀ &gt; 0 and
    /// c = 12γ₀² + J²g_eff². θ = arctan(√(c/b² − 1)) above
    /// the EP (c > b² ⟺ J²g_eff² > 4γ₀²), zero at the EP, and NaN below.</summary>
    public double EpF95Angle(double gamma0, double q, double gEff)
    {
        RequirePositiveGammaZero(gamma0);
        RequireValidQ(q);
        RequirePositiveFiniteGEff(gEff);
        if (q == QEp(gEff)) return 0.0;
        double j = q * gamma0;
        double b = EpAnchorB(gamma0);
        double c = 12.0 * gamma0 * gamma0 + j * j * gEff * gEff;
        return F95.ThetaGeneral(c, b);
    }

    /// <summary>At J = 0 the two decay rates are the absorption rungs ⟨n_XY⟩ = 1 and 3 (2γ₀ and
    /// 6γ₀), whose midpoint is the anchor b = 4γ₀. Computed from the Liouvillian roots.</summary>
    public bool UncoupledRatesAreRungsOneAndThreeAroundTheAnchor(double gamma0)
    {
        var (plus, minus) = EpLiouvillianRoots(gamma0, 0.0, 1.0);
        double slow = -plus.Real, fast = -minus.Real;
        return Math.Abs(slow - 2.0 * gamma0 * 1) < 1e-14 * gamma0
            && Math.Abs(fast - 2.0 * gamma0 * 3) < 1e-14 * gamma0
            && Math.Abs((slow + fast) / 2.0 - EpAnchorB(gamma0)) < 1e-14 * gamma0;
    }

    /// <summary>Numerical check of the algebraic identity: at and above the EP the clock Rotation
    /// expression equals the F95 expression within a fixed tolerance.</summary>
    public bool EpClockAngleEqualsF95Angle(double gamma0, double q, double gEff)
    {
        double a = EpClockAngle(gamma0, q, gEff);
        double f = EpF95Angle(gamma0, q, gEff);
        if (double.IsNaN(a) || double.IsNaN(f)) return false;
        return Math.Abs(a - f) < 1e-12;
    }

    public override string DisplayName =>
        "Two distinct positive-b quadratic applications of the F95 root-angle coordinate";

    public override string Summary =>
        $"two distinct positive-b quadratic applications: z_rec at b=½ and F86 z_decay=−λ at b=4γ₀ " +
        $"(Q_EP=2/g_eff). Within the F86 polynomial, its clock and F95-coordinate evaluations agree " +
        $"numerically (sample g_eff=4/3, Q=2.5: {EpClockAngleEqualsF95Angle(1.0, 2.5, 4.0 / 3.0)}). " +
        $"This is a coordinate comparison, not an object identity ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("the cusp (the TransitionBridge), b = ½",
                summary: "dimensionless recurrence z_rec²−z_rec+CΨ, double root at CΨ=¼=(½)²; θ=arctan(√(4CΨ−1)); at ¼ the recurrence roots meet. This does not name a Liouvillian eigenvalue.");
            yield return new InspectableNode("the genuine F86 toy 2×2 EP, b = 4γ₀",
                summary: "F86 2-level in z_decay=−λ: z_decay²−8γ₀z_decay+(12γ₀²+J²g_eff²), double root (EP) at Q_EP=2/g_eff; b=4γ₀>0, the absorption rung ⟨n_XY⟩=2 between the uncoupled rates 2γ₀ and 6γ₀; above it the F95 angle and clock Rotation are the same real-arithmetic expression.");
            yield return new InspectableNode("the coordinate comparison (algebraic; numerically checked)",
                summary: $"within F86, clock angle = F95-coordinate angle (g_eff=0.8,Q=3: {EpClockAngleEqualsF95Angle(1.0, 3.0, 0.8)}; g_eff=4/3,Q=2: {EpClockAngleEqualsF95Angle(1.0, 2.0, 4.0 / 3.0)}). Separately, the recurrence coordinate at CΨ=1/3 is {CuspAngle(1.0 / 3.0) * 180.0 / Math.PI:F1}° (= 30°). Different variables and anchors remain different objects.");
            yield return new InspectableNode("how they connect, and how not",
                summary: "reused coordinate: the F95 angle at a positive-b quadratic zero. Distinct objects: two polynomials, two anchors b (½ vs 4γ₀), and two variables (dimensionless z_rec vs decay z_decay=−λ). No hidden identity or shared dynamics.");
            yield return new InspectableNode("FRAGILE_BRIDGE boundary",
                summary: "FRAGILE_BRIDGE threshold, at generic couplings and two qubits per chain: a defective EP2 on the real γ axis of a separate Σγ=0 gain-loss system, where, in the first popcount block of its Liouvillian to go unstable, two eigenvalues on the imaginary axis meet at zero decay and form a 2×2 Jordan block (at five exact couplings the threshold is zero, with no EP). It is not the F86 toy 2×2 quadratic, whose double root sits at the positive decay z_decay=4γ₀, and it is not licensed by this Tier-1 siblinghood.");
            yield return F95;
        }
    }

    private static void RequirePositiveGammaZero(double gamma0)
    {
        if (!double.IsFinite(gamma0) || gamma0 <= 0.0)
            throw new ArgumentOutOfRangeException(
                nameof(gamma0),
                "gamma0 must be positive for the finite positive-b F95 embedding");
    }

    private static void RequirePositiveFiniteGEff(double gEff)
    {
        if (!double.IsFinite(gEff) || gEff <= 0.0)
            throw new ArgumentOutOfRangeException(
                nameof(gEff),
                "gEff must be finite and positive");
    }

    private static void RequireValidQ(double q)
    {
        if (!double.IsFinite(q) || q < 0.0)
            throw new ArgumentOutOfRangeException(
                nameof(q),
                "q must be finite and non-negative for the declared J>=0 coordinate");
    }
}
