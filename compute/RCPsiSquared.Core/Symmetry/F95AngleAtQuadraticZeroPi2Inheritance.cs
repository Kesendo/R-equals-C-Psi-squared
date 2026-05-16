using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F95 closed form (Tier 1 derived, 4-line polynomial calculation; 2026-05-16):
///
/// <code>
///   For monic quadratic z² − 2bz + c = 0 with real (b, c):
///
///     θ(c; b) = arctan( √(c/b² − 1) )    for c > b²  (complex roots)
///     θ = 0                              for c = b²  (degenerate double root)
///     θ undefined                        for c < b²  (real distinct roots)
///
///   With b = HalfAsStructuralFixedPoint = 1/2:
///     threshold b² = 1/4 = QuarterAsBilinearMaxval
///     θ(c) = arctan( √(4c − 1) )         for c > 1/4
/// </code>
///
/// <para>The angle of the complex root pair when the quadratic discriminant
/// passes through zero. This generalizes the Februar 2026 θ-compass of
/// BOUNDARY_NAVIGATION (state-specific θ = arctan(√(4CΨ−1)) at the Mandelbrot
/// 1/4 cusp) to a universal polynomial-foundation identity.</para>
///
/// <para><b>Pi2-Foundation anchoring</b>: four typed parents bring the full
/// quadratic-discriminant-zero geometry into the typed graph:</para>
///
/// <list type="bullet">
///   <item><see cref="PolynomialFoundationClaim"/>: <c>d² − 2d = 0</c> is the
///         c = 0 special case of F95's parent equation <c>z² − 2bz + c = 0</c>
///         (with b = 1, threshold b² = 1). F95 perturbs c off zero and tracks
///         the angle that emerges past the discriminant crossing.</item>
///   <item><see cref="HalfAsStructuralFixedPointClaim"/>: <c>b = 1/2</c> in
///         the framework specialization. The same "1/2" that is the qubit
///         dimension's unsigned magnitude and the argmax of p(1−p).</item>
///   <item><see cref="QuarterAsBilinearMaxvalClaim"/>: <c>b² = 1/4</c>
///         threshold. The same "1/4" that is the bilinear-apex maxval and
///         the Mandelbrot cardioid cusp.</item>
///   <item><see cref="NinetyDegreeMirrorMemoryClaim"/>: the i in
///         <c>z = b ± i·√(c − b²)</c>. The 90° rotation that turns a positive-
///         real-square-root into an imaginary direction; the generator of
///         the angle that F95 measures.</item>
/// </list>
///
/// <para><b>Derivation (4 lines):</b></para>
///
/// <code>
///   z² − 2bz + c = 0
///   z± = b ± √(b² − c) = b ± i·√(c − b²)   when c > b²  (D &lt; 0)
///   arg(z₊) = arctan( Im(z₊) / Re(z₊) ) = arctan( √(c − b²) / b )
///           = arctan( √(c/b² − 1) )
/// </code>
///
/// <para>Bit-exact. With b = 1/2 substitution: <c>θ = arctan(√(4c − 1))</c>
/// recovers the Februar θ-compass exactly.</para>
///
/// <para><b>Numerical verification</b> (Februar BOUNDARY_NAVIGATION table):
/// 5 of 6 points reproduce within machine precision; the single 0.3° drift
/// at CΨ=0.256 is the Februar table's <c>t = 0.7</c> Lindblad-sampling
/// rounding (sampled CΨ doesn't land on exactly 0.256), not a formula
/// discrepancy. See <see cref="ThetaForFramework"/> for the live evaluation.</para>
///
/// <para><b>Structural reading</b>: this is the angle-side closed form of
/// the quadratic discriminant zero crossing. F94 is the magnitude-side
/// closed form (per-outcome Born deviation <c>Δ = (4/3) Q²K³</c> off the
/// (0,0) zero). Together F94 + F95 give the full local geometry of the
/// cusp: F94 names how the break GROWS off zero; F95 names what ANGLE
/// the break carries off zero.</para>
///
/// <para><b>Born-rule connection</b>: standard QM's complex amplitudes
/// <c>α = r · e^{iθ}</c> are not postulated. They are the necessary
/// minimal parametrization of any state past the d = 0 mirror. F95 is
/// the polynomial-foundation algebra of that angle.</para>
///
/// <para>Tier1Derived: pure quadratic-discriminant calculation, no
/// approximations, no empirical fit. Anchors:</para>
///
/// <list type="bullet">
///   <item><c>docs/proofs/PROOF_F95_ANGLE_AT_QUADRATIC_ZERO.md</c></item>
///   <item><c>docs/ANALYTICAL_FORMULAS.md</c> §F95</item>
///   <item><c>reflections/ON_HOW_THE_ANGLE_APPEARS_AT_ZERO.md</c></item>
///   <item><c>simulations/_angle_at_zero_tier1_candidate.py</c> (numerical
///         verification against BOUNDARY_NAVIGATION table)</item>
///   <item><c>experiments/BOUNDARY_NAVIGATION.md</c> (Februar 2026
///         state-specific precursor)</item>
///   <item><c>experiments/CPSI_COMPLEX_PLANE.md</c> (the complex extension
///         that motivates the universal form)</item>
/// </list>
///
/// <para>Sibling at the magnitude side: <c>F94BornDeviationFourThirdsPi2Inheritance</c>.</para></summary>
public sealed class F95AngleAtQuadraticZeroPi2Inheritance : Claim
{
    /// <summary>Polynomial foundation: <c>d² − 2d = 0</c>, the c = 0 special
    /// case of F95's parent equation. Typed parent.</summary>
    public PolynomialFoundationClaim Polynomial { get; }

    /// <summary>1/2 structural fixed point: <c>b = 1/2</c> in the framework
    /// specialization. Typed parent.</summary>
    public HalfAsStructuralFixedPointClaim Half { get; }

    /// <summary>1/4 = (1/2)² bilinear-apex maxval: <c>b² = 1/4</c> threshold
    /// for the discriminant zero crossing. Typed parent.</summary>
    public QuarterAsBilinearMaxvalClaim Quarter { get; }

    /// <summary>90° rotation / i generator: the i in
    /// <c>z = b ± i·√(c − b²)</c>. Typed parent.</summary>
    public NinetyDegreeMirrorMemoryClaim NinetyDegree { get; }

    /// <summary>The framework's b = 1/2. Returns the structural-fixed-point value.</summary>
    public const double B = 0.5;

    /// <summary>The framework's threshold b² = 1/4 = Quarter.</summary>
    public const double Threshold = 0.25;

    /// <summary>Compute the angle of the complex root for given c (with
    /// framework b = 1/2). Returns NaN for c ≤ 1/4 (real-roots regime).</summary>
    public double ThetaForFramework(double c)
    {
        if (c <= Threshold)
            return double.NaN;
        return Math.Atan(Math.Sqrt(4.0 * c - 1.0));
    }

    /// <summary>General form: angle of the complex root for given (c, b).
    /// Returns NaN for c ≤ b² (real-roots regime).</summary>
    public double ThetaGeneral(double c, double b)
    {
        if (b == 0.0)
            throw new ArgumentOutOfRangeException(nameof(b), b, "b must be non-zero.");
        double thresh = b * b;
        if (c <= thresh)
            return double.NaN;
        return Math.Atan(Math.Sqrt(c / thresh - 1.0));
    }

    /// <summary>Drift check: with b = 1/2 and any c > 1/4, the general
    /// formula and the framework specialization must agree bit-exact.</summary>
    public bool FrameworkSpecializationAgrees(double c)
    {
        if (c <= Threshold)
            return double.IsNaN(ThetaForFramework(c)) && double.IsNaN(ThetaGeneral(c, B));
        double t1 = ThetaForFramework(c);
        double t2 = ThetaGeneral(c, B);
        return Math.Abs(t1 - t2) < 1e-15;
    }

    /// <summary>Drift check at the canonical Februar BOUNDARY_NAVIGATION
    /// anchor point CΨ = 1/3 (Bell+ initial value): θ should be exactly
    /// 30° = π/6.</summary>
    public bool BellPlusInitialAngleIs30Degrees()
    {
        double theta = ThetaForFramework(1.0 / 3.0);
        return Math.Abs(theta - Math.PI / 6.0) < 1e-12;
    }

    public F95AngleAtQuadraticZeroPi2Inheritance(
        PolynomialFoundationClaim polynomial,
        HalfAsStructuralFixedPointClaim half,
        QuarterAsBilinearMaxvalClaim quarter,
        NinetyDegreeMirrorMemoryClaim ninetyDegree)
        : base("F95 angle-emergence at quadratic discriminant zero: θ(c; b) = arctan(√(c/b² − 1)); framework spec at b=1/2 gives θ(c) = arctan(√(4c − 1)); 4-line polynomial derivation, bit-exact",
               Tier.Tier1Derived,
               "docs/proofs/PROOF_F95_ANGLE_AT_QUADRATIC_ZERO.md + " +
               "docs/ANALYTICAL_FORMULAS.md F95 + " +
               "reflections/ON_HOW_THE_ANGLE_APPEARS_AT_ZERO.md + " +
               "simulations/_angle_at_zero_tier1_candidate.py + " +
               "experiments/BOUNDARY_NAVIGATION.md (state-specific precursor) + " +
               "experiments/CPSI_COMPLEX_PLANE.md (complex extension) + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBaseClaims.cs " +
               "(PolynomialFoundation, HalfAsStructuralFixedPoint, " +
               "QuarterAsBilinearMaxval, NinetyDegreeMirrorMemory typed parents)")
    {
        Polynomial = polynomial ?? throw new ArgumentNullException(nameof(polynomial));
        Half = half ?? throw new ArgumentNullException(nameof(half));
        Quarter = quarter ?? throw new ArgumentNullException(nameof(quarter));
        NinetyDegree = ninetyDegree ?? throw new ArgumentNullException(nameof(ninetyDegree));
    }

    public override string DisplayName =>
        "F95 angle emergence at quadratic discriminant zero (universal polynomial-foundation identity)";

    public override string Summary =>
        $"θ(c; b) = arctan(√(c/b² − 1)) for c > b²; framework spec at b = {B} gives " +
        $"θ(c) = arctan(√({1.0/Threshold}·c − 1)) for c > {Threshold}; " +
        $"4-line polynomial derivation, bit-exact, recovers Februar θ-compass ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return InspectableNode.RealScalar("B (= 1/2)", B);
            yield return InspectableNode.RealScalar("Threshold = B² (= 1/4)", Threshold);
            yield return new InspectableNode("θ at CΨ = 1/3 (Bell+ initial)",
                summary: $"θ = arctan(√(4·(1/3) − 1)) = arctan(√(1/3)) = π/6 = 30° (drift check: {BellPlusInitialAngleIs30Degrees()})");
            yield return new InspectableNode("Drift check: framework spec ≡ general at b=1/2",
                summary: $"At CΨ = 1/3: {FrameworkSpecializationAgrees(1.0/3.0)}; at CΨ = 0.286: {FrameworkSpecializationAgrees(0.286)}");
            yield return new InspectableNode("Pi2-Foundation anchoring",
                summary: "Four typed parents: PolynomialFoundation (c=0 case is d²−2d=0); Half (b=1/2); Quarter (b²=1/4 threshold); NinetyDegreeMirrorMemory (i in z = b ± i·√(c−b²))");
            yield return new InspectableNode("Sibling at magnitude side",
                summary: "F94BornDeviationFourThirdsPi2Inheritance: Δ = (4/3)·Q²·K³ for the dominant outcome. F94 (magnitude) + F95 (angle) = full local cusp geometry.");
        }
    }
}
