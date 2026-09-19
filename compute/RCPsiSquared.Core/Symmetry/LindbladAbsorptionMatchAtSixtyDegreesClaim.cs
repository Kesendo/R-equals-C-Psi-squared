using System.Numerics;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>
/// A local algebraic cross-reading of one named two-level model. If that model is written as
/// <c>lambda_± = -gamma0 ± iJ</c>, <c>gamma0 &gt; 0</c>, and <c>Q=J/gamma0</c>, then
/// <c>|lambda_±|/gamma0 = sqrt(1+Q²)</c>. At <c>Q=sqrt(3)</c> the ratio is 2 and the
/// positive decay roots <c>z=-lambda</c> obey
/// <c>z²-2gamma0*z+(gamma0²+J²)=0</c>. Thus <c>b=gamma0&gt;0</c> and the F95
/// principal angle is <c>atan(sqrt(3)) = pi/3</c>.
///
/// <para>The value 2 in this calculation is a literal magnitude ratio of this model. It does not
/// inherit from the absorption theorem, identify a universal one-disagreement decay rate, or make
/// the angle a descendant of a canonical-angle or Pi2 ladder. The sole typed parent is F95, whose
/// quadratic-angle formula is actually used.</para>
/// </summary>
public sealed class LindbladAbsorptionMatchAtSixtyDegreesClaim : Claim
{
    /// <summary>The positive solution of <c>sqrt(1+Q²)=2</c>.</summary>
    public static readonly double QValue = Math.Sqrt(3.0);

    /// <summary>The corresponding angle in degrees: <c>atan(sqrt(3))=60</c>.</summary>
    public const double CanonicalAngleDegrees = 60.0;

    /// <summary>Genuine parent: the positive-b F95 quadratic-angle law.</summary>
    public F95AngleAtQuadraticZeroPi2Inheritance F95 { get; }

    public LindbladAbsorptionMatchAtSixtyDegreesClaim(
        F95AngleAtQuadraticZeroPi2Inheritance f95)
        : base(
            "Named two-level cross-reading: Q=sqrt(3) gives |lambda|/gamma0=2 and the F95 angle pi/3",
            Tier.Tier1Derived,
            "compute/RCPsiSquared.Core/Symmetry/F95AngleAtQuadraticZeroPi2Inheritance.cs + " +
            "docs/proofs/PROOF_F95_ANGLE_AT_QUADRATIC_ZERO.md + docs/Q_REGIME_ANCHORS.md")
    {
        F95 = f95 ?? throw new ArgumentNullException(nameof(f95));
    }

    /// <summary>Builds the reading with a fresh, parentless F95 claim.</summary>
    public static LindbladAbsorptionMatchAtSixtyDegreesClaim Build() =>
        new(new F95AngleAtQuadraticZeroPi2Inheritance());

    public static LindbladAbsorptionMatchAtSixtyDegreesClaim Shared { get; } = Build();

    /// <summary>The named Liouvillian pair at Q=sqrt(3).</summary>
    public (Complex Plus, Complex Minus) LiouvillianRootsAtQSqrt3(double gammaZero)
    {
        RequirePositiveGammaZero(gammaZero);
        double j = QValue * gammaZero;
        return (new Complex(-gammaZero, j), new Complex(-gammaZero, -j));
    }

    /// <summary>The positive-decay pair constructed explicitly as z=-lambda.</summary>
    public (Complex FromLambdaPlus, Complex FromLambdaMinus) PositiveDecayRootsAtQSqrt3(
        double gammaZero)
    {
        var lambda = LiouvillianRootsAtQSqrt3(gammaZero);
        return (-lambda.Plus, -lambda.Minus);
    }

    /// <summary>Positive quadratic anchor b=gamma0 for the decay roots.</summary>
    public double PositiveDecayAnchorB(double gammaZero)
    {
        RequirePositiveGammaZero(gammaZero);
        return gammaZero;
    }

    /// <summary>|lambda|² at Q=sqrt(3), obtained from the constructed pair.</summary>
    public double LindbladMagnitudeSquaredAtQSqrt3(double gammaZero)
    {
        var lambda = LiouvillianRootsAtQSqrt3(gammaZero);
        return lambda.Plus.Magnitude * lambda.Plus.Magnitude;
    }

    /// <summary>Floating reconstruction of the exact ratio 2. gamma0=0 is
    /// outside the finite positive-b F95 domain and is rejected.</summary>
    public double LindbladMagnitudeOverGamma0Computed(double gammaZero)
    {
        var lambda = LiouvillianRootsAtQSqrt3(gammaZero);
        return lambda.Plus.Magnitude / gammaZero;
    }

    /// <summary>F95 principal angle reconstructed from the positive-decay
    /// quadratic b=gamma0, c=gamma0²+J².</summary>
    public double F95AngleAtQSqrt3Degrees(double gammaZero)
    {
        double b = PositiveDecayAnchorB(gammaZero);
        double j = QValue * gammaZero;
        double c = gammaZero * gammaZero + j * j;
        return AnchorConstants.RadiansToDegrees(F95.ThetaGeneral(c, b));
    }

    public override string DisplayName =>
        "Named two-level Q=sqrt(3) cross-reading (magnitude ratio 2; F95 angle pi/3)";

    public override string Summary =>
        $"For lambda_±=-gamma0±iJ in the named two-level model, Q=sqrt(3)≈{QValue:F4} gives " +
        $"|lambda_±|/gamma0=2; its z=-lambda roots have positive b=gamma0 and F95 angle " +
        $"{CanonicalAngleDegrees} degrees. The equality is " +
        "local algebra, not an absorption-rate identity or canonical-angle ancestry.";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return InspectableNode.RealScalar("Q (= sqrt(3))", QValue, "F6");
            yield return InspectableNode.RealScalar(
                "|lambda_±|/gamma0 in the named model", LindbladMagnitudeOverGamma0Computed(1.0), "F11");
            yield return InspectableNode.RealScalar(
                "F95 angle from z=-lambda (degrees)", F95AngleAtQSqrt3Degrees(1.0), "F8");
            yield return new InspectableNode(
                "scope seam",
                summary: "The number 2 is computed directly from sqrt(1+Q²). No AbsorptionTheorem, Niven, Pi2-ladder, or triple-axis parent is asserted.");
            yield return F95;
        }
    }

    private static void RequirePositiveGammaZero(double gammaZero)
    {
        if (!double.IsFinite(gammaZero) || gammaZero <= 0.0)
            throw new ArgumentOutOfRangeException(
                nameof(gammaZero),
                "gammaZero must be positive for the finite positive-b F95 ratio");
    }
}
