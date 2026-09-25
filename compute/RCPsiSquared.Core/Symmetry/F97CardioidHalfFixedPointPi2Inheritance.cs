using System.Numerics;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>
/// F97 is the period-one boundary of z ↦ z² + c. The marginal selected fixed point satisfies
/// |2z*| = 1, hence z*(φ) = e^(iφ)/2, and the boundary is c(φ) = e^(iφ)/2 − e^(2iφ)/4, which is
/// the bilinear form z*(1 − z*) of the selected root. That form is the typed parent's object:
/// <see cref="QuarterAsBilinearMaxvalClaim"/> (over real p, max p(1 − p) = ¼ at p = ½), reached here
/// on the real axis at φ = 0, the cusp c = ¼ where the two fixed points meet; around the complex
/// cardioid |c| runs from ¼ at the cusp to ¾ at φ = π, so the maximum is the real-axis statement.
/// The inheritance runs through this bilinear form only. F95's root angle applies instead to
/// real c > b² and finite b>0; it is a sibling coordinate that shares the value ¼ at b = ½, not an
/// ancestor. The identity is about the iteration's parameter plane; a plotted hardware CΨ trace is
/// not an orbit of z² + c.
/// </summary>
public sealed class F97CardioidHalfFixedPointPi2Inheritance : Claim
{
    public const double B = 0.5;
    public const double Threshold = B * B;

    /// <summary>Typed parent: the bilinear form p(1 − p) and its maximum ¼ at p = ½, which the
    /// cardioid's c = z*(1 − z*) takes at the cusp φ = 0.</summary>
    public QuarterAsBilinearMaxvalClaim Quarter { get; }

    public F97CardioidHalfFixedPointPi2Inheritance(QuarterAsBilinearMaxvalClaim quarter)
        : base(
            "F97 period-one cardioid: z*(φ) = e^(iφ)/2, c(φ) = e^(iφ)/2 − e^(2iφ)/4 = z*(1 − z*), |2z*| = 1",
            Tier.Tier1Derived,
            "docs/proofs/PROOF_F97_CARDIOID_HALF_FIXED_POINT.md + " +
            "docs/ANALYTICAL_FORMULAS.md F97 + " +
            "simulations/cardioid_parametrization_tier1.py + " +
            "experiments/CPSI_COMPLEX_PLANE.md")
    {
        Quarter = quarter ?? throw new ArgumentNullException(nameof(quarter));
    }

    public Complex CardioidFixedPoint(double phi) =>
        B * Complex.Exp(Complex.ImaginaryOne * phi);

    /// <summary>The closed form c(φ) = b·e^(iφ) − b²·e^(2iφ), computed from the exponential,
    /// independently of the fixed point.</summary>
    public Complex CardioidC(double phi)
    {
        Complex e = Complex.Exp(Complex.ImaginaryOne * phi);
        return B * e - Threshold * (e * e);
    }

    /// <summary>The bilinear form z*(1 − z*) of the selected root: the second route to c(φ).</summary>
    public Complex BilinearFormOfFixedPoint(double phi)
    {
        Complex z = CardioidFixedPoint(phi);
        return z * (Complex.One - z);
    }

    /// <summary>Both fixed points of z ↦ z² + c, the roots of z² − z + c = 0.</summary>
    public static (Complex Minus, Complex Plus) FixedPoints(Complex c)
    {
        Complex root = Complex.Sqrt(Complex.One - 4.0 * c);
        return ((Complex.One - root) / 2.0, (Complex.One + root) / 2.0);
    }

    public Complex OtherFixedPoint(double phi) => Complex.One - CardioidFixedPoint(phi);

    public static double MultiplierMagnitude(Complex z) => Complex.Abs(2.0 * z);

    public double FixedPointMagnitude(double phi) => Complex.Abs(CardioidFixedPoint(phi));

    public double FixedPointMagnitudeSquared(double phi)
    {
        double magnitude = FixedPointMagnitude(phi);
        return magnitude * magnitude;
    }

    public double FixedPointArgument(double phi)
    {
        Complex z = CardioidFixedPoint(phi);
        return Math.Atan2(z.Imaginary, z.Real);
    }

    public bool SelectedMultiplierIsMarginal(double phi) =>
        Math.Abs(MultiplierMagnitude(CardioidFixedPoint(phi)) - 1.0) < 1e-14;

    /// <summary>Two routes to c(φ) agree (exponential closed form against z*(1 − z*)), and the
    /// selected root solves the quadratic at that c.</summary>
    public bool AlgebraicIdentityHolds(double phi) => AlgebraicIdentityHolds(phi, CardioidC);

    /// <summary>The same two-route gate for any candidate closed form, so that a wrong closed
    /// form can be put through the gate the claim itself uses.</summary>
    public bool AlgebraicIdentityHolds(double phi, Func<double, Complex> closedForm)
    {
        Complex c = closedForm(phi);
        Complex z = CardioidFixedPoint(phi);
        return Complex.Abs(c - BilinearFormOfFixedPoint(phi)) < 1e-14
            && Complex.Abs(z * z - z + c) < 1e-14;
    }

    /// <summary>The quadratic formula at c(φ) returns the selected root as one of its two fixed
    /// points: the parametrization is recovered, not assumed. Near the cusp √(1 − 4c) is
    /// ill-conditioned and the root error grows like ε/|1 − 2z*|, so the gate is on the error times
    /// |1 − 2z*| against a few ε. The model holds while φ² ≫ ε (φ ≳ 1e-7); closer to the cusp
    /// 1 − 4c rounds to zero and the gate stops discriminating.</summary>
    public bool QuadraticFormulaRecoversSelectedRoot(double phi)
    {
        var (minus, plus) = FixedPoints(CardioidC(phi));
        Complex z = CardioidFixedPoint(phi);
        double error = Math.Min(Complex.Abs(minus - z), Complex.Abs(plus - z));
        return error * Complex.Abs(Complex.One - 2.0 * z) < 1e-14;
    }

    public bool MagnitudeInvariantAroundCardioid() =>
        CanonicalAngles().All(phi => Math.Abs(FixedPointMagnitude(phi) - B) < 1e-14);

    public bool SquaredMagnitudeInvariantAroundCardioid() =>
        CanonicalAngles().All(phi => Math.Abs(FixedPointMagnitudeSquared(phi) - Threshold) < 1e-14);

    /// <summary>At the cusp the bilinear form takes Quarter's value: c(0) = ¼ at z* = ½, exactly.</summary>
    public bool CuspTakesQuartersValue()
    {
        Complex c = CardioidC(0.0);
        return c == new Complex(Threshold, 0.0) && CardioidFixedPoint(0.0) == new Complex(B, 0.0);
    }

    public bool TailAtMinusThreeQuarters()
    {
        Complex c = CardioidC(Math.PI);
        return Math.Abs(c.Real + 0.75) < 1e-14 && Math.Abs(c.Imaginary) < 1e-14;
    }

    /// <summary>At c = 0 the quadratic formula gives two distinct fixed points, 0 (multiplier 0,
    /// superattracting) and 1 (multiplier 2, repelling).</summary>
    public static bool OriginFixedPointsAreDistinct()
    {
        var (minus, plus) = FixedPoints(Complex.Zero);
        return minus == Complex.Zero && plus == Complex.One
            && MultiplierMagnitude(minus) == 0.0 && MultiplierMagnitude(plus) == 2.0;
    }

    private static double[] CanonicalAngles() =>
        [0.0, Math.PI / 3, Math.PI / 2, Math.PI, 4 * Math.PI / 3, 5 * Math.PI / 3];

    public override string DisplayName => "F97 period-one cardioid c(φ) = z*(1 − z*), |z*| = ½";

    public override string Summary =>
        "z*(φ) = e^(iφ)/2 follows from |2z*| = 1; c(φ) = e^(iφ)/2 − e^(2iφ)/4 = z*(1 − z*), " +
        "the bilinear form whose real-axis maximum ¼ (Quarter) is the cusp; only the selected root is marginal in general.";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return InspectableNode.RealScalar("selected fixed-point magnitude", B);
            yield return InspectableNode.RealScalar("cusp c(0) = Quarter's maximum", Threshold);
            Complex selected = CardioidFixedPoint(Math.PI / 2);
            Complex other = OtherFixedPoint(Math.PI / 2);
            yield return new InspectableNode(
                "φ = π/2 two-root control",
                summary: $"selected z = {selected}, |2z| = {MultiplierMagnitude(selected):G6}; " +
                         $"other z = {other}, |2z| = {MultiplierMagnitude(other):G6} = √5.");
            yield return new InspectableNode(
                "Object boundary",
                summary: "The cardioid is a stability boundary in the iteration's parameter plane. A CΨ_com overlay is a coordinate comparison, not an orbit claim.");
        }
    }
}
