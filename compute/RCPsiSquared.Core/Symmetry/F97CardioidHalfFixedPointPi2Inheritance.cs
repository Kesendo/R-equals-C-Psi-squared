using System.Numerics;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>
/// F97 is the period-one boundary of z -> z^2+c. The marginal selected fixed
/// point satisfies |2z*|=1, hence z*(phi)=exp(i phi)/2 and
/// c(phi)=z*(1-z*). This exact parameter-plane identity is not evidence that
/// a plotted hardware CPsi_com trace is a Mandelbrot orbit. F95 applies instead
/// to real c > b² and finite b>0; its root angle and this complex-c locus are
/// sibling views, not executable ancestry.
/// </summary>
public sealed class F97CardioidHalfFixedPointPi2Inheritance : Claim
{
    public const double B = 0.5;
    public const double Threshold = B * B;

    public Complex CardioidFixedPoint(double phi) =>
        B * Complex.Exp(Complex.ImaginaryOne * phi);

    public Complex CardioidC(double phi)
    {
        Complex z = CardioidFixedPoint(phi);
        return z * (Complex.One - z);
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

    public bool AlgebraicIdentityHolds(double phi)
    {
        Complex c = CardioidC(phi);
        Complex z = CardioidFixedPoint(phi);
        return Complex.Abs(z * z - z + c) < 1e-14;
    }

    public bool MagnitudeInvariantAroundCardioid() =>
        CanonicalAngles().All(phi => Math.Abs(FixedPointMagnitude(phi) - B) < 1e-14);

    public bool SquaredMagnitudeInvariantAroundCardioid() =>
        CanonicalAngles().All(phi => Math.Abs(FixedPointMagnitudeSquared(phi) - Threshold) < 1e-14);

    public bool CuspAgreesWithF95Threshold()
    {
        Complex c = CardioidC(0.0);
        return Math.Abs(c.Real - Threshold) < 1e-15 && Math.Abs(c.Imaginary) < 1e-15;
    }

    public bool TailAtMinusThreeQuarters()
    {
        Complex c = CardioidC(Math.PI);
        return Math.Abs(c.Real + 0.75) < 1e-14 && Math.Abs(c.Imaginary) < 1e-14;
    }

    public bool OriginRootsAreDistinct()
    {
        Complex z0 = Complex.Zero;
        Complex z1 = Complex.One;
        return z0 != z1 && z0 * z0 - z0 == Complex.Zero && z1 * z1 - z1 == Complex.Zero;
    }

    private static double[] CanonicalAngles() =>
        [0.0, Math.PI / 3, Math.PI / 2, Math.PI, 4 * Math.PI / 3, 5 * Math.PI / 3];

    public F97CardioidHalfFixedPointPi2Inheritance()
        : base(
            "F97 period-one cardioid: z*=exp(i phi)/2, c=z*(1-z*), |2z*|=1",
            Tier.Tier1Derived,
            "docs/proofs/PROOF_F97_CARDIOID_HALF_FIXED_POINT.md + " +
            "docs/ANALYTICAL_FORMULAS.md F97 + " +
            "simulations/cardioid_parametrization_tier1.py")
    {
    }

    public override string DisplayName => "F97 parentless period-one cardioid";

    public override string Summary =>
        "z*(phi)=exp(i phi)/2 follows from |2z*|=1; c(phi)=z*(1-z*). " +
        "Only the selected root is marginal in general; numerical evaluations are machine-precision checks.";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return InspectableNode.RealScalar("selected fixed-point magnitude", B);
            yield return InspectableNode.RealScalar("real cusp c(0)", Threshold);
            Complex selected = CardioidFixedPoint(Math.PI / 2);
            Complex other = OtherFixedPoint(Math.PI / 2);
            yield return new InspectableNode(
                "phi=pi/2 two-root control",
                summary: $"selected z={selected}, |2z|={MultiplierMagnitude(selected):G6}; " +
                         $"other z={other}, |2z|={MultiplierMagnitude(other):G6}=sqrt(5).");
            yield return new InspectableNode(
                "Object boundary",
                summary: "The cardioid is an iteration-parameter stability boundary. A CPsi_com overlay is a coordinate comparison, not an orbit claim.");
        }
    }
}
