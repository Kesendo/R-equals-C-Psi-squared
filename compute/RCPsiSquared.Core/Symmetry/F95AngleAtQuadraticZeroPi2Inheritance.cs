using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>
/// F95 is the principal root angle of z^2-2bz+c=0 for real c and finite b>0.
/// For c>b^2 the upper root is b+i sqrt(c-b^2), so
/// theta=atan(sqrt(c/b^2-1)); theta is zero at c=b^2 and undefined below it.
/// At b=1/2 the discriminant-zero locus is the single real point c=1/4.
/// A circle |c|=1/4 is a different radial object. This coordinate identity
/// does not derive superposition, the Born rule, or a physical phase change.
/// </summary>
public sealed class F95AngleAtQuadraticZeroPi2Inheritance : Claim
{
    public const double B = 0.5;
    public const double Threshold = B * B;

    /// <summary>Returns zero at c = 1/4 for the framework specialization.</summary>
    public double ThetaForFramework(double c) => ThetaGeneral(c, B);

    public double ThetaGeneral(double c, double b)
    {
        if (!double.IsFinite(b) || b <= 0.0)
            throw new ArgumentOutOfRangeException(nameof(b), b, "b must be finite and > 0.");
        if (!double.IsFinite(c))
            throw new ArgumentOutOfRangeException(nameof(c), c, "c must be finite.");
        double thresh = b * b;
        if (c < thresh)
            return double.NaN;
        return Math.Atan(Math.Sqrt(c / thresh - 1.0));
    }

    public bool FrameworkSpecializationAgrees(double c)
    {
        double general = ThetaGeneral(c, B);
        double specialized;
        if (c < Threshold)
            specialized = double.NaN;
        else
            specialized = Math.Atan(Math.Sqrt(4.0 * c - 1.0));
        return double.IsNaN(general)
            ? double.IsNaN(specialized)
            : Math.Abs(general - specialized) < 1e-15;
    }

    public bool BellPlusInitialAngleIs30Degrees() =>
        Math.Abs(ThetaForFramework(1.0 / 3.0) - Math.PI / 6.0) < 1e-12;

    public F95AngleAtQuadraticZeroPi2Inheritance()
        : base(
            "F95 positive-b quadratic root angle theta(c;b)=atan(sqrt(c/b^2-1))",
            Tier.Tier1Derived,
            "docs/proofs/PROOF_F95_ANGLE_AT_QUADRATIC_ZERO.md + " +
            "docs/ANALYTICAL_FORMULAS.md F95 + " +
            "simulations/angle_at_zero_tier1_candidate.py + " +
            "experiments/BOUNDARY_NAVIGATION.md")
    {
    }

    public override string DisplayName => "F95 positive-b quadratic root angle";

    public override string Summary =>
        $"theta(c;b)=atan(sqrt(c/b^2-1)) for finite b > 0 and c>b^2, theta=0 at c=b^2, " +
        $"and the real angle is undefined below the boundary; at b={B}, " +
        $"the double root is the single real point c={Threshold}; parentless.";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return InspectableNode.RealScalar("b (= 1/2 specialization)", B);
            yield return InspectableNode.RealScalar("b^2 (= 1/4 double-root point)", Threshold);
            yield return new InspectableNode(
                "Finite readout",
                summary: $"At c=0.286, theta={ThetaForFramework(0.286):G6}. This is a finite coordinate readout, not a canonical crossing.");
            yield return new InspectableNode(
                "Object boundary",
                summary: "The radial set |c|=1/4 is separate from the period-one recurrence locus. " +
                         "The formula alone implies no Born, superposition, or hardware mechanism.");
        }
    }
}
