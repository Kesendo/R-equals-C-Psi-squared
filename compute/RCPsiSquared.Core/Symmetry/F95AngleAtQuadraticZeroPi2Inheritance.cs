using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>
/// F95 is the principal root angle of z² − 2bz + c = 0 for real c and finite b &gt; 0.
/// For c &gt; b² the upper root is b + i√(c − b²), so θ = arctan(√(c/b² − 1)); θ is zero at
/// c = b², the double root, and undefined below it. At b = ½ the discriminant-zero locus is the
/// single real point c = ¼ and θ = arctan(√(4c − 1)), the compass F15. A circle |c| = ¼ is a
/// different radial object. This coordinate identity does not derive superposition, the Born
/// rule, or a physical phase change.
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
            "F95 positive-b quadratic root angle θ(c; b) = arctan(√(c/b² − 1))",
            Tier.Tier1Derived,
            "docs/proofs/PROOF_F95_ANGLE_AT_QUADRATIC_ZERO.md + " +
            "docs/ANALYTICAL_FORMULAS.md F95 + " +
            "simulations/angle_at_zero_tier1_candidate.py + " +
            "experiments/BOUNDARY_NAVIGATION.md")
    {
    }

    public override string DisplayName => "F95 root angle θ(c; b) at the quadratic's double root";

    public override string Summary =>
        $"θ(c; b) = arctan(√(c/b² − 1)) for finite b > 0 and c > b², θ = 0 at the double root c = b², " +
        $"undefined below it; at b = ½ the double root is the single real point c = ¼.";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return InspectableNode.RealScalar("b (the ½ specialization)", B);
            yield return InspectableNode.RealScalar("b² (the ¼ double-root point)", Threshold);
            yield return new InspectableNode(
                "Bell⁺ start",
                summary: $"At c = ⅓, θ = {ThetaForFramework(1.0 / 3.0) * 180.0 / Math.PI:G6}° (the Bell⁺ initial point).");
            yield return new InspectableNode(
                "Object boundary",
                summary: "The radial set |c| = ¼ is separate from the period-one recurrence locus. " +
                         "The formula alone implies no Born, superposition, or hardware mechanism.");
        }
    }
}
