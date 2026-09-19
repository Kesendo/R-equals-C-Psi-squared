using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>
/// F96 is the setup-specific subdominant table for the same named N=4 ring as
/// F94. Its slopes are derived independently from its own Dyson and unitary
/// elements: M3=-4, U2=3/4 and M5=-20, U4=3/2. Ring/K4 equality is a finite
/// comparison; the chain value -4/3 for |10> refutes topology universality.
/// </summary>
public sealed class F96BornSubdominantSlopesPi2Inheritance : Claim
{
    public const int N = 4;
    public const int M3_SingleFlipped = -4;
    public const int U2_SingleFlipped_TimesFour = 3;
    public const int M5_DoubleFlipped = -20;
    public const int U4_DoubleFlipped_TimesTwo = 3;

    public double SlopeSingleFlipped =>
        (double)M3_SingleFlipped / (3.0 * U2_SingleFlipped_TimesFour / 4.0);

    public double SlopeDoubleFlipped =>
        (double)M5_DoubleFlipped / (5.0 * U4_DoubleFlipped_TimesTwo / 2.0);

    public double DeltaSingleFlipped(double K)
    {
        if (!double.IsFinite(K) || K < 0.0)
            throw new ArgumentOutOfRangeException(nameof(K), K, "K must be finite and >= 0.");
        return SlopeSingleFlipped * K;
    }

    public double DeltaDoubleFlipped(double K)
    {
        if (!double.IsFinite(K) || K < 0.0)
            throw new ArgumentOutOfRangeException(nameof(K), K, "K must be finite and >= 0.");
        return SlopeDoubleFlipped * K;
    }

    /// <summary>A numerical comparison inside the named ring table, not ancestry.</summary>
    public bool SingleFlipSlopeEqualsMinusFourThirdsSquared() =>
        Math.Abs(SlopeSingleFlipped + (4.0 / 3.0) * (4.0 / 3.0)) < 1e-15;

    /// <summary>A numerical comparison inside the named ring table, not ancestry.</summary>
    public bool DoubleFlipSlopeEqualsMinusTwoTimesFourThirds() =>
        Math.Abs(SlopeDoubleFlipped + 2.0 * (4.0 / 3.0)) < 1e-15;

    public static int[] AbsoluteThirdOrderDiagonal() => [8, -4, -4, 0];

    public static bool AbsoluteThirdOrderDiagonalSumsToZero() =>
        AbsoluteThirdOrderDiagonal().Sum() == 0;

    public F96BornSubdominantSlopesPi2Inheritance()
        : base(
            "F96 named N=4-ring slopes (-16/9,-16/9,-8/3) from M/U elements",
            Tier.Tier1Derived,
            "docs/proofs/PROOF_F96_BORN_SUBDOMINANT_SLOPES.md + " +
            "docs/ANALYTICAL_FORMULAS.md F96 + " +
            "simulations/born_rule_subdominant_dyson.py")
    {
    }

    public override string DisplayName => "F96 named N=4-ring subdominant slopes";

    public override string Summary =>
        $"pair (0,2): (-16/9,-16/9,-8/3)=({SlopeSingleFlipped:G17}," +
        $"{SlopeSingleFlipped:G17},{SlopeDoubleFlipped:G17}) from independent M/U elements; parentless.";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return InspectableNode.RealScalar("single-flip slope (-16/9)", SlopeSingleFlipped);
            yield return InspectableNode.RealScalar("double-flip slope (-8/3)", SlopeDoubleFlipped);
            yield return new InspectableNode(
                "Exact inputs",
                summary: $"single: M3={M3_SingleFlipped}, U2={U2_SingleFlipped_TimesFour}/4; " +
                         $"double: M5={M5_DoubleFlipped}, U4={U4_DoubleFlipped_TimesTwo}/2.");
            yield return new InspectableNode(
                "Trace control and topology fence",
                summary: "The absolute third-order diagonal (8,-4,-4,0) sums to zero. " +
                         "Relative deviations are not summed as probabilities. Chain |10> gives -4/3.");
        }
    }
}
