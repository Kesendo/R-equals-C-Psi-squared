using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>
/// F94 is the leading dominant-outcome Born-deviation term for one named
/// calculation: the N=4 Heisenberg ring, initial state |0+0+>, retained pair
/// (0,2), outcome |00>, and local Z dephasing.
///
/// The exact Dyson/Pauli reduction gives sym3 = 8, hence
/// Delta = (8/3!) J^2 gamma t^3 = (4/3) Q^2 K^3. The higher-order remainder
/// is not assigned a universal monomial. The coefficient is owned by this
/// calculation; same-number dyadic anchors are comparisons, not ancestry.
/// </summary>
public sealed class F94BornDeviationFourThirdsPi2Inheritance : Claim
{
    public const int N = 4;
    public const int RetainedSiteA = 0;
    public const int RetainedSiteB = 2;
    public const int Sym3PartialTraceInteger = 8;
    public const int TaylorThreeFactorial = 6;
    public const int SurvivingDysonDiagrams = 32;
    public const int RawPauliPerDiagram = 4;
    public const int CellA_Ord1XX_AdjKeptSide = 8;
    public const int CellB_Ord2XX_SelfOrAdjKeptSide = 16;
    public const int CellC_Ord2YY_Self = 8;

    public double Coefficient =>
        (double)Sym3PartialTraceInteger / TaylorThreeFactorial;

    public double DeltaDominant(double Q, double K)
    {
        if (!double.IsFinite(Q) || Q < 0.0)
            throw new ArgumentOutOfRangeException(nameof(Q), Q, "Q must be finite and >= 0.");
        if (!double.IsFinite(K) || K < 0.0)
            throw new ArgumentOutOfRangeException(nameof(K), K, "K must be finite and >= 0.");
        return Coefficient * Q * Q * K * K * K;
    }

    public double DeltaP_Dominant(double J, double gamma, double t)
    {
        if (!double.IsFinite(J) || J < 0.0)
            throw new ArgumentOutOfRangeException(nameof(J), J, "J must be finite and >= 0.");
        if (!double.IsFinite(gamma) || gamma < 0.0)
            throw new ArgumentOutOfRangeException(nameof(gamma), gamma, "gamma must be finite and >= 0.");
        if (!double.IsFinite(t) || t < 0.0)
            throw new ArgumentOutOfRangeException(nameof(t), t, "t must be finite and >= 0.");
        return Coefficient * J * J * gamma * t * t * t;
    }

    /// <summary>The leading-order ratio 1 + Delta, not a Born-rule derivation.</summary>
    public double C_DominantOutcome(double Q, double K) =>
        1.0 + DeltaDominant(Q, K);

    public bool CoefficientAgreesWithSym3() =>
        Coefficient == (double)Sym3PartialTraceInteger / TaylorThreeFactorial;

    public bool CellCountsSumToSurvivingDiagrams() =>
        CellA_Ord1XX_AdjKeptSide + CellB_Ord2XX_SelfOrAdjKeptSide + CellC_Ord2YY_Self
        == SurvivingDysonDiagrams;

    public bool StructuralDecompositionRecoversSym3() =>
        (long)SurvivingDysonDiagrams * RawPauliPerDiagram
        == (long)Sym3PartialTraceInteger * 16;

    public F94BornDeviationFourThirdsPi2Inheritance()
        : base(
            "F94 named N=4-ring leading Born deviation: Delta_|00> = (8/6) Q^2 K^3 = (4/3) Q^2 K^3",
            Tier.Tier1Derived,
            "docs/proofs/PROOF_F94_BORN_DOMINANT_FOUR_THIRDS.md + " +
            "docs/ANALYTICAL_FORMULAS.md F94 + " +
            "simulations/born_rule_tier1_derivation.py + " +
            "simulations/born_rule_delta_dominant_coefficient.py")
    {
    }

    public override string DisplayName =>
        "F94 named N=4-ring dominant-outcome coefficient 8/6 = 4/3";

    public override string Summary =>
        $"|0+0+>, pair (0,2), |00>: leading Delta = ({Coefficient:G17}) Q^2 K^3; " +
        "exact analytic sym3 integer 8 over Taylor 3!, with unspecified higher order; parentless.";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return InspectableNode.RealScalar("Coefficient (= 8/6 = 4/3)", Coefficient);
            yield return InspectableNode.RealScalar("Sym3PartialTraceInteger", Sym3PartialTraceInteger);
            yield return InspectableNode.RealScalar("TaylorThreeFactorial", TaylorThreeFactorial);
            yield return new InspectableNode(
                "Named setup and scope",
                summary: "N=4 Heisenberg ring; |0+0+>; retained pair (0,2); |00>; local Z dephasing. " +
                         "The coefficient is exact; linked floating producers are corroborations.");
            yield return new InspectableNode(
                "Structural decomposition",
                summary: $"{CellA_Ord1XX_AdjKeptSide}+{CellB_Ord2XX_SelfOrAdjKeptSide}+" +
                         $"{CellC_Ord2YY_Self}={SurvivingDysonDiagrams} diagrams and " +
                         $"{SurvivingDysonDiagrams}*{RawPauliPerDiagram}/16={Sym3PartialTraceInteger}.");
        }
    }
}
