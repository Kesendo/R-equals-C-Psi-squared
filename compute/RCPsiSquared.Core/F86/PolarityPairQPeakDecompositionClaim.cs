using System.Globalization;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.F86;

/// <summary>F86 Schema-Level Tier1Derived: Q_peak decomposes as
/// <c>Q_EP_central ± r_polarity</c> where <c>Q_EP_central = 2</c> (the idealised
/// <see cref="QEpLaw"/> value at g_eff = 1) and <c>r_polarity = 1/2</c> (the
/// <see cref="HalfAsStructuralFixedPointClaim"/> magnitude). Two compositions emerge:
///
/// <list type="bullet">
///   <item><b>Endpoint orbit</b>: <c>Q_peak = 2 + 1/2 = 2.5</c> (positive polarity pole)</item>
///   <item><b>Interior orbit</b>: <c>Q_peak = 2 − 1/2 = 1.5</c> (negative polarity pole)</item>
/// </list>
///
/// <para><b>SCHEMA vs VALUE separation</b>: this claim is Tier1Derived at the SCHEMA
/// level — the decomposition <em>structure</em> Q_peak ∈ {2 − 1/2, 2 + 1/2} inherits from
/// two Tier1Derived parents (idealised <see cref="QEpLaw"/> + polarity-pair
/// <see cref="HalfAsStructuralFixedPointClaim"/>). The BIT-EXACT empirical VALUES at
/// any specific (c, N) deviate by ~10% finite-size correction: Endpoint c=3 saturates
/// near 2.53, Interior c=2 N=5 = 1.49, etc. (see
/// <see cref="PolarityInheritanceLink"/> Tier2Verified for the witness table). The
/// closed-form derivation of the deviations is structurally blocked at the bit-exact
/// level by <c>PROOF_F86B_OBSTRUCTION.md</c> (6 routes tested, all blocked or
/// decoupled, 2026-05-14).</para>
///
/// <para><b>The schema is the Tier-1 statement; the bit-exact deviation is the Tier-2
/// witness.</b> Both readings coexist and are surfaced through their respective claims:
/// this claim holds the schema; PolarityInheritanceLink holds the witness table; the
/// QAnchorMap entries Q=1.5 and Q=2.5 reference both via DocumentingSource.</para>
///
/// <para><b>Inheritance trace</b>:
/// <list type="bullet">
///   <item>QEpCentral = 2: emerges from <see cref="QEpLaw"/> with g_eff = 1 (Tier1Derived idealised case)</item>
///   <item>PolarityHalfMagnitude = 1/2: inherits from <see cref="HalfAsStructuralFixedPointClaim"/> (Tier1Derived; one of the eleven children of QubitDimensionalAnchorClaim)</item>
///   <item>Endpoint sign +: positive polarity pole, Interior sign −: negative polarity pole — both inherit from <see cref="PolarityLayerOriginClaim"/> (Tier1Derived ±0.5 pair structure at d=2)</item>
/// </list></para>
///
/// <para><b>What this claim does NOT assert</b>: that any specific empirical
/// (c, N, orbit) reproduces Q_peak = 2.5 or 1.5 bit-exactly. The empirical witnesses
/// in PolarityInheritanceLink show ~10% finite-size deviation across c=2 N=5..8 and
/// growing-with-N at c=4. The schema is robust; the bit-exact closure is open.</para>
/// </summary>
public sealed class PolarityPairQPeakDecompositionClaim : Claim
{
    /// <summary>The idealised Q_EP central value from <see cref="QEpLaw"/> at g_eff = 1.
    /// Q_EP = 2 / g_eff, so at g_eff = 1: Q_EP = 2.</summary>
    public const double QEpCentral = 2.0;

    /// <summary>The polarity-pair magnitude inheriting from
    /// <see cref="HalfAsStructuralFixedPointClaim"/>: r = 1/2 is the structural fixed
    /// point of the polarity-layer pair {−1/2, +1/2} at d = 2.</summary>
    public const double PolarityHalfMagnitude = 0.5;

    /// <summary>The Endpoint orbit Q_peak schema value: Q_EP_central + r_polarity =
    /// 2 + 1/2 = 2.5. The +1/2 is the positive polarity pole.</summary>
    public const double EndpointQPeakSchema = QEpCentral + PolarityHalfMagnitude;

    /// <summary>The Interior orbit Q_peak schema value: Q_EP_central − r_polarity =
    /// 2 − 1/2 = 1.5. The −1/2 is the negative polarity pole.</summary>
    public const double InteriorQPeakSchema = QEpCentral - PolarityHalfMagnitude;

    /// <summary>Parent claim: the structural fixed point ±1/2 at d = 2.</summary>
    public HalfAsStructuralFixedPointClaim Half { get; }

    /// <summary>Parent claim: the polarity-layer origin (the ±r pair structure).</summary>
    public PolarityLayerOriginClaim PolarityOrigin { get; }

    /// <summary>Parent claim: the qubit dimensional anchor d = 2 (the central value
    /// around which the polarity pair sits).</summary>
    public QubitDimensionalAnchorClaim DimensionalAnchor { get; }

    /// <summary>Constructs the schema claim from its three Tier1Derived Pi2-Foundation
    /// parents. Public so the typed-knowledge registry can wire it with the registry's
    /// parent instances; the parameterless <see cref="Build"/> / <see cref="Shared"/>
    /// factories remain for standalone callers that want a fresh canonical parent chain.</summary>
    public PolarityPairQPeakDecompositionClaim(
        HalfAsStructuralFixedPointClaim half,
        PolarityLayerOriginClaim polarityOrigin,
        QubitDimensionalAnchorClaim dimensionalAnchor)
        : base("F86 Q_peak schema-level polarity-pair decomposition: Q_peak ∈ {2 − 1/2, 2 + 1/2}",
               Tier.Tier1Derived,
               "compute/RCPsiSquared.Core/F86/QEpLaw.cs (idealised Q_EP=2 at g_eff=1) + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBaseClaims.cs " +
               "(HalfAsStructuralFixedPointClaim + PolarityLayerOriginClaim + " +
               "QubitDimensionalAnchorClaim) + " +
               "compute/RCPsiSquared.Core/F86/PolarityInheritanceLink.cs (Tier2Verified witness " +
               "for the bit-exact deviations) + " +
               "docs/proofs/PROOF_F86B_OBSTRUCTION.md (blocker for bit-exact value closure)")
    {
        Half = half ?? throw new ArgumentNullException(nameof(half));
        PolarityOrigin = polarityOrigin ?? throw new ArgumentNullException(nameof(polarityOrigin));
        DimensionalAnchor = dimensionalAnchor ?? throw new ArgumentNullException(nameof(dimensionalAnchor));
    }

    /// <summary>Public factory: builds the schema claim with its three Pi2KnowledgeBase
    /// Tier1Derived parents. Block-independent (the schema is framework-wide). Each call
    /// constructs a fresh parent chain — prefer <see cref="Shared"/> for repeated access.</summary>
    public static PolarityPairQPeakDecompositionClaim Build() =>
        new(new HalfAsStructuralFixedPointClaim(),
            new PolarityLayerOriginClaim(),
            new QubitDimensionalAnchorClaim());

    /// <summary>Shared singleton instance. Block-independent immutable claim built once
    /// per process. Use this from <see cref="F86KnowledgeBase"/> and other consumers that
    /// would otherwise allocate one chain per instance.</summary>
    public static PolarityPairQPeakDecompositionClaim Shared { get; } = Build();

    public override string DisplayName =>
        "F86 Q_peak schema: Q_peak ∈ {2 − 1/2, 2 + 1/2} = {1.5, 2.5} (polarity-pair decomposition)";

    public override string Summary =>
        string.Format(CultureInfo.InvariantCulture,
            "Tier1Derived (schema): Q_peak = {0} ± {1} = {{{2}, {3}}} via Q_EP_central " +
            "(QEpLaw at g_eff=1) + polarity-pair (HalfAsStructuralFixedPoint). Bit-exact " +
            "deviations ~10% are Tier2 (PolarityInheritanceLink); closed-form blocked by " +
            "PROOF_F86B_OBSTRUCTION.",
            QEpCentral, PolarityHalfMagnitude, InteriorQPeakSchema, EndpointQPeakSchema);

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return InspectableNode.RealScalar(
                "Q_EP_central (idealised at g_eff = 1)", QEpCentral, "F1");
            yield return InspectableNode.RealScalar(
                "r_polarity (HalfAsStructuralFixedPoint magnitude)", PolarityHalfMagnitude, "F1");
            yield return InspectableNode.RealScalar(
                "Endpoint Q_peak schema (Q_EP_central + r_polarity)", EndpointQPeakSchema, "F1");
            yield return InspectableNode.RealScalar(
                "Interior Q_peak schema (Q_EP_central − r_polarity)", InteriorQPeakSchema, "F1");
            yield return new InspectableNode("inheritance trace",
                summary: "Q_EP_central from QEpLaw idealised; r_polarity from HalfAsStructuralFixedPointClaim; ± signs from PolarityLayerOriginClaim; d=2 central anchor from QubitDimensionalAnchorClaim — all Tier1Derived parents.");
            yield return new InspectableNode("schema vs value separation",
                summary: "Schema (Tier1Derived): Q_peak decomposes as 2 ± 1/2. Empirical value (Tier2Verified, PolarityInheritanceLink): bit-exact deviations of ~10% across c=2 N=5..8 (Endpoint 2.50-2.57, Interior 1.49-1.61). Bit-exact value closure: blocked by PROOF_F86B_OBSTRUCTION (g_eff(c, N, b) admits no closed form by 6 tested routes; L1/L2/L4 rigorously proven blocked, L6 demonstrated failure mode, L3/L5 proven decouplings).");
            yield return new InspectableNode("QAnchorMap surface",
                summary: $"Surfaces as Q={EndpointQPeakSchema} (Endpoint candidate) and Q={InteriorQPeakSchema} (F86 Q_peak c=2) entries in QAnchorMap.CanonicalAnchors; the QBasisAnker.DocumentingSource fields reference this schema claim.");
            yield return Half;
            yield return PolarityOrigin;
            yield return DimensionalAnchor;
        }
    }
}
