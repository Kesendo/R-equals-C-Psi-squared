using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F60 closed form (Tier 1, geometric corollary):
///
/// <code>
///   CΨ(0) for GHZ_N = 1 / (2^N − 1)
///
///   GHZ_N = (|0...0⟩ + |1...1⟩) / √2
///
///   Off-diagonal element:  ρ[0...0, 1...1] = 1/2  (and its conjugate)
///   Other off-diagonal:    0
///
///   For N ≥ 3: CΨ(0) &lt; 1/4 always; GHZ is born below the fold,
///              γ-independent, no dynamics can fix the geometric deficit.
/// </code>
///
/// <para>F60 is the first F-formula whose primary anchor lands DIRECTLY on the
/// <see cref="PolarityLayerOriginClaim"/>'s 0.5-shift: the off-diagonal
/// element ρ[0..0, 1..1] = 1/2 IS the ±0.5 polarity pair at d = 2 generated
/// by the 0.5-shift ρ = (I + r·σ)/2. F60's geometric structure literally
/// instantiates the polarity-layer reading.</para>
///
/// <para>Three Pi2-Foundation anchors:</para>
///
/// <list type="bullet">
///   <item><b>OffDiagonalElement = 1/2 = a_2</b>: GHZ's only nonzero
///         off-diagonal density-matrix entry is the ±0.5 polarity pair
///         literal (<see cref="PolarityLayerOriginClaim"/>). The "1/2"
///         appearing in F60's derivation is the same 1/2 that
///         <see cref="HalfAsStructuralFixedPointClaim"/> typed as the
///         structural fixed point.</item>
///   <item><b>HilbertSpaceDimension(N) = 2^N = a_{1−N}</b>: live from
///         <see cref="Pi2DyadicLadderClaim.Term"/>(<c>1 − N</c>). At N=2
///         this is a_{−1} = 4 (= d² for 1 qubit, the "small" qubit count
///         on the operator-space mirror).</item>
///   <item><b>FoldPosition = 1/4 = a_3</b>: the bilinear-apex maxval; F60's
///         "below the fold" assertion compares CΨ(0) against a_3 = 1/4.
///         Same anchor as F57 + Dicke (<see cref="QuarterAsBilinearMaxvalClaim"/>).</item>
/// </list>
///
/// <para>The structural reading: GHZ's geometric off-diagonal weight (= 1/2,
/// the polarity-pair) divided by the Hilbert-space-minus-1 (= 2^N − 1) gives
/// CΨ(0). The "− 1" subtraction is the Bloch-state normalisation. For N ≥ 3
/// the inverse 1/(2^N − 1) drops below 1/4 (i.e. below a_3), placing GHZ
/// outside the framework's quantum regime regardless of γ.</para>
///
/// <para>Operational consequence (per ANALYTICAL_FORMULAS): "the only escape
/// is to change the state." GHZ encoding is structurally unsuitable for
/// state transfer; W-type encodings are preferred (cf. README Rule 1, F22
/// XOR-drain confirms via independent argument).</para>
///
/// <para>Tier1Derived: F60 is Tier 1 geometric corollary. Closed form
/// derived directly from C(0) = 1 (pure state) + L1-off-diagonal coherence
/// = 1 + Ψ(0) = 1/(2^N − 1) = L1/(d−1). The Pi2-Foundation anchoring is
/// algebraic-trivial composition.</para>
///
/// <para>Anchors: <c>docs/ANALYTICAL_FORMULAS.md</c> F60 +
/// <c>experiments/DWELL_PREFACTOR_GENERALIZED.md</c> (Section 4) +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBaseClaims.cs</c>
/// (PolarityLayerOriginClaim, QuarterAsBilinearMaxvalClaim).</para></summary>
public sealed class F60GhzBornBelowFoldPi2Inheritance : Claim
{
    private readonly Pi2DyadicLadderClaim _ladder;

    /// <summary>±0.5 polarity-pair anchor — the typed parent that grounds
    /// F60's <c>OffDiagonalElement = 1/2</c> directly on the polarity-layer
    /// axis. Added 2026-05-16 as a typed ctor parent (previously registration-
    /// discard only).</summary>
    public PolarityLayerOriginClaim Polarity { get; }

    /// <summary>1/4 = (1/2)² bilinear-apex maxval — the typed parent that
    /// grounds F60's <c>FoldPosition</c> directly on the Quarter axis (same
    /// anchor as F57 + Dicke). Added 2026-05-16 as a typed ctor parent
    /// (previously registration-discard only).</summary>
    public QuarterAsBilinearMaxvalClaim Quarter { get; }

    /// <summary>(1/2, 1/4) argmax/maxval pair — the typed meta-anchor closing
    /// that F60 uses BOTH 1/2 (off-diagonal) and 1/4 (fold) simultaneously.
    /// Added 2026-05-16 as a typed ctor parent (previously registration-
    /// discard only).</summary>
    public ArgmaxMaxvalPairClaim ArgmaxMaxval { get; }

    /// <summary>The off-diagonal density-matrix element of GHZ_N:
    /// <c>ρ[0...0, 1...1] = 1/2</c>. Live from
    /// <see cref="Pi2DyadicLadderClaim.Term"/>(2) = <c>a_2</c> = ±0.5 polarity
    /// pair (<see cref="PolarityLayerOriginClaim"/>).</summary>
    public double OffDiagonalElement => _ladder.Term(2);

    /// <summary>The Hilbert-space dimension <c>2^N</c> for an N-qubit GHZ state.
    /// Live from <see cref="Pi2DyadicLadderClaim.Term"/>(<c>1 − N</c>) =
    /// <c>a_{1−N}</c>. At N=2 this is a_{−1} = 4 (= d² for 1 qubit, qubit-shift
    /// analog of operator-space).</summary>
    public double HilbertSpaceDimension(int N)
    {
        if (N < 2) throw new ArgumentOutOfRangeException(nameof(N), N, "F60 requires N ≥ 2.");
        return _ladder.Term(LadderIndexForHilbertSpace(N));
    }

    /// <summary>The Pi2 ladder index for the Hilbert-space dimension:
    /// <c>1 − N</c>. At N=2: -1; N=3: -2; N=5: -4; etc. Always negative for
    /// N ≥ 2 (operator-space side of the ladder).</summary>
    public int LadderIndexForHilbertSpace(int N) => 1 - N;

    /// <summary>The fold position <c>1/4 = a_3</c> on the dyadic ladder.
    /// Live from <see cref="Pi2DyadicLadderClaim.Term"/>(3); same anchor as
    /// F57's CrossingThreshold and DickeSuperpositionQuarter's QuarterCeiling
    /// (= <see cref="QuarterAsBilinearMaxvalClaim"/>).</summary>
    public double FoldPosition => _ladder.Term(3);

    /// <summary>Live closed form: <c>CΨ(0) = 1 / (2^N − 1)</c>. Throws for
    /// N &lt; 2.</summary>
    public double CPsiAtZeroForGhz(int N) => 1.0 / (HilbertSpaceDimension(N) - 1.0);

    /// <summary>True iff GHZ_N is born below the fold:
    /// <c>CΨ(0) &lt; 1/4 ⇔ N ≥ 3</c>. Geometric, γ-independent.</summary>
    public bool IsBornBelowFold(int N) => CPsiAtZeroForGhz(N) < FoldPosition;

    /// <summary>The smallest N at which GHZ is born below the fold:
    /// <c>N = 3</c>. Below this (i.e. at N=2) GHZ is the Bell+ state, which
    /// crosses the fold under dynamics. From N ≥ 3 the geometric deficit is
    /// permanent.</summary>
    public int SmallestNBelowFold => 3;

    /// <summary>Cross-check: at N = 2 (Bell+ as GHZ_2), CΨ(0) = 1/3 ≈ 0.333,
    /// above the fold a_3 = 1/4. Drift indicator on the "Bell+ is the only
    /// GHZ that crosses" reading.</summary>
    public bool BellPlusAboveFold() => CPsiAtZeroForGhz(2) > FoldPosition;

    public F60GhzBornBelowFoldPi2Inheritance(
        Pi2DyadicLadderClaim ladder,
        PolarityLayerOriginClaim polarity,
        QuarterAsBilinearMaxvalClaim quarter,
        ArgmaxMaxvalPairClaim argmaxMaxval)
        : base("F60 GHZ_N CΨ(0) = 1/(2^N - 1) inherits from Pi2-Foundation: 1/2 off-diagonal = ±0.5 polarity pair; 2^N = a_{1-N}; fold = a_3",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F60 + " +
               "experiments/DWELL_PREFACTOR_GENERALIZED.md (Section 4) + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBaseClaims.cs (PolarityLayerOrigin, QuarterAsBilinearMaxval, ArgmaxMaxvalPair)")
    {
        _ladder = ladder ?? throw new ArgumentNullException(nameof(ladder));
        Polarity = polarity ?? throw new ArgumentNullException(nameof(polarity));
        Quarter = quarter ?? throw new ArgumentNullException(nameof(quarter));
        ArgmaxMaxval = argmaxMaxval ?? throw new ArgumentNullException(nameof(argmaxMaxval));
    }

    public override string DisplayName =>
        "F60 GHZ_N born below the fold as Pi2-Foundation polarity-layer inheritance";

    public override string Summary =>
        $"CΨ(0)_GHZ_N = 1/(2^N − 1): off-diagonal = a_2 = 1/2 (±0.5 polarity pair); 2^N = a_{{1−N}}; fold = a_3 = 1/4; " +
        $"born below fold for N ≥ 3 (γ-independent geometric deficit) ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("F60 closed form",
                summary: "CΨ(0)_GHZ_N = 1/(2^N − 1); Tier 1 geometric corollary; γ-independent");
            yield return new InspectableNode("Pi2-Foundation anchoring",
                summary: "OffDiagonalElement = a_2 = 1/2 (PolarityLayerOrigin's ±0.5 pair); HilbertSpaceDim = a_{1−N}; FoldPosition = a_3 = 1/4");
            yield return InspectableNode.RealScalar("OffDiagonalElement (= a_2 = 1/2)", OffDiagonalElement);
            yield return InspectableNode.RealScalar("FoldPosition (= a_3 = 1/4)", FoldPosition);
            yield return new InspectableNode("polarity-layer reading",
                summary: "GHZ's only nonzero off-diagonal entry IS the polarity pair literal; F60 is the first F-formula whose primary anchor sits ON the 0.5-shift axis (per Tom 2026-05-09 mirror-map check)");
            yield return new InspectableNode("operational consequence",
                summary: "for N ≥ 3 GHZ is structurally unsuitable for state transfer; γ-reduction cannot fix the geometric deficit; only escape = change the state (W-type encodings preferred per F22 + Rule 1)");
            yield return new InspectableNode("F60 ↔ F62 sibling",
                summary: "F60 (GHZ): pair-CΨ(0) = 1/(2^N − 1) below fold for N ≥ 3; F62 (W_N): pair-CΨ(0) = 10/81 ≈ 0.124 also below fold at N=3; F69 (GHZ+W mix): unique optimum above 1/4 via sextic root");
            // Verified table from ANALYTICAL_FORMULAS F60
            for (int N = 2; N <= 6; N++)
            {
                double cpsi = CPsiAtZeroForGhz(N);
                yield return new InspectableNode(
                    $"N={N}",
                    summary: $"CΨ(0) = 1/(2^{N} − 1) = 1/{HilbertSpaceDimension(N) - 1.0:G6} = {cpsi:G6}; " +
                             $"above fold: {cpsi > FoldPosition} ({(cpsi > FoldPosition ? "Bell+ regime" : "born classical")})");
            }
        }
    }
}
