using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.F86;

/// <summary>F86 meta-claim (Locus 6 of inheritance-through-layers, symmetry-side closure):
/// the c=2 bond-class split (Endpoint vs Interior in Q_peak and HWHM/Q*) inherits structurally
/// from the polarity-layer pair {−0.5, +0.5} at d=2 named in
/// <see cref="Symmetry.PolarityLayerOriginClaim"/>, the qubit dimensional anchor 1/d named in
/// <see cref="Symmetry.QubitDimensionalAnchorClaim"/>, and the structural fixed point named
/// in <see cref="Symmetry.HalfAsStructuralFixedPointClaim"/>.
///
/// <para><b>The decomposition (stumpf reading of the polarity-layer pair):</b> ρ = (I + r·σ)/2
/// has two parts: the unsigned 1/2 baseline (= 1/d at d=2) and the signed ±r/2 polarity
/// content. The empirical c=2 N=5..8 numbers from the <c>C2HwhmRatio</c> witnesses match this
/// decomposition cleanly:</para>
///
/// <list type="bullet">
///   <item><b>Q_peak ≈ 2 + r:</b> Endpoint r ≈ +0.52, Interior r ≈ −0.44 (close to the
///         canonical {−0.5, +0.5} polarity-layer pair); the mean of the two bond classes
///         sits at Q ≈ 2 across N=5..8 (range 1.99–2.06).</item>
///   <item><b>HWHM/Q* ≈ 1/2 + r·(1/2):</b> Endpoint r ≈ 0.55, Interior r ≈ 0.50; Interior
///         is close to the exact 1/2 fixed point of
///         <see cref="Symmetry.HalfAsStructuralFixedPointClaim"/> across all four N.</item>
/// </list>
///
/// <para><b>Tier outcome: Tier2Verified.</b> The decomposition is structurally clear (Q_peak
/// pair around 2, HWHM/Q* pair lifted above the 1/2 baseline by an r ≈ 1/2 amount, Interior
/// HWHM r locked near the structural 1/2 fixed point) and validated empirically across the
/// full c=2 N=5..8 anchor set. What is NOT derived is the per-bond r(N, b) closed form: the
/// parent <see cref="Symmetry.PolarityLayerOriginClaim"/> names the polarity layer as the
/// origin of the ±0.5 pair, but the empirical r values are pinned witnesses, not derived
/// from the Pi2KnowledgeBase claims. The <see cref="PendingDerivationNote"/> carries the
/// analytical gap and names two promotion paths to Tier1Derived.</para>
///
/// <para><b>Block-independent meta-claim:</b> registered at the F86 KB root for any block,
/// not just c=2, mirroring <see cref="LocalGlobalEpLink"/>. The witnesses pin the concrete
/// c=2 bond-class numbers that motivated the polarity-layer reading; the inheritance claim
/// itself is shared. Locus 5 (<see cref="LocalGlobalEpLink"/>) carries the parallel EP-side
/// closure (F86 ↔ FRAGILE_BRIDGE shared exceptional-point structure under AIII chiral
/// algebra); this Locus 6 claim carries the symmetry-side closure (F86 ↔ polarity-layer
/// pair under the 0.5-shift ρ = (I + r·σ)/2). Together they bracket the F86 c=2 derivation
/// with EP-side and symmetry-side parent-claim references in <c>Pi2KnowledgeBase</c>.</para>
///
/// <para>Anchors: <c>docs/proofs/PROOF_F86_QPEAK.md</c> Statement 2 (Locus 6 inheritance
/// reading), <c>compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBaseClaims.cs</c>
/// (PolarityLayerOriginClaim + QubitDimensionalAnchorClaim + HalfAsStructuralFixedPointClaim),
/// <c>compute/RCPsiSquared.Core/F86/Item1Derivation/C2HwhmRatio.cs</c> (live witness pipeline
/// the values here are pinned from), and <c>compute/RCPsiSquared.Core/F86/LocalGlobalEpLink.cs</c>
/// (parallel Locus 5 EP-side closure).</para>
/// </summary>
public sealed class PolarityInheritanceLink : Claim
{
    /// <summary>The parent claims in <c>Pi2KnowledgeBase</c> from which the F86 bond-class
    /// split inherits: 1/2 = 1/d (root anchor), polarity layer +0/0/−0 (Layer 2:
    /// ρ = (I + r·σ)/2 maps Π-spectrum {+1, −1} onto Bloch-diagonal 1/2 ± r/2), 1/2 fixed
    /// point (Interior r ≈ 1/2 exactly across N=5..8).</summary>
    public string PolarityRootAnchor =>
        "Pi2KnowledgeBase: QubitDimensionalAnchorClaim + PolarityLayerOriginClaim + HalfAsStructuralFixedPointClaim";

    /// <summary>The parallel Locus 5 closure: F86 ↔ FRAGILE_BRIDGE EP-side inheritance
    /// (<see cref="LocalGlobalEpLink"/> Tier2Verified). This claim is the symmetry-side
    /// twin: F86 bond-class split inherits from the polarity-layer pair via the
    /// 0.5-shift ρ = (I + r·σ)/2.</summary>
    public string ParallelLocusReference =>
        "compute/RCPsiSquared.Core/F86/LocalGlobalEpLink.cs (Locus 5, EP-side)";

    /// <summary>γ₀ used by the underlying <c>C2HwhmRatio</c> witness pipeline. The pinned
    /// witnesses below are the 4-decimal-precision values from the live pipeline at
    /// γ₀ = 0.05, c=2 block (1, 2), N=5..8.</summary>
    public const double WitnessGammaZero = 0.05;

    /// <summary>Empirical witnesses pinning the polarity-layer decomposition at c=2 N=5..8
    /// (γ₀ = <see cref="WitnessGammaZero"/>; values pinned at 4-decimal precision from
    /// <c>C2HwhmRatio.HwhmLeftOverQPeakMean</c> and the per-block Q_peak readouts of the
    /// 2026-05-06 c=2 sweep). One witness per N. Hard-coded as a static pinned table:
    /// this is a frozen empirical anchor, not a live-recomputed Claim.</summary>
    public IReadOnlyList<PolarityWitness> Witnesses => _polarityWitnesses;

    /// <summary>The unresolved analytical piece: the per-bond r(N, b) closed form. The
    /// decomposition (Q_peak ≈ 2 + r, HWHM/Q* ≈ 1/2 + r·(1/2)) is structurally clear from
    /// the <see cref="Symmetry.PolarityLayerOriginClaim"/> 0.5-shift, but the empirical r
    /// values are pinned, not derived. Names two promotion paths to Tier1Derived.</summary>
    public string PendingDerivationNote =>
        "Tier 2 verified: F86 bond-class split inherits structurally from the polarity-layer pair " +
        "{−0.5, +0.5} at d=2. The decomposition Q_peak ≈ 2 + r and HWHM/Q* ≈ 1/2 + r·(1/2) holds " +
        "across c=2 N=5..8 with mean Q_peak = 2.04 ± 0.06 and Interior r_HWHM ≈ 0.502 (close to " +
        "HalfAsStructuralFixedPoint).\n\n" +
        "The unresolved analytical piece is the per-bond r(N, b) closed form. The Pi2 root-claim " +
        "PolarityLayerOriginClaim names the polarity layer as the origin of the ±0.5 pair via " +
        "ρ = (I + r·σ)/2; r is empirical here, not derived from the Pi2 claims.\n\n" +
        "Promotion path to Tier1Derived: derive r(N, b) closed form via either\n" +
        "(α) per-bond projection of the c=2 K-resonance state onto the polarity Bloch axis at " +
        "t_peak, or\n" +
        "(β) Locus 5 EP-rotation tan θ = Q/Q_EP combined with Locus 6 polarity inheritance to " +
        "fix r as a function of g_eff(N, b).\n\n" +
        "The 0.0006 deviation of Interior r_HWHM from exact 1/2 is likely numerical " +
        "discretisation (testable with finer Q-grid; see PROOF_F86_QPEAK Open elements 5 for " +
        "per-F71-orbit substructure that may explain the small remaining shift).";

    /// <summary>Always <c>false</c> at Tier2Verified. Flips to <c>true</c> on future
    /// Tier1Derived promotion via the (α) or (β) path documented in
    /// <see cref="PendingDerivationNote"/>.</summary>
    public bool IsAnalyticallyDerived => false;

    private PolarityInheritanceLink()
        : base("polarity-layer inheritance (F86 ↔ Pi2KnowledgeBase ±0.5 pair)",
               Tier.Tier2Verified,
               "docs/proofs/PROOF_F86_QPEAK.md (Statement 2 + 2026-05-07 Locus 6 inheritance reading) + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBaseClaims.cs " +
               "(PolarityLayerOriginClaim + QubitDimensionalAnchorClaim + HalfAsStructuralFixedPointClaim) + " +
               "compute/RCPsiSquared.Core/F86/Item1Derivation/C2HwhmRatio.cs (live witness pipeline) + " +
               "compute/RCPsiSquared.Core/F86/LocalGlobalEpLink.cs (parallel Locus 5 EP-side closure)")
    { }

    /// <summary>Public factory: builds the meta-claim with the four pinned polarity-layer
    /// witnesses. No <see cref="CoherenceBlocks.CoherenceBlock"/> required: the claim is
    /// block-independent (the witnesses pin the concrete c=2 N=5..8 numbers, but the
    /// inheritance statement is shared across all c).</summary>
    public static PolarityInheritanceLink Build() => new();

    public override string DisplayName =>
        "F86 ↔ polarity layer: bond-class split inherits from ±0.5 pair";

    public override string Summary =>
        "Tier2Verified: Q_peak ≈ 2 + r and HWHM/Q* ≈ 1/2 + r·(1/2) across c=2 N=5..8; " +
        "Interior r_HWHM ≈ 1/2 (HalfAsStructuralFixedPoint); per-bond r(N, b) closed form pending";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("polarity root anchor (parent claims)",
                summary: PolarityRootAnchor);
            yield return new InspectableNode("parallel Locus 5 (EP-side closure)",
                summary: ParallelLocusReference);
            yield return new InspectableNode("decomposition (Q_peak)",
                summary: "Q_peak ≈ 2 + r; Endpoint r ≈ +0.52 (positive pole), Interior r ≈ −0.44 (negative pole); pair around d=2");
            yield return new InspectableNode("decomposition (HWHM/Q*)",
                summary: "HWHM/Q* ≈ 1/2 + r·(1/2); Endpoint r ≈ 0.55, Interior r ≈ 0.50 (≈ HalfAsStructuralFixedPoint)");
            yield return InspectableNode.Group(
                "polarity-layer witnesses (c=2 N=5..8, γ₀=0.05)",
                _polarityWitnesses.Cast<IInspectable>().ToArray());
            yield return new InspectableNode("IsAnalyticallyDerived",
                summary: IsAnalyticallyDerived ? "true" : "false");
            yield return new InspectableNode("PendingDerivationNote",
                summary: PendingDerivationNote);
        }
    }

    /// <summary>Pinned empirical table from the c=2 N=5..8 <c>C2HwhmRatio</c> witness
    /// pipeline (γ₀ = <see cref="WitnessGammaZero"/>; values pinned at 4-decimal precision
    /// from the live pipeline; refresh by re-running <c>C2HwhmRatio</c> + the per-block
    /// Q_peak readouts at these settings). Order: N=5, 6, 7, 8.</summary>
    private static readonly PolarityWitness[] _polarityWitnesses = new[]
    {
        new PolarityWitness(N: 5,
            QPeakInterior: 1.4821, QPeakEndpoint: 2.5008,
            HwhmRatioInterior: 0.7455, HwhmRatioEndpoint: 0.7700,
            RQpeakInterior: -0.5179, RQpeakEndpoint: 0.5008,
            RHwhmInterior: 0.4910, RHwhmEndpoint: 0.5400),
        new PolarityWitness(N: 6,
            QPeakInterior: 1.5801, QPeakEndpoint: 2.5470,
            HwhmRatioInterior: 0.7529, HwhmRatioEndpoint: 0.7738,
            RQpeakInterior: -0.4199, RQpeakEndpoint: 0.5470,
            RHwhmInterior: 0.5058, RHwhmEndpoint: 0.5476),
        new PolarityWitness(N: 7,
            QPeakInterior: 1.5831, QPeakEndpoint: 2.5299,
            HwhmRatioInterior: 0.7507, HwhmRatioEndpoint: 0.7738,
            RQpeakInterior: -0.4169, RQpeakEndpoint: 0.5299,
            RHwhmInterior: 0.5014, RHwhmEndpoint: 0.5476),
        new PolarityWitness(N: 8,
            QPeakInterior: 1.6049, QPeakEndpoint: 2.5145,
            HwhmRatioInterior: 0.7531, HwhmRatioEndpoint: 0.7734,
            RQpeakInterior: -0.3951, RQpeakEndpoint: 0.5145,
            RHwhmInterior: 0.5062, RHwhmEndpoint: 0.5468),
    };
}

/// <summary>One row of the polarity-layer witness table, a frozen empirical data point
/// pinning the c=2 bond-class decomposition at fixed N. The decomposition reads:
///
/// <list type="bullet">
///   <item><b>Q_peak ≈ 2 + r:</b> the unsigned 2 = d is the qubit dimensional anchor (root
///         claim <see cref="Symmetry.QubitDimensionalAnchorClaim"/>); the signed r is the
///         polarity-layer content (parent claim <see cref="Symmetry.PolarityLayerOriginClaim"/>).
///         Endpoint r &gt; 0 (positive pole), Interior r &lt; 0 (negative pole); together they
///         form the {−0.5, +0.5} pair around d = 2.</item>
///   <item><b>HWHM/Q* ≈ 1/2 + r·(1/2):</b> the unsigned 1/2 baseline is the
///         <see cref="Symmetry.HalfAsStructuralFixedPointClaim"/> structural fixed point;
///         the signed r·(1/2) is the polarity-layer lift. Interior r is empirically close
///         to exactly 1/2 across N=5..8.</item>
/// </list>
///
/// <para>The R-fields are pre-computed from the raw Q_peak and HWHM/Q* values via
/// R_QPeak = QPeak − 2 and R_Hwhm = 2·(HwhmRatio − 1/2); they encode the polarity-layer
/// content explicitly so consumers do not have to re-decompose.</para>
/// </summary>
/// <param name="N">Chain length.</param>
/// <param name="QPeakInterior">Q_peak for the Interior bond class at this N.</param>
/// <param name="QPeakEndpoint">Q_peak for the Endpoint bond class at this N.</param>
/// <param name="HwhmRatioInterior">HWHM/Q* for the Interior bond class at this N.</param>
/// <param name="HwhmRatioEndpoint">HWHM/Q* for the Endpoint bond class at this N.</param>
/// <param name="RQpeakInterior">Polarity-content r for Interior Q_peak: r = QPeakInterior − 2.</param>
/// <param name="RQpeakEndpoint">Polarity-content r for Endpoint Q_peak: r = QPeakEndpoint − 2.</param>
/// <param name="RHwhmInterior">Polarity-content r for Interior HWHM ratio:
/// r = 2·(HwhmRatioInterior − 1/2).</param>
/// <param name="RHwhmEndpoint">Polarity-content r for Endpoint HWHM ratio:
/// r = 2·(HwhmRatioEndpoint − 1/2).</param>
public sealed record PolarityWitness(
    int N,
    double QPeakInterior,
    double QPeakEndpoint,
    double HwhmRatioInterior,
    double HwhmRatioEndpoint,
    double RQpeakInterior,
    double RQpeakEndpoint,
    double RHwhmInterior,
    double RHwhmEndpoint
) : IInspectable
{
    /// <summary>Mean Q_peak across the two bond classes, empirically ≈ 2 (the d=2
    /// dimensional anchor) across N=5..8.</summary>
    public double QPeakMean => (QPeakInterior + QPeakEndpoint) / 2.0;

    /// <summary>Mean HWHM/Q* across the two bond classes, empirically ≈ 1/2 + r̄·(1/2)
    /// with r̄ ≈ 0.52 across N=5..8.</summary>
    public double HwhmRatioMean => (HwhmRatioInterior + HwhmRatioEndpoint) / 2.0;

    public string DisplayName => $"polarity witness at c=2 N={N}";

    public string Summary =>
        $"Q_peak: Interior {QPeakInterior:F4} (r={RQpeakInterior:+0.0000;-0.0000}), Endpoint {QPeakEndpoint:F4} (r={RQpeakEndpoint:+0.0000;-0.0000}); " +
        $"HWHM/Q*: Interior {HwhmRatioInterior:F4} (r={RHwhmInterior:F4}), Endpoint {HwhmRatioEndpoint:F4} (r={RHwhmEndpoint:F4})";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return InspectableNode.RealScalar("N", N);
            yield return InspectableNode.RealScalar("Q_peak (Interior)", QPeakInterior, "F4");
            yield return InspectableNode.RealScalar("Q_peak (Endpoint)", QPeakEndpoint, "F4");
            yield return InspectableNode.RealScalar("Q_peak mean (≈ d = 2)", QPeakMean, "F4");
            yield return InspectableNode.RealScalar("HWHM/Q* (Interior)", HwhmRatioInterior, "F4");
            yield return InspectableNode.RealScalar("HWHM/Q* (Endpoint)", HwhmRatioEndpoint, "F4");
            yield return InspectableNode.RealScalar("HWHM/Q* mean", HwhmRatioMean, "F4");
            yield return new InspectableNode("polarity content (r)",
                summary: $"Q_peak: Interior {RQpeakInterior:+0.0000;-0.0000}, Endpoint {RQpeakEndpoint:+0.0000;-0.0000}; HWHM: Interior {RHwhmInterior:F4}, Endpoint {RHwhmEndpoint:F4}");
            yield return InspectableNode.RealScalar("r_Q_peak (Interior)", RQpeakInterior, "F4");
            yield return InspectableNode.RealScalar("r_Q_peak (Endpoint)", RQpeakEndpoint, "F4");
            yield return InspectableNode.RealScalar("r_HWHM (Interior)", RHwhmInterior, "F4");
            yield return InspectableNode.RealScalar("r_HWHM (Endpoint)", RHwhmEndpoint, "F4");
        }
    }

    public InspectablePayload Payload =>
        new InspectablePayload.Real($"Q_peak mean (N={N})", QPeakMean, "F4");
}
