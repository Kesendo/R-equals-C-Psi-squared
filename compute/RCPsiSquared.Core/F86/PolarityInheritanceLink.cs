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
/// <para><b>Tier outcome: Tier2Verified.</b> The decomposition is structurally clear
/// (Q_peak pair around 2, HWHM/Q* pair lifted above the 1/2 baseline by an r ≈ 1/2 amount,
/// Interior HWHM r locked near the structural 1/2 fixed point) and validated empirically
/// across the full c=2 N=5..8 anchor set. The 2026-05-07 Direction (α) attempt extracted
/// an algebraic composition-reading of r_Q via two existing Tier 1 facts. The 2026-05-08
/// code review reclassified Direction (α)'s "no sign split in the polarity-Bloch projection"
/// finding as <b>structurally tautological under the uniform-J 4-mode reduction</b> rather
/// than an empirical falsification: the 4-mode L_eff(Q) is bond-summed by design (see
/// <see cref="Decomposition.FourModeEffective.LEffAtQ"/>), so its eigenstates cannot
/// carry bond-class information no matter the projection axis. The bond-class signature
/// must enter through dL/dJ_b per-bond V_b in the K-resonance, not through L_eff spectrum.</para>
///
/// <para>The composition-reading itself stands:</para>
///
/// <code>
///   r_Q(N, b) = BareDoubledPtfXPeak · Q_EP(N, b) − 2 = 4.39382 / g_eff(N, b) − 2
/// </code>
///
/// where <see cref="Item1Derivation.C2HwhmRatio.BareDoubledPtfXPeak"/> = 2.196910 (universal
/// in C2HwhmRatio, Tier 1 derived) and Q_EP(N, b) = 2/g_eff(N, b) (F86 Statement 1, Tier 1
/// derived). This is mathematically tautological if g_eff is defined via Q_peak inversion;
/// the genuine Tier-1 content is the universality of BareDoubledPtfXPeak, so the bond-class
/// signature must live entirely in g_eff(N, b). The polarity-pair signature {−0.5, +0.5} of
/// r_Q is then the polarity-layer fingerprint operationalised via the universal-shape
/// constant. See <see cref="ClosedFormCompositionNote"/>. What remains open is the closed-form
/// g_eff(N, b) per bond class. Direction (α)-test g_eff_E ≈ σ_0(N)·√(3/8) matches Δ ≤ 0.01
/// for N ≥ 6 but does NOT pin at tolerance 0.005. This inherits A3's
/// <see cref="Item1Derivation.C2InterChannelAnalytical"/> Tier 2 obstruction at the
/// |u_0⟩, |v_0⟩ closed-form level. The <see cref="PendingDerivationNote"/> carries the
/// refined next directions (α'), (β'), (γ').</para>
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

    /// <summary>Composition reading from the 2026-05-07 Direction (α) attempt: the per-bond
    /// r_Q(N, b) reduces to a composition of two existing Tier 1 facts in the c=2 chain:
    ///
    /// <code>
    ///   r_Q(N, b) = Q_peak(N, b) − 2 = BareDoubledPtfXPeak · Q_EP(N, b) − 2
    ///             = (BareDoubledPtfXPeak · 2 / g_eff(N, b)) − 2
    ///             = 4.39382 / g_eff(N, b) − 2
    /// </code>
    ///
    /// where:
    /// <list type="bullet">
    ///   <item><see cref="Item1Derivation.C2HwhmRatio.BareDoubledPtfXPeak"/> = 2.196910 is the
    ///         universal post-EP Q_peak/Q_EP ratio (Tier 1 derived in C2HwhmRatio).</item>
    ///   <item>Q_EP(N, b) = 2 / g_eff(N, b) per F86 Statement 1 (Tier 1 derived).</item>
    /// </list>
    ///
    /// <para><b>Tier reading:</b> the composition itself is mathematically a tautology when
    /// g_eff(N, b) is defined as 4.39382 / Q_peak(N, b); that just inverts Q_peak. The
    /// genuine structural content is that <see cref="Item1Derivation.C2HwhmRatio.BareDoubledPtfXPeak"/>
    /// is universal (g_eff-invariant in dimensionless x = Q/Q_EP), so the bond-class signature
    /// must live entirely in g_eff(N, b). Empirically per-class g_eff is N-stable:
    /// g_eff_Endpoint ≈ 1.74 (range 1.725..1.757 across N=5..8) and g_eff_Interior ≈ 2.81
    /// (range 2.738..2.964). For N ≥ 6 the Endpoint g_eff approximately matches σ_0(N)·√(3/8)
    /// (Δ ≤ 0.01), but tolerance is not bit-exact; for Interior, σ_0(N) alone is closer for
    /// N ≥ 6 (Δ ≤ 0.05) but Δ exceeds 0.10 at N=5 and N=8. Closed-form g_eff(N, b) thus inherits
    /// A3's <see cref="Item1Derivation.C2InterChannelAnalytical"/> Tier 2 obstruction at the
    /// |u_0⟩, |v_0⟩ closed-form level,the per-class g_eff witness numbers are pinned but not
    /// derived from σ_0(N) at tolerance 0.005.</para></summary>
    public string ClosedFormCompositionNote =>
        "r_Q(N, b) = BareDoubledPtfXPeak · Q_EP(N, b) − 2 = 4.39382 / g_eff(N, b) − 2 " +
        "(composition of C2HwhmRatio.BareDoubledPtfXPeak and F86 Statement 1's Q_EP = 2/g_eff). " +
        "The polarity-pair signature {−0.5, +0.5} of r_Q is the polarity-layer fingerprint " +
        "operationalised via BareDoubledPtfXPeak. Empirically per-class g_eff is N-stable: " +
        "g_eff_E ≈ 1.74, g_eff_I ≈ 2.81 across N=5..8. Direction-(α) test g_eff_E ≈ σ_0(N)·√(3/8) " +
        "matches at Δ ≤ 0.01 for N ≥ 6 but does not pin at tolerance 0.005. Closed-form " +
        "g_eff(N, b) inherits C2InterChannelAnalytical's Tier 2 obstruction.";

    /// <summary>Tier 2 empirical witness (sub-pattern): the asymptotic value of
    /// Q_peak_Endpoint + Q_peak_Interior across c=2 N=5..8 (pinned ≈ 4.12 from N=6..8 with a
    /// ~0.13 deviation at N=5). Equivalently: <c>1/g_eff_E + 1/g_eff_I</c> approaches ≈ 0.937
    /// for N ≥ 6. This is an empirical sub-pattern observed in the 2026-05-07 Direction (α)
    /// attempt; no closed form pinned (tested 2/π, (N+1)/(N+3), 1−1/√(N+1), 1−1/N, none match).
    /// Pinned as a witness for the next-direction work, consistent with the parent class's
    /// Tier2Verified discipline.</summary>
    public const double EmpiricalSumQPeakAsymptote = 4.12;

    /// <summary>The unresolved analytical piece: the per-bond r(N, b) closed form. The
    /// composition r_Q = BareDoubledPtfXPeak · Q_EP − 2 is mathematically tautological
    /// (g_eff is defined by Q_peak inversion in the empirical pipeline); the genuine open piece
    /// is the closed-form g_eff(N, b) which inherits A3's |u_0⟩, |v_0⟩ Tier 2 obstruction.</summary>
    public string PendingDerivationNote =>
        "Tier 2 verified: F86 bond-class split inherits structurally from the polarity-layer " +
        "pair {−0.5, +0.5} at d=2. The decomposition Q_peak ≈ 2 + r and HWHM/Q* ≈ 1/2 + r·(1/2) " +
        "holds across c=2 N=5..8 with mean Q_peak = 2.04 ± 0.06 and Interior r_HWHM ≈ 0.502 " +
        "(close to HalfAsStructuralFixedPoint).\n\n" +
        "## Composition-reading (2026-05-07 Direction (α) attempt)\n\n" +
        "r_Q(N, b) = BareDoubledPtfXPeak · Q_EP(N, b) − 2 = 4.39382 / g_eff(N, b) − 2.\n" +
        "Composition of two Tier 1 facts: BareDoubledPtfXPeak (= 2.196910 universal in " +
        "C2HwhmRatio) and Q_EP = 2/g_eff (F86 Statement 1). Caveat: mathematically tautological " +
        "if g_eff is defined via Q_peak inversion. The genuine Tier-1 content is the universality " +
        "of BareDoubledPtfXPeak,bond-class signature lives entirely in g_eff(N, b). The " +
        "polarity-pair signature {−0.5, +0.5} of r_Q is then the polarity-layer fingerprint " +
        "operationalised via the universal-shape constant. See ClosedFormCompositionNote.\n\n" +
        "## What was tried this session (Direction (α))\n\n" +
        "Per-bond polarity-Bloch projection of ρ_K(Q_peak, t_peak) onto channel-uniform " +
        "polarity axes (c_1 ± c_3)/√2 in the 4-mode basis B = [c_1, c_3, u_0, v_0].\n" +
        "Result: r_polarity(N, b) ≈ +0.39 → +0.21 across N=5..8 (decays with N), POSITIVE " +
        "for both bond classes (no sign split). Does not match r_Q (= ±0.5) or r_H (= ≈ +0.5) " +
        "at tolerance 0.005.\n\n" +
        "## Structural reframing (2026-05-08 code review)\n\n" +
        "Direction (α) is **structurally tautological under the uniform-J 4-mode reduction**, " +
        "not an empirical falsification. The 4-mode L_eff(Q) = D_eff + Q·γ₀·MhTotalEff is " +
        "**bond-summed by design** (FourModeEffective.LEffAtQ). Per-bond information enters " +
        "ONLY through the per-bond MhPerBondEff[b] used in the Duhamel K_b(Q, t) evaluation, " +
        "NOT through the L_eff spectrum or eigenstates. Any attempt to extract the bond-class " +
        "signature from the 4-mode L_eff EIGENSTATES alone (Direction (α)'s polarity-Bloch " +
        "projection at t_peak) is therefore guaranteed to be bond-class-blind by construction. " +
        "What looked like 'empirical falsification' is actually the design constraint that " +
        "L_eff is the same operator for every bond at uniform J. The bond-class signature " +
        "lives in the K-resonance via per-bond V_b, not in the L_eff spectrum.\n\n" +
        "Direction (α)-test for g_eff_E: hypothesis g_eff_E ≈ σ_0(N)·√(3/8). Δ across N=5..8: " +
        "0.063, 0.008, 0.004, 0.008. Matches at Δ ≤ 0.01 for N ≥ 6 but does NOT pin at " +
        "tolerance 0.005 for all N. σ_0(N) alone for Interior g_eff misses Δ > 0.10 at N=5.\n\n" +
        "## Refined next directions\n\n" +
        "(α') Full block-L per-bond Q_peak derivation from the K_b Duhamel formula (not 4-mode), " +
        "yielding the closed-form g_eff(N, b) directly. Joins direction (b'') of C2HwhmRatio.\n" +
        "(β') Locus 5 EP-rotation × Locus 6 polarity inheritance: tan θ = Q/Q_EP combined with " +
        "the polarity-pair shift to give r as a function of g_eff(N, b). Same residual.\n" +
        "(γ') Asymptotic structural constant: 1/g_eff_E + 1/g_eff_I ≈ 0.937 across N=6..8 " +
        "(empirical, EmpiricalSumQPeakAsymptote = 4.12 = 4.394·0.937). Test: does the " +
        "harmonic-mean structural reduce to the F86KB EP-rotation via a sum-rule? Witnessed " +
        "but not derived.\n\n" +
        "The 0.0006 deviation of Interior r_HWHM from exact 1/2 is likely numerical " +
        "discretisation (testable with finer Q-grid; see PROOF_F86_QPEAK Open elements 5 for " +
        "per-F71-orbit substructure that may explain the small remaining shift).";

    /// <summary>Tier 2 reflection: a Tier 1 sub-result (the closed-form composition for r_Q
    /// via BareDoubledPtfXPeak · Q_EP) DID land in the Direction (α) attempt, but does not by
    /// itself fully resolve r(N, b) because g_eff(N, b) is left open at the |u_0⟩, |v_0⟩
    /// closed-form level (A3 inheritance). Flips to <c>true</c> on future closure of
    /// g_eff(N, b) per the (α') / (β') / (γ') paths documented in
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
        "Tier2Verified: r_Q(N, b) = BareDoubledPtfXPeak · Q_EP(N, b) − 2 (composition-reading from " +
        "C2HwhmRatio.BareDoubledPtfXPeak + F86 Statement 1's Q_EP = 2/g_eff); g_eff(N, b) closed " +
        "form open (Direction (α)-test g_eff_E ≈ σ_0·√(3/8) matches Δ ≤ 0.01 for N ≥ 6 but not " +
        "tol 0.005; inherits A3 Tier 2 obstruction)";

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
            yield return new InspectableNode("Tier 1 sub-result: r_Q closed-form composition",
                summary: ClosedFormCompositionNote);
            yield return InspectableNode.RealScalar(
                "EmpiricalSumQPeakAsymptote (witness)", EmpiricalSumQPeakAsymptote, "F2");
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
    // Refresh: re-run C2HwhmRatio.HwhmLeftOverQPeakMean(BondClass.{Interior,Endpoint}) and the
    // per-block Q_peak readouts at γ₀=0.05 for each N=5..8 to update the pinned table below.
    private static readonly PolarityWitness[] _polarityWitnesses = new[]
    {
        new PolarityWitness(N: 5, QPeakInterior: 1.4821, QPeakEndpoint: 2.5008, HwhmRatioInterior: 0.7455, HwhmRatioEndpoint: 0.7700),
        new PolarityWitness(N: 6, QPeakInterior: 1.5801, QPeakEndpoint: 2.5470, HwhmRatioInterior: 0.7529, HwhmRatioEndpoint: 0.7738),
        new PolarityWitness(N: 7, QPeakInterior: 1.5831, QPeakEndpoint: 2.5299, HwhmRatioInterior: 0.7507, HwhmRatioEndpoint: 0.7738),
        new PolarityWitness(N: 8, QPeakInterior: 1.6049, QPeakEndpoint: 2.5145, HwhmRatioInterior: 0.7531, HwhmRatioEndpoint: 0.7734),
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
/// <para>The R-fields are exposed as computed properties from the raw Q_peak and HWHM/Q*
/// values via R_QPeak = QPeak − 2 and R_Hwhm = 2·(HwhmRatio − 1/2); they encode the
/// polarity-layer content explicitly so consumers do not have to re-decompose.</para>
/// </summary>
/// <param name="N">Chain length.</param>
/// <param name="QPeakInterior">Q_peak for the Interior bond class at this N.</param>
/// <param name="QPeakEndpoint">Q_peak for the Endpoint bond class at this N.</param>
/// <param name="HwhmRatioInterior">HWHM/Q* for the Interior bond class at this N.</param>
/// <param name="HwhmRatioEndpoint">HWHM/Q* for the Endpoint bond class at this N.</param>
public sealed record PolarityWitness(
    int N,
    double QPeakInterior,
    double QPeakEndpoint,
    double HwhmRatioInterior,
    double HwhmRatioEndpoint
) : IInspectable
{
    /// <summary>Polarity-content r for Interior Q_peak: r = QPeakInterior − 2.</summary>
    public double RQpeakInterior => QPeakInterior - 2.0;

    /// <summary>Polarity-content r for Endpoint Q_peak: r = QPeakEndpoint − 2.</summary>
    public double RQpeakEndpoint => QPeakEndpoint - 2.0;

    /// <summary>Polarity-content r for Interior HWHM ratio: r = 2·(HwhmRatioInterior − 1/2).</summary>
    public double RHwhmInterior => 2.0 * (HwhmRatioInterior - 0.5);

    /// <summary>Polarity-content r for Endpoint HWHM ratio: r = 2·(HwhmRatioEndpoint − 1/2).</summary>
    public double RHwhmEndpoint => 2.0 * (HwhmRatioEndpoint - 0.5);

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
