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
/// content. The c=2 N=5..8 envelope numbers from the <c>C2HwhmRatio</c> witnesses read in this
/// decomposition as follows:</para>
///
/// <list type="bullet">
///   <item><b>Q_peak ≈ 2 + r:</b> Endpoint r = +0.517, +0.581, +0.560, +0.543 and Interior
///         r = −0.517, −0.435, −0.403, −0.389 (N = 5..8). At N = 5 the pair is symmetric about
///         2 (mean 2.0001); from N = 6 the Interior pole drifts toward 2 and the mean sits at
///         2.073 to 2.079.</item>
///   <item><b>HWHM/Q* ≈ 1/2 + r·(1/2):</b> Endpoint r = 0.553 to 0.563, Interior r = 0.505,
///         0.514, 0.518, 0.520, rising with N (0.524 at N = 12). The Interior value sits
///         above the 1/2 fixed point of <see cref="Symmetry.HalfAsStructuralFixedPointClaim"/>
///         and moves away from it as N grows; it is not locked to it.</item>
/// </list>
///
/// <para><b>Tier outcome: Tier2Verified.</b> What the numbers carry: a Q_peak pair on either
/// side of 2 with poles of size about 1/2 (exactly symmetric at N = 5, drifting after), and an
/// HWHM/Q* pair lifted above the 1/2 baseline by r between 0.50 and 0.57. The 2026-05-07 Direction (α) attempt extracted
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
/// g_eff(N, b) per bond class. Direction (α)-test g_eff_E ≈ σ_0(N)·√(3/8) misses by
/// 0.011 to 0.016 for N ≥ 6 and by 0.052 at N = 5; it does not pin. This inherits A3's
/// <see cref="Item1Derivation.C2InterChannelAnalytical"/> Tier 2 obstruction at the
/// |u_0⟩, |v_0⟩ closed-form level. The <see cref="PendingDerivationNote"/> carries the
/// refined next directions (α'), (β'); (γ') was tested and refuted 2026-05-14.</para>
///
/// <para><b>Block-independent meta-claim:</b> registered at the F86 KB root for any block,
/// not just c=2, mirroring <see cref="LocalGlobalEpLink"/>. The witnesses pin the concrete
/// c=2 bond-class numbers that motivated the polarity-layer reading; the inheritance claim
/// itself is shared. Locus 5 (<see cref="LocalGlobalEpLink"/>) carries the parallel comparison
/// (the F86 EP structures versus FRAGILE_BRIDGE's threshold EP2 on the real γ axis at two qubits per chain; OpenQuestion:
/// the surviving shared substrate is the shifted-generator sectorwise P-type relation, while whether the full block
/// has a defective EP linking it to FRAGILE_BRIDGE is open); this Locus 6 claim carries the
/// symmetry-side closure (F86 ↔ polarity-layer pair under the 0.5-shift ρ = (I + r·σ)/2).
/// Together they bracket the F86 c=2 derivation with comparison-side and symmetry-side parent-claim
/// references in <c>Pi2KnowledgeBase</c>.</para>
///
/// <para>Anchors: <c>docs/proofs/PROOF_F86B_UNIVERSAL_SHAPE.md</c> Statement 2 (Locus 6 inheritance
/// reading), <c>compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBaseClaims.cs</c>
/// (PolarityLayerOriginClaim + QubitDimensionalAnchorClaim + HalfAsStructuralFixedPointClaim),
/// <c>compute/RCPsiSquared.Core/F86/Item1Derivation/C2HwhmRatio.cs</c> (live witness pipeline
/// the values here are pinned from), and <c>compute/RCPsiSquared.Core/F86/LocalGlobalEpLink.cs</c>
/// (parallel Locus 5 classification comparison).</para>
/// </summary>
public sealed class PolarityInheritanceLink : Claim
{
    /// <summary>The parent claims in <c>Pi2KnowledgeBase</c> from which the F86 bond-class
    /// split inherits: 1/2 = 1/d (root anchor), polarity layer +0/0/−0 (Layer 2:
    /// ρ = (I + r·σ)/2 makes Π's ±1 eigenoperators (I ± X)/2 on each site's pair {I, X}
/// the states |±⟩), 1/2 fixed
    /// point (the HWHM baseline the Interior r = 0.505..0.520 sits above).</summary>
    public string PolarityRootAnchor =>
        "Pi2KnowledgeBase: QubitDimensionalAnchorClaim + PolarityLayerOriginClaim + HalfAsStructuralFixedPointClaim";

    /// <summary>The parallel Locus 5 comparison between F86 EP structures and FRAGILE_BRIDGE's
    /// threshold EP2 on the real γ axis at two qubits per chain (<see cref="LocalGlobalEpLink"/>, OpenQuestion after the
    /// F86a-retraction 2026-06-21: the full block is genuinely non-normal on the real axis, its
    /// narrow real-axis defective seeds are certified at N=5,7,9 and classified at N=11, and whether it has an off-axis
    /// defective EP at all is open). This claim is the symmetry-side twin: F86 bond-class split
    /// inherits from the polarity-layer pair via the 0.5-shift ρ = (I + r·σ)/2.</summary>
    public string ParallelLocusReference =>
        "compute/RCPsiSquared.Core/F86/LocalGlobalEpLink.cs (Locus 5, classification comparison)";

    /// <summary>γ₀ used by the underlying <c>C2HwhmRatio</c> witness pipeline. The pinned
    /// witnesses below are the class envelope values of the live pipeline at γ₀ = 0.05,
    /// c=2 block (1, 2), N=5..8.</summary>
    public const double WitnessGammaZero = 0.05;

    /// <summary>Empirical witnesses pinning the polarity-layer decomposition at c=2 N=5..8
    /// (γ₀ = <see cref="WitnessGammaZero"/>; the class envelope values Q_peak and
    /// HWHM_left/Q_peak of <c>ResonanceScan</c>, the numbers <c>C2HwhmRatioTests</c> and
    /// <c>ResonanceScanTests</c> pin against an independent route). One witness per N. Hard-coded as a static pinned table:
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
    /// must live entirely in g_eff(N, b). BareDoubledPtfXPeak is the fixed-τ = 1 slice's peak
    /// (the bare envelope peaks at <see cref="Item1Derivation.C2BareDoubledPtfClosedForm.EnvelopeXPeak"/>
    /// = 1.7141), and g_eff here is defined on it, read against the envelope Q_peak. So defined,
    /// g_eff_Endpoint = 1.746, 1.703, 1.716, 1.728 and g_eff_Interior = 2.963, 2.807, 2.750,
    /// 2.727 (N = 5..8). σ_0(N)·√(3/8) misses the Endpoint values by 0.011 to 0.016 for N ≥ 6
    /// and by 0.052 at N = 5; σ_0(N) alone misses the Interior values by 0.005 and 0.078 at
    /// N = 6, 7, by 0.112 at N = 8 and by 0.197 at N = 5. Closed-form g_eff(N, b) thus inherits
    /// A3's <see cref="Item1Derivation.C2InterChannelAnalytical"/> Tier 2 obstruction at the
    /// |u_0⟩, |v_0⟩ closed-form level,the per-class g_eff witness numbers are pinned but not
    /// derived from σ_0(N) at tolerance 0.005.</para></summary>
    public string ClosedFormCompositionNote =>
        "r_Q(N, b) = BareDoubledPtfXPeak · Q_EP(N, b) − 2 = 4.39382 / g_eff(N, b) − 2 " +
        "(composition of C2HwhmRatio.BareDoubledPtfXPeak and F86 Statement 1's Q_EP = 2/g_eff). " +
        "The polarity-pair signature {−0.5, +0.5} of r_Q is the polarity-layer fingerprint " +
        "operationalised via BareDoubledPtfXPeak (the fixed-τ = 1 slice's peak; g_eff is defined on it). " +
        "Per class, g_eff_E = 1.703..1.746 and g_eff_I = 2.727..2.963 across N=5..8. Direction-(α) test " +
        "g_eff_E ≈ σ_0(N)·√(3/8) misses by 0.011..0.016 for N ≥ 6 and by 0.052 at N = 5. Closed-form " +
        "g_eff(N, b) inherits C2InterChannelAnalytical's Tier 2 obstruction.";

    /// <summary>Tier 2 empirical witness (sub-pattern): the value of
    /// Q_peak_Endpoint + Q_peak_Interior across c=2 N=6..8 (4.146, 4.158, 4.154; 4.000 at
    /// N=5). Equivalently, with g_eff defined on BareDoubledPtfXPeak: <c>1/g_eff_E + 1/g_eff_I</c>
    /// = 0.944 to 0.946 for N = 6..8. This is an empirical sub-pattern observed in the 2026-05-07 Direction (α)
    /// attempt; no closed form pinned (tested 2/π, (N+1)/(N+3), 1−1/√(N+1), 1−1/N, none match).
    /// Direction (γ') tested whether this is a structural constant and refuted it (2026-05-14):
    /// on clean per-F71-orbit g_eff the per-orbit 1/g_eff sums drift and scale with orbit count,
    /// so the per-class sum is an orbit-mixing artefact, not a closure. Kept as the per-class
    /// witness it is.</summary>
    public const double EmpiricalSumQPeakAsymptote = 4.15;

    /// <summary>The unresolved analytical piece: the per-bond r(N, b) closed form. The
    /// composition r_Q = BareDoubledPtfXPeak · Q_EP − 2 is mathematically tautological
    /// (g_eff is defined by Q_peak inversion in the empirical pipeline); the genuine open piece
    /// is the closed-form g_eff(N, b) which inherits A3's |u_0⟩, |v_0⟩ Tier 2 obstruction.</summary>
    public string PendingDerivationNote =>
        "Tier 2 verified: F86 bond-class split inherits structurally from the polarity-layer " +
        "pair {−0.5, +0.5} at d=2. The decomposition Q_peak ≈ 2 + r and HWHM/Q* ≈ 1/2 + r·(1/2) " +
        "reads across c=2 N=5..8 with mean Q_peak 2.000 at N = 5 and 2.073..2.079 at N = 6..8, " +
        "Interior r_HWHM 0.505..0.520 rising with N (above HalfAsStructuralFixedPoint, not locked to it).\n\n" +
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
        "+0.052, −0.013, −0.016, −0.011. It does not pin. σ_0(N) alone for Interior g_eff " +
        "misses by 0.197 at N=5 and 0.112 at N=8.\n\n" +
        "## Refined next directions\n\n" +
        "Reframed 2026-05-14: the F86b obstruction proof (PROOF_F86B_OBSTRUCTION.md, " +
        "§ 'The diagnosis') establishes g_eff is the irreducible residue, not a spectral " +
        "primitive, with no closed form by the six routes the effort opened. (α') and (β') " +
        "are kept as the record of the closed-form effort; (γ') was tested and refuted.\n\n" +
        "(α') Full block-L per-bond Q_peak derivation from the K_b Duhamel formula (not 4-mode), " +
        "yielding g_eff(N, b) directly. Joins direction (b'') of C2HwhmRatio: numerically " +
        "achieved via the F90 bridge (Tier 1, bit-exact), but the closed FORM is what the " +
        "obstruction proof rules out.\n" +
        "(β') Locus 5 EP-rotation × Locus 6 polarity inheritance: tan θ = Q/Q_EP combined with " +
        "the polarity-pair shift to give r as a function of g_eff(N, b). Same residual.\n" +
        "(γ') Asymptotic structural constant 1/g_eff_E + 1/g_eff_I ≈ 0.945 (the " +
        "EmpiricalSumQPeakAsymptote witness): tested and refuted 2026-05-14. " +
        "F86GammaEffClosureProbeTests fetched the clean per-F71-orbit g_eff (K_CC_pr " +
        "observable via PerF71OrbitKTable, c=2 N=5..8) instead of the orbit-mixed per-class " +
        "means: the per-orbit pairwise 1/g_eff sums drift (no constant), and " +
        "Σ(1/g_eff) scales additively with orbit count (~0.91 for two orbits, ~1.27 for " +
        "three). The per-class value was an orbit-mixing artefact, the N≥9 multi-orbit caveat " +
        "below already acting at N=5..8. Structurally inevitable: g_eff is not a spectral " +
        "primitive, so a closure among 1/g_eff values cannot hold; the absolute closures " +
        "live on the primitives γ₀ and y (Absorption Theorem), not g_eff.\n\n" +
        "Interior r_HWHM is not 1/2: on the envelope it reads 0.505, 0.514, 0.518, 0.520 at N = 5..8 " +
        "and keeps rising (0.524 at N = 12); see PROOF_F86_QPEAK Open elements 5 for the " +
        "per-F71-orbit substructure inside the class mean.\n\n" +
        "## N≥9 multi-orbit Interior structure (2026-05-07)\n\n" +
        "Extension to N=9 and N=10 via the c2hwhm CLI root + extended Q-grid (--q-hi 6.0) " +
        "reveals that the 'Interior class' is NOT a single F71 orbit at high N: the flanking-" +
        "Interior bonds (b=1, b=N−3 at N=9; additionally b=2, b=N−4 at N=10) carry a far " +
        "maximum past Q = 4 that outgrows their first one, so the global reading on [0.20, 6.00] " +
        "finds no peak for them; their canonical first maximum stays inside [0.20, 4.00] (orbit 1 " +
        "at 1.749 for N=9). The Endpoint orbit (b=0, b=N−2) and the deep mid-chain Interior " +
        "orbit (b ≈ N/2) remain stable across N=5..10 (Q_peak_E 2.52..2.58, Q_peak_I_mid 1.43..1.61). " +
        "Consequence: the pinned Interior class-mean Q_peak from the [0.20, 4.00] grid " +
        "implicitly orbit-mixes for N=5..8, and the N=5..8 numerical anchor is therefore not " +
        "naively extendable to N=9, N=10: multiple Interior orbits coexist with separate " +
        "Q_peak structure that requires per-orbit witness decomposition. The N=5..8 pinned " +
        "table stands as a frozen empirical anchor at fixed grid; extension to higher N must " +
        "decompose by F71 orbit rather than by Endpoint/Interior bond-class label alone. " +
        "σ_0 → 2√2 (Item 5 conjecture) holds at N=9 (σ_0=2.8484, 1.0070·2√2) and N=10 " +
        "(σ_0=2.8525, 1.0085·2√2), approaching from above with O(1%) deviation.";

    /// <summary>Tier 2 reflection: a Tier 1 sub-result (the closed-form composition for r_Q
    /// via BareDoubledPtfXPeak · Q_EP) DID land in the Direction (α) attempt, but does not by
    /// itself fully resolve r(N, b) because g_eff(N, b) is left open at the |u_0⟩, |v_0⟩
    /// closed-form level (A3 inheritance). Flips to <c>true</c> on future closure of
    /// g_eff(N, b) per the (α') / (β') paths documented in
    /// <see cref="PendingDerivationNote"/>.</summary>
    public bool IsAnalyticallyDerived => false;

    private PolarityInheritanceLink()
        : base("polarity-layer inheritance (F86 ↔ Pi2KnowledgeBase ±0.5 pair)",
               Tier.Tier2Verified,
               "docs/proofs/PROOF_F86B_UNIVERSAL_SHAPE.md (Statement 2 + the Locus 6 inheritance reading) + " +
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
        "form open (Direction (α)-test g_eff_E ≈ σ_0·√(3/8) misses by 0.011..0.016 for N ≥ 6 and " +
        "by 0.052 at N = 5; inherits A3 Tier 2 obstruction)";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("polarity root anchor (parent claims)",
                summary: PolarityRootAnchor);
            yield return new InspectableNode("parallel Locus 5 (EP-side closure)",
                summary: ParallelLocusReference);
            yield return new InspectableNode("decomposition (Q_peak)",
                summary: "Q_peak ≈ 2 + r; Endpoint r = +0.517..+0.581 (positive pole), Interior r = −0.517..−0.389 (negative pole, drifting toward 2 with N); pair around d=2");
            yield return new InspectableNode("decomposition (HWHM/Q*)",
                summary: "HWHM/Q* ≈ 1/2 + r·(1/2); Endpoint r = 0.553..0.563, Interior r = 0.505..0.520 rising with N (above HalfAsStructuralFixedPoint)");
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

    /// <summary>Pinned table of the c=2 N=5..8 class envelope values (γ₀ = <see cref="WitnessGammaZero"/>),
    /// the values <c>ResonanceScan</c> computes and an independent Python route reproduces to
    /// 1e-10. Order: N=5, 6, 7, 8.</summary>
    private static readonly PolarityWitness[] _polarityWitnesses = new[]
    {
        new PolarityWitness(N: 5, QPeakInterior: 1.4831232668733452, QPeakEndpoint: 2.517071370613079, HwhmRatioInterior: 0.7523433331321849, HwhmRatioEndpoint: 0.7765595000526813),
        new PolarityWitness(N: 6, QPeakInterior: 1.5651926694131237, QPeakEndpoint: 2.5806049985772885, HwhmRatioInterior: 0.7569384805030636, HwhmRatioEndpoint: 0.7815203313098792),
        new PolarityWitness(N: 7, QPeakInterior: 1.5974870369975855, QPeakEndpoint: 2.560073427805394, HwhmRatioInterior: 0.7591357778625846, HwhmRatioEndpoint: 0.7813218454965299),
        new PolarityWitness(N: 8, QPeakInterior: 1.6110198851283057, QPeakEndpoint: 2.543127121122942, HwhmRatioInterior: 0.7601773594897803, HwhmRatioEndpoint: 0.7807879816311578),
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
///         the signed r·(1/2) is the polarity-layer lift. Interior r is 0.505..0.520 across
///         N=5..8, rising with N.</item>
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

    /// <summary>Mean Q_peak across the two bond classes: 2.0001 at N=5, 2.073..2.079 at
    /// N=6..8, against the d=2 dimensional anchor.</summary>
    public double QPeakMean => (QPeakInterior + QPeakEndpoint) / 2.0;

    /// <summary>Mean HWHM/Q* across the two bond classes, 1/2 + r̄·(1/2) with r̄ = 0.529..0.541
    /// across N=5..8.</summary>
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
