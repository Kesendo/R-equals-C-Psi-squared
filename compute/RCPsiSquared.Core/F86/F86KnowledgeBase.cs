using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.Decomposition;
using RCPsiSquared.Core.F86.Item1Derivation;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Resonance;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.F86;

/// <summary>"Everything we know about F86" as a typed OOP knowledge graph attached to a
/// <see cref="CoherenceBlock"/>. The root <see cref="IInspectable"/> for <c>--root f86</c>.
///
/// <para>Below the block node and the named Q-anchor map, the claims are organised into
/// five groups; the c=2-only members appear only when <c>Block.C == 2</c>:</para>
/// <list type="bullet">
///   <item><b>Tier 1 (derived)</b>: <see cref="TPeakLaw"/>, <see cref="QEpLaw"/>,
///         <see cref="TwoLevelEpModel"/> (pre/at/post-EP traversal and the higher-k
///         hierarchy), <see cref="ChiralAiiiClassification"/>, <see cref="F71MirrorInvariance"/>,
///         <see cref="LEffMirrorAxisClaim"/>, <see cref="PolarityPairQPeakDecompositionClaim"/>;
///         c=2-only <see cref="F90F86C2BridgeIdentity"/> and the F86e identity
///         <see cref="SigmaZeroCommutatorNormClaim"/> (σ_0 = ‖[Π_HD1, M_H]‖).</item>
///   <item><b>Tier 1 (candidate)</b>: <see cref="UniversalShapePrediction"/> for Interior
///         and Endpoint with <see cref="UniversalShapeWitness"/> data across c=2..4 N=5..8,
///         <see cref="ShapeFunctionWitnesses"/>, <see cref="DressedModeWeightClaim"/>,
///         <see cref="F86HwhmClosedFormClaim"/>; c=2-only <see cref="C2UniversalShapeDerivation"/>.</item>
///   <item><b>Tier 2</b> (verified and empirical): <see cref="SigmaZeroChromaticityScaling"/>,
///         the per-block and per-bond Q_peak tables, <see cref="PerF71OrbitObservation"/>,
///         <see cref="PolarityInheritanceLink"/>, <see cref="IbmBlockCpsiHardwareTable"/>;
///         c=2-only the live orbit-K and full-block σ-anatomy tables.</item>
///   <item><b>Retracted</b>: the refuted F86 closed forms / mechanisms
///         (<see cref="RetractedClaim.Standard"/>): the two 2026-05-02 Q_peak csc(...)
///         conjectures plus the 2026-06-21 F86a real-axis-EP mechanism (no eigenvalue
///         coalescence on the real Q axis; the Petermann factor is genuine non-normality,
///         not an EP artifact).</item>
///   <item><b>Open questions</b>: <see cref="F86OpenQuestions.Standard"/>, the items still
///         missing for full Tier-1 promotion, plus <see cref="LocalGlobalEpLink"/> (demoted to
///         <see cref="Knowledge.Tier.OpenQuestion"/> by the F86a-retraction 2026-06-21:
///         genuine non-normality on the real axis but no real-axis EP; whether the full block
///         has an off-axis defective EP at all is open), followed by the 4-mode insufficiency
///         note.</item>
/// </list>
/// </summary>
public sealed class F86KnowledgeBase : IInspectable
{
    public CoherenceBlock Block { get; }
    public WitnessCache WitnessCache { get; }
    public TPeakLaw TPeak { get; }
    public InterChannelSvd? Svd { get; }
    public QEpLaw? QEp { get; }
    public IReadOnlyList<TwoLevelEpModel> EpTraversal { get; }
    public IReadOnlyList<TwoLevelEpModel> HigherKLevels { get; }
    public UniversalShapePrediction InteriorShape { get; }
    public UniversalShapePrediction EndpointShape { get; }
    public ShapeFunctionWitnesses InteriorShapeFunction { get; }
    public ShapeFunctionWitnesses EndpointShapeFunction { get; }
    public IReadOnlyList<PerBlockQPeakClaim> PerBlockQPeaks { get; }
    public PerBondQPeakWitnessTable EndpointPerBondTable { get; }
    public PerBondQPeakWitnessTable InteriorPerBondTable { get; }
    public PerF71OrbitObservation PerOrbitSubstructure { get; }

    /// <summary>The canonical Q-anchors on the Q = J/γ₀ axis as typed records with
    /// per-anchor (Q, J@γ₀=0.05, band, role, tier, source) metadata. See
    /// <see cref="QAnchorMap"/> for the full list. Block-independent (the anchors are
    /// framework-wide, not per-block); constructed eagerly mirroring
    /// <see cref="PerBlockQPeaks"/>.</summary>
    public QAnchorMap QAnchors { get; }

    private readonly Lazy<InspectableNode> _qAnchorsChildrenGroup;
    public DressedModeWeightClaim DressedModeWeight { get; }
    public ChiralAiiiClassification AlgebraicClass { get; }
    public F71MirrorInvariance F71Mirror { get; }
    public SigmaZeroChromaticityScaling Sigma0Scaling { get; }
    /// <summary>F86 c=2 HWHM_left/Q_peak prediction per BondSubClass (Tier 1 candidate):
    /// bare floor 0.671535 IS derived (`C2BareDoubledPtfClosedForm`); per-sub-class
    /// (alpha, beta) are fitted via polyfit on N=5..8 F90-bridge anchors. Block-independent
    /// meta-claim.</summary>
    public F86HwhmClosedFormClaim F86HwhmClosedForm { get; }
    /// <summary>The c=2-specific top-level Claim wrapping the Stage D2
    /// <see cref="Item1Derivation.C2HwhmRatio"/>: empirical anchor + directional
    /// Endpoint &gt; Interior split, both pinned for the c=2 stratum (witness anchor:
    /// c=2 N=5..8 in the canonical sweep). Non-null iff <c>Block.C == 2</c>.
    /// Lazily built on first access — the underlying full-block Q-scan
    /// (153×21 eigendecompositions) is paid only by consumers that read the property.</summary>
    public C2UniversalShapeDerivation? C2UniversalShape => _c2UniversalShape.Value;

    private readonly Lazy<C2UniversalShapeDerivation?> _c2UniversalShape;

    /// <summary>Block-independent meta-claim: the F86 ↔ FRAGILE_BRIDGE EP-side relation.
    /// <see cref="Knowledge.Tier.OpenQuestion"/> after the F86a-retraction 2026-06-21
    /// (demoted from Tier2Verified): the full Σγ=N·γ₀ block is genuinely strongly non-normal
    /// on the real Q axis but has NO eigenvalue coalescence there (eigenvalues simple), the
    /// Petermann factor large but finite (artifact-free ‖P‖); whether it has an off-axis
    /// defective EP at all, and thus the "same EP" link to FRAGILE_BRIDGE, is open. The
    /// surviving shared substrate is the AIII chiral algebra. Exposed at the KB root for any
    /// block, not just c=2. Lazy: a single static-data Claim with no compute cost.</summary>
    public LocalGlobalEpLink LocalGlobalEpLink => _localGlobalEpLink.Value;

    private readonly Lazy<LocalGlobalEpLink> _localGlobalEpLink;

    /// <summary>Block-independent meta-claim (Locus 6, symmetry-side closure): the F86 c=2
    /// bond-class split (Q_peak 2 ± 0.5, HWHM/Q* 1/2 + r·(1/2)) inherits structurally from
    /// the polarity-layer pair {−0.5, +0.5} at d=2 named in <c>Pi2KnowledgeBase</c>
    /// (PolarityLayerOriginClaim + QubitDimensionalAnchorClaim + HalfAsStructuralFixedPointClaim).
    /// <see cref="Knowledge.Tier.Tier2Verified"/>: the decomposition is empirically pinned
    /// across c=2 N=5..8; per-bond r(N, b) closed form is the documented gap. Companion
    /// to <see cref="LocalGlobalEpLink"/> (Locus 5, EP-side closure): together they bracket
    /// the F86 c=2 derivation with EP-side and symmetry-side parent-claim references.
    /// Exposed at the KB root for any block, not just c=2; the inheritance statement is
    /// shared across all c. Lazy: a single static-data Claim with no compute cost.</summary>
    public PolarityInheritanceLink PolarityInheritanceLink => _polarityInheritanceLink.Value;

    private readonly Lazy<PolarityInheritanceLink> _polarityInheritanceLink;

    /// <summary>F86 schema-level meta-claim (Tier1Derived): Q_peak ∈ {2 − 1/2, 2 + 1/2} =
    /// {1.5, 2.5} via composition of two Tier1Derived parents: idealised <see cref="QEpLaw"/>
    /// at g_eff = 1 (Q_EP_central = 2) and <see cref="Symmetry.HalfAsStructuralFixedPointClaim"/>
    /// (polarity magnitude r = 1/2). Companion to <see cref="PolarityInheritanceLink"/>
    /// (Tier2Verified, holds the bit-exact witness table with ~10% finite-size deviation):
    /// the schema is Tier1Derived, the bit-exact value is Tier2Verified, the closed form
    /// for the deviations is blocked by PROOF_F86B_OBSTRUCTION (six tested routes).
    /// Block-independent; resolves to a process-wide
    /// <see cref="PolarityPairQPeakDecompositionClaim.Shared"/> singleton.</summary>
    public PolarityPairQPeakDecompositionClaim PolarityPairQPeakDecomposition =>
        PolarityPairQPeakDecompositionClaim.Shared;

    /// <summary>F86 meta-claim (Tier1Derived): in the 2-level reduction L_eff, the real
    /// part −4γ₀ = −2γ₀·2 is the mirror axis of the channel pair (−2γ₀, −6γ₀); the EP is
    /// the coalescence onto it; g_eff lives in the imaginary part as the branch's relative
    /// clock. Pure 2×2 algebra, block-independent. See <see cref="LEffMirrorAxisClaim"/>.</summary>
    public LEffMirrorAxisClaim LEffMirrorAxis => _lEffMirrorAxis.Value;

    private readonly Lazy<LEffMirrorAxisClaim> _lEffMirrorAxis;

    /// <summary>F71-orbit-grouped live K-resonance witness table for c=2
    /// (Tier2Verified). Live counterpart to <see cref="PerOrbitSubstructure"/>
    /// (frozen 9-case sweep). Non-null iff <c>Block.C == 2</c>. Lazily built on
    /// first access — the underlying full-block Q-scan inside
    /// <see cref="PerF71OrbitKTable.Build"/> runs only once and only on read.</summary>
    public PerF71OrbitKTable? OrbitKTable => _orbitKTable.Value;

    /// <summary>F89 per-Bloch-mode σ_n extraction via R†·S·R diagonal at the c=2 stratum
    /// (Tier 2 verified). Bit-exact match against
    /// <see cref="F89UnifiedFaClosedFormClaim.Sigma"/> for path-3..9 (path-7 closed form
    /// derived via this anatomy; boundary case nBlock = N = k+1, no bare site).
    /// Non-null iff <c>Block.C == 2</c>. Lazily built on first access,
    /// the full-block eigendecomposition runs only once and only on read.</summary>
    public C2FullBlockSigmaAnatomy? FullBlockSigmaAnatomy => _fullBlockSigmaAnatomy.Value;

    private readonly Lazy<C2FullBlockSigmaAnatomy?> _fullBlockSigmaAnatomy;

    private readonly Lazy<PerF71OrbitKTable?> _orbitKTable;

    /// <summary>Block-independent meta-claim (Tier2Verified): the IBM 2026-04-26
    /// framework_snapshots viewed through Theorem 2's universal C_block ≤ 1/4 ceiling
    /// (PROOF_BLOCK_CPSI_QUARTER). Pinned 32-row witness table over (Aer / Marrakesh /
    /// Kingston / Fez) × (heisenberg / truly / soft / hard) × ((0,1) / (1,2)) blocks of
    /// the (q0, q2) reduced ρ on the N=3 |+−+⟩ chain at t=0.8, J=1.0. Independent of
    /// <see cref="Block"/>; the lens applies at every c-block.</summary>
    public IbmBlockCpsiHardwareTable IbmBlockCpsiHardwareTable => _ibmBlockCpsiHardwareTable.Value;

    private readonly Lazy<IbmBlockCpsiHardwareTable> _ibmBlockCpsiHardwareTable;

    /// <summary>F90 c=2 ↔ F89 bridge identity (Tier1Derived): F86 c=2 K_b(Q,t)
    /// IS F89 path-(N−1) (SE,DE) per-bond Hellmann-Feynman, modulo the J convention J_F89 = J_F86/2
    /// (operator-exact: ‖L_F86(J) − L_F89(J/2)‖ = 0 at N=5..8).
    /// Bit-exact verified at 20/22 bonds across N=5..8 including orbit escapes. Resolves
    /// Direction (b'') (full block-L derivation, NOT 4-mode) numerically Tier-1; closed-form
    /// HWHM_left/Q_peak per bond class via F89's AT-locked F_a/F_b structure remains the
    /// analytical target. Block-independent meta-claim — the bridge identity holds for any
    /// c=2 block N ≥ 3. See <c>docs/proofs/PROOF_F90_F86C2_BRIDGE.md</c>.</summary>
    public F90F86C2BridgeIdentity F89BridgeIdentity => _f89BridgeIdentity.Value;

    private readonly Lazy<F90F86C2BridgeIdentity> _f89BridgeIdentity;

    /// <summary>F86e c=2 meta-claim (Tier1Derived): the inter-channel SVD-top singular value
    /// σ_0 IS the operator norm of the commutator [Π_HD1, M_H]. On the c=2 block HD ∈ {1, 3}
    /// only, so Π_HD1 + Π_HD3 = I and V_inter = Π_HD1·M_H·(I−Π_HD1); the lemma
    /// ‖P·M·(1−P)‖ = ‖[P, M]‖ then gives σ_0 = ‖[Π_HD1, M_H]‖. Non-null iff <c>Block.C == 2</c>
    /// — the identity is c=2-specific (at c ≥ 3 the HD spectrum has more than two values).
    /// Lazily built on first access: the underlying InterChannelSvd plus the square-projector
    /// commutator norm run only once and only on read. See
    /// <see cref="SigmaZeroCommutatorNormClaim"/>.</summary>
    public SigmaZeroCommutatorNormClaim? SigmaZeroCommutatorNorm => _sigmaZeroCommutatorNorm.Value;

    private readonly Lazy<SigmaZeroCommutatorNormClaim?> _sigmaZeroCommutatorNorm;

    public IReadOnlyList<RetractedClaim> Retracted { get; }
    public IReadOnlyList<OpenQuestion> OpenQuestions { get; }
    public InspectableNode FourModeInsufficiencyNote { get; }

    public F86KnowledgeBase(CoherenceBlock block, WitnessCache? witnessCache = null)
    {
        Block = block;
        WitnessCache = witnessCache ?? new WitnessCache();
        TPeak = new TPeakLaw(block.GammaZero);

        // Q_EP is meaningful only when chromaticity ≥ 2 (we need HD=1 and HD=3 channels for
        // the inter-channel SVD that produces the natural g_eff).
        if (block.C >= 2)
        {
            Svd = WitnessCache.GetOrComputeSvd(block.C, block.N, block.GammaZero);
            double gEff = Svd.Sigma0;
            QEp = new QEpLaw(gEff);
            // Traverse pre/at/post EP at Q ∈ {0.5·Q_EP, Q_EP, 1.5·Q_EP} for slowest-pair k=1 inspection.
            double qEp = QEp.Value;
            EpTraversal = new[]
            {
                TwoLevelEpModel.AtQ(block.GammaZero, qEp * 0.5, gEff, k: 1),
                TwoLevelEpModel.AtQ(block.GammaZero, qEp,         gEff, k: 1),
                TwoLevelEpModel.AtQ(block.GammaZero, qEp * 1.5,   gEff, k: 1),
            };
            // Higher-k EP hierarchy: one model per k ∈ {1, …, c−1} at the EP, showing the
            // 1/(4γ₀·k) decay-time spectrum.
            var higher = new TwoLevelEpModel[block.C - 1];
            for (int k = 1; k <= block.C - 1; k++)
                higher[k - 1] = TwoLevelEpModel.AtQ(block.GammaZero, qEp, gEff, k: k);
            HigherKLevels = higher;
        }
        else
        {
            EpTraversal = Array.Empty<TwoLevelEpModel>();
            HigherKLevels = Array.Empty<TwoLevelEpModel>();
        }

        // Witnesses are self-computing via the shared cache; the prediction's expected
        // ratio is the empirical anchor (mean across the historical step_f/g sweeps),
        // not derived live from the witnesses themselves to avoid triggering ~18 scans
        // every time someone reads the prediction's Summary.
        InteriorShape = new UniversalShapePrediction(
            BondClass.Interior,
            expectedRatio: 0.756,
            tolerance: 0.005,
            witnesses: BuildInteriorWitnesses(block.GammaZero, WitnessCache));
        EndpointShape = new UniversalShapePrediction(
            BondClass.Endpoint,
            expectedRatio: 0.770,
            tolerance: 0.005,
            witnesses: BuildEndpointWitnesses(block.GammaZero, WitnessCache));

        InteriorShapeFunction = ShapeFunctionWitnesses.BuildInterior(block.GammaZero, WitnessCache);
        EndpointShapeFunction = ShapeFunctionWitnesses.BuildEndpoint(block.GammaZero, WitnessCache);

        PerBlockQPeaks = PerBlockQPeakClaim.Standard;
        EndpointPerBondTable = PerBondQPeakWitnessTable.BuildEndpoint(block.GammaZero, WitnessCache);
        InteriorPerBondTable = PerBondQPeakWitnessTable.BuildInterior(block.GammaZero, WitnessCache);
        PerOrbitSubstructure = new PerF71OrbitObservation();
        QAnchors = new QAnchorMap();
        _qAnchorsChildrenGroup = new Lazy<InspectableNode>(() =>
            InspectableNode.Group("named Q-anchors (Q = J/γ₀ axis)",
                QAnchors.Anchors.Select(a =>
                    (IInspectable)new InspectableNode(
                        displayName: $"Q = {a.Q:F3}",
                        summary: $"{a.Role}, J(γ₀=0.05) = {a.JAtGamma0Point05:G4}, θ = {a.ThetaDegrees():F1}°, {a.Tier.Label()}"))
                    .ToArray()));

        DressedModeWeight = new DressedModeWeightClaim();
        AlgebraicClass = new ChiralAiiiClassification();
        F71Mirror = new F71MirrorInvariance();
        Sigma0Scaling = new SigmaZeroChromaticityScaling(block.GammaZero, cache: WitnessCache);
        F86HwhmClosedForm = new F86HwhmClosedFormClaim();

        // c=2-specific top-level synthesis of Stages A–D (Stage E1 integration). Only built
        // for c=2 blocks; null otherwise — c≥3 strata fall under OpenQuestions Item 4'
        // (multi-k extension to c≥3, out of scope for the c=2 derivation plan). Lazy so the
        // full-block Q-scan inside C2UniversalShapeDerivation.Build runs only on first access.
        // Shared WitnessCache so OrbitKTable below reuses the C2HwhmRatio result.
        _c2UniversalShape = new Lazy<C2UniversalShapeDerivation?>(() =>
            block.C == 2 ? C2UniversalShapeDerivation.Build(block, WitnessCache) : null);

        // Block-independent meta-claim wiring the F86 local EP at real Q_EP to
        // FRAGILE_BRIDGE's complex-γ-plane EP under shared AIII chiral algebra.
        // Available for any block — the algebraic statement is shared across all c;
        // the pinned witnesses are the c=2 N=5..8 Petermann-K sweep.
        _localGlobalEpLink = new Lazy<LocalGlobalEpLink>(() => LocalGlobalEpLink.Build());

        // Block-independent meta-claim (Locus 6, symmetry-side closure): F86 c=2 bond-class
        // split inherits from the polarity-layer pair {−0.5, +0.5} at d=2 in Pi2KnowledgeBase.
        // Companion to LocalGlobalEpLink (Locus 5, EP-side); the inheritance statement is
        // shared across all c, the pinned witnesses pin the c=2 N=5..8 numbers.
        _polarityInheritanceLink = new Lazy<PolarityInheritanceLink>(() => PolarityInheritanceLink.Build());

        // Block-independent meta-claim (Tier1Derived): the L_eff mirror axis −4γ₀ = −2γ₀·2,
        // the channel pair as a mirror pair about it, the EP as the coalescence, g_eff as the
        // branch's relative clock in Im(λ). Pure 2×2 algebra; cheap to build, lazy for
        // consistency with the sibling meta-claims.
        _lEffMirrorAxis = new Lazy<LEffMirrorAxisClaim>(() => LEffMirrorAxisClaim.Build());

        // Live F71-orbit K-resonance witness table for c=2. Lazy so the C2HwhmRatio.Build
        // cost is paid only on first read; shared WitnessCache reuses any C2HwhmRatio that
        // C2UniversalShape already built.
        _orbitKTable = new Lazy<PerF71OrbitKTable?>(() =>
            block.C == 2 ? PerF71OrbitKTable.Build(block, WitnessCache) : null);

        // F89 per-Bloch-mode sigma anatomy (Tier2Verified): extracts σ_n via R†·S·R diagonal
        // for all F_a modes of the c=2 uniform-J block-L. Bit-exact vs F89UnifiedFaClosedFormClaim
        // for path-3..7 (path-7 closed form derived from this anatomy).
        // Only meaningful at c=2 (throws otherwise); null for all other c.
        _fullBlockSigmaAnatomy = new Lazy<C2FullBlockSigmaAnatomy?>(() =>
            block.C == 2 ? C2FullBlockSigmaAnatomy.Build(block) : null);

        // Block-independent Tier-2-Verified table: IBM 2026-04-26 framework_snapshots
        // through Theorem 2's C_block lens. Static-data Claim, no compute cost; lazy
        // construction only on first read. The Quarter parent is the universal 1/4
        // ceiling against which every witness row is asserted.
        _ibmBlockCpsiHardwareTable = new Lazy<IbmBlockCpsiHardwareTable>(
            () => new IbmBlockCpsiHardwareTable(new Symmetry.QuarterAsBilinearMaxvalClaim()));

        // F90 bridge identity (Tier1Derived): F86 c=2 K_b ↔ F89 path-(N−1)
        // per-bond Hellmann-Feynman. Block-independent — the algebraic statement is
        // c=2-stratum-wide and N-independent (verified bit-exact at N=5..8). Cheap to build
        // (constructor only stores typed parent edges); lazy for consistency with sibling
        // meta-claims. The F89 parent chain is built inline since F86KnowledgeBase doesn't
        // receive Pi2/F89 instances; one fresh allocation per F86KB.
        _f89BridgeIdentity = new Lazy<F90F86C2BridgeIdentity>(() =>
        {
            var f89 = new F89TopologyOrbitClosure(new Pi2DyadicLadderClaim());
            var atLock = new F89PathKAtLockMechanismClaim(f89);
            return new F90F86C2BridgeIdentity(f89, atLock);
        });

        // F86e σ_0 = ‖[Π_HD1, M_H]‖ identity (Tier1Derived): c=2-specific (needs
        // Π_HD1 + Π_HD3 = I). Lazy so the InterChannelSvd + square-projector commutator
        // norm are paid only on first read; null for all c ≠ 2 blocks.
        _sigmaZeroCommutatorNorm = new Lazy<SigmaZeroCommutatorNormClaim?>(() =>
            block.C == 2 ? SigmaZeroCommutatorNormClaim.Build(block) : null);

        Retracted = RetractedClaim.Standard;
        OpenQuestions = F86OpenQuestions.Standard;

        FourModeInsufficiencyNote = new InspectableNode(
            "4-mode minimal effective insufficient",
            summary: "Interior HWHM/Q ≈ 0.74 partially preserved, Q_peak shifted ~2× and Endpoint goes off-grid; more modes needed for full Tier 1 (PROOF_F86_QPEAK Item 1')");
    }

    /// <summary>Compare measured Interior + Endpoint peaks against the predictions; returns
    /// inspectable match objects for both bond classes.</summary>
    public IReadOnlyList<PredictionMatch> CompareTo(KCurve measured) => new[]
    {
        InteriorShape.CompareTo(measured.Peak(BondClass.Interior)),
        EndpointShape.CompareTo(measured.Peak(BondClass.Endpoint)),
    };

    public string DisplayName =>
        $"F86 knowledge base (c={Block.C}, N={Block.N}, n={Block.LowerPopcount}, γ₀={Block.GammaZero:G3})";

    public string Summary =>
        $"t_peak={TPeak.Value:G4}" +
        (QEp is not null ? $", Q_EP={QEp.Value:G4}" : "") +
        $", universal: Interior {InteriorShape.ExpectedHwhmOverQPeak:F3}, Endpoint {EndpointShape.ExpectedHwhmOverQPeak:F3}, " +
        $"{Retracted.Count} retracted";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return new InspectableNode("Block (CoherenceBlock)",
                summary: $"N={Block.N}, n={Block.LowerPopcount}, c={Block.C}, γ₀={Block.GammaZero:G3}");

            yield return _qAnchorsChildrenGroup.Value;

            yield return InspectableNode.Group("Tier 1 (derived)",
                CollectTier1Derived().ToArray());

            yield return InspectableNode.Group("Tier 1 (candidate)",
                CollectTier1Candidate().ToArray());

            yield return InspectableNode.Group("Tier 2",
                CollectTier2().ToArray());

            yield return InspectableNode.Group("retracted",
                Retracted.Cast<IInspectable>().ToArray());

            yield return InspectableNode.Group("open questions",
                OpenQuestions.Cast<IInspectable>()
                    // LocalGlobalEpLink (OpenQuestion since the F86a-retraction 2026-06-21):
                    // whether the full Σγ=N·γ₀ block has an off-axis defective EP at all.
                    .Append(LocalGlobalEpLink)
                    .ToArray());

            yield return FourModeInsufficiencyNote;
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;

    private IEnumerable<IInspectable> CollectTier1Derived()
    {
        yield return TPeak;
        if (QEp is not null) yield return QEp;
        if (EpTraversal.Count > 0)
            yield return InspectableNode.Group("2-level EP traversal (pre / at / post, k=1)",
                EpTraversal.Cast<IInspectable>().ToArray());
        if (HigherKLevels.Count > 1)
            yield return InspectableNode.Group("higher-k EP hierarchy (decay times 1/(4γ₀·k))",
                HigherKLevels.Cast<IInspectable>().ToArray());
        yield return AlgebraicClass;
        yield return F71Mirror;
        yield return LEffMirrorAxis;
        yield return PolarityPairQPeakDecomposition;
        // F90 bridge identity is c=2-stratum specific (only the c=2 K_b ↔ F89 path-(N−1)
        // identity is verified bit-exact). Surface it under Tier 1 derived for c=2 blocks.
        if (Block.C == 2) yield return F89BridgeIdentity;
        // σ_0 = ‖[Π_HD1, M_H]‖ is c=2-specific (needs Π_HD1 + Π_HD3 = I). Non-null only at c=2.
        if (SigmaZeroCommutatorNorm is not null) yield return SigmaZeroCommutatorNorm;
    }

    private IEnumerable<IInspectable> CollectTier1Candidate()
    {
        yield return InteriorShape;
        yield return EndpointShape;
        yield return InteriorShapeFunction;
        yield return EndpointShapeFunction;
        yield return DressedModeWeight;
        yield return F86HwhmClosedForm;
        // c=2 top-level synthesis from Stages A–D: only present for c=2 blocks.
        if (C2UniversalShape is not null)
            yield return C2UniversalShape;
    }

    private IEnumerable<IInspectable> CollectTier2()
    {
        if (OrbitKTable is not null) yield return OrbitKTable;
        if (FullBlockSigmaAnatomy is not null) yield return FullBlockSigmaAnatomy;

        yield return Sigma0Scaling;
        yield return InspectableNode.Group("per-block Q_peak (Q_SCALE convention)",
            PerBlockQPeaks.Cast<IInspectable>().ToArray());
        yield return EndpointPerBondTable;
        yield return InteriorPerBondTable;
        yield return PerOrbitSubstructure;
        yield return PolarityInheritanceLink;
        yield return IbmBlockCpsiHardwareTable;
        // LocalGlobalEpLink is NOT yielded here: it was demoted Tier2Verified → OpenQuestion
        // by the F86a-retraction 2026-06-21, so it now displays under the "open questions"
        // group (see Children) to keep the tier-group labels honest (guarded by
        // F86KnowledgeBase_TierGroups_HoldOnlyClaimsOfMatchingTier).
    }

    /// <summary>Build the Interior witness list — each witness computes its HWHM/Q from the
    /// cache on demand. Use <paramref name="cache"/> to share across multiple knowledge bases.</summary>
    public static IReadOnlyList<UniversalShapeWitness> BuildInteriorWitnesses(
        double gammaZero = 0.05, WitnessCache? cache = null) =>
        F86StandardLocations.Full
            .Select(loc => new UniversalShapeWitness(loc.C, loc.N, gammaZero, BondClass.Interior, cache))
            .ToArray();

    /// <summary>Build the Endpoint witness list — same self-computing pattern.</summary>
    public static IReadOnlyList<UniversalShapeWitness> BuildEndpointWitnesses(
        double gammaZero = 0.05, WitnessCache? cache = null) =>
        F86StandardLocations.EndpointDefault
            .Select(loc => new UniversalShapeWitness(loc.C, loc.N, gammaZero, BondClass.Endpoint, cache))
            .ToArray();
}
