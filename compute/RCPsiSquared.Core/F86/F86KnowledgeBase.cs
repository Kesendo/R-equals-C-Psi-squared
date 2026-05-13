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
/// <para>What is in here:</para>
/// <list type="bullet">
///   <item>Tier-1 derived: <see cref="TPeakLaw"/> (universal), <see cref="QEpLaw"/> at the
///         block's natural g_eff (= σ_0 from <see cref="InterChannelSvd"/>),
///         <see cref="TwoLevelEpModel"/> traversed at multiple Q values to show pre/at/post-EP.</item>
///   <item>Tier-1 candidate: <see cref="UniversalShapePrediction"/> for Interior + Endpoint,
///         each with empirical <see cref="UniversalShapeWitness"/> data points across c=2..4
///         N=5..8. For c=2 specifically, <see cref="C2UniversalShapeDerivation"/> wraps
///         the Stage D2 <c>C2HwhmRatio</c> with empirical anchor + directional Endpoint &gt;
///         Interior split derived.</item>
///   <item>Retracted: the two refuted closed forms (<see cref="RetractedClaim.Standard"/>) —
///         the PTF-lesson reminder.</item>
///   <item>4-mode insufficiency note: 2026-05-02 finding that the minimal effective fails
///         to reproduce the universal shape numerically — promoted from open-question to
///         answered-no.</item>
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
    public DressedModeWeightClaim DressedModeWeight { get; }
    public ChiralAiiiClassification AlgebraicClass { get; }
    public F71MirrorInvariance F71Mirror { get; }
    public SigmaZeroChromaticityScaling Sigma0Scaling { get; }
    /// <summary>F86 c=2 HWHM_left/Q_peak closed-form prediction per BondSubClass
    /// (Tier 1 derived; 2026-05-13). Closes Item 1' from PROOF_F90_F86C2_BRIDGE.md via the
    /// bare 4-mode floor + linear lift alpha*g_eff + beta composition. Block-independent
    /// meta-claim; the per-sub-class fit is anchored to c=2 N=5..8 F90-bridge data.</summary>
    public F86HwhmClosedFormClaim F86HwhmClosedForm { get; }
    /// <summary>The c=2-specific top-level Claim wrapping the Stage D2
    /// <see cref="Item1Derivation.C2HwhmRatio"/>: empirical anchor + directional
    /// Endpoint &gt; Interior split, both pinned for the c=2 stratum (witness anchor:
    /// c=2 N=5..8 in the canonical sweep). Non-null iff <c>Block.C == 2</c>.
    /// Lazily built on first access — the underlying full-block Q-scan
    /// (153×21 eigendecompositions) is paid only by consumers that read the property.</summary>
    public C2UniversalShapeDerivation? C2UniversalShape => _c2UniversalShape.Value;

    private readonly Lazy<C2UniversalShapeDerivation?> _c2UniversalShape;

    /// <summary>Block-independent meta-claim: the F86 local EP at real Q_EP and the
    /// FRAGILE_BRIDGE global EP in the complex-γ plane are the same exceptional-point
    /// structure under shared AIII chiral algebra, validated empirically by the c=2
    /// N=5..8 Petermann-K sweep (2026-05-06). <see cref="Knowledge.Tier.Tier2Verified"/>:
    /// shared algebra + real-axis hit empirically pinned; complex-γ analytic continuation
    /// is the documented gap. Exposed at the KB root for any block, not just c=2 — the
    /// algebraic statement is shared across all c. Lazy: a single static-data Claim with
    /// no compute cost.</summary>
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

    /// <summary>F71-orbit-grouped live K-resonance witness table for c=2
    /// (Tier2Verified). Live counterpart to <see cref="PerOrbitSubstructure"/>
    /// (frozen 9-case sweep). Non-null iff <c>Block.C == 2</c>. Lazily built on
    /// first access — the underlying full-block Q-scan inside
    /// <see cref="PerF71OrbitKTable.Build"/> runs only once and only on read.</summary>
    public PerF71OrbitKTable? OrbitKTable => _orbitKTable.Value;

    /// <summary>F89 per-Bloch-mode σ_n extraction via R†·S·R diagonal at the c=2 stratum
    /// (Tier 2 verified). Bit-exact match against
    /// <see cref="F89UnifiedFaClosedFormClaim.Sigma"/> for path-3..6 (boundary case
    /// nBlock = N = k+1, no bare site); provides path-7 data at c=2 N=9 (no
    /// analytical oracle). <c>null</c> for c ≠ 2 blocks.</summary>
    public C2FullBlockSigmaAnatomy? FullBlockSigmaAnatomy { get; }

    private readonly Lazy<PerF71OrbitKTable?> _orbitKTable;

    /// <summary>Block-independent meta-claim (Tier2Verified): the IBM 2026-04-26
    /// framework_snapshots viewed through Theorem 2's universal C_block ≤ 1/4 ceiling
    /// (PROOF_BLOCK_CPSI_QUARTER). Pinned 32-row witness table over (Aer / Marrakesh /
    /// Kingston / Fez) × (heisenberg / truly / soft / hard) × ((0,1) / (1,2)) blocks of
    /// the (q0, q2) reduced ρ on the N=3 |+−+⟩ chain at t=0.8, J=1.0. Independent of
    /// <see cref="Block"/>; the lens applies at every c-block.</summary>
    public IbmBlockCpsiHardwareTable IbmBlockCpsiHardwareTable => _ibmBlockCpsiHardwareTable.Value;

    private readonly Lazy<IbmBlockCpsiHardwareTable> _ibmBlockCpsiHardwareTable;

    /// <summary>F90 c=2 ↔ F89 bridge identity (Tier1Derived 2026-05-11): F86 c=2 K_b(Q,t)
    /// IS F89 path-(N−1) (SE,DE) per-bond Hellmann-Feynman, modulo F89-J = 2·F86-J convention.
    /// Bit-exact verified at 20/22 bonds across N=5..8 including orbit escapes. Resolves
    /// Direction (b'') (full block-L derivation, NOT 4-mode) numerically Tier-1; closed-form
    /// HWHM_left/Q_peak per bond class via F89's AT-locked F_a/F_b structure remains the
    /// analytical target. Block-independent meta-claim — the bridge identity holds for any
    /// c=2 block N ≥ 3. See <c>docs/proofs/PROOF_F90_F86C2_BRIDGE.md</c>.</summary>
    public F90F86C2BridgeIdentity F89BridgeIdentity => _f89BridgeIdentity.Value;

    private readonly Lazy<F90F86C2BridgeIdentity> _f89BridgeIdentity;

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

        // Live F71-orbit K-resonance witness table for c=2. Lazy so the C2HwhmRatio.Build
        // cost is paid only on first read; shared WitnessCache reuses any C2HwhmRatio that
        // C2UniversalShape already built.
        _orbitKTable = new Lazy<PerF71OrbitKTable?>(() =>
            block.C == 2 ? PerF71OrbitKTable.Build(block, WitnessCache) : null);

        // F89 per-Bloch-mode sigma anatomy (Tier2Verified): extracts σ_n via R†·S·R diagonal
        // for all F_a modes of the c=2 uniform-J block-L. Bit-exact vs F89UnifiedFaClosedFormClaim
        // for path-3..6; lifts to path-7 (c=2 N=9) where no analytical oracle exists.
        // Only meaningful at c=2 (throws otherwise); null for all other c.
        FullBlockSigmaAnatomy = block.C == 2 ? C2FullBlockSigmaAnatomy.Build(block) : null;

        // Block-independent Tier-2-Verified table: IBM 2026-04-26 framework_snapshots
        // through Theorem 2's C_block lens. Static-data Claim, no compute cost; lazy
        // construction only on first read.
        _ibmBlockCpsiHardwareTable = new Lazy<IbmBlockCpsiHardwareTable>(() => new IbmBlockCpsiHardwareTable());

        // F90 bridge identity (Tier1Derived 2026-05-11): F86 c=2 K_b ↔ F89 path-(N−1)
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

            yield return InspectableNode.Group("Tier 1 (derived)",
                CollectTier1Derived().ToArray());

            yield return InspectableNode.Group("Tier 1 (candidate)",
                CollectTier1Candidate().ToArray());

            yield return InspectableNode.Group("Tier 2 (empirical)",
                CollectTier2Empirical().ToArray());

            yield return InspectableNode.Group("retracted",
                Retracted.Cast<IInspectable>().ToArray());

            yield return InspectableNode.Group("open questions",
                OpenQuestions.Cast<IInspectable>().ToArray());

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
        yield return DressedModeWeight;
        yield return AlgebraicClass;
        yield return F71Mirror;
        yield return F86HwhmClosedForm;
        // F90 bridge identity is c=2-stratum specific (only the c=2 K_b ↔ F89 path-(N−1)
        // identity is verified bit-exact). Surface it under Tier 1 derived for c=2 blocks.
        if (Block.C == 2) yield return F89BridgeIdentity;
    }

    private IEnumerable<IInspectable> CollectTier1Candidate()
    {
        yield return InteriorShape;
        yield return EndpointShape;
        yield return InteriorShapeFunction;
        yield return EndpointShapeFunction;
        yield return Sigma0Scaling;
        // c=2 top-level synthesis from Stages A–D: only present for c=2 blocks.
        if (C2UniversalShape is not null)
            yield return C2UniversalShape;
    }

    private IEnumerable<IInspectable> CollectTier2Empirical()
    {
        if (OrbitKTable is not null) yield return OrbitKTable;
        if (FullBlockSigmaAnatomy is not null) yield return FullBlockSigmaAnatomy;

        yield return InspectableNode.Group("per-block Q_peak (Q_SCALE convention)",
            PerBlockQPeaks.Cast<IInspectable>().ToArray());
        yield return EndpointPerBondTable;
        yield return InteriorPerBondTable;
        yield return PerOrbitSubstructure;
        yield return LocalGlobalEpLink;
        yield return PolarityInheritanceLink;
        yield return IbmBlockCpsiHardwareTable;
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
