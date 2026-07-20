using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.F86;

/// <summary>F86 open question (was an F86a meta-claim, demoted 2026-06-21): whether the
/// full Σγ = N·γ₀ (n,n+1)-coherence block has a defective exceptional point OFF the real
/// Q axis AT ALL, and — if so — whether it is the same EP structure as the firmly-established
/// genuine EPs (the toy 2×2 reduction and the SEPARATE Σγ = 0 gain-loss system,
/// FRAGILE_BRIDGE, complex-γ, K=403).
///
/// <para><b>Tier outcome: OpenQuestion</b> (demoted from Tier2Verified 2026-06-21 by the
/// F86a-retraction review; the retraction itself corrected in part 2026-07-07, the CORRECTED
/// para below). The prior reading — "F86's local EP at real Q_EP is a real-axis hit of the
/// same EP FRAGILE_BRIDGE detects at K=403" — does NOT survive an artifact-free
/// re-verification from below: at every point of the Petermann sweep's Q grid the eigenvalues
/// are simple (nearest-neighbour gap ~0.25–0.35), and the sweep's peak magnitudes are grid
/// artifacts. The real-axis defective structure the block DOES carry (the F89 seeds, CORRECTED
/// para below) sits between the grid points, at scattered q*, not at Q_EP: nothing is "hit"
/// at the resonance peak.</para>
///
/// <para><b>What is genuine (the phenomenon stays):</b> the block IS strongly NON-NORMAL on
/// the real axis (cond(V) = 48.7 / 50.9 / 268.5 at N = 5 / 6 / 7 — the opposite of a diabolic,
/// i.e. normal, crossing). This is NOT an <c>eig</c> artifact and NOT a degenerate-eigenspace
/// effect: the artifact-free Riesz spectral-projector norm ‖P‖ reproduces the large Petermann
/// factor on a SIMPLE, isolated eigenvalue at Re ≈ −4γ₀ (N=5: ‖P‖ = 19.4 = √375, gap ~0.25).
/// The Petermann factor is large but FINITE at every grid point sampled: a near-EP shadow,
/// whose demonstrated caster is the F89 real-axis seed nearby (at N=5 the seed at carrier
/// Q ≈ 1.286 sits ~1.9e-3 from the spike location Q ≈ 1.288; CORRECTED para below).</para>
///
/// <para><b>What is retracted (the numbers go):</b> (i) reading the real-axis Petermann sweep
/// as a defective EP on the real axis (there is no coalescence there), and (ii) the peak
/// MAGNITUDES and the law built on them — K swings 2–4× over ΔQ = 1e-3 (grid-sensitive), so
/// "6× above FRAGILE_BRIDGE", the specific peak K=2384.7, and the "within-parity monotonic
/// growth law" + the odd-dominates-even parity-asymmetry reading are grid artifacts and are
/// dropped. The four <see cref="PetermannSpikeWitness"/> rows are RETAINED as a cautionary
/// record of genuine non-normality near each Q_peak, NOT as EP evidence and NOT as a magnitude
/// law.</para>
///
/// <para><b>The open question:</b> whether the full Σγ = N·γ₀ block has an off-axis defective
/// EP at all is OPEN. A search on 2026-06-21 found the nearest complex-Q eigenvalue
/// coalescences of the full block are themselves DIABOLIC (‖P‖ = 1, departure-from-normality
/// = 0), not defective — so even off the real axis the same-EP-structure connection to the
/// toy 2×2 / FRAGILE_BRIDGE is not established.</para>
///
/// <para><b>Scope (do not conflate with F89's octic — reconciled 2026-06-27):</b> this OPEN
/// question concerns the FULL N≥5 (1,2)/c=2 block (dim 50–224, the c=2 sweep above ran N=5–8)
/// and its FRAGILE_BRIDGE connection specifically. It is NOT F89's N=4 octic (1,2)-FACTOR
/// (dim 8), whose off-axis defective branch points ARE established (Tier-1): the squarefree
/// P₁₀ zeros of the λ-discriminant are defective √-branch transpositions generating S₈, while
/// the real-axis q_EP≈0.659 is the diabolic (3q⁴+q²−1) double-zero (silent) — see
/// <c>F89OcticMonodromyClaim</c> / <c>inspect --root galoismonodromy</c> (G2: "the silent
/// diabolic vs the braiding EPs, side by side"). The "nearest coalescences are diabolic"
/// finding here is consistent with that: those are the abundant free-fermion sum-coincidences
/// (codim-1 by additivity, <c>reference_nonhermitian_diabolic_codimension</c>), not the rare
/// S_d-generating defective branch points. The defective/diabolic split is read off the
/// discriminant's squarefree/square factorization — the same object as the Galois parity split.
/// This claim's own real-axis Petermann sweep read a non-normal but SIMPLE eigenvalue at Q_peak
/// (the shadow of an off-axis EP); the defective-vs-diabolic rule (single-particle coalescence vs
/// frequency-difference coincidence) and the cross-block seam live in <c>EpCharacterWitness</c>
/// (siblings: <c>CoherenceHorizonClaim</c>, <c>F89Path3OcticEpClaim</c>).</para>
///
/// <para><b>CORRECTED 2026-07-07 (the real-axis point is NOT uniformly simple):</b> the "SIMPLE
/// eigenvalue at Q_peak" reading above was Q_peak-local and grid-coarse. F89's exact nullity count
/// (r(0+) - r(inf) = N - 1) proves a real-to-complex transition on this SAME full (1,2) block at every
/// odd N, and the Kato simple-zero lemma makes it a real-axis DEFECTIVE EP (a Jordan block) at every seed
/// tested, census-confirmed to N=11 (beta-exotic-scoped to Puiseux p ~ 0.5; only the codim-2 beta-exotic
/// genericity stays open for all N). Those seeds sit at scattered q* the Petermann sweep's dQ ~ 0.029 grid (121 pts over
/// [0.5, 4]) never bisected onto: a defective sqrt-EP splits its pair by ~sqrt|q - q*|, visible only
/// within a window |q - q*| &lt; ~1e-3, some 20-30x narrower than the grid step (shown from below in
/// <c>F86aSeedMaskingTests</c>: at the N=9 seed q*=0.849011 the pair is a Jordan block, gap ~2.7e-3, but
/// ~0.35 apart one grid step away). So the block IS defective on the real axis at the seeds; the F89
/// grid-robust count-change detector, not this Petermann sweep, is the instrument that sees them. What
/// remains genuinely open is the DISTINCT off-real-axis complex-Q EP structure. See
/// <c>SeedHolonomyClaim</c>, <c>RealSeedCensusTests</c>, and PROOF_F86A_EP_MECHANISM section The real-axis EP.</para>
///
/// <para>Block-independent meta-claim: registered at the F86 KB root for any block, not
/// just c=2 (similar to <see cref="ChiralAiiiClassification"/>). The witnesses pin the
/// concrete c=2 sweep that motivated the original connection; the open question itself is
/// shared across all c.</para>
///
/// <para>Anchors: <c>docs/proofs/PROOF_F86A_EP_MECHANISM.md</c> (F86 local instance, Statement 1),
/// <c>hypotheses/FRAGILE_BRIDGE.md</c> (separate Σγ=0 gain-loss instance, K=403), and
/// <c>experiments/PT_SYMMETRY_ANALYSIS.md</c> (AIII chiral classification, the shared
/// algebraic substrate; its Petermann reading carries the same two-half verdict).</para>
/// </summary>
public sealed class LocalGlobalEpLink : Claim
{
    /// <summary>Reference parameters for the empirical Petermann-K sweep. The
    /// <see cref="Witnesses"/> below are pinned to these values; a future refresh
    /// of the witnesses requires re-running
    /// <c>F86PetermannProbe.Probe_PetermannFineGrid_C2_VsN</c> at these settings.
    /// γ₀=0.05, c=2 block (1,2), Q-grid <see cref="SweepQPoints"/> points uniform
    /// on [<see cref="SweepQMin"/>, <see cref="SweepQMax"/>].</summary>
    public const double SweepGammaZero = 0.05;
    /// <inheritdoc cref="SweepGammaZero"/>
    public const double SweepQMin = 0.50;
    /// <inheritdoc cref="SweepGammaZero"/>
    public const double SweepQMax = 4.00;
    /// <inheritdoc cref="SweepGammaZero"/>
    public const int SweepQPoints = 121;

    /// <summary>FRAGILE_BRIDGE Petermann K=403 ballpark — the near-singularity reference
    /// value of the SEPARATE Σγ = 0 gain-loss system (complex-γ) from
    /// <c>hypotheses/FRAGILE_BRIDGE.md</c>. Kept as a documentary reference to one of the
    /// two firmly-established genuine EPs; the prior "the real-axis K-sweep sits ~6× above
    /// this" comparison was retracted 2026-06-21 (the magnitudes are grid-sensitive; the
    /// block's real-axis defective seeds sit off the sweep grid, see the class summary).</summary>
    public const double FragileBridgeKReference = 403.0;

    /// <summary>Same-sign-imaginary 2×2 algebra, the shared algebraic object:
    /// L_eff − (trace/2)·I = [[−Δ/2, +iJ·g_eff], [+iJ·g_eff, +Δ/2]] (PROOF_F86A_EP_MECHANISM
    /// Statement 1). AIII chiral classification (PT_SYMMETRY_ANALYSIS): linear Π, Π⁴=I,
    /// anti-commutation {Π, L_c}=0.</summary>
    public string SharedAlgebra => "same-sign-imaginary 2×2; AIII chiral";

    /// <summary>F86 local instance (Σγ = N·γ₀ ≠ 0, real Q_EP = 2/g_eff).</summary>
    public string LocalInstanceAnchor => "F86 Statement 1, PROOF_F86A_EP_MECHANISM.md";

    /// <summary>FRAGILE_BRIDGE global instance (Σγ = 0, complex-γ-plane EP, Petermann K=403).</summary>
    public string GlobalInstanceAnchor => "hypotheses/FRAGILE_BRIDGE.md";

    /// <summary>Cautionary record from the c=2 N=5..8 Petermann-K sweep
    /// (pinned to <see cref="SweepGammaZero"/>, <see cref="SweepQMin"/>,
    /// <see cref="SweepQMax"/>, <see cref="SweepQPoints"/>; c=2 block (1,2);
    /// 2026-05-06). One row per N. RETAINED to document genuine non-normality near each
    /// Q_peak (large but finite Petermann factor), NOT as EP evidence: the peak magnitudes
    /// are grid-sensitive (K swings 2–4× over ΔQ = 1e-3) and were retracted 2026-06-21, so
    /// these numbers do not support a defective-EP reading or a growth law. Hard-coded as a
    /// static pinned table — a frozen record, not a live-recomputed Claim.</summary>
    public IReadOnlyList<PetermannSpikeWitness> Witnesses => _petermannWitnesses;

    /// <summary>The 2026-06-21 F86a-retraction correction and the resulting open question.
    /// Records what was retracted (the real-axis defective-EP reading and the grid-sensitive
    /// magnitudes), what survives (genuine non-normality, artifact-free via the Riesz
    /// projector norm), and the open question (whether the full block has an off-axis
    /// defective EP at all).</summary>
    public string PendingDerivationNote =>
        "Correction (2026-06-21, independently re-verified from below, artifact-free): the " +
        "full Σγ = N·γ₀ (n,n+1) block has NO eigenvalue coalescence on the real Q axis — its " +
        "eigenvalues stay SIMPLE there (nearest-neighbour gap ~0.25–0.35), so there is no " +
        "defective EP on the real axis. The block IS genuinely strongly NON-NORMAL on the " +
        "real axis (cond(V) = 48.7 / 50.9 / 268.5 at N = 5 / 6 / 7); this is NOT an eig " +
        "artifact: the artifact-free Riesz spectral-projector norm ‖P‖ reproduces the large " +
        "Petermann factor on a SIMPLE, isolated eigenvalue at Re ≈ −4γ₀ (N=5: ‖P‖ = 19.4 = " +
        "√375, gap ~0.25). Adopting PT_SYMMETRY_ANALYSIS verbatim: no real-axis EP; the " +
        "Petermann factor is large but FINITE on the real axis, signalling a nearby EP in " +
        "the complex parameter plane. RETRACTED: (i) reading the sweep as a defective EP on " +
        "the real axis, and (ii) the peak magnitudes and the law on them — K swings 2–4× over " +
        "ΔQ = 1e-3, so '6×', the peak K=2384.7, and the within-parity monotonic growth / " +
        "odd-dominates-even parity-asymmetry readings are grid artifacts. The firmly-" +
        "established genuine defective EPs are ONLY the toy 2×2 reduction (PROOF_F86A_EP_MECHANISM " +
        "Statement 1) and the SEPARATE Σγ = 0 gain-loss system (FRAGILE_BRIDGE, complex-γ, " +
        "K=403). OPEN QUESTION: whether the full Σγ = N·γ₀ block has an off-axis defective EP " +
        "AT ALL — a 2026-06-21 search found the nearest complex-Q coalescences of the full " +
        "block are themselves DIABOLIC (‖P‖ = 1, departure-from-normality = 0), not " +
        "defective. The surviving shared algebraic substrate (AIII chiral classification, " +
        "linear Π with Π⁴=I and {Π, L_c}=0) is unchanged; what is unestablished is the " +
        "defective-EP connection between the full block and the two genuine EPs. " +
        "SCOPE (2026-06-27 reconciliation): this concerns the FULL N≥5 (1,2) block ↔ " +
        "FRAGILE_BRIDGE specifically, NOT F89's N=4 octic (1,2)-factor, whose off-axis " +
        "defective branch points ARE established (F89OcticMonodromyClaim, the squarefree-P₁₀ " +
        "simple-zeros = transpositions generating S₈; its real-axis q_EP≈0.659 is the diabolic " +
        "(3q⁴+q²−1) double-zero). The nearest-diabolic coalescences here are the abundant " +
        "free-fermion sum-coincidences (codim-1 by additivity), not those rare S_d-generators. " +
        "FURTHER CORRECTED (2026-07-07): the 'NO real-axis defective EP' conclusion above is " +
        "superseded. F89's exact nullity count (r(0⁺)−r(∞)=N−1) proves a real-to-complex transition " +
        "on this full (1,2) block at every odd N, and the Kato simple-zero lemma makes it defective " +
        "(a Jordan block) at every seed tested, census-confirmed to N=11 (β-exotic-scoped to Puiseux " +
        "p≈0.5, only the codim-2 β-exotic genericity open for all N). The 2026-06-21 sit-and-" +
        "characterize scan missed it because a defective √-EP splits its pair by ~√|q−q*|, visible " +
        "only within |q−q*|≲1e-3, while the scan's ΔQ≈0.029 grid is 20-30× coarser and never sat " +
        "inside that window (shown from below, F86aSeedMaskingTests). What survives: the grid-" +
        "artifact call on the Petermann magnitudes; what stays genuinely open is the DISTINCT off-" +
        "real-axis complex-Q EP. See PROOF_F86A_EP_MECHANISM §The real-axis EP, SeedHolonomyClaim, " +
        "RealSeedCensusTests.";

    private LocalGlobalEpLink()
        : base("local–global EP link (F86 ↔ FRAGILE_BRIDGE) — OPEN",
               Tier.OpenQuestion,
               "docs/proofs/PROOF_F86A_EP_MECHANISM.md (F86 Statement 1, toy 2×2 EP) + " +
               "hypotheses/FRAGILE_BRIDGE.md (separate Σγ=0 gain-loss EP, K=403) + " +
               "experiments/PT_SYMMETRY_ANALYSIS.md (AIII chiral, the shared algebraic substrate) + " +
               "docs/superpowers/syntheses/2026-05-06-petermann-sweep-c2.md (cautionary non-normality record)")
    { }

    /// <summary>Public factory: builds the meta-claim with the four pinned Petermann
    /// spike witnesses. No <see cref="CoherenceBlocks.CoherenceBlock"/> required — the
    /// claim is block-independent (the witnesses pin the concrete c=2 sweep, but the
    /// algebraic statement is shared across all c).</summary>
    public static LocalGlobalEpLink Build() => new();

    public override string DisplayName =>
        "F86 ↔ FRAGILE_BRIDGE: shared EP under AIII chiral algebra";

    public override string Summary =>
        "OpenQuestion (F86a-retraction 2026-06-21, CORRECTED 2026-07-07): the full Σγ=N·γ₀ block " +
        "is genuinely non-normal on the real Q axis, and the Petermann magnitudes were grid-" +
        "sensitive artifacts (dropped). But the retraction's 'no real-axis defective EP' was " +
        "itself an over-correction: F89 proves a real-to-complex transition on this (1,2) block at " +
        "every odd N (nullity count), defective (a Jordan block, Kato) at every seed tested, census " +
        "to N=11 (β-exotic genericity still open for all N) — which the 2026-06-21 scan's ΔQ≈0.029 " +
        "grid missed, being 20-30× coarser than the √-EP window ~1e-3 (F86aSeedMaskingTests). What " +
        "stays OPEN besides that genericity is the DISTINCT off-real-axis complex-Q EP. Genuine EPs still separate: the " +
        "toy 2×2 and the Σγ=0 gain-loss system (FRAGILE_BRIDGE, K=403).";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("shared algebra", summary: SharedAlgebra);
            yield return new InspectableNode("local instance anchor",
                summary: LocalInstanceAnchor);
            yield return new InspectableNode("global instance anchor",
                summary: GlobalInstanceAnchor);
            yield return InspectableNode.Group(
                "Petermann-K cautionary record (c=2 N=5..8 real-axis sweep; genuine non-normality, magnitudes grid-sensitive, NOT EP evidence)",
                _petermannWitnesses.Cast<IInspectable>().ToArray());
            yield return new InspectableNode("correction + open question (F86a-retraction)",
                summary: PendingDerivationNote);
        }
    }

    /// <summary>Pinned cautionary table from the 2026-05-06 c=2 N=5..8 Petermann-K sweep.
    /// Pinned to <see cref="SweepGammaZero"/>, <see cref="SweepQMin"/>,
    /// <see cref="SweepQMax"/>, <see cref="SweepQPoints"/>; c=2 block (1,2).
    /// Order: N=5, 6, 7, 8. RETAINED as a record of genuine non-normality near each Q_peak;
    /// the peak magnitudes are grid-sensitive (K swings 2–4× over ΔQ = 1e-3, retracted
    /// 2026-06-21) and do NOT constitute EP evidence, a magnitude law, or a parity asymmetry.</summary>
    // Pinned data from F86PetermannProbe.Probe_PetermannFineGrid_C2_VsN; magnitudes are grid-sensitive (see class summary).
    private static readonly PetermannSpikeWitness[] _petermannWitnesses = new[]
    {
        new PetermannSpikeWitness(N: 5, BlockDim: 50,
            MaxKGlobal: 1333.6, ArgMaxQ: 1.288,
            MaxKInterior: 1333.6, MaxKEndpoint:  252.5),
        new PetermannSpikeWitness(N: 6, BlockDim: 90,
            MaxKGlobal:  337.9, ArgMaxQ: 0.938,
            MaxKInterior:  337.9, MaxKEndpoint:  131.7),
        new PetermannSpikeWitness(N: 7, BlockDim: 147,
            MaxKGlobal: 2384.7, ArgMaxQ: 1.842,
            MaxKInterior: 2384.7, MaxKEndpoint: 1341.2),
        new PetermannSpikeWitness(N: 8, BlockDim: 224,
            MaxKGlobal:  795.4, ArgMaxQ: 2.046,
            MaxKInterior:  580.0, MaxKEndpoint:  795.4),
    };
}

/// <summary>One row of the c=2 Petermann-K sweep — a frozen, cautionary data point
/// recording genuine non-normality near Q_peak. It is NOT EP evidence: the full block has
/// no eigenvalue coalescence on the real axis, and the peak magnitude is grid-sensitive
/// (retracted 2026-06-21; see <see cref="LocalGlobalEpLink"/>).
///
/// <para>Hosts the per-N spike data: the global Petermann factor maximum across the
/// Q-grid, the Q at which it peaked, and the per-bond-class breakdown. The
/// <see cref="Parity"/> label (derived from <paramref name="N"/>) is kept only as a
/// descriptive tag; the earlier "odd-N peaks dominate even-N by 2–4×" parity-asymmetry
/// reading was a grid artifact and is retracted (the magnitudes are unreliable).</para>
///
/// <para>Block dimension is the dimension of the c=2 block (n=1 → n=2 coherence
/// subspace) — verifies the sweep ran on the correct block.</para>
/// </summary>
/// <param name="N">Chain length.</param>
/// <param name="BlockDim">Dimension of the c=2 block at this N.</param>
/// <param name="MaxKGlobal">max Petermann K across the Q-grid (max over both bond
/// classes — the global maximum).</param>
/// <param name="ArgMaxQ">Q value at which <paramref name="MaxKGlobal"/> peaks.</param>
/// <param name="MaxKInterior">max Petermann K restricted to Interior bonds.</param>
/// <param name="MaxKEndpoint">max Petermann K restricted to Endpoint bonds.</param>
public sealed record PetermannSpikeWitness(
    int N,
    int BlockDim,
    double MaxKGlobal,
    double ArgMaxQ,
    double MaxKInterior,
    double MaxKEndpoint
) : IInspectable
{
    /// <summary>Parity class derived from <see cref="N"/>: "odd" or "even". Descriptive tag
    /// only; the earlier "odd-N peaks dominate even-N by 2–4×" reading was a grid artifact
    /// (retracted 2026-06-21 — the peak magnitudes are unreliable).</summary>
    public string Parity => N % 2 == 0 ? "even" : "odd";

    public string DisplayName => $"Petermann-K spike at c=2 N={N} ({Parity})";

    public string Summary => $"max K = {MaxKGlobal:F1} at Q = {ArgMaxQ:F3}";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return InspectableNode.RealScalar("N", N);
            yield return new InspectableNode("parity", summary: Parity);
            yield return InspectableNode.RealScalar("block dim (c=2)", BlockDim);
            yield return InspectableNode.RealScalar("max K (global)", MaxKGlobal, "F1");
            yield return InspectableNode.RealScalar("argmax Q", ArgMaxQ, "F3");
            yield return InspectableNode.RealScalar("max K (Interior)", MaxKInterior, "F1");
            yield return InspectableNode.RealScalar("max K (Endpoint)", MaxKEndpoint, "F1");
        }
    }

    public InspectablePayload Payload =>
        new InspectablePayload.Real($"max K (N={N}, {Parity})", MaxKGlobal, "F1");
}
