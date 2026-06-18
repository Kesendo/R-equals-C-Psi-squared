using System.Globalization;
using RCPsiSquared.Core.BlockSpectrum;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Numerics;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The chain/ring/star survivor trichotomy as one swept artifact: per (topology, N, Q, γ-profile)
/// the survivor sector + freeze-route + Δn + darkness + birth-canal deviation. Pure assembly of the
/// existing route-detectors (horizon / starseam / ceiling / surface); no new linear algebra, no new Claim.
/// Breadcrumbs StarFrozenSeamClaim / CoherenceHorizonClaim / HandoverFloorClaim / SecondClockRegimeClaim /
/// VacuumBlockReductionClaim. C# twin of simulations/birth_canal_junction_nature.py + star_frozen_seam.py.
/// Spec: docs/superpowers/specs/2026-06-18-trichotomy-witness-design.md.</summary>
public sealed class TrichotomyWitness : IInspectable
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    // Tolerances (reuse the shared birth-canal tolerance later; do not invent literals).
    private const double ImTol = 1e-6;          // matches StarFrozenSeamWitness frozen/oscillating split
    // The commutant ceiling 2γ·CommutantDarkest is the high-Q ASYMPTOTE; at a finite carbon Q the frozen
    // star (p,p) rate approaches it from below (relErr 1.6% at Q=8, 0.25% at Q=20, → 0 — gate-measured),
    // so the match is RELATIVE to the rate, |rate − ceilingRate| / ceilingRate ≤ tol. 3e-2 catches the
    // star's finite-Q approach at Q=8 yet rejects the ring interior, whose (p,p) survivor sits nowhere near
    // its commutant ceiling (relErr 64–85% at N=5 — a level crossing, not a commutant ceiling). The window
    // 1.6% < tol < 64% is wide; 3e-2 sits well inside it.
    private const double CommutantRelTol = 3e-2; // RELATIVE match Rate ≈ 2γ·CommutantDarkest (high-Q asymptote)

    // TWO physical sweeps, two reads (the convention defect of a single Classify, resolved 2026-06-18):
    //  - Route: the CARBON un-freeze trichotomy (Q = J/γ, uniform γ = 1/Q). Frozen (p,p) interior below
    //    Q*, oscillating (0,1) band edge above. Read by ClassifyUnfreeze.
    //  - SeamKind: the ABSOLUTE Δn-seam (fixed γ-profile, the global slowest's drift). Read by ClassifySeam.
    public enum Route { UnfreezingSeEp, FrozenLevelCrossing, FrozenCommutant }
    public enum SeamKind { Sterile, OddDrift, Junction }

    public readonly record struct SurvivorReading(
        int PCol, int PRow, int Dn, Route Route,
        double ImMax, double NXy, double Rate, double Rigidity, double Deviation);

    public readonly record struct SeamReading(int DnLo, int DnHi, double Deviation, SeamKind Kind);

    public TrichotomyWitness(int n = 6, double q = 1.5)
    {
        if (n < 3 || n > 8) throw new ArgumentOutOfRangeException(nameof(n), n, "N in 3..8");
        if (q <= 0) throw new ArgumentOutOfRangeException(nameof(q), q, "Q > 0");
        N = n; Q = q;
    }

    public int N { get; }
    public double Q { get; }

    private const double RateTieTol = 1e-9; // an F1-conjugate (p,p)<->(N-p,N-p) rate-tie is exact in physics;
                                            // independent eigensolves split it by ~ULPs. A new block must beat
                                            // the incumbent by more than this to win, so the canonical (smallest-p)
                                            // representative of a tied floor orbit is reported, not a noise pick.

    /// <summary>The survivor sector + its rate, scanning the candidate low-light blocks
    /// ((p,p) for p=1..N-1 plus the (0,1) band edge) with <see cref="SectorReductionWitness.SectorSlowest"/>
    /// on the absolute convention (H = Q·H_unit, per-site γ-profile). One scale for Rate, the route gates,
    /// and Deviation (a later task). The reported (PCol,PRow) of a Δn=1 survivor is the canonical (0,1)
    /// representative of the band-edge floor orbit; for a (p,p) survivor it is the smallest-p representative
    /// of the F1-conjugate ((p,p)~(N-p,N-p)) rate-tied orbit. The route keys on Δn + rate, not the literal
    /// sector. Ties (within <see cref="RateTieTol"/>) keep the incumbent so the canonical representative is
    /// stable against eigensolve noise.</summary>
    public static (int PCol, int PRow, double Rate) SurvivorSector(
        TopologyKind topo, int n, double q, IReadOnlyList<double> gammaProfile)
    {
        int bestPc = 0, bestPr = 1;
        double best = SectorReductionWitness.SectorSlowest(n, q, gammaProfile, 0, 1, topo);
        for (int p = 1; p <= n - 1; p++)
        {
            double r = SectorReductionWitness.SectorSlowest(n, q, gammaProfile, p, p, topo);
            if (r < best - RateTieTol) { best = r; bestPc = p; bestPr = p; }
        }
        return (bestPc, bestPr, best);
    }

    private const double KernelTol = 1e-7;

    /// <summary>H = Q·H_unit = ChainSystem(N, J=2Q, XY).BuildHamiltonian() (off-diag Q): the SAME absolute
    /// convention SectorReductionWitness.QHUnit uses, so the block here is the survivor block whose rate
    /// <see cref="SurvivorSector"/> computed. Copied (not shared) to keep this witness self-contained.</summary>
    private static ComplexMatrix QHUnit(int n, double q, TopologyKind topology) =>
        new ChainSystem(n, 2.0 * q, 1.0, HamiltonianType.XY, topology).BuildHamiltonian();

    /// <summary>The flat indices of the (pCol, pRow) sector. Copied verbatim from
    /// SectorReductionWitness.SectorFlat (do NOT edit the shipped witness to expose it).</summary>
    private static int[] SectorFlat(int n, int pCol, int pRow)
    {
        var decomp = JointPopcountSectorBuilder.Build(n);
        var s = decomp.SectorRanges.First(r => r.PCol == pCol && r.PRow == pRow);
        var flat = new int[s.Size];
        for (int k = 0; k < s.Size; k++) flat[k] = decomp.Permutation[s.Offset + k];
        return flat;
    }

    /// <summary>The sterile/canal split: the rate-drift of the GLOBAL slowest mode between Q_lo and Q_hi.
    /// It SWITCHES sectors (that switch is the junction; a fixed sector cannot see it — R2). For chain N≤6 use
    /// the validated <see cref="PostEpFlowField.BirthCanalDeviation"/> (it computes the global slowest at its
    /// own Q_lo=1.5/Q_hi=1000 on the full 4^N Liouvillian); for the sector-reduced path (N>6, or ring/star,
    /// which PostEpFlowField's FlowTopology cannot represent) recompute the global min over the candidate
    /// sectors at each Q via <see cref="SurvivorSector"/> and difference. A sector-PINNED (0,1) deviation
    /// gives 0.304 for canal-N5 — WRONG; the correct global value is 0.085.</summary>
    public static double Deviation(TopologyKind topo, int n, IReadOnlyList<double> gammaProfile)
    {
        const double qLo = PostEpFlowField.BirthCanalProbeQLow;   // 1.5
        const double qHi = PostEpFlowField.BirthCanalProbeQHigh;  // 1000.0
        if (topo == TopologyKind.Chain && n <= 6)
        {
            // Same construction as BirthCanalSurfaceWitness.ReadPoint: the probe-Q grid, a 2-point
            // τ-grid (PostEpFlowField needs ≥2), the literal per-site profile, default FlowTopology.Chain.
            var flow = new PostEpFlowField(n,
                new[] { qLo, qHi }, new[] { 0.0, 1.0 }, gammaProfile);
            return flow.BirthCanalDeviation;
        }
        double GlobalSlowest(double q) => SurvivorSector(topo, n, q, gammaProfile).Rate;
        return GlobalSlowest(qHi) - GlobalSlowest(qLo);
    }

    // ===================== the two-read Classify ===================================================

    /// <summary>The carbon-convention coordinates of <see cref="IncompletenessSurvivorWitness.Survivor"/>:
    /// off-diag hopping Qh=0.5 (the carbon J), uniform γ = 1/Q. <see cref="SurvivorSector"/> and
    /// <see cref="SectorSlowest"/> are ABSOLUTE (H = 2q·H_unit, off-diag q), so to read the SAME block
    /// Survivor read, pass q'=Qh=0.5 and the uniform 1/Q profile. The carbon-mapping gate
    /// (CarbonBlock_RateMatchesSurvivor) pins this: CarbonSlowestRate == Survivor.Gap to 9 digits.</summary>
    private const double CarbonQh = 0.5; // = IncompletenessSurvivorWitness.Qh (its private carbon J builder)

    private static IReadOnlyList<double> CarbonProfile(int n, double Q) =>
        Enumerable.Repeat(1.0 / Q, n).ToList();

    /// <summary>The slowest non-kernel rate of the carbon (pc,pr) block — the SAME block, built the SAME
    /// way, as <see cref="IncompletenessSurvivorWitness.Survivor"/> (off-diag Qh=0.5, uniform γ=1/Q). Used
    /// by the carbon-mapping gate to certify the |Im|/rigidity block is the canonical Survivor block.</summary>
    public static double CarbonSlowestRate(TopologyKind topo, int n, double Q, int pc, int pr) =>
        SectorReductionWitness.SectorSlowest(n, CarbonQh, CarbonProfile(n, Q), pc, pr, topo);

    /// <summary>|Im| / phase-rigidity of the SINGLE slowest non-kernel mode of the carbon (pc,pr) survivor
    /// block (q'=Qh=0.5, uniform γ=1/Q — bit-identical to <see cref="IncompletenessSurvivorWitness.Survivor"/>'s
    /// block; pinned by <see cref="CarbonSlowestRate"/> == Survivor.Gap). The un-freeze read
    /// is about the SURVIVOR itself: the single slowest mode (not a 2-mode pair, which would catch a band-edge
    /// {0,2} EP). Below Q* the chain (1,1) survivor is the
    /// overdamped-REAL SE-EP (|Im|≈0, frozen) while a far faster pair in the same block oscillates — taking
    /// the pair would mislabel the frozen survivor. This matches the Python verifier's slowest() (the single
    /// slowest mode's |Im|), the discriminator birth_canal_junction_nature.py reads.</summary>
    public static (double ImMax, double Rigidity) CarbonImAndRigidity(
        TopologyKind topo, int n, double Q, int pc, int pr)
    {
        var h = QHUnit(n, CarbonQh, topo);
        int[] flat = SectorFlat(n, pc, pr);
        var block = PerBlockLiouvillianBuilder.BuildBlockZ(h, CarbonProfile(n, Q), flat);
        var slow = PhaseRigidity.Compute(block)
            .Where(m => m.Lambda.Magnitude > KernelTol)   // drop the kernel (steady states)
            .OrderByDescending(m => m.Lambda.Real)         // slowest = largest (least negative) Re
            .FirstOrDefault();
        if (slow.Right is null) return (0.0, 1.0);
        return (Math.Abs(slow.Lambda.Imaginary), slow.Rigidity);
    }

    /// <summary>The carbon-convention un-freeze read (Q = J/γ, uniform γ = 1/Q): frozen (p,p) interior below
    /// Q*, oscillating (0,1) band edge above. Drives RouteSweep + ThresholdLadder. The survivor sector is the
    /// canonical IncompletenessSurvivorWitness.Survivor (J=1, γ=1/Q, validated bit-for-bit vs the full 4^N);
    /// |Im|/rigidity are read off the SAME carbon block.</summary>
    public static SurvivorReading ClassifyUnfreeze(TopologyKind topo, int n, double Q)
    {
        var (rate, pc, pr, nXy) = IncompletenessSurvivorWitness.Survivor(n, Q, topo); // (Gap, PCol, PRow, NXy)
        int dn = Math.Abs(pc - pr);
        var (imMax, rigidity) = CarbonImAndRigidity(topo, n, Q, pc, pr);
        Route route;
        if (dn == 1) route = Route.UnfreezingSeEp;       // the (0,1) band edge above Q* (oscillating)
        else
        {
            // Dn==0 frozen interior — the route is the topology's MECHANISM (not the instantaneous rate):
            //   star  → FrozenCommutant: the flat-band hub's survivor IS the [H,A]=0 commutant coherence at
            //           EVERY Q (that is why it is frozen); its rate only approaches the ceiling g2=4/(N−1)
            //           asymptotically, so the mechanism is commutant throughout, not just near the ceiling.
            //   chain → UnfreezingSeEp: the dispersive-band SE-EP (frozen below Q*, un-freezes above).
            //   ring  → FrozenLevelCrossing: the wrap-bond (2,2) level crossing — EXCEPT the N=4 ring (2,2),
            //           which sits ON the K_4/ring-4 commutant ceiling (g2=1) and routes FrozenCommutant,
            //           caught by the relative rate-match to the high-Q commutant asymptote.
            double? commutant = StructuralCeilingWitness.CommutantDarkest(TopoString(topo), n, pc, pr);
            double gamma = 1.0 / Q;
            double relErr = commutant.HasValue
                ? Math.Abs(rate - 2.0 * gamma * commutant.Value) / (2.0 * gamma * commutant.Value)
                : double.PositiveInfinity;
            route = topo switch
            {
                TopologyKind.Chain => Route.UnfreezingSeEp,
                TopologyKind.Star => Route.FrozenCommutant,
                TopologyKind.Ring => relErr <= CommutantRelTol ? Route.FrozenCommutant : Route.FrozenLevelCrossing,
                _ => Route.FrozenLevelCrossing,
            };
        }
        // Deviation (last field) is the SEAM read's quantity; the carbon un-freeze read has none → 0.0.
        return new SurvivorReading(pc, pr, dn, route, imMax, nXy, rate, rigidity, 0.0);
    }

    /// <summary>The absolute-convention sterile/odd-drift/junction read (fixed γ-profile, the global slowest's
    /// drift between Q_lo and Q_hi). Drives DeltaNSeam. Reproduces birth_canal_junction_nature.py.</summary>
    public static SeamReading ClassifySeam(TopologyKind topo, int n, IReadOnlyList<double> gammaProfile)
    {
        const double qLo = 1.5, qHi = 1000.0;
        int dnLo = AbsDn(SurvivorSector(topo, n, qLo, gammaProfile));
        int dnHi = AbsDn(SurvivorSector(topo, n, qHi, gammaProfile));
        double dev = Deviation(topo, n, gammaProfile);
        SeamKind kind = Math.Abs(dev) < PostEpFlowField.BirthCanalTolerance ? SeamKind.Sterile
                      : dnLo != dnHi ? SeamKind.Junction      // the interior overtakes at low Q (Δn flips 0→1)
                      : SeamKind.OddDrift;                      // the (0,1) edge drifts, no Δn switch
        return new SeamReading(dnLo, dnHi, dev, kind);
    }

    private static int AbsDn((int PCol, int PRow, double Rate) s) => Math.Abs(s.PCol - s.PRow);

    /// <summary>The lowercase topology name <see cref="StructuralCeilingWitness.CommutantDarkest"/> keys on.</summary>
    private static string TopoString(TopologyKind topo) => topo switch
    {
        TopologyKind.Chain => "chain",
        TopologyKind.Ring => "ring",
        TopologyKind.Star => "star",
        _ => throw new ArgumentOutOfRangeException(nameof(topo), topo, "unsupported topology"),
    };

    public string DisplayName => $"TrichotomyWitness (N={N}, Q={Q.ToString("0.###", Inv)})";
    public string Summary => "the chain/ring/star survivor trichotomy as one sweep; see --root starseam / horizon / surface";

    // ===================== the IInspectable tree (the four slices) ==================================

    private static readonly double[] QGrid = { 1.0, 1.5, 2.0, 3.0, 6.0, 12.0, 25.0, 50.0 };
    private static readonly TopologyKind[] Trio = { TopologyKind.Chain, TopologyKind.Ring, TopologyKind.Star };
    private static IReadOnlyList<double> Uniform(int n, double gamma) => Enumerable.Repeat(gamma, n).ToList();

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return TheRouteSweep();
            yield return TheThresholdLadder();
            yield return TheDeltaNSeam();
            yield return TheVocabulary();
        }
    }

    /// <summary>1. RouteSweep — the CARBON un-freeze read across Q, per topology (the trichotomy visible in one
    /// place). Chain flips frozen→oscillating at Q*(N), ring at Q_h, star never (N≥5). Drives
    /// <see cref="ClassifyUnfreeze"/>. The per-Q <see cref="ClassifyUnfreeze"/> work is LAZY: each topology
    /// node's rows are a <c>yield</c>-based iterator (<see cref="RouteRows"/>), so the eigensolves run only
    /// when that topology node is enumerated by the renderer — a shallow render fires none of them.</summary>
    private InspectableNode TheRouteSweep()
    {
        var topos = new List<IInspectable>();
        foreach (var topo in Trio)
        {
            var t = topo; // capture per-iteration for the lazy iterator
            topos.Add(new InspectableNode(TopoString(t),
                summary: $"survivor + route across Q (carbon, N={N})", children: RouteRows(t)));
        }
        return new InspectableNode("the route sweep (carbon): survivor + freeze-route across Q",
            summary: "chain flips frozen→oscillating at Q*(N), ring at Q_h, star never (N≥5)", children: topos);
    }

    /// <summary>The per-Q rows of one topology in <see cref="TheRouteSweep"/>, computed lazily: the
    /// <see cref="ClassifyUnfreeze"/> eigensolve for each Q runs only as the renderer pulls the next row,
    /// so an un-enumerated topology node costs nothing.</summary>
    private IEnumerable<IInspectable> RouteRows(TopologyKind topo)
    {
        foreach (double q in QGrid)
        {
            var r = ClassifyUnfreeze(topo, N, q);
            string frozen = r.ImMax < ImTol ? "frozen" : "oscillating";
            yield return new InspectableNode($"Q={q.ToString("0.##", Inv)}",
                summary: $"({r.PCol},{r.PRow}) Δn={r.Dn} | {r.Route} | {frozen} " +
                         $"|Im|={r.ImMax.ToString("0.##e+0", Inv)} | r={r.Rigidity.ToString("0.###", Inv)} " +
                         $"| ⟨n_XY⟩={r.NXy.ToString("0.###", Inv)}");
        }
    }

    /// <summary>2. ThresholdLadder — the three thresholds over N: chain Q*(N), ring Q_h≈0.29N, star g2=4/(N−1).
    /// N=4 star is the lone (2,2)/K₄ outlier (g2=4/3>1, un-freezes). The per-N rows are LAZY
    /// (<see cref="LadderRows"/> is a <c>yield</c> iterator): the <see cref="IncompletenessSurvivorWitness.HandoverQ"/>
    /// and <see cref="StructuralCeilingWitness.CommutantDarkest"/> compute for each N runs only when this node's
    /// children are enumerated, so a shallow render that stops at this slice header pays nothing.</summary>
    private InspectableNode TheThresholdLadder() => new(
        "the threshold ladder over N",
        summary: "chain Q*(N) / ring Q_h≈0.29N / star g2=4/(N−1); N=4 star outlier", children: LadderRows());

    /// <summary>The per-N rows of <see cref="TheThresholdLadder"/>, computed lazily: each N's HandoverQ +
    /// CommutantDarkest runs only as the renderer pulls the next row.</summary>
    private IEnumerable<IInspectable> LadderRows()
    {
        for (int n = 4; n <= 8; n++)
        {
            double chainQStar = IncompletenessSurvivorWitness.HandoverQ(n, TopologyKind.Chain);
            double ringQh = IncompletenessSurvivorWitness.HandoverQ(n, TopologyKind.Ring);
            double starCeil = StructuralCeilingWitness.CommutantDarkest("star", n, 1, 1) ?? double.NaN;
            string starVerdict = starCeil <= 1.0 ? "frozen (g2≤1)" : "UN-FREEZES (g2>1; (0,1) edge at √(N−1)·J)";
            yield return new InspectableNode($"N={n}",
                summary: $"chain Q*={chainQStar.ToString("0.###", Inv)} | ring Q_h={ringQh.ToString("0.###", Inv)} | " +
                         $"star g2=4/(N−1)={starCeil.ToString("0.###", Inv)} → {starVerdict}");
        }
    }

    /// <summary>3. DeltaNSeam — the ABSOLUTE seam read (uniform/canal/deep-edge), reproducing
    /// simulations/birth_canal_junction_nature.py. NOTE: for chain N=6 this calls the slow PostEpFlowField
    /// path (full 4^N, ~2 min) — it is LAZY (<see cref="SeamRows"/> is a <c>yield</c> iterator), so the
    /// <see cref="ClassifySeam"/> work runs ONLY when this node's children are enumerated by the renderer,
    /// not when this method is called. A shallow render (this slice header only) fires none of it. Drives
    /// <see cref="ClassifySeam"/>.</summary>
    private InspectableNode TheDeltaNSeam() => new(
        "the Δn seam (absolute): sterile / odd-drift / junction",
        summary: "junction ⟹ birth canal, not conversely (birth_canal_junction_nature.py); slow node (~2 min)",
        children: SeamRows());

    /// <summary>The seam-case rows of <see cref="TheDeltaNSeam"/>, computed lazily: the slow per-case
    /// <see cref="ClassifySeam"/> (full 4^N PostEpFlowField for chain N=6, ~2 min total) runs only as the
    /// renderer pulls the next row, so the parent <see cref="Children"/> getter can list this slice header
    /// without paying for the seam compute.</summary>
    private IEnumerable<IInspectable> SeamRows()
    {
        // The EXACT profiles from simulations/birth_canal_junction_nature.py.
        var cases = new (string Name, int N, double[] Profile)[]
        {
            ("uniform N=5",   5, Uniform(5, 0.5).ToArray()),
            ("canal N=5",     5, new[] { 0.25, 1.5, 1.5, 1.5, 0.25 }),
            ("deep-edge N=6", 6, new[] { 0.25, 1.375, 1.375, 1.375, 1.375, 0.25 }),
        };
        foreach (var (name, n, profile) in cases)
        {
            var s = ClassifySeam(TopologyKind.Chain, n, profile);
            yield return new InspectableNode(name,
                summary: $"{s.Kind} | Deviation={s.Deviation.ToString("0.####", Inv)} | Δn {s.DnLo}→{s.DnHi}");
        }
    }

    /// <summary>4. Vocabulary — prose + cross-links: the two reads and where they live.</summary>
    private static InspectableNode TheVocabulary() => new(
        "the vocabulary: rate_slow = min over Δn-sorted joint-popcount sectors, two reads",
        summary: "CARBON un-freeze read (Q=J/γ): UnfreezingSeEp (chain) / FrozenLevelCrossing (ring) / " +
                 "FrozenCommutant (star). ABSOLUTE seam read (fixed γ, vary profile): Sterile / OddDrift / Junction. " +
                 "See --root horizon / starseam / ceiling / surface / survivor / secondclock; docs " +
                 "THE_STAR_FROZEN_SEAM.md, STERILE_BIRTHCANAL_AND_THE_JUNCTION.md.");

    public InspectablePayload Payload => InspectablePayload.Empty;
}
