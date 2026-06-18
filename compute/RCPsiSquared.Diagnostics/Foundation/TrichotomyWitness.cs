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
    private const double RigTol = 1e-2;          // matches the horizon EP detector
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

    /// <summary>max|Im| of the slowest non-kernel pair (frozen ≈ 0 vs oscillating) and the min Petermann
    /// phase rigidity (→ 0 at an EP) of the (pCol,pRow) survivor block, on the SAME absolute convention as
    /// <see cref="SurvivorSector"/> (H = 2q·H_unit, the same per-site γ-profile, the same SectorFlat indices),
    /// so the |Im|/rigidity belongs to the block whose rate SurvivorSector reported. Builds the block with the
    /// SectorReductionWitness recipe (no 4^N). The star (1,1) commutant is frozen at every Q (|Im|≈0); the
    /// chain (0,1) band edge oscillates above Q*(N).</summary>
    public static (double ImMax, double Rigidity) ImAndRigidity(
        TopologyKind topo, int n, double q, IReadOnlyList<double> gammaProfile, int pCol, int pRow)
    {
        var h = QHUnit(n, q, topo);
        int[] flat = SectorFlat(n, pCol, pRow);
        var block = PerBlockLiouvillianBuilder.BuildBlockZ(h, gammaProfile, flat);
        // Mirror SectorReductionWitness.cs:126: PhaseRigidity.Compute(block) takes the BuildBlockZ result
        // (a Matrix<Complex> = ComplexMatrix) verbatim; no conversion/adapter.
        var modes = PhaseRigidity.Compute(block)
            .Where(m => m.Lambda.Magnitude > KernelTol)   // drop the kernel (steady states)
            .OrderByDescending(m => m.Lambda.Real)         // slowest = largest (least negative) Re
            .ToList();
        if (modes.Count == 0) return (0.0, 1.0);
        double imMax = modes.Take(2).Max(m => Math.Abs(m.Lambda.Imaginary));
        double rigidity = modes.Take(2).Min(m => m.Rigidity);
        return (imMax, rigidity);
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

    // ===================== TASK 5: the two-read Classify ===========================================

    /// <summary>The carbon-convention coordinates of <see cref="IncompletenessSurvivorWitness.Survivor"/>:
    /// off-diag hopping Qh=0.5 (the carbon J), uniform γ = 1/Q. <see cref="ImAndRigidity"/> and
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
    /// block; pinned by <see cref="CarbonSlowestRate"/> == Survivor.Gap). Unlike the absolute
    /// <see cref="ImAndRigidity"/> (which reads a 2-mode pair for the band-edge {0,2} EP), the un-freeze read
    /// is about the SURVIVOR itself: the single slowest mode. Below Q* the chain (1,1) survivor is the
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
            double? commutant = StructuralCeilingWitness.CommutantDarkest(TopoString(topo), n, pc, pr);
            double gamma = 1.0 / Q;
            double ceilingRate = commutant.HasValue ? 2.0 * gamma * commutant.Value : double.PositiveInfinity;
            // RELATIVE match to the high-Q commutant asymptote: the star (p,p) survivor sits ON the ceiling
            // (relErr → 0); the ring's interior level crossing sits far below it (relErr O(1)).
            double relErr = double.IsPositiveInfinity(ceilingRate) ? double.PositiveInfinity
                          : Math.Abs(rate - ceilingRate) / ceilingRate;
            if (relErr <= CommutantRelTol)
                route = Route.FrozenCommutant;            // star N≥5 (flat-band non-dispersive hub)
            else
                route = topo == TopologyKind.Chain ? Route.UnfreezingSeEp : Route.FrozenLevelCrossing; // R1 topo-key
        }
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
    public IEnumerable<IInspectable> Children { get { yield break; } } // filled in a later task
    public InspectablePayload Payload => InspectablePayload.Empty;
}
