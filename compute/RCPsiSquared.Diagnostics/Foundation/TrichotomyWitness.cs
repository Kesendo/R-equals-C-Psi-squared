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
    private const double CommutantRelTol = 1e-3; // relative match Rate ≈ 2γ·CommutantDarkest

    public enum Route { UnfreezingSeEp, FrozenLevelCrossing, FrozenCommutant, OddDriftBandEdge, SterileBandEdge }
    // Five labels, four physical routes: the (0,1) band edge carries two (Sterile/OddDrift, split by |Deviation|).

    public readonly record struct SurvivorReading(
        int PCol, int PRow, int Dn, Route Route,
        double ImMax, double NXy, double Rate, double Rigidity, double Deviation);

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

    public string DisplayName => $"TrichotomyWitness (N={N}, Q={Q.ToString("0.###", Inv)})";
    public string Summary => "the chain/ring/star survivor trichotomy as one sweep; see --root starseam / horizon / surface";
    public IEnumerable<IInspectable> Children { get { yield break; } } // filled in a later task
    public InspectablePayload Payload => InspectablePayload.Empty;
}
