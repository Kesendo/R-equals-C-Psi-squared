using System.Globalization;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Inspection;

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

    public string DisplayName => $"TrichotomyWitness (N={N}, Q={Q.ToString("0.###", Inv)})";
    public string Summary => "the chain/ring/star survivor trichotomy as one sweep; see --root starseam / horizon / surface";
    public IEnumerable<IInspectable> Children { get { yield break; } } // filled in a later task
    public InspectablePayload Payload => InspectablePayload.Empty;
}
