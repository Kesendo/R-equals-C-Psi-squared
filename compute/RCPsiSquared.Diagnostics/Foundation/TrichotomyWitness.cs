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

    public string DisplayName => $"TrichotomyWitness (N={N}, Q={Q.ToString("0.###", Inv)})";
    public string Summary => "the chain/ring/star survivor trichotomy as one sweep; see --root starseam / horizon / surface";
    public IEnumerable<IInspectable> Children { get { yield break; } } // filled in a later task
    public InspectablePayload Payload => InspectablePayload.Empty;
}
