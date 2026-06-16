using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>The second clock's regime map (Tier 1 candidate, "the stitch"): the {0,2}/half-filling
/// coherence (the SECOND clock, the one that competes with the band-edge survivor for "slowest mode")
/// is ONE phenomenon across the arcs, and its fate is selected by the single-particle band through TWO
/// knobs. This is the seam that joins three previously separate readings into one node.
///
/// <para><b>The two knobs.</b> The first clock, the band-edge coherence |vac⟩⟨ψ_1| (n_XY = 1), always
/// sits on the −2γ floor (g2 = 1) for every graph at every Q (the Absorption Theorem). A SECOND mode,
/// the {0,2}/half-filling coherence, competes with it; what it does is read off the band:
/// <list type="bullet">
/// <item><b>Knob 1, DEGENERACY m (the high-Q fate).</b> A dark-enough degenerate single-particle manifold
/// pulls the second clock BELOW the −2γ floor: a STRUCTURAL CEILING g2 &lt; 1 at all Q (the
/// <see cref="StructuralCeilingClaim"/>: complete K_N = 4/N, star N≥6 = 4/(N−1); both = 4/(m+1) for the
/// symmetric manifold, &lt; 1 iff m ≥ 4). Otherwise the second clock reaches the floor.</item>
/// <item><b>Knob 2, DISPERSION (the low-Q character once it reaches the floor).</b> On a DISPERSIVE band
/// (a real cosine spectrum: the chain m=1, the disordered chain, the ring, star N=4) the two clocks
/// COALESCE — a square-root EP at a finite Q*(N), the COHERENCE HORIZON (<see cref="CoherenceHorizonClaim"/>):
/// below Q* a real overdamped mode is the gap, above it the band edge protects, g2 → 1 SHARPLY. On a
/// flat / marginal band (star N=5, m=3, 4/(m+1)=1) there is only ASYMPTOTIC protection (g2 = 1 − c/Q², no
/// sharp horizon — the "star no-horizon" reading).</item>
/// </list></para>
///
/// <para><b>The stitch.</b> CoherenceHorizon (the EP regime), StructuralCeiling (the CEILING regime,
/// 4/(m+1) the bridge), and star-no-horizon (the GRADUAL regime) are not three findings: they are one
/// mode whose regime = map(degeneracy, dispersion). The N=4 anomalies (ring-4 GRADUAL, complete-4 CEILING)
/// are the (2,2) half-filling sector, the same N=3/N=4 specials the rest of the project keeps meeting.</para>
///
/// <para>Tier1Candidate: gate-verified as a 2D regime map (15/15 across chain / disordered-chain / ring /
/// star / complete × N=4,5,6, simulations/second_clock_regime_axis.py) and live (the N=4 frame is a
/// self-validating full-Liouvillian gate, inspect --root secondclock). It is a candidate, not derived,
/// because it is capped by its weaker parent <see cref="CoherenceHorizonClaim"/> (Tier1Candidate, its own
/// half-filling V-Effect seam still open); the other parent <see cref="StructuralCeilingClaim"/> is
/// Tier1Derived. The two parents ARE the two regimes the map joins.</para>
///
/// <para>Live witness: <c>inspect --root secondclock</c>
/// (<c>compute/RCPsiSquared.Diagnostics/Foundation/SecondClockRegimeWitness.cs</c>).</para></summary>
public sealed class SecondClockRegimeClaim : Claim
{
    /// <summary>Parent: the EP regime. On a dispersive band the second clock coalesces with the band edge
    /// at the finite coherence horizon Q*(N) (a √-EP). The weaker parent (Tier1Candidate); caps this child.</summary>
    public CoherenceHorizonClaim Horizon { get; }

    /// <summary>Parent: the CEILING regime. On a degenerate band the second clock is the darkest commutant
    /// coherence g2 = 4/(m+1) below the floor (Tier1Derived). The 4/(m+1) form is the bridge between the
    /// regimes (it equals 1 exactly at the marginal m=3 star, where the regime turns GRADUAL).</summary>
    public StructuralCeilingClaim Ceiling { get; }

    /// <summary>The symmetric-manifold high-Q ceiling 4/(m+1) for a degeneracy-m single-particle band: the
    /// bridge formula. Complete K_N (m=N−1) gives 4/N; star (m=N−2 for the leaf manifold) gives 4/(N−1); it
    /// drops below the floor (a structural ceiling) iff m ≥ 4, equals 1 (marginal, GRADUAL) at m = 3.</summary>
    public static double SymmetricManifoldCeiling(int m) => 4.0 / (m + 1);

    public SecondClockRegimeClaim(CoherenceHorizonClaim horizon, StructuralCeilingClaim ceiling)
        : base("The second clock (the {0,2}/half-filling coherence) is ONE mode whose regime is selected by " +
               "the single-particle band through two knobs. Knob 1 (DEGENERACY m): a dark-enough degenerate " +
               "manifold pulls it below the −2γ floor = a structural ceiling g2=4/(m+1)<1 (m≥4: complete K_N=4/N, " +
               "star N≥6=4/(N−1)); else it reaches the floor. Knob 2 (DISPERSION): on a dispersive band the two " +
               "clocks coalesce at a √-EP coherence horizon Q*(N) (sharp protection); on a flat/marginal band " +
               "(m=3, 4/(m+1)=1) only asymptotic protection (1−c/Q², no sharp horizon). The stitch: the coherence " +
               "horizon (EP), the structural ceiling (CEILING, 4/(m+1) the bridge), and star-no-horizon (GRADUAL) " +
               "are one mode, regime = map(degeneracy, dispersion). N=4 anomalies (ring-4 GRADUAL, complete-4 " +
               "CEILING) are the (2,2) half-filling sector. Gate-verified 15/15 over chain/disordered/ring/star/" +
               "complete × N=4,5,6.",
               Tier.Tier1Candidate,
               "simulations/second_clock_regime_axis.py + " +
               "compute/RCPsiSquared.Diagnostics/Foundation/SecondClockRegimeWitness.cs (SecondClockRegimeWitness, inspect --root secondclock)")
    {
        Horizon = horizon ?? throw new ArgumentNullException(nameof(horizon));
        Ceiling = ceiling ?? throw new ArgumentNullException(nameof(ceiling));
    }

    public override string DisplayName =>
        "The second clock's regime map: one mode, regime = map(band degeneracy, dispersion)";

    public override string Summary =>
        $"the {{0,2}}/half-filling coherence is ONE mode; knob 1 (degeneracy m) → high-Q ceiling 4/(m+1) " +
        $"(below the floor iff m≥4); knob 2 (dispersion) → low-Q √-EP horizon (dispersive) vs gradual (flat). " +
        $"Stitches CoherenceHorizon (EP) + StructuralCeiling (CEILING) + star-no-horizon (GRADUAL) " +
        $"({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("knob 1: degeneracy m → the high-Q fate (4/(m+1))",
                summary: "a dark-enough degenerate single-particle manifold pulls the second clock below the −2γ " +
                         "floor: g2 = 4/(m+1), a structural ceiling iff m ≥ 4 (complete K_N: m=N−1, g2=4/N; star: the " +
                         "leaf manifold, g2=4/(N−1)). At m=3 it equals 1 exactly (marginal); m≤2 reaches the floor.");
            yield return new InspectableNode("knob 2: dispersion → the low-Q character (EP vs gradual)",
                summary: "once the second clock reaches the floor, the band's dispersion sets HOW: a dispersive band " +
                         "(real cosine spectrum: chain, disordered chain, ring, star N=4) makes the two clocks coalesce " +
                         "at a sharp √-EP coherence horizon Q*(N); a flat/marginal band (star N=5) gives only asymptotic " +
                         "1−c/Q² protection (no sharp horizon).");
            yield return new InspectableNode("the N=4 anomalies = the (2,2) half-filling sector",
                summary: "ring-4 (GRADUAL) and complete-4 (CEILING, 2−2/√3) live in the (2,2) two-excitation sector — " +
                         "the same N=3/N=4 specials the project keeps meeting (the n3_special_cases arc).");
            yield return new InspectableNode("the stitch (why this is one node)",
                summary: "CoherenceHorizonClaim (EP), StructuralCeilingClaim (CEILING, 4/(m+1) the bridge), and " +
                         "star-no-horizon (GRADUAL) are not three findings but one mode whose regime = " +
                         "map(degeneracy, dispersion). Everything is one object, seen at max zoom on different facets.");
            yield return Horizon;   // typed parent edge (the EP regime, Tier1Candidate, caps this child)
            yield return Ceiling;   // typed parent edge (the CEILING regime, Tier1Derived)
        }
    }

    public static SecondClockRegimeClaim Build() =>
        new SecondClockRegimeClaim(CoherenceHorizonClaim.Build(), StructuralCeilingClaim.Build());

    public static SecondClockRegimeClaim Shared { get; } = Build();
}
