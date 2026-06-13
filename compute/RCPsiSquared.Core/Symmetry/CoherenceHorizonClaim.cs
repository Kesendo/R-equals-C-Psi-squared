using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>The coherence horizon Q*(N) (Tier 1 candidate). Q*(N) is the threshold where, sweeping
/// Q = J/γ downward, the slowest non-zero Liouvillian mode stops oscillating (the coherence hand
/// freezes): Q*(2)=1, Q*(3)=√2, Q*(4)≈1.8785, Q*(5)≈2.3722. The witness verifies these match, bit
/// for bit, the carbon Frost-Hückel coherent↔incoherent threshold (√2 / 1.879 / 2.372 at N=3/4/5)
/// under the label swap J ↔ |β|: the XY chain's coherence horizon IS the Frost-Hückel threshold of
/// the same polyene (the cross-substrate identity). N=2 (Q*=1) is the exceptional point itself, the
/// base rung the carbon polyene layer (N≥3) cannot reach; the quantum side supplies it.
///
/// <para>What is VERIFIED N=2..5: the ladder Q*(N) and its carbon identity, the band-edge
/// coincidence Q*(N) = 2cos(π/(N+1)) at N=2,3 only, and (corrected 2026-06-13, phase-rigidity
/// review) the freezing mechanism. The mode that COALESCES at Q*(N) is the {0,2}-coherence
/// (population / antisymmetric-coherence block, n_diff histogram {0: ½, 2: ½}, ⟨n_diff⟩ = 1) at
/// ALL N=2..5, a genuine square-root EP (phase rigidity r → 0, Im ∝ √(Q−Q*)); there is NO sector
/// bifurcation at N=4. The band edge 2cos(π/(N+1)) is the co-located SURVIVOR (the |vac⟩⟨ψ_k|
/// coherence hand, Uhr 1, γ-protected, r ≈ 1) sharing the gap Re = −2γ only because the Absorption
/// Theorem pins both (both ⟨n_diff⟩ = 1). So Q*(N) is at once a {0,2}-coherence EP (the erasure
/// point, which climbs the ladder) and a band-edge crossing (Uhr 1 survives). What is OPEN: the
/// closed form = the discriminant-zero of the {0,2}-coherence block (which grows with N). Hence
/// Tier1Candidate: ladder, carbon identity and EP mechanism verified, closed form the open
/// remainder.</para>
///
/// <para>Live witness: <c>inspect --root horizon</c>
/// (<c>compute/RCPsiSquared.Diagnostics/Foundation/CoherenceHorizonWitness.cs</c>).</para>
///
/// <para>Typed parents: <see cref="ClockHandLadderClaim"/> (the horizon IS the clock's Q-floor made
/// exact: where the coherence hand finally freezes) and <see cref="F2bXyChainSpectrumPi2Inheritance"/>
/// (the band edge 2cos(π/(N+1)) the horizon coincides with at N=2,3).</para>
///
/// <para>Anchors: <c>docs/ANALYTICAL_FORMULAS.md</c> F2b corollary "Coherence horizon Q*(N)" +
/// <c>docs/carbon/FROST_CIRCLE_AS_THE_CLOCK_FACE.md</c> (the carbon-layer twin) +
/// <c>compute/RCPsiSquared.Diagnostics/Foundation/CoherenceHorizonWitness.cs</c> (the live lab,
/// <c>inspect --root horizon</c>).</para></summary>
public sealed class CoherenceHorizonClaim : Claim
{
    /// <summary>Parent: the two clocks. The coherence horizon is the clock's Q-floor made exact,
    /// the Q below which the coherence hand finally freezes.</summary>
    public ClockHandLadderClaim Horizon { get; }

    /// <summary>Parent: the F2b band edge 2J·cos(π/(N+1)) the horizon coincides with at N=2,3.</summary>
    public F2bXyChainSpectrumPi2Inheritance BandEdge { get; }

    public CoherenceHorizonClaim(
        ClockHandLadderClaim horizon,
        F2bXyChainSpectrumPi2Inheritance bandEdge)
        : base("The coherence horizon Q*(N) = 1 / √2 / 1.8785 / 2.3722 for N=2..5, verified equal to the carbon " +
               "Frost-Hückel coherent↔incoherent threshold (√2 / 1.879 / 2.372 at N=3/4/5) under the label swap " +
               "J ↔ |β| (the cross-substrate identity); N=2 (Q*=1) is the exceptional point, the base rung the " +
               "carbon polyene layer (N≥3) cannot reach. The mode that coalesces at Q*(N) is the {0,2}-coherence " +
               "(population/antisymmetric block) at ALL N=2..5, a genuine √-EP (phase rigidity r → 0); the band edge " +
               "2cos(π/(N+1)) is the co-located γ-protected survivor, not the freezer (no bifurcation at N=4). The " +
               "closed form of Q*(N) is OPEN: the discriminant-zero of the {0,2}-coherence block.",
               Tier.Tier1Candidate,
               "docs/ANALYTICAL_FORMULAS.md F2b corollary \"Coherence horizon Q*(N)\" + " +
               "docs/carbon/FROST_CIRCLE_AS_THE_CLOCK_FACE.md (the carbon-layer twin) + " +
               "compute/RCPsiSquared.Diagnostics/Foundation/CoherenceHorizonWitness.cs (CoherenceHorizonWitness, inspect --root horizon)")
    {
        Horizon = horizon ?? throw new ArgumentNullException(nameof(horizon));
        BandEdge = bandEdge ?? throw new ArgumentNullException(nameof(bandEdge));
    }

    public override string DisplayName => "The coherence horizon Q*(N): the carbon coherent↔incoherent threshold made exact";

    public override string Summary =>
        $"Q*(N) = 1/√2/1.8785/2.3722 (N=2..5) = the carbon Frost-Hückel coherent↔incoherent threshold under J ↔ |β|; " +
        $"N=2 (Q*=1) is the EP base; the {{0,2}}-coherence coalesces (√-EP) at every N, the band edge is the co-located " +
        $"survivor; closed form OPEN (the {{0,2}}-block discriminant) ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("the ladder = carbon (verified N=2..5)",
                summary: "Q*(N) = 1 / √2 / 1.8785 / 2.3722 for N=2..5, bisected on Symphony.Clock.Omega at J=1, equals " +
                         "the carbon Frost-Hückel coherent↔incoherent threshold (√2 / 1.879 / 2.372 at N=3/4/5) under the " +
                         "label swap J ↔ |β|: the XY chain's coherence horizon IS the Frost-Hückel threshold of the same " +
                         "polyene. Live: inspect --root horizon.");
            yield return new InspectableNode("the EP base (N=2)",
                summary: "Q*(2) = 1 is the exceptional point itself (γ=J), the base rung the carbon polyene layer (N≥3) " +
                         "cannot reach (a polyene needs ≥3 sites); the quantum side supplies it.");
            yield return new InspectableNode("the band-edge coincidence (N=2,3)",
                summary: "Q*(N) = 2cos(π/(N+1)) at N=2,3 ONLY: 1 = 2cos60°, √2 = 2cos45°. A low-N accident, departing at " +
                         "N≥4 (Q*(4)=1.8785 ≠ φ=1.618); this is why the √2 looked exact at N=3 and the rest awaited a clean form.");
            yield return new InspectableNode("the EP mechanism: the {0,2}-coherence coalesces at every N",
                summary: "corrected 2026-06-13 (phase-rigidity review, superseding an earlier \"sector bifurcation at N=4\" " +
                         "reading that was an argmax-Re / Im-tracking artifact). The mode that coalesces at Q*(N) is the " +
                         "{0,2}-coherence (population/antisymmetric block, n_diff histogram {0: ½, 2: ½}, ⟨n_diff⟩=1) at ALL " +
                         "N=2..5, a genuine square-root EP (phase rigidity r → 0; r at Q* = 0.000/0.015/0.026 at N=3/4/5, " +
                         "Im ∝ √(Q−Q*)). NO bifurcation at N=4. The band edge 2cos(π/(N+1)) (Im = φ at N=4, √3 at N=5) is the " +
                         "co-located SURVIVOR, the γ-protected |vac⟩⟨ψ_k| coherence hand (Uhr 1, r ≈ 1), sharing the gap " +
                         "Re = −2γ only because the Absorption Theorem pins both (both ⟨n_diff⟩=1). So Q*(N) is at once a " +
                         "{0,2}-coherence EP (the erasure point, climbs the ladder) and a band-edge crossing (Uhr 1 survives). " +
                         "What is OPEN: the closed form = the discriminant-zero of the {0,2}-coherence block (grows with N). " +
                         "Probe: simulations/_coherence_horizon_closed_form.py (phase-rigidity sweep to be added).");
            yield return Horizon;   // typed parent edge
            yield return BandEdge;  // typed parent edge
        }
    }

    public static CoherenceHorizonClaim Build()
    {
        // Reuse one ladder within the F2b band-edge subtree; the ClockHandLadder subtree builds its
        // own (standalone Build()/Shared only; the live registry shares both by type via b.Get).
        var ladder = new Pi2DyadicLadderClaim();
        var qubit = new QubitDimensionalAnchorClaim();

        // F2b band-edge parent: F2b(ladder, f65); f65(ladder, f66); f66(ladder, qubit). Same wiring
        // as ClockHandLadderClaim.Build(), reusing the one shared ladder.
        var f66 = new F66PoleModesPi2Inheritance(ladder, qubit);
        var f65 = new F65XxChainSpectrumPi2Inheritance(ladder, f66);
        var bandEdge = new F2bXyChainSpectrumPi2Inheritance(ladder, f65);

        // The two-clocks parent is a fully wired sibling claim; reuse its own Build().
        var horizon = ClockHandLadderClaim.Build();

        return new CoherenceHorizonClaim(horizon, bandEdge);
    }

    public static CoherenceHorizonClaim Shared { get; } = Build();
}
