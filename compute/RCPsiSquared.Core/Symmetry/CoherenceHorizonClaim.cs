using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>The coherence horizon Q*(N) (Tier 1 candidate).
///
/// <para><b>Where this sits</b> (returning-reader compass, so the EP-character seams stop tripping us):
/// this is the <b>(1,1) SE-SE coherence block</b>'s real-axis <b>DEFECTIVE</b> √-EP — a single-particle
/// ε-pair driven together by γ (a Jordan block, at the rung midpoint Re = −2γ). It is a DIFFERENT pair
/// from the same block's <b>−4γ DIABOLIC</b> point (a sum-coincidence ε_a−ε_b = ε_c−ε_d, semisimple,
/// present already at γ=0). One free-fermion-kernel rule covers both: <i>defective</i> = a γ-driven
/// single-particle coalescence (Jordan, √-branch, at a rung midpoint); <i>diabolic</i> = a
/// frequency-difference coincidence (semisimple, number-theoretic, at γ=0). The artifact-free test is
/// <c>EpCharacterWitness</c> / <c>inspect --root epcharacter</c>. Convention: Q = J/γ in the carrier /
/// (J/2)(XX+YY) normalization (the F89 octic's q = Q/2; see <c>docs/GLOSSARY.md</c> "The coupling ratio
/// q and Q"). Siblings on the same EP-character seam: <c>F89Path3OcticEpClaim</c> (the (1,2) octic's
/// on-axis diabolic), <c>experiments/F89_MONODROMY_MIRROR.md</c> (the (1,2) octic's off-axis defective
/// branch points), <c>F86/LocalGlobalEpLink</c> (the N≥5 (1,2) block, non-normal-simple on axis).</para>
///
/// Q*(N) is the threshold where, sweeping
/// Q = J/γ downward, the slowest non-zero Liouvillian mode stops oscillating (the coherence hand
/// freezes): Q*(2)=1, Q*(3)=√2, Q*(4)≈1.8787, Q*(5)≈2.3737. The witness verifies these match, bit
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
/// point, which climbs the ladder) and a band-edge crossing (Uhr 1 survives). The closed form
/// (resolved 2026-06-13, single-excitation reduction): Q*(N) reduces 4^N → N² (the coalescing mode
/// is single-excitation, so Q*(N) is the EP of the Haken-Strobl Liouvillian); at N=2,3 the pair are
/// roots of λ²+4γλ+c·J²=0 (c=4, 2), giving Q* = 2/√c = 1, √2 exactly; at N≥4 the pair is collectively
/// dressed (no clean 2×2), the exact condition transcendental, Q* a diffusive ~linear growth with asymptotic
/// slope DERIVED = 2/π (2026-06-15, PROOF_COHERENCE_HORIZON_SLOPE: the resummed coherence-range ladder gives
/// λ²+8γλ+4J²q², EP at Q*=2/q_min → 2N/π; the ring sibling is exactly half, 1/π).
/// Tier1Candidate — NOT for the gap-dominance (now PROVEN, docs/proofs/PROOF_CHAIN_GAP_DOMINANCE.md, so the
/// parent <see cref="ClockHandLadderClaim"/> graduated to Tier1Derived), but for its OWN open piece: the
/// ring 2-excitation (2,2)/(N−2,N−2) doublet V-Effect seam (NOT half-filling: at N=6 the (2,2) block's
/// handover Q_h matches the full-L survivor, the half-filling (3,3) block does not; the sector gate in
/// docs/proofs/PROOF_RING_HANDOVER_SLOPE.md). That proof derives the seam's slope Q_h → N·√3/(2π) (the
/// darkness-1 sibling of Q*) but its review is pending; this claim graduates when it is reviewed. Concrete
/// carbon-ring probe: benzene C₆ Q* = 1.609 (the benzene node, simulations/carbon/benzene_two_clocks.py).</para>
///
/// <para>Live witnesses: <c>inspect --root horizon</c>
/// (<c>compute/RCPsiSquared.Diagnostics/Foundation/CoherenceHorizonWitness.cs</c>, the ladder + the eig
/// phase-rigidity reading) and <c>inspect --root epcharacter</c>
/// (<c>compute/RCPsiSquared.Diagnostics/Foundation/EpCharacterWitness.cs</c>, the ARTIFACT-FREE EP-
/// character verdict: the non-eig sibling of <see cref="Numerics.PhaseRigidity"/> that confirms the
/// EP is genuinely DEFECTIVE — Riesz ‖P‖, departure-from-normality ≈ 4, geometric mult 1 &lt; algebraic
/// mult 2 — gate-first on a known defective + diabolic toy).</para>
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
        : base("The coherence horizon Q*(N) = 1 / √2 / 1.8787 / 2.3737 for N=2..5, verified equal to the carbon " +
               "Frost-Hückel coherent↔incoherent threshold (√2 / 1.879 / 2.372 at N=3/4/5) under the label swap " +
               "J ↔ |β| (the cross-substrate identity); N=2 (Q*=1) is the exceptional point, the base rung the " +
               "carbon polyene layer (N≥3) cannot reach. The mode that coalesces at Q*(N) is the {0,2}-coherence " +
               "at ALL N=2..5, a genuine √-EP; the band edge is the co-located γ-protected survivor (no bifurcation " +
               "at N=4). Closed form (2026-06-13): Q*(N) reduces to the single-excitation (Haken-Strobl) Liouvillian; " +
               "N=2,3 the pair are roots of λ²+4γλ+c·J²=0 (c=4,2) so Q*=2/√c=1,√2 exactly; N≥4 collectively dressed, " +
               "the exact condition transcendental; the asymptotic slope is DERIVED = 2/π " +
               "(PROOF_COHERENCE_HORIZON_SLOPE: the resummed coherence-ladder dispersion λ²+8γλ+4J²q², EP at Q*=2/q_min → 2N/π).",
               Tier.Tier1Candidate,
               "docs/ANALYTICAL_FORMULAS.md F2b corollary \"Coherence horizon Q*(N)\" + " +
               "docs/carbon/FROST_CIRCLE_AS_THE_CLOCK_FACE.md (the carbon-layer twin) + " +
               "compute/RCPsiSquared.Diagnostics/Foundation/CoherenceHorizonWitness.cs (CoherenceHorizonWitness, inspect --root horizon) + " +
               "docs/proofs/PROOF_COHERENCE_HORIZON_SLOPE.md (the asymptotic slope 2/π) + simulations/coherence_horizon_slope.py")
    {
        Horizon = horizon ?? throw new ArgumentNullException(nameof(horizon));
        BandEdge = bandEdge ?? throw new ArgumentNullException(nameof(bandEdge));
    }

    public override string DisplayName => "The coherence horizon Q*(N): the carbon coherent↔incoherent threshold made exact";

    public override string Summary =>
        $"Q*(N) = 1/√2/1.8787/2.3737 (N=2..5) = the carbon Frost-Hückel coherent↔incoherent threshold under J ↔ |β|; " +
        $"N=2 (Q*=1) is the EP base; the {{0,2}}-coherence coalesces (√-EP) at every N. Closed form (2026-06-13): Q*(N) is " +
        $"the single-excitation (Haken-Strobl) EP, λ²+4γλ+c·J²=0 at N=2,3 (Q*=2/√c=1,√2), transcendental N≥4; " +
        $"asymptotic slope 2/π (derived) ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("the ladder = carbon (verified N=2..5)",
                summary: "Q*(N) = 1 / √2 / 1.8787 / 2.3737 for N=2..5, bisected on Symphony.Clock.Omega at J=1, equals " +
                         "the carbon Frost-Hückel coherent↔incoherent threshold (√2 / 1.879 / 2.372 at N=3/4/5) under the " +
                         "label swap J ↔ |β|: the XY chain's coherence horizon IS the Frost-Hückel threshold of the same " +
                         "polyene. Live: inspect --root horizon.");
            yield return new InspectableNode("the EP base (N=2)",
                summary: "Q*(2) = 1 is the exceptional point itself (γ=J), the base rung the carbon polyene layer (N≥3) " +
                         "cannot reach (a polyene needs ≥3 sites); the quantum side supplies it.");
            yield return new InspectableNode("the band-edge coincidence (N=2,3)",
                summary: "Q*(N) = 2cos(π/(N+1)) at N=2,3 ONLY: 1 = 2cos60°, √2 = 2cos45°. A low-N accident, departing at " +
                         "N≥4 (Q*(4)=1.8787 ≠ φ=1.618); this is why the √2 looked exact at N=3 (the clean 2×2) and the rest is transcendental.");
            yield return new InspectableNode("the EP mechanism: the {0,2}-coherence coalesces at every N",
                summary: "corrected 2026-06-13 (phase-rigidity review, superseding an earlier \"sector bifurcation at N=4\" " +
                         "reading that was an argmax-Re / Im-tracking artifact). The mode that coalesces at Q*(N) is the " +
                         "{0,2}-coherence (population/antisymmetric block, n_diff histogram {0: ½, 2: ½}, ⟨n_diff⟩=1) at ALL " +
                         "N=2..5, a genuine square-root EP. CONFIRMED ARTIFACT-FREE (gate-first on known defective + diabolic " +
                         "toys, diabolic control at γ=0): the coalescing pair has departure-from-normality ≈ 4 as the split→0, " +
                         "geometric mult 1 < algebraic mult 2, Schur off-diagonal = the departure (the Jordan coupling), and the " +
                         "two eigenvectors merge (|cos|→1). The eig-phase-rigidity r → 0 (the older read; r at Q* ≈ 0.0/0.015/0.026 " +
                         "at N=3/4/5) is grid-sensitive (Kmax swings orders of magnitude under ΔQ=1e-3, sitting on the coalescence) " +
                         "and is the F86a-misfire-prone family, so it CORROBORATES the EP here but is not load-bearing; the " +
                         "defective-Jordan measure above is. Im ∝ √(Q−Q*). NO bifurcation at N=4. The band edge 2cos(π/(N+1)) " +
                         "(Im = φ at N=4, √3 at N=5) is the " +
                         "co-located SURVIVOR, the γ-protected |vac⟩⟨ψ_k| coherence hand (Uhr 1, r ≈ 1), sharing the gap " +
                         "Re = −2γ only because the Absorption Theorem pins both (both ⟨n_diff⟩=1). So Q*(N) is at once a " +
                         "{0,2}-coherence EP (the erasure point, climbs the ladder) and a band-edge crossing (Uhr 1 survives). " +
                         "Probe: simulations/coherence_horizon_se_block.py + simulations/review_coherence_horizon_ep.py. " +
                         "Live artifact-free witness (the non-eig sibling of PhaseRigidity, gate-first): " +
                         "inspect --root epcharacter (EpCharacterWitness), which recomputes the Riesz ‖P‖ / " +
                         "departure-from-normality / geo-vs-alg verdict at inspect time and reads DEFECTIVE at every N=2..5.");
            yield return new InspectableNode("the closed form: the single-excitation (Haken-Strobl) reduction",
                summary: "resolved 2026-06-13 (Approach A). Q*(N) reduces 4^N → N²: the coalescing mode is single-" +
                         "excitation, so Q*(N) is the EP of the single-excitation (Haken-Strobl) Liouvillian, the N-site " +
                         "dephased tight-binding chain (validated bit-for-bit as a sub-spectrum of the full L). At N=2,3 the " +
                         "coalescing pair are the roots of λ²+4γλ+c·J²=0 with c constant (sum=−4γ, product=c·J² are γ-" +
                         "independent identities; c=4, 2), so Q*=2/√c = 1, √2 EXACTLY — the structural form of the low-N " +
                         "accident (the whole clean 2×2, not just the value, exists only at N=2,3). At N≥4 the pair is " +
                         "collectively dressed (the trace departs from −4γ), no clean 2×2, the exact condition transcendental: " +
                         "a diffusive long-wavelength critical damping with asymptotic slope DERIVED = 2/π (the resummed coherence-" +
                         "range ladder dispersion λ²+8γλ+4J²q², EP at Q*=2/q_min → 2N/π; the nearest-neighbour telegrapher truncation " +
                         "gives the wrong √2/π; canonical Q*(4)=1.87874, Q*(5)=2.37367). Verifiers coherence_horizon_se_block.py + " +
                         "coherence_horizon_slope.py (self-validating); proof PROOF_COHERENCE_HORIZON_SLOPE.md.");
            yield return new InspectableNode("the benzene ring instance + the V-Effect seam (ring-specific)",
                summary: "the ring carbon case (2026-06-13, verifier simulations/carbon/benzene_two_clocks.py). " +
                         "Benzene C₆'s own coherence horizon (Uhr 2, the single-excitation {0,2}-EP via the same " +
                         "Haken-Strobl reduction) is Q* = 1.609, transcendental like N≥4 and below every open polyene " +
                         "(the closed ring, beating at the full 2|β| band edge, holds coherence to a lower Q). But " +
                         "benzene (even N) is where the 2-excitation V-Effect seam opens: the mode that " +
                         "overtakes the band-edge beat in the FULL 4^6 Liouvillian is a DOUBLE-excitation coherence " +
                         "(filling sector (2,2)/(4,4)), so the full-L handover (~1.95) SPLITS from the clean SE-EP Uhr 2 " +
                         "(1.609). For the open chains the Absorption Theorem co-locates the two at Re = −2γ (that IS the " +
                         "ladder); benzene breaks it. The split is RING-SPECIFIC: the open even-N chain N=6 overtaker " +
                         "spreads across all fillings at its own SE-EP, so the double-excitation seam is a feature of the " +
                         "closed ring's interior 2-excitation doublet, not even N alone. Aromaticity (Hückel 4n+2) " +
                         "is NOT the discriminant (tested C4/C8 in simulations/carbon/aromatic_ring_v_effect.py): the " +
                         "seam is RING-UNIVERSAL, the off-centre (2,2)/(N−2,N−2) doublet (not half-filling; the " +
                         "sector gate in PROOF_RING_HANDOVER_SLOPE.md), a sibling of the " +
                         "incompleteness V-Effect (docs/HIERARCHY_OF_INCOMPLETENESS.md). The 4n anti-aromatics C4 and C8 " +
                         "do NOT group (C4 is a small-ring anomaly, its seam dominates even at weak dephasing). Carbon " +
                         "doc: FROST_CIRCLE_AS_THE_CLOCK_FACE.md.");
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
