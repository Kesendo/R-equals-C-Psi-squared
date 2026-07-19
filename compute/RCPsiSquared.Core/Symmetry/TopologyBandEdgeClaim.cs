using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>The topology-general XY band edge (Tier 1 derived). The single-excitation band edge (the
/// largest |E_k| of the vac↔single-excitation coherences) = J × the hopping graph's adjacency spectral
/// radius ρ: chain 2cos(π/(N+1)), star √(N−1), ring 2 — generalizing the chain-only F2b band edge to any
/// topology. This is the Im/L_H side; the Re=−2γ floor (the band-edge coherence has n_XY=1) is the
/// Absorption Theorem, cited not re-derived. Whether the clock/seam reads that band edge as ω_mem is a
/// topology-specific regime (witnessed, not separately claimed): a Q-horizon for chain (all N), star
/// (N≤5), and odd rings; a structural ceiling for star N≥6 (strict gap saturates below 2γ); a
/// co-occupied-floor mismatch for ring N=4 (a (2,2) mode at the floor with Im=2√2·J > band edge). The
/// quantitative Re-side ceiling (the closed forms g2(K_N)=4/N, star √(N−1)→4/(N−1), K_4=2−2/√3, derived
/// from the commutant rep structure) is the companion <see cref="StructuralCeilingClaim"/> (F122).
///
/// <para>Map note (why the magnitude ρ suffices): the band edge uses only ρ = |E_max|. For BIPARTITE hopping
/// graphs (chain, star, even rings) the single-particle spectrum is ±symmetric under the chiral K-mirror
/// (<see cref="ChiralKClaim"/>: K·H·K = −H ⟹ E ↔ −E, the sibling of F1 acting on H), so ρ = E_max = −E_min
/// and the band edge carries a mirror partner at −ρ; non-bipartite graphs (complete K_N, odd rings) have no
/// such ±symmetry. The ring-4 (2,2) co-occupier 2√2·J is the anti-periodic even-sector image of this same
/// chiral structure (PROOF_STRUCTURAL_CEILING §4).</para>
///
/// <para>Tier: Tier1Derived (2026-06-16). The band-edge LAW (SE block of XY = J·adjacency) is exact; the
/// Re=−2γ floor is the Tier1Derived <see cref="AbsorptionTheoremClaim"/>; and the chain gap-dominance that its
/// parent <see cref="ClockHandLadderClaim"/> had left open is now PROVEN (docs/proofs/PROOF_CHAIN_GAP_DOMINANCE.md:
/// the exact-(−2γ) modes are free fermions c_k^(†)·f(N_tot) at ±E_k, max E1), so ClockHandLadder graduated to
/// Tier1Derived and this inherited cap is lifted. The topology-specific gap-dominance map (below) is
/// witnessed, and independently understood: the star has no coherence horizon (flat band,
/// PROOF_STRUCTURAL_CEILING §7), and complete/star structurally ceiling (g2=4/N, 4/(N−1); F122 /
/// <see cref="StructuralCeilingClaim"/>). The sibling <c>CoherenceHorizonClaim</c> stays Tier1Candidate for
/// its OWN open piece (the ring 2-excitation (2,2)/(N−2,N−2) doublet V-Effect seam, pending the
/// PROOF_RING_HANDOVER_SLOPE review), not the gap-dominance.</para>
///
/// <para>Typed parents: <see cref="ClockHandLadderClaim"/> (the chain instance it generalizes — its
/// ω_mem = 2J·cos(π/(N+1)) is the chain ρ), <see cref="AbsorptionTheoremClaim"/> (the Re=−2γ floor).</para></summary>
public sealed class TopologyBandEdgeClaim : Claim
{
    /// <summary>Parent: the chain instance this generalizes (its ω_mem = 2J·cos(π/(N+1)) is the chain ρ).</summary>
    public ClockHandLadderClaim ChainInstance { get; }

    /// <summary>Parent: the Absorption Theorem — the band-edge coherence has n_XY=1, so Re=−2γ on any graph
    /// (the floor the map's gap-dominance is read against). Cited, not re-derived.</summary>
    public AbsorptionTheoremClaim Absorption { get; }

    public TopologyBandEdgeClaim(ClockHandLadderClaim chainInstance, AbsorptionTheoremClaim absorption)
        : base("Topology band edge: the XY single-excitation band edge = J × the hopping graph's adjacency " +
               "spectral radius ρ (chain 2cos(π/(N+1)), star √(N−1), ring 2), generalizing the chain-only F2b " +
               "band edge to any topology (the Im/L_H side; the Re=−2γ floor is the Absorption Theorem, n_XY=1). " +
               "Gap-dominance (whether ω_mem reads the band edge) is topology-specific and witnessed: Q-horizon " +
               "for chain (all N) / star (N≤5) / odd rings; structural ceiling for star N≥6; co-occupied floor " +
               "for ring N=4 (a (2,2) mode at −2γ with Im=2√2·J > band edge).",
               Tier.Tier1Derived,
               "docs/proofs/PROOF_ABSORPTION_THEOREM.md (§4.3 gap=2γ, §4.5 n_XY=Hamming — the −2γ floor) + " +
               "docs/proofs/PROOF_CHAIN_GAP_DOMINANCE.md (the chain gap-dominance proof that lifts the cap) + " +
               "compute/RCPsiSquared.Diagnostics/Foundation/TopologyBandEdgeWitness.cs (inspect --root bandedge) + " +
               "simulations/topology_band_edge_review.py (the gate-first verifier, full Q-sweep to Q=1000)")
    {
        ChainInstance = chainInstance ?? throw new ArgumentNullException(nameof(chainInstance));
        Absorption = absorption ?? throw new ArgumentNullException(nameof(absorption));
    }

    public override string DisplayName => "Topology band edge: J × the hopping-graph adjacency spectral radius";

    public override string Summary =>
        $"band edge = J·ρ (chain 2cos(π/(N+1)), star √(N−1), ring 2); gap-dominance topology-specific " +
        $"(chain/star≤5/odd-ring horizon, star≥6 ceiling, ring-4 co-occupied floor) ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("the law (Im side)",
                summary: "band edge = J × adjacency spectral radius; chain 2cos(π/(N+1)), star √(N−1), ring 2 " +
                         "(SE block of XY = J × adjacency). Generalizes F2b (chain) to any topology.");
            yield return new InspectableNode("the floor (Re side) = the Absorption Theorem",
                summary: "the band-edge coherence |vac⟩⟨ψ_k| has n_XY=1, so it sits at Re=−2γ for any topology " +
                         "(PROOF_ABSORPTION_THEOREM.md §4.3/§4.5). Cited, not re-derived.");
            yield return new InspectableNode("the gap-dominance map (witnessed)",
                summary: "Q-horizon: chain (all N), star (N≤5), odd rings. Structural ceiling: star N≥6 " +
                         "(strict gap saturates g2=4/(N−1)·2γ, N-dependent: N=6→0.80). Co-occupied floor: ring N=4 (a (2,2) mode at −2γ, " +
                         "Im=2√2·J > band edge). Live: inspect --root bandedge.");
            yield return new InspectableNode("ring gap-dominance (the dihedral lock, characterized 2026-06-17)",
                summary: "PROOF_RING_GAP_DOMINANCE.md / ring_gap_dominance.py (gate-first N=3..6): max|Im| over the " +
                         "exact-(−2γ) ring modes = 2J = J·ρ (the periodic band top = the k=0 uniform single-excitation " +
                         "mode fixed by C_N), reached via the (0,1) sector for general N — for all N EXCEPT N=4, the " +
                         "lone exception where the half-filling (2,2) {0,2} √-EP reaches 2√2·J > 2J (the same (2,2) " +
                         "sector as K_4's ceiling; the ring analogue of the chain's N=3 special, but ABOVE the band " +
                         "top not below). The ring sibling of chain gap-dominance; both give max|Im| = J·ρ. " +
                         "Free-fermion COMPLETENESS gate-verified 2026-06-20 (ring_gap_completeness.py): the n_XY=1 " +
                         "family spans the exact-(−2γ) subspace (V_1 dim = full-L dim at N=5,6; chain sanity 32/50/72), " +
                         "so nothing exceeds 2J; dihedral lock carried to N=7 in V_1; N=4 is the lone exception (full-L " +
                         "dim > V_1, the {0,2} extras). Residual all-N step = the chain's, split by wrap-bond parity.");
            yield return ChainInstance;   // typed parent edge
            yield return Absorption;      // typed parent edge
        }
    }

    public static TopologyBandEdgeClaim Build()
    {
        var clock = ClockHandLadderClaim.Build();
        return new TopologyBandEdgeClaim(clock, clock.Absorption);
    }

    public static TopologyBandEdgeClaim Shared { get; } = Build();
}
