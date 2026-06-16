using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>The topology-general XY band edge (Tier 1 candidate). The single-excitation band edge (the
/// largest |E_k| of the vac↔single-excitation coherences) = J × the hopping graph's adjacency spectral
/// radius ρ: chain 2cos(π/(N+1)), star √(N−1), ring 2 — generalizing the chain-only F2b band edge to any
/// topology. This is the Im/L_H side; the Re=−2γ floor (the band-edge coherence has n_XY=1) is the
/// Absorption Theorem, cited not re-derived. Whether the clock/seam reads that band edge as ω_mem is a
/// topology-specific regime (witnessed, not separately claimed): a Q-horizon for chain (all N), star
/// (N≤5), and odd rings; a structural ceiling for star N≥6 (strict gap saturates below 2γ); a
/// co-occupied-floor mismatch for ring N=4 (a (2,2) mode at the floor with Im=2√2·J > band edge).
///
/// <para>Tier: Tier1Candidate, not Derived. The band-edge LAW (SE block of XY = J·adjacency) is exact, and
/// the Re=−2γ floor rests on the Tier1Derived <see cref="AbsorptionTheoremClaim"/>. But the typed parent edge
/// to <see cref="ClockHandLadderClaim"/> (the chain instance this generalizes) makes the claim rest on a
/// Tier1Candidate, and the registry's tier-inheritance invariant forbids a child stronger than its weakest
/// parent. This is honest: the gap-dominance "map" reuses exactly the chain's open gap-dominance proof (the
/// general proof that the max frequency in the protected n_XY=1 subspace is the band edge — still open in
/// ClockHandLadderClaim). The sibling <c>CoherenceHorizonClaim</c>, the other ClockHandLadder child, is
/// Tier1Candidate for the same reason. (The design spec said Tier1Derived; the tier-inheritance gate caught
/// the overlooked candidate parent — see commit message.)</para>
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
               Tier.Tier1Candidate,
               "docs/superpowers/specs/2026-06-16-topology-band-edge-design.md + " +
               "docs/proofs/PROOF_ABSORPTION_THEOREM.md (§4.3 gap=2γ, §4.5 n_XY=Hamming — the −2γ floor) + " +
               "compute/RCPsiSquared.Diagnostics/Foundation/TopologyBandEdgeWitness.cs (inspect --root bandedge) + " +
               "simulations/_topology_band_edge_review.py (the gate-first verifier, full Q-sweep to Q=1000)")
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
                         "(strict gap saturates ≈0.80·2γ). Co-occupied floor: ring N=4 (a (2,2) mode at −2γ, " +
                         "Im=2√2·J > band edge). Live: inspect --root bandedge.");
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
