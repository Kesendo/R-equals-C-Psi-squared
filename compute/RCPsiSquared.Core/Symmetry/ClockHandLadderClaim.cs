using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>The two clocks (Tier 1 derived: the gap-dominance lemma is now proven via free fermions,
/// docs/proofs/PROOF_CHAIN_GAP_DOMINANCE.md). The Symphony clock node's coherence hand ω_mem is, for N≥3, the F2b band edge 2J·cos(π/(N+1))
/// (√2 / φ / √3 at N=3/4/5), and it is γ-protected: the |vac⟩⟨ψ_k| modes are simultaneous eigenoperators
/// of L_D (eigenvalue −2γ, the Absorption Theorem) and L_H (eigenvalue iE_k, F2b), so γ moves the Takt hand
/// (Gap = 2γ) but not the coherence hand. For N=2 the hand is γ-pulled to 2√(J²−γ²) (a different block,
/// exactly derivable) and stops at the exceptional point Q=1. The γ-protection for N≥3 is the carrier-
/// blindness of the two-tempo certification read from the other side.
///
/// <para>What is DERIVED: the γ-protection mechanism (two Tier-1 results compose), the N=2 closed form, and
/// — as of 2026-06-16 — the gap-dominance lemma itself (docs/proofs/PROOF_CHAIN_GAP_DOMINANCE.md, verifier
/// simulations/chain_gap_dominance.py). The off-gap argument's first half was already closed: the Q-floor
/// below which the band edge stops being the gap IS the Coherence Horizon Q*(N), derived → 2N/π
/// (PROOF_COHERENCE_HORIZON_SLOPE). The remaining gate — that the MAX|Im| among the exact-(−2γ) modes is E1 —
/// is now proven: although L_H leaks n_XY=1→3 (the dephased chain is interacting), ON the exact-(−2γ)
/// subspace L_D=−2γ is a scalar so L=L_H−2γ acts freely; via Jordan-Wigner the modes c_k^(†)·f(N_tot)
/// (a single fermion dressed by a function of the total number) are exactly that subspace, oscillating at the
/// single-particle energies ±E_k≤E1 (they SPAN it for N≥4; N=3 adds a {0,2} √-EP family at
/// √(E1²−(2γ)²)<E1). So max|Im|=E1: the band edge is gap-dominant. Hence Tier1Derived. (Scope: the chain;
/// Jordan-Wigner is 1D. The N=3 extra family is one more entry in the chain's list of N=3 accidents.)</para>
///
/// <para>Live witness: <c>inspect --root clock</c>
/// (<c>compute/RCPsiSquared.Diagnostics/Foundation/ClockHandLadderWitness.cs</c>).</para>
///
/// <para>Typed parents: <see cref="F2bXyChainSpectrumPi2Inheritance"/> (the band edge the hand is protected
/// to), <see cref="AbsorptionTheoremClaim"/> (the −2γ that does the protecting),
/// <see cref="UniversalCarrierClaim"/> (carrier-blindness).</para>
///
/// <para>Anchors: <c>docs/ANALYTICAL_FORMULAS.md</c> F2b corollary +
/// <c>reflections/ON_HOW_GAMMA_BECAME_THE_TICK.md</c> + <c>reflections/ON_TWO_TIMES.md</c> +
/// <c>reflections/ON_HOW_THE_CARRIER_SHOWS_ITSELF.md</c> +
/// <c>compute/RCPsiSquared.Diagnostics/Foundation/ClockHandLadderWitness.cs</c> (the live lab).
/// The carbon-layer twin of the same ladder is <c>docs/carbon/FROST_CIRCLE_AS_THE_CLOCK_FACE.md</c>,
/// a parked seam.</para></summary>
public sealed class ClockHandLadderClaim : Claim
{
    /// <summary>Parent: the F2b band edge 2J·cos(π/(N+1)) the coherence hand is γ-protected to (N≥3).</summary>
    public F2bXyChainSpectrumPi2Inheritance BandEdge { get; }

    /// <summary>Parent: the Absorption Theorem's −2γ, the eigenvalue of L_D that does the protecting.</summary>
    public AbsorptionTheoremClaim Absorption { get; }

    /// <summary>Parent: γ₀ as Universal Carrier; the γ-protection IS its carrier-blindness from the other side.</summary>
    public UniversalCarrierClaim Carrier { get; }

    public ClockHandLadderClaim(
        F2bXyChainSpectrumPi2Inheritance bandEdge,
        AbsorptionTheoremClaim absorption,
        UniversalCarrierClaim carrier)
        : base("The two clocks: ω_mem (coherence hand) = γ-protected F2b band edge 2J·cos(π/(N+1)) for N≥3 " +
               "(√2/φ/√3, gap-dominance PROVEN via free fermions, PROOF_CHAIN_GAP_DOMINANCE); γ-pulled 2√(J²−γ²) " +
               "for N=2, stopping at the exceptional point Q=1. Takt hand = 2γ. θ = arctan(Q·cos(π/(N+1))) [N≥3] / arctan(√(Q²−1)) [N=2]",
               Tier.Tier1Derived,
               "docs/proofs/PROOF_CHAIN_GAP_DOMINANCE.md (the gap-dominance proof) + " +
               "docs/ANALYTICAL_FORMULAS.md F2b corollary + " +
               "reflections/ON_HOW_GAMMA_BECAME_THE_TICK.md + " +
               "reflections/ON_TWO_TIMES.md + " +
               "reflections/ON_HOW_THE_CARRIER_SHOWS_ITSELF.md + " +
               "compute/RCPsiSquared.Diagnostics/Foundation/ClockHandLadderWitness.cs (ClockHandLadderWitness, inspect --root clock) + " +
               "simulations/chain_gap_dominance.py (the gap-dominance gate-first verifier)")
    {
        BandEdge = bandEdge ?? throw new ArgumentNullException(nameof(bandEdge));
        Absorption = absorption ?? throw new ArgumentNullException(nameof(absorption));
        Carrier = carrier ?? throw new ArgumentNullException(nameof(carrier));
    }

    public override string DisplayName => "The two clocks: the γ-protected band-edge ladder and the N=2 exceptional point";

    public override string Summary =>
        $"ω_mem = 2J·cos(π/(N+1)) γ-protected (N≥3, √2/φ/√3, gap-dominance PROVEN, PROOF_CHAIN_GAP_DOMINANCE); " +
        $"2√(J²−γ²) γ-pulled (N=2, stops at Q=1); Takt = 2γ ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("the ladder (N≥3, γ-protected)",
                summary: "ω_mem = 2J·cos(π/(N+1)) = F2b band edge (√2/φ/√3 at N=3/4/5), held fixed under a γ-sweep " +
                         "because |vac⟩⟨ψ_k| modes are simultaneous L_D/L_H eigenoperators. Gap-dominance PROVEN " +
                         "(PROOF_CHAIN_GAP_DOMINANCE: the exact-(−2γ) modes are free fermions c_k^(†)·f(N_tot) at ±E_k, " +
                         "max E1), in the H-competitive regime above the Coherence Horizon Q*(N) (below it a slower real " +
                         "mode takes the gap and the band edge ceases to be the coherence hand). Live: inspect --root clock.");
            yield return new InspectableNode("the exceptional point (N=2)",
                summary: "ω_mem = 2√(J²−γ²), a different block (not the band edge), exactly derivable; stops at Q=1 (γ=J). " +
                         "docs/ANALYTICAL_FORMULAS.md F2b corollary.");
            yield return new InspectableNode("carrier-blindness = γ-protection",
                summary: "the N≥3 γ-protection IS the two-tempo certification's carrier-blindness read from the other " +
                         "side: in the H-competitive regime ω_mem/Q does not depend on γ, so the carrier γ₀ is invisible " +
                         "from inside.");
            yield return BandEdge;    // typed parent edge
            yield return Absorption;  // typed parent edge
            yield return Carrier;     // typed parent edge
        }
    }

    public static ClockHandLadderClaim Build()
    {
        // Shared Pi2-Foundation roots (all parameterless), reused across every parent so the graph is consistent.
        var ladder = new Pi2DyadicLadderClaim();
        var polynomial = new PolynomialFoundationClaim();
        var qubit = new QubitDimensionalAnchorClaim();

        // F2b band-edge parent: F2b(ladder, f65); f65(ladder, f66); f66(ladder, qubit).
        var f66 = new F66PoleModesPi2Inheritance(ladder, qubit);
        var f65 = new F65XxChainSpectrumPi2Inheritance(ladder, f66);
        var bandEdge = new F2bXyChainSpectrumPi2Inheritance(ladder, f65);

        // One shared AbsorptionTheoremClaim feeds both the Absorption parent edge and the UniversalCarrier parent.
        var absorption = new AbsorptionTheoremClaim(ladder);

        // UniversalCarrier(absorption, ladder, discriminant); discriminant(polynomial, qubit, ladder).
        var discriminant = new PolynomialDiscriminantAnchorClaim(polynomial, qubit, ladder);
        var carrier = new UniversalCarrierClaim(absorption, ladder, discriminant);

        return new ClockHandLadderClaim(bandEdge, absorption, carrier);
    }

    public static ClockHandLadderClaim Shared { get; } = Build();
}
