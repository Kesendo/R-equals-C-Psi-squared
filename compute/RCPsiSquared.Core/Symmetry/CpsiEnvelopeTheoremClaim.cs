using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>The CΨ Envelope Theorem in the typed graph: for any 2-qubit state (N=2) under any
/// Hamiltonian + local Z-dephasing the local maxima of CΨ form a strictly non-increasing sequence
/// (¼ is the absorbing boundary) — proven Tier-1 (PROOF_MONOTONICITY_CPSI Part 5). The N≥3 FULL-state
/// envelope is OPEN: the earlier "verified N=3-5" was channel-monotonicity + GHZ/W subsystem-PAIRS
/// (2-qubit reduced densities staying below ¼), NOT a full-state envelope test. At N≥4 strong coupling
/// (J≫γ) the full-state envelope GENUINELY RISES (the live <c>EnvelopeTheoremWitness</c>): the internal
/// J-coupling acts as the Part-6 coherence injector (Corollary 3, coupled resonators), internalized;
/// N=3 holds in the same regime. The reduced carrier-pair CΨ is likewise outside the N=2 theorem — its
/// beat envelope rises (the freedom).
///
/// <para>Two typed parent edges: <see cref="F25CPsiBellPlusPi2Inheritance"/> (the 2-qubit Bell+ closed
/// form the proof's base case rests on) and <see cref="QuarterAsBilinearMaxvalClaim"/> (¼ as the
/// bilinear maxval, the absorbing boundary).</para>
///
/// <para>Anchors: <c>docs/proofs/PROOF_MONOTONICITY_CPSI.md</c> (Part 5 the 2-qubit theorem, Part 6 +
/// Corollary 3 the coherence-injection mechanism) +
/// <c>compute/RCPsiSquared.Diagnostics/Foundation/QuarterEnvelope.cs</c> +
/// <c>compute/RCPsiSquared.Diagnostics/Foundation/EnvelopeTheoremWitness.cs</c>.</para></summary>
public sealed class CpsiEnvelopeTheoremClaim : Claim
{
    /// <summary>Parent: the 2-qubit Bell+ closed form CΨ(t)=f(1+f²)/6 the proof's base case is built on.</summary>
    public F25CPsiBellPlusPi2Inheritance F25 { get; }

    /// <summary>Parent: ¼ as the bilinear maxval, the absorbing boundary the envelope settles into.</summary>
    public QuarterAsBilinearMaxvalClaim Quarter { get; }

    public CpsiEnvelopeTheoremClaim(
        F25CPsiBellPlusPi2Inheritance f25,
        QuarterAsBilinearMaxvalClaim quarter)
        : base("For any 2-qubit state (N=2) under any Hamiltonian + local Z-dephasing the local maxima of " +
               "CΨ form a strictly non-increasing sequence (the Envelope Theorem); ¼ is the absorbing " +
               "boundary. Proven Tier-1 for N=2 (Part 5). The N≥3 full-state envelope is OPEN: the earlier " +
               "'verified N=3-5' was channel-monotonicity + GHZ/W subsystem-pairs (2-qubit reduced), NOT a " +
               "full-state envelope test; at N≥4 strong coupling it GENUINELY RISES (witness: the internal " +
               "J-coupling is the Part-6 coherence injector; N=3 holds). The reduced-subsystem CΨ likewise " +
               "rises (the freedom). The N≥3 rise BOUNDARY is now charted (EnvelopeBoundaryTests): a pure " +
               "(N, Q=J/γ) observable (the J- and γ-sweeps collapse, certified bit-identical), with an N≥4 " +
               "floor (N=3 never rises, Q_c(3)=∞) and, above it, a threshold Q_c(N) that climbs with N " +
               "(Q_c(4)≈27, Q_c(5)≈45).",
               Tier.Tier1Derived,
               "docs/proofs/PROOF_MONOTONICITY_CPSI.md + " +
               "compute/RCPsiSquared.Diagnostics/Foundation/QuarterEnvelope.cs + " +
               "compute/RCPsiSquared.Diagnostics/Foundation/EnvelopeTheoremWitness.cs + " +
               "compute/RCPsiSquared.Diagnostics.Tests/Foundation/EnvelopeBoundaryTests.cs + " +
               "experiments/ENVELOPE_RISE_BOUNDARY.md")
    {
        F25 = f25 ?? throw new ArgumentNullException(nameof(f25));
        Quarter = quarter ?? throw new ArgumentNullException(nameof(quarter));
    }

    /// <summary>Builds the claim with fresh parent chains (mirroring ApproachFamilyCarrierClaim.Build()):
    /// one shared QuarterAsBilinearMaxvalClaim feeds both F25 and this claim's ¼-boundary edge.</summary>
    public static CpsiEnvelopeTheoremClaim Build()
    {
        var ladder = new Pi2DyadicLadderClaim();
        var quarter = new QuarterAsBilinearMaxvalClaim();
        var f25 = new F25CPsiBellPlusPi2Inheritance(ladder, quarter);
        return new CpsiEnvelopeTheoremClaim(f25, quarter);
    }

    /// <summary>Shared singleton; a structural statement, block-independent.</summary>
    public static CpsiEnvelopeTheoremClaim Shared { get; } = Build();

    public override string DisplayName =>
        "CΨ Envelope Theorem (2-qubit peaks non-increasing — proven Tier-1 for N=2; N≥3 OPEN, RISES at N≥4 strong coupling; the reduced pair escapes: the freedom)";

    public override string Summary =>
        "for any 2-qubit state the local maxima of CΨ form a non-increasing sequence under any H + local " +
        "Z-dephasing (¼ absorbing; proven Tier-1 for N=2). The N≥3 full-state envelope is OPEN and GENUINELY " +
        "RISES at N≥4 strong coupling (internal J = the Part-6 injector; N=3 holds); the reduced carrier-pair " +
        "also rises (the freedom). The boundary is charted (an N≥4 floor and a Q_c(N) contour climbing with " +
        $"N): see EnvelopeTheoremWitness / EnvelopeBoundaryTests ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("the theorem (2-qubit, proven) + the N≥3 boundary (open)",
                summary: "for any 2-qubit state under any Hamiltonian + local Z-dephasing the local maxima of CΨ " +
                         "form a strictly non-increasing sequence (¼ absorbing) — proven Tier-1 for N=2 " +
                         "(PROOF_MONOTONICITY_CPSI Part 5). The N≥3 full-state envelope is OPEN: the 'verified " +
                         "N=3-5' was channel-monotonicity + GHZ/W subsystem-pairs (2-qubit reduced), not a " +
                         "full-state test; at N≥4 strong coupling it RISES (witness; internal J = the Part-6 " +
                         "injector, Corollary 3; N=3 holds).");
            yield return new InspectableNode("the boundary (the freedom)",
                summary: "reduced subsystems escape the theorem: the carrier-pair CΨ is a reduced open " +
                         "subsystem with no such guarantee, and its beat envelope can genuinely RISE (beating). " +
                         "The live EnvelopeTheoremWitness shows global = 0 rises at N=3 vs the full state rising " +
                         "at N≥4 (and the reduced pair rising).");
            yield return F25;
            yield return Quarter;
        }
    }
}
