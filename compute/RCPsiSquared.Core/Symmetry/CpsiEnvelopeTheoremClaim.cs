using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>The CΨ Envelope Theorem in the typed graph: under any Hamiltonian + local Z-dephasing the
/// local maxima of the FULL-state CΨ form a strictly non-increasing sequence (¼ is the absorbing
/// boundary). Proven for 2-qubit states, Tier-2 verified N=3-5 (PROOF_MONOTONICITY_CPSI). The reduced
/// carrier-pair CΨ is explicitly OUTSIDE the theorem — its beat envelope can genuinely RISE (the
/// freedom), which the live <c>EnvelopeTheoremWitness</c> demonstrates.
///
/// <para>Two typed parent edges: <see cref="F25CPsiBellPlusPi2Inheritance"/> (the 2-qubit Bell+ closed
/// form the proof's base case rests on) and <see cref="QuarterAsBilinearMaxvalClaim"/> (¼ as the
/// bilinear maxval, the absorbing boundary).</para>
///
/// <para>Anchors: <c>docs/proofs/PROOF_MONOTONICITY_CPSI.md</c> +
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
        : base("Under any Hamiltonian + local Z-dephasing the local maxima of the full-state CΨ form a " +
               "strictly non-increasing sequence (the Envelope Theorem); ¼ is the absorbing boundary. " +
               "Proven for 2-qubit states, Tier-2 verified N=3-5. The reduced-subsystem CΨ is explicitly " +
               "NOT covered — its beat envelope can rise (the freedom).",
               Tier.Tier1Derived,
               "docs/proofs/PROOF_MONOTONICITY_CPSI.md + " +
               "compute/RCPsiSquared.Diagnostics/Foundation/QuarterEnvelope.cs + " +
               "compute/RCPsiSquared.Diagnostics/Foundation/EnvelopeTheoremWitness.cs")
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
        "CΨ Envelope Theorem (full-state peaks non-increasing — proven N=2 / verified N≥3; the reduced pair escapes: the freedom)";

    public override string Summary =>
        "the local maxima of the full-state CΨ form a non-increasing sequence under any H + local Z-dephasing " +
        "(¼ absorbing; proven 2-qubit, verified N=3-5); the reduced carrier-pair CΨ is NOT covered and its beat " +
        $"envelope can rise (the freedom) — see EnvelopeTheoremWitness ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("the theorem (global, full state)",
                summary: "under any Hamiltonian + local Z-dephasing the local maxima of CΨ form a strictly " +
                         "non-increasing sequence; ¼ is the absorbing boundary. Proven for 2-qubit states " +
                         "(PROOF_MONOTONICITY_CPSI), Tier-2 verified N=3-5.");
            yield return new InspectableNode("the boundary (the freedom)",
                summary: "reduced subsystems escape the theorem: the carrier-pair CΨ is a reduced open " +
                         "subsystem with no such guarantee, and its beat envelope can genuinely RISE (beating). " +
                         "The live EnvelopeTheoremWitness shows global = 0 rises vs the reduced pair rising.");
            yield return F25;
            yield return Quarter;
        }
    }
}
