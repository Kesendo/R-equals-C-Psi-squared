using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>The stable typed home of the repaired CΨ-envelope boundary. The historical Envelope
/// package is false in its universal pointwise, forward-invariant/absorbing, and local-control or
/// trajectory readings. Exact local Markov counterexamples establish those negative statements.
/// They do not settle the narrower autonomous sequence-of-local-maxima question: the autonomous N=2
/// decreasing-successive-maxima statement remains unproved by the old argument and is not called
/// refuted here.
///
/// <para>What survives positively is deliberately smaller: the named Bell+ channel formulae (with
/// F25's pure-Z closed form as the direct typed parent), one-quarter as the algebraic bilinear maxval,
/// and the conditional convergence implication. If a continuous trajectory converges to ρ* with
/// CΨ(ρ*) &lt; 1/4, then a start above 1/4 eventually stays below it. A finite N/Q/K rise atlas is a
/// reading, not a theorem or an absence proof.</para>
///
/// <para>Anchors: <c>docs/proofs/PROOF_MONOTONICITY_CPSI.md</c> (named channel formulae, failed
/// historical proof step, and exact counterexamples) +
/// <c>compute/RCPsiSquared.Diagnostics/Foundation/QuarterEnvelope.cs</c> +
/// <c>compute/RCPsiSquared.Diagnostics/Foundation/EnvelopeTheoremWitness.cs</c>.</para></summary>
public sealed class CpsiEnvelopeTheoremClaim : Claim
{
    /// <summary>Direct parent: F25's named Bell+/local-Z closed form.</summary>
    public F25CPsiBellPlusPi2Inheritance F25 { get; }

    /// <summary>Direct parent: one-quarter as an algebraic bilinear maxval, not a dynamical absorber.</summary>
    public QuarterAsBilinearMaxvalClaim Quarter { get; }

    public CpsiEnvelopeTheoremClaim(
        F25CPsiBellPlusPi2Inheritance f25,
        QuarterAsBilinearMaxvalClaim quarter)
        : base("The historical Envelope package is false in its universal pointwise-monotonicity, forward-invariant/absorbing-boundary, and local-control/trajectory readings. The " +
               "autonomous N=2 decreasing-successive-maxima statement remains unproved by the old " +
               "argument; the exact counterexamples do not themselves refute that narrower statement. " +
               "Retained: named Bell+ channel formulae, the algebraic quarter maxval, and the conditional " +
               "convergence implication (continuous rho(t) -> rho*, Cpsi(rho*) < 1/4, start above => " +
               "eventually stays below). Finite N/Q/K atlas rows are numerical readings only.",
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

    /// <summary>Build with one shared algebraic-quarter parent across F25 and the direct edge.</summary>
    public static CpsiEnvelopeTheoremClaim Build()
    {
        var ladder = new Pi2DyadicLadderClaim();
        var quarter = new QuarterAsBilinearMaxvalClaim();
        var f25 = new F25CPsiBellPlusPi2Inheritance(ladder, quarter);
        return new CpsiEnvelopeTheoremClaim(f25, quarter);
    }

    public static CpsiEnvelopeTheoremClaim Shared { get; } = Build();

    public override string DisplayName =>
        "CΨ envelope boundary (historical theorem name retained; universal package false)";

    public override string Summary =>
        "the universal pointwise/absorbing/local-control Envelope package is false; the old proof does " +
        "not establish the autonomous N=2 peak-sequence claim. Named Bell+ formulae, the algebraic " +
        "quarter, and a conditional convergence implication survive. The live object is a finite rise " +
        $"atlas, not a theorem verdict ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("the repaired negative boundary",
                summary: "exact examples give CΨ'(0)=+1/6 under fixed local Z dephasing, show a local " +
                         "basis change 0 -> 1/3, flip a derivative with an active Z pulse, and cross " +
                         "upward from below under one fixed local Markov semigroup. They refute the " +
                         "historical package but do not settle finite local-maxima sequences.");
            yield return new InspectableNode("the surviving conditional statement",
                summary: "continuity plus rho(t) -> rho* with CΨ(rho*) < 1/4 implies that any trajectory " +
                         "starting above 1/4 eventually stays below; locality or Markovianity alone does not.");
            yield return F25;
            yield return Quarter;
        }
    }
}
