using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>The typed ARGUMENT of INCOMPLETENESS_PROOF (2026-07-02): dephasing noise cannot
/// originate WITHIN the d(d−2)=0 ontology, so it must come from OUTSIDE. The dimension conclusion
/// d∈{0,2} is already carried by the parents (<see cref="PolynomialFoundationClaim"/> is literally
/// d²−2d=0 with d=1 excluded; <see cref="QubitDimensionalAnchorClaim"/> is the d=2 anchor); this
/// claim types the FULL five-candidate elimination that turns "d∈{0,2}" into "the noise is external".
///
/// <para>The five internal candidates for the noise origin, and why each fails: (1) internal /
/// self-generated — a CONSTRAINT, not an elimination ([Π², L]=0 decouples the parity sectors, typed
/// as F63, leaving underdetermination); (2) single-qubit decay — ELIMINATED (γ_eff=0, 0-of-16
/// palindromic pairs); (3) many-qubit bath — ELIMINATED (N instances of (2) + infinite regress);
/// (4) nothing, d=0 — ELIMINATED by definition (no properties); (5) something other than a qubit or
/// nothing — ELIMINATED by the algebra (d²−2d=0 has only the roots {0,2}). No internal source
/// remains, so the noise is external: the incompleteness the proof names (the V-Effect one dimension
/// up).</para>
///
/// <para>Live witness (light: recomputes Candidate 5's dimension algebra, surfaces the rest):
/// <c>compute/RCPsiSquared.Diagnostics/Foundation/NoiseOriginExclusionWitness.cs</c>,
/// <c>inspect --root noise-origin</c>. The heavier Candidate-2 process-tomography compute is
/// deferred (see the proof + simulations/failed_third.py).</para></summary>
public sealed class NoiseOriginExclusionClaim : Claim
{
    /// <summary>Typed parent: the minimum-memory polynomial d²−2d=0 (Candidate 5's algebra).</summary>
    public PolynomialFoundationClaim Polynomial { get; }

    /// <summary>Typed parent: the d=2 qubit dimensional anchor (the surviving internal dimension).</summary>
    public QubitDimensionalAnchorClaim Dimension { get; }

    public NoiseOriginExclusionClaim(PolynomialFoundationClaim polynomial, QubitDimensionalAnchorClaim dimension)
        : base("Noise-origin exclusion: no internal source in the d(d−2)=0 ontology generates the dephasing; it comes from OUTSIDE (the 5-candidate elimination)",
               Tier.Tier1Derived,
               "docs/proofs/INCOMPLETENESS_PROOF.md")
    {
        Polynomial = polynomial ?? throw new ArgumentNullException(nameof(polynomial));
        Dimension = dimension ?? throw new ArgumentNullException(nameof(dimension));
    }

    public override string DisplayName =>
        "Noise-origin exclusion: no internal source generates the dephasing (5-candidate elimination) ⟹ it is external";

    public override string Summary =>
        "the five candidates for an INTERNAL noise origin all fail inside the d(d−2)=0 ontology: (1) self-generated " +
        "is only a constraint ([Π²,L]=0, F63); (2) single-qubit decay is eliminated (γ_eff=0, 0/16 palindromic); " +
        "(3) many-qubit bath is eliminated (N×(2) + infinite regress); (4) nothing (d=0) is eliminated by definition; " +
        "(5) any other dimension is eliminated by the algebra (d²−2d=0 has only roots {0,2}). So the dephasing noise " +
        "must originate OUTSIDE the system — the incompleteness (V-Effect one dimension up). Dimension algebra typed " +
        "as the parents; this node types the full elimination (live: NoiseOriginExclusionWitness).";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return Polynomial;
            yield return Dimension;
            yield return new InspectableNode("the 5-candidate elimination",
                summary: "(1) internal/self-generated: [Π²,L]=0 decouples the parity sectors (F63) — a CONSTRAINT, " +
                         "the elimination is carried by (2)-(3). (2) single-qubit decay: ELIMINATED, γ_eff=0, process " +
                         "tomography 0/16 palindromic (ref Z-dephasing 16/16). (3) many-qubit bath: ELIMINATED, N " +
                         "instances of (2), infinite regress. (4) nothing (d=0): ELIMINATED by definition (no " +
                         "properties). (5) any other dimension: ELIMINATED by the algebra, d²−2d=0 ⟹ d∈{0,2}, d=1 and " +
                         "d≥3 excluded.");
            yield return new InspectableNode("the conclusion: external noise",
                summary: "no internal source survives, so the dephasing noise originates OUTSIDE the d(d−2)=0 " +
                         "ontology. This is the incompleteness the proof names: the noise the system cannot make from " +
                         "within is the V-Effect one dimension up (INCOMPLETENESS_PROOF.md).");
            yield return new InspectableNode("live witness",
                summary: "NoiseOriginExclusionWitness (inspect --root noise-origin) recomputes Candidate 5 live " +
                         "(d²−2d across d, roots {0,2}, d=1/d≥3 excluded) and enumerates all five candidates; the " +
                         "heavier Candidate-2 process-tomography compute is deferred in this light reading.");
        }
    }
}
