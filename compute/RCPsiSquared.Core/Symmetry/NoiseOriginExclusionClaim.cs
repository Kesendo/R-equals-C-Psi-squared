using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>The typed ARGUMENT of INCOMPLETENESS_PROOF, and what is left of it after 2026-08-29.
/// It read: dephasing noise cannot originate WITHIN the d(d−2)=0 ontology, so it must come from
/// OUTSIDE. That elimination does NOT hold. What holds is the trace identity: the commutator part
/// of any generator is traceless and each jump costs the trace r(|tr F|² − d·tr(F†F)) ≤ 0
/// (Cauchy-Schwarz on ⟨I,F⟩, equality only for F ∝ I, which dissipates nothing), so FOR A
/// COMPLETELY POSITIVE generator trace(L) = 0 iff the dissipator vanishes iff the system is closed.
/// Equivalently: open iff some eigenvalue is off the imaginary axis. The PALINDROME is an
/// instrument rather than a premise: where the spectrum pairs, its centre is trace(L)/dim, so a
/// measured pairing lets the trace be read off a spectrum whose generator was never in hand.
/// The word delivered is OPEN, not EXTERNAL.
/// The dimension conclusion
/// d∈{0,2} is already carried by the parents (<see cref="PolynomialFoundationClaim"/> is literally
/// d²−2d=0 with d=1 excluded; <see cref="QubitDimensionalAnchorClaim"/> is the d=2 anchor); this
/// claim carried the five-candidate elimination that was to turn "d∈{0,2}" into "the noise is
/// external"; it now carries what survived the re-measurement, which is the trace identity and two
/// open candidates.
///
/// <para>The five internal candidates for the noise origin, and why each fails: (1) internal /
/// self-generated — a CONSTRAINT, not an elimination ([Π², L]=0 decouples the parity sectors, typed
/// as F63, leaving underdetermination); (2) single-qubit decay — OPEN since 2026-08-29: its test
/// reads a two-qubit MARGINAL over a still-coupled spectator, which fails to pair equally for
/// external noise (0.177), internal noise (0.116-0.149) and NO NOISE AT ALL (0.094), while the
/// three-qubit generator it builds is an exact palindrome at F137's centre; and its "internal"
/// source is a Lindblad jump, an external bath assumed rather than tested; (3) many-qubit
/// bath — the inheritance half falls with (2), the regress stands as a shift of the question;
/// (4) nothing, d=0 — ELIMINATED by definition (no properties); (5) something other than a qubit or
/// nothing — ELIMINATED by the algebra (d²−2d=0 has only the roots {0,2}). (4) and (5) are
/// definitional: they say what cannot EXIST, never that an existing qubit is not the source. So no
/// internal source is eliminated and the origin question is open; what the framework does deliver is
/// openness, and the reason it cannot deliver more is that a Lindblad dissipator assumes the
/// environment while the unitary alternative has no decay to explain.</para>
///
/// <para>Live witness (light: recomputes Candidate 5's dimension algebra, surfaces the rest):
/// <c>compute/RCPsiSquared.Diagnostics/Foundation/NoiseOriginExclusionWitness.cs</c>,
/// <c>inspect --root noise-origin</c>. Candidate 2's evidence is re-measured in
/// <c>simulations/incompleteness_candidate2_evidence.py</c> (29 gates); the original run it
/// replaces is <c>simulations/failed_third.py</c>, left as it was run.</para></summary>
public sealed class NoiseOriginExclusionClaim : Claim
{
    /// <summary>Typed parent: the minimum-memory polynomial d²−2d=0 (Candidate 5's algebra).</summary>
    public PolynomialFoundationClaim Polynomial { get; }

    /// <summary>Typed parent: the d=2 qubit dimensional anchor (the surviving internal dimension).</summary>
    public QubitDimensionalAnchorClaim Dimension { get; }

    public NoiseOriginExclusionClaim(PolynomialFoundationClaim polynomial, QubitDimensionalAnchorClaim dimension)
        : base("Openness from the trace: a completely positive generator is closed iff trace(L) = 0, so a palindrome centred away from zero certifies an OPEN system (the noise ORIGIN is a separate, open question)",
               Tier.Tier1Derived,
               "docs/proofs/INCOMPLETENESS_PROOF.md")
    {
        Polynomial = polynomial ?? throw new ArgumentNullException(nameof(polynomial));
        Dimension = dimension ?? throw new ArgumentNullException(nameof(dimension));
    }

    public override string DisplayName =>
        "Openness from the trace, exactly; the noise ORIGIN is open and the 5-candidate elimination does not close it";

    public override string Summary =>
        "what is EXACT: the commutator part of any generator is traceless and each jump costs the trace " +
        "r(|tr F|^2 - d*tr(F^dag F)) <= 0, so a generator with all rates non-negative has trace(L) = 0 iff " +
        "every jump with a nonzero rate is a multiple of I iff the dissipator vanishes iff the system is " +
        "closed; equivalently, open iff some eigenvalue lies off the imaginary axis. The palindrome is the INSTRUMENT that makes trace(L)/dim readable off a " +
        "measured spectrum, not a premise. That is OPEN, not EXTERNAL. What is NOT settled: (1) self-generated is only a constraint ([Π²,L]=0, " +
        "F63); (2) single-qubit decay is OPEN, its test scores a marginal over a coupled spectator and " +
        "returns the same verdict for external, internal and absent noise; (3) the bath's inheritance " +
        "from (2) falls with it, only the regress survives; (4) nothing (d=0) and (5) any other " +
        "dimension are eliminated, but definitionally, saying what cannot exist rather than clearing an " +
        "existing qubit. Dimension algebra typed as the parents (live: NoiseOriginExclusionWitness).";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return Polynomial;
            yield return Dimension;
            yield return new InspectableNode("the 5 candidates, and what each is worth now",
                summary: "(1) internal/self-generated: [Π²,L]=0 decouples the parity sectors (F63) — a CONSTRAINT, " +
                         "never an elimination. (2) single-qubit decay: OPEN. Its instrument reads the two-qubit " +
                         "marginal over a coupled spectator, whose pairing residual is 0.094 with NO noise, 0.177 " +
                         "with EXTERNAL dephasing and 0.116-0.149 with the internal mechanisms, so it separates " +
                         "nothing; the published 0/16 came from a seeded centre search that an exactly palindromic " +
                         "spectrum also fails. The decay does cost the neighbours purity (1.000 → 0.471/0.375/0.335). " +
                         "(3) many-qubit bath: the half that inherited (2) is struck; the infinite regress survives " +
                         "as a shift of the question. (4) nothing (d=0): ELIMINATED by definition (no properties). " +
                         "(5) any other dimension: ELIMINATED by the algebra, d²−2d=0 ⟹ d∈{0,2}, d=1 and d≥3 " +
                         "excluded. (4) and (5) are definitional: they say what cannot exist, not that an existing " +
                         "qubit is not the source.");
            yield return new InspectableNode("the conclusion: an OPEN system, and a question the framework cannot pose",
                summary: "the elimination does not reach 'external'. What is exact is the trace identity: centre = 0 " +
                         "iff the dissipator vanishes iff the system is closed, so a palindrome centred away from " +
                         "zero certifies a non-unitary generator. Between that and 'external' sits the reason the " +
                         "gap did not close: every internal source in the proof is written as a Lindblad " +
                         "dissipator, which assumes the environment, and the unitary alternative is a closed system " +
                         "whose palindrome sits at zero. The framework cannot express a candidate origin without " +
                         "already granting an outside; that is the incompleteness, sharper than the elimination was " +
                         "reaching for (INCOMPLETENESS_PROOF.md section 3).");
            yield return new InspectableNode("live witness",
                summary: "NoiseOriginExclusionWitness (inspect --root noise-origin) recomputes Candidate 5 live " +
                         "(d²−2d across d, roots {0,2}, d=1/d≥3 excluded) and enumerates all five candidates; the " +
                         "heavier Candidate-2 process-tomography compute is deferred in this light reading.");
        }
    }
}
