using System.Globalization;
using System.Linq;
using RCPsiSquared.Core.Inspection;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The live witness for the INCOMPLETENESS_PROOF argument and its 2026-08-29 standing. The
/// argument was that dephasing noise cannot originate WITHIN the d(d−2)=0 ontology, so it must come
/// from OUTSIDE; the elimination that carried it does not hold, and what remains is the trace
/// identity, which gives OPEN rather than EXTERNAL
/// (typed as <see cref="RCPsiSquared.Core.Symmetry.NoiseOriginExclusionClaim"/>). The
/// dimension conclusion d∈{0,2} is typed as the parents
/// (PolynomialFoundationClaim, QubitDimensionalAnchorClaim); this witness enumerates all five
/// internal candidates and recomputes the one that is elementary arithmetic (Candidate 5).
///
/// <para>This is the LIGHT reading: Candidate 5 (the framework's own algebra, d²−2d = 0 ⟹ d∈{0,2},
/// d=1 and d≥3 excluded) is RECOMPUTED live; Candidate 1 (the [Π², L] = 0 sector decoupling, a
/// constraint not an elimination) is carried by the already-typed F63; and Candidates 2-4 (single-
/// qubit decay, now OPEN, the many-qubit bath's infinite regress, and d=0's
/// property-lessness) are surfaced from the proof, their heavier process-tomography compute
/// deliberately deferred (see the proof + <c>simulations/failed_third.py</c>). Nothing is claimed as
/// recomputed that is not.</para>
///
/// <para>Anchors: <c>docs/proofs/INCOMPLETENESS_PROOF.md</c> +
/// <c>simulations/bootstrap_test.py</c> + <c>simulations/failed_third.py</c>.</para></summary>
public sealed class NoiseOriginExclusionWitness : IInspectable
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    public const int DefaultMaxD = 5;
    public int MaxD { get; }

    /// <summary>The minimum-memory polynomial d²−2d = d(d−2). Memory needs at least two values; this
    /// polynomial's roots {0, 2} are the only dimensions the palindromic mirror condition balances
    /// (d=0 the nothing-axis, d=2 the qubit). d=1 and d≥3 give a nonzero value and are excluded.</summary>
    public static int MinMemoryPolynomial(int d) => d * d - 2 * d;

    public sealed record DimensionRow(int D, int PolyValue)
    {
        public bool IsRoot => PolyValue == 0;
    }

    /// <summary>d²−2d evaluated across d = 0..MaxD.</summary>
    public IReadOnlyList<DimensionRow> Dimensions { get; }

    /// <summary>The dimensions the framework's algebra allows: the roots of d²−2d, i.e. {0, 2}.</summary>
    public IReadOnlyList<int> AllowedDimensions { get; }

    /// <summary>One of the five internal noise-origin candidates the proof rules out. <see cref="Recomputed"/>
    /// is true only for the candidate this light witness computes live (Candidate 5).</summary>
    public sealed record Candidate(int Index, string Name, string Verdict, bool Recomputed);

    public IReadOnlyList<Candidate> Candidates { get; }

    public NoiseOriginExclusionWitness(int maxD = DefaultMaxD)
    {
        if (maxD < 2)
            throw new ArgumentOutOfRangeException(nameof(maxD), maxD, "maxD must be ≥ 2 so the qubit root d=2 is in range.");
        MaxD = maxD;

        var dims = new List<DimensionRow>(maxD + 1);
        for (int d = 0; d <= maxD; d++) dims.Add(new DimensionRow(d, MinMemoryPolynomial(d)));
        Dimensions = dims;
        AllowedDimensions = dims.Where(r => r.IsRoot).Select(r => r.D).ToList();

        Candidates = new[]
        {
            new Candidate(1, "internal / self-generated",
                "CONSTRAINT, not elimination: [Π², L] = 0 exactly decouples the parity sectors (typed as F63), " +
                "so the bootstrap yields UNDERdetermination. It used to say the elimination was carried by " +
                "candidates 2-3; since 2026-08-29 neither carries it.",
                Recomputed: false),
            new Candidate(2, "single-qubit decay",
                "OPEN since 2026-08-29, previously ELIMINATED. Its test scores the two-qubit MARGINAL over a " +
                "spectator that is still coupled, and that marginal fails to pair for reasons unrelated to origin: " +
                "residual 0.094 with NO noise anywhere, 0.177 with EXTERNAL dephasing, 0.116-0.149 with the " +
                "internal mechanisms, and 2.2e-16 only when the spectator is fully decoupled. The published 0/16 " +
                "came from a seeded centre search that an exactly palindromic spectrum also fails, on a spectrum " +
                "shifted by a double division by d. Two further readings, gamma_eff = 0 and a non-Markovian " +
                "signature, were guard else-branches. The decay does cost the neighbours purity, " +
                "1.000 → 0.471/0.375/0.335. And the three-qubit generator itself is an EXACT palindrome at F137's " +
                "centre in every mechanism. See simulations/incompleteness_candidate2_evidence.py (29 gates).",
                Recomputed: false),
            new Candidate(3, "many-qubit bath",
                "HALF STANDING. The half that inherited candidate 2 (N instances of an eliminated source) falls " +
                "with it. The regress survives on its own: a bath that is to be the ORIGIN needs a reason for its " +
                "own dynamics, and both ways of giving it one move the question outward rather than answering it, " +
                "as a Lindblad dissipator (an environment assumed) or as a closed system (traceless, palindrome " +
                "at zero, no decay).",
                Recomputed: false),
            new Candidate(4, "nothing (d=0)",
                "ELIMINATED by definition: d=0 has no Hilbert space, no operators, no dynamics, no properties; " +
                "an entity with no properties cannot generate anything, noise included.",
                Recomputed: false),
            new Candidate(5, "something other than a qubit or nothing",
                "ELIMINATED by the framework's own algebra: the minimum-memory polynomial d²−2d = 0 has exactly the " +
                "roots {0, 2}; no other dimension balances the palindromic mirror condition (d=1 gives −1, d≥3 gives >0).",
                Recomputed: true),
        };
    }

    private static string Sign(int v) => v > 0 ? ">0" : v < 0 ? "<0" : "=0";

    public string DisplayName => $"NoiseOriginExclusionWitness (d²−2d over d=0..{MaxD}; the 5 noise-origin candidates and what each is worth)";

    public string Summary
    {
        get
        {
            string allowed = "{" + string.Join(", ", AllowedDimensions) + "}";
            int eliminated = Candidates.Count(c => c.Verdict.StartsWith("ELIMINATED"));
            int open = Candidates.Count(c => c.Verdict.StartsWith("OPEN"));
            return $"the minimum-memory polynomial d²−2d = 0 allows exactly d ∈ {allowed} (d=0 nothing, d=2 qubit); " +
                   $"d=1 and d≥3 are excluded, and that half is live arithmetic. The noise-origin elimination built " +
                   $"on top of it is NOT: of the five internal candidates, {eliminated} are eliminated and both of " +
                   $"those definitionally, {open} is open, one is a structural constraint and one keeps only its " +
                   $"regress, so 'the noise is external' does not follow. What does follow, exactly, is that a " +
                   $"completely positive generator is closed iff its trace vanishes, so an open system is " +
                   $"exactly one with an eigenvalue off the imaginary axis; where the spectrum pairs, its " +
                   $"centre trace(L)/dim makes that readable without the generator in hand. " +
                   $"Candidate 5 (the dimension algebra) is recomputed live here; candidate 1 is the typed F63 " +
                   $"[Π², L]=0; candidates 2-4 are surfaced from the proof.";
        }
    }

    public IEnumerable<IInspectable> Children
    {
        get
        {
            // The live dimension algebra (Candidate 5): d²−2d across d, roots {0,2}, others excluded.
            yield return new InspectableNode(
                displayName: $"the dimension algebra (Candidate 5, live): d²−2d = 0 ⟹ d ∈ {{{string.Join(", ", AllowedDimensions)}}}",
                summary: "d²−2d = d(d−2) across d: " +
                         string.Join("; ", Dimensions.Select(r =>
                             $"d={r.D}→{r.PolyValue} {Sign(r.PolyValue)}{(r.IsRoot ? " (allowed)" : " (excluded)")}")) +
                         ". Only d=0 (nothing) and d=2 (qubit) are roots; d=1 and d≥3 are excluded, so nothing " +
                         "'other than a qubit or nothing' can carry the palindromic mirror condition.",
                provenance: NodeProvenance.Live);

            // One node per candidate; only Candidate 5 is recomputed (Live), the rest surfaced (Stored).
            foreach (var c in Candidates)
                yield return new InspectableNode(
                    displayName: $"Candidate {c.Index}: {c.Name}{(c.Recomputed ? " (live)" : "")}",
                    summary: c.Verdict,
                    provenance: c.Recomputed ? NodeProvenance.Live : NodeProvenance.Stored);

            // The conclusion: OPEN by the trace, and the origin question left open.
            yield return new InspectableNode(
                displayName: "conclusion: an OPEN system, and an origin the formalism cannot ask about",
                summary: "the elimination does not reach OUTSIDE: (1) is a constraint, (2) and (3) have no evidence, " +
                         "and (4) and (5) say what cannot EXIST rather than clearing an existing qubit. What is exact " +
                         "is the trace: the commutator part contributes nothing and each jump costs " +
                         "r(|tr F|^2 - d*tr(F^dag F)) <= 0, so a generator with all rates non-negative is closed exactly " +
                         "when its trace vanishes, i.e. the system is open exactly when some eigenvalue lies off the " +
                         "imaginary axis. Where the spectrum pairs, its centre trace(L)/dim makes that readable " +
                         "without the generator in hand. Between OPEN and EXTERNAL sits the reason the gap does not " +
                         "close: an internal source can only be written here as a Lindblad dissipator, which already " +
                         "is a coupling to an environment (INCOMPLETENESS_PROOF.md section 3).",
                provenance: NodeProvenance.Live);
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
