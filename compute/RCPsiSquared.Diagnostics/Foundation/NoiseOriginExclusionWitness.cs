using System.Globalization;
using System.Linq;
using RCPsiSquared.Core.Inspection;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The live witness for the INCOMPLETENESS_PROOF argument (2026-07-02), the noise-origin
/// 5-candidate elimination: dephasing noise cannot originate WITHIN the d(d−2)=0 ontology, so it
/// must come from OUTSIDE (typed as <see cref="RCPsiSquared.Core.Symmetry.NoiseOriginExclusionClaim"/>). The
/// dimension conclusion d∈{0,2} is typed as the parents
/// (PolynomialFoundationClaim, QubitDimensionalAnchorClaim); this witness enumerates all five
/// internal candidates and recomputes the one that is elementary arithmetic (Candidate 5).
///
/// <para>This is the LIGHT reading: Candidate 5 (the framework's own algebra, d²−2d = 0 ⟹ d∈{0,2},
/// d=1 and d≥3 excluded) is RECOMPUTED live; Candidate 1 (the [Π², L] = 0 sector decoupling, a
/// constraint not an elimination) is carried by the already-typed F63; and Candidates 2-4 (single-
/// qubit decay γ_eff = 0 / 0-of-16 palindromic, the many-qubit bath's infinite regress, and d=0's
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
                "so the bootstrap yields UNDERdetermination; the actual elimination is carried by candidates 2-3.",
                Recomputed: false),
            new Candidate(2, "single-qubit decay",
                "ELIMINATED: γ_eff = 0 for all four mechanisms; process tomography gives 0/16 palindromic pairs " +
                "(reference standard Z-dephasing gives 16/16). Heavier compute (3-qubit evolution + partial trace) " +
                "deferred in this light witness; see the proof + simulations/failed_third.py.",
                Recomputed: false),
            new Candidate(3, "many-qubit bath",
                "ELIMINATED: a bath of N qubits is N instances of candidate 2; if one finite qubit cannot generate " +
                "the required noise, a collection cannot either, and the bath would need its own source → infinite regress.",
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

    public string DisplayName => $"NoiseOriginExclusionWitness (d²−2d over d=0..{MaxD}, the 5-candidate elimination)";

    public string Summary
    {
        get
        {
            string allowed = "{" + string.Join(", ", AllowedDimensions) + "}";
            int eliminated = Candidates.Count(c => c.Verdict.StartsWith("ELIMINATED"));
            return $"the minimum-memory polynomial d²−2d = 0 allows exactly d ∈ {allowed} (d=0 nothing, d=2 qubit); " +
                   $"d=1 and d≥3 are excluded. All five internal noise-origin candidates fail inside the d(d−2)=0 " +
                   $"ontology ({eliminated} eliminated, 1 a constraint), so dephasing noise must originate OUTSIDE it. " +
                   $"Candidate 5 (the dimension algebra) is recomputed live here; candidate 1 is the typed F63 " +
                   $"[Π², L]=0; candidates 2-4 are surfaced from the proof (heavier compute deferred).";
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

            // The conclusion: every internal origin fails ⟹ the noise comes from outside.
            yield return new InspectableNode(
                displayName: "conclusion: the noise comes from OUTSIDE d(d−2)=0",
                summary: "every candidate internal to the d(d−2)=0 ontology is either eliminated (2, 3, 4, 5) or only " +
                         "a constraint (1); no internal source generates the required dephasing. Therefore the noise " +
                         "originates OUTSIDE the system's own algebra — the incompleteness the proof names " +
                         "(INCOMPLETENESS_PROOF.md). The V-Effect one dimension up.",
                provenance: NodeProvenance.Live);
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
