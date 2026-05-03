using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Resonance;

namespace RCPsiSquared.Core.F86;

/// <summary>F86-specific query helpers on top of the generic <see cref="InspectionQuery"/>
/// layer. These are the named questions we actually ask of the F86 knowledge base —
/// "what's at tier X", "give me the witness at (c, N)", "compare measured to predicted".
///
/// <para>Each method is a typed query, not a tree-walk-and-string-match. Use these when
/// asking F86 questions; fall back to <c>kb.Walk()...</c> for ad-hoc exploration.</para>
/// </summary>
public static class F86Query
{
    /// <summary>All <see cref="F86Claim"/> objects in the tree (recursive).</summary>
    public static IEnumerable<F86Claim> AllClaims(this F86KnowledgeBase kb) =>
        kb.Walk().OfType<F86Claim>();

    /// <summary>All claims at a specific tier.</summary>
    public static IEnumerable<F86Claim> ClaimsAtTier(this F86KnowledgeBase kb, Tier tier) =>
        kb.AllClaims().Where(c => c.Tier == tier);

    /// <summary>Counts of claims grouped by tier — useful for the "tier inventory" summary.</summary>
    public static IReadOnlyDictionary<Tier, int> CountByTier(this F86KnowledgeBase kb) =>
        kb.AllClaims()
            .GroupBy(c => c.Tier)
            .ToDictionary(g => g.Key, g => g.Count());

    /// <summary>Unique anchor pointers referenced by any claim in this knowledge base.</summary>
    public static IReadOnlySet<string> AnchorsReferenced(this F86KnowledgeBase kb) =>
        kb.AllClaims().Select(c => c.Anchor).ToHashSet();

    /// <summary>The Interior universal-shape witness at (c, N), or null if not in the
    /// witness table.</summary>
    public static UniversalShapeWitness? InteriorWitnessAt(this F86KnowledgeBase kb, int c, int N) =>
        kb.InteriorShape.Witnesses.FirstOrDefault(w => w.Chromaticity == c && w.N == N);

    /// <summary>The Endpoint universal-shape witness at (c, N), or null.</summary>
    public static UniversalShapeWitness? EndpointWitnessAt(this F86KnowledgeBase kb, int c, int N) =>
        kb.EndpointShape.Witnesses.FirstOrDefault(w => w.Chromaticity == c && w.N == N);

    /// <summary>Both bond-class witnesses at (c, N).</summary>
    public static (UniversalShapeWitness? Interior, UniversalShapeWitness? Endpoint) WitnessesAt(
        this F86KnowledgeBase kb, int c, int N) =>
        (kb.InteriorWitnessAt(c, N), kb.EndpointWitnessAt(c, N));

    /// <summary>Per-block Q_peak (Q_SCALE) for chromaticity c, or null if not in the standard
    /// table.</summary>
    public static PerBlockQPeakClaim? PerBlockQPeak(this F86KnowledgeBase kb, int c) =>
        kb.PerBlockQPeaks.FirstOrDefault(p => p.Chromaticity == c);

    /// <summary>Per-bond Q_peak fine-grid value at (c, N, BondClass), or null if not in the
    /// witness table.</summary>
    public static PerBondQPeakWitness? PerBondQPeak(this F86KnowledgeBase kb, int c, int N, BondClass cls)
    {
        var table = cls == BondClass.Endpoint ? kb.EndpointPerBondTable : kb.InteriorPerBondTable;
        return table.Witnesses.FirstOrDefault(w => w.Chromaticity == c && w.N == N);
    }

    /// <summary>γ₀ extraction protocol: given measured peak J* on a chromaticity-c block,
    /// what is γ₀? Uses the Q_SCALE per-block Q_peak; null if the chromaticity has no
    /// per-block claim.</summary>
    public static double? ExtractGammaZero(this F86KnowledgeBase kb, double measuredJStar, int c) =>
        kb.PerBlockQPeak(c)?.ExtractGammaZero(measuredJStar);

    /// <summary>Compare a measured KCurve against the universal-shape predictions, returning
    /// per-class match results. Convenience wrapper over <see cref="F86KnowledgeBase.CompareTo"/>.
    /// </summary>
    public static IReadOnlyList<PredictionMatch> ComparePredictions(this F86KnowledgeBase kb,
        KCurve measured) => kb.CompareTo(measured);

    /// <summary>All retracted claims — the "what we've ruled out" view.</summary>
    public static IReadOnlyList<RetractedClaim> Retractions(this F86KnowledgeBase kb) =>
        kb.Retracted;

    /// <summary>All open theoretical items — the "what's still missing" view.</summary>
    public static IReadOnlyList<OpenQuestion> OpenItems(this F86KnowledgeBase kb) =>
        kb.OpenQuestions;

    /// <summary>One-line tier inventory string: "T1d=5, T1c=5, T2e=3, retracted=2, open=3".</summary>
    public static string TierInventoryLine(this F86KnowledgeBase kb)
    {
        var counts = kb.CountByTier();
        return string.Join(", ", new[]
        {
            (Tier.Tier1Derived, "T1d"),
            (Tier.Tier1Candidate, "T1c"),
            (Tier.Tier2Empirical, "T2e"),
            (Tier.Tier2Verified, "T2v"),
            (Tier.OpenQuestion, "open"),
            (Tier.Retracted, "retracted"),
        }.Where(t => counts.TryGetValue(t.Item1, out var count) && count > 0)
         .Select(t => $"{t.Item2}={counts[t.Item1]}"));
    }
}
