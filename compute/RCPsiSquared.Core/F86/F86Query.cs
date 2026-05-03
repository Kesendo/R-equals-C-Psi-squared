using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Resonance;

namespace RCPsiSquared.Core.F86;

/// <summary>F86-specific query helpers: typed lookups into the F86 knowledge base
/// (witnesses at (c, N), γ₀ extraction, measured-vs-predicted comparison, etc.).
///
/// <para>Cross-cutting tier/anchor queries (<c>AllClaims</c>, <c>ClaimsAtTier</c>,
/// <c>CountByTier</c>, <c>AnchorsReferenced</c>, <c>TierInventoryLine</c>) live in
/// <see cref="KnowledgeBaseQuery"/> and apply to every F-theorem KB.</para>
/// </summary>
public static class F86Query
{
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
}
