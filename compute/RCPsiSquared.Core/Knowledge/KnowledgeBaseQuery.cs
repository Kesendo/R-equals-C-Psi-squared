using RCPsiSquared.Core.Inspection;

namespace RCPsiSquared.Core.Knowledge;

/// <summary>Generic query layer over any <see cref="IInspectable"/> tree containing
/// <see cref="Claim"/> nodes. Use these methods on any F-theorem knowledge base (F86, F71,
/// future F-roots) for tier inventory, anchor enumeration, and tier-filtered claim lookup.
///
/// <para>F-theorem-specific helpers (e.g. <c>InteriorWitnessAt</c>, <c>ExtractGammaZero</c>)
/// continue to live next to their KB; this layer covers the cross-cutting tier/anchor
/// queries shared by every KB.</para>
/// </summary>
public static class KnowledgeBaseQuery
{
    /// <summary>All <see cref="Claim"/> objects in the tree (recursive depth-first).</summary>
    public static IEnumerable<Claim> AllClaims(this IInspectable root) =>
        root.Walk().OfType<Claim>();

    /// <summary>All claims at a specific <see cref="Tier"/>.</summary>
    public static IEnumerable<Claim> ClaimsAtTier(this IInspectable root, Tier tier) =>
        root.AllClaims().Where(c => c.Tier == tier);

    /// <summary>Counts of claims grouped by tier; useful for tier-inventory summaries.</summary>
    public static IReadOnlyDictionary<Tier, int> CountByTier(this IInspectable root) =>
        root.AllClaims().GroupBy(c => c.Tier).ToDictionary(g => g.Key, g => g.Count());

    /// <summary>Unique anchor pointers referenced by any claim in this tree.</summary>
    public static IReadOnlySet<string> AnchorsReferenced(this IInspectable root) =>
        root.AllClaims().Select(c => c.Anchor).ToHashSet();

    /// <summary>One-line tier inventory string, e.g. <c>"T1d=5, T1c=5, T2e=3, retracted=2, open=3"</c>.
    /// Empty tiers are skipped.</summary>
    public static string TierInventoryLine(this IInspectable root)
    {
        var counts = root.CountByTier();
        return string.Join(", ", TierShortCodes
            .Where(t => counts.TryGetValue(t.Tier, out var count) && count > 0)
            .Select(t => $"{t.Code}={counts[t.Tier]}"));
    }

    private static readonly (Tier Tier, string Code)[] TierShortCodes =
    {
        (Tier.Tier1Derived, "T1d"),
        (Tier.Tier1Candidate, "T1c"),
        (Tier.Tier2Empirical, "T2e"),
        (Tier.Tier2Verified, "T2v"),
        (Tier.OpenQuestion, "open"),
        (Tier.Retracted, "retracted"),
    };
}
