using System.Text;
using KnowledgeTier = RCPsiSquared.Core.Knowledge.Tier;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Orchestration.Cli;

/// <summary>CLI inspector consumer. Renders deterministic ASCII output for the supported
/// <see cref="KnowledgeQuery"/> shapes. Output is keyed and memoised; the registry is
/// immutable, so the rendered string is cacheable as long as both registry and query are
/// reference-equal.</summary>
public sealed class KnowledgeCli
{
    private readonly ClaimRegistry _registry;
    private readonly Dictionary<KnowledgeQuery, string> _cache = new();

    public KnowledgeCli(ClaimRegistry registry)
    {
        _registry = registry;
    }

    public string Render(KnowledgeQuery query)
    {
        if (_cache.TryGetValue(query, out var cached)) return cached;

        var output = query switch
        {
            KnowledgeQuery.Tier t => RenderTier(t.Level),
            KnowledgeQuery.Ancestors a => RenderAncestors(a.Of),
            KnowledgeQuery.Descendants d => RenderDescendants(d.Of),
            KnowledgeQuery.All => RenderAll(),
            _ => throw new ArgumentOutOfRangeException(nameof(query)),
        };
        _cache[query] = output;
        return output;
    }

    private string RenderTier(KnowledgeTier level) => RenderList(
        header: $"Tier {level.Label()}",
        emptyText: "(empty)",
        claims: _registry.AllOfTier(level),
        formatLine: c => $"  - {c.GetType().Name}: {c.DisplayName}");

    private string RenderAncestors(Type childType) => RenderList(
        header: $"Ancestors of {childType.Name}",
        emptyText: "(no ancestors)",
        claims: _registry.AncestorsOf(childType),
        formatLine: c => $"  - {c.GetType().Name}: {c.DisplayName} [{c.Tier.Label()}]");

    private string RenderDescendants(Type parentType) => RenderList(
        header: $"Descendants of {parentType.Name}",
        emptyText: "(no descendants)",
        claims: _registry.DescendantsOf(parentType),
        formatLine: c => $"  - {c.GetType().Name}: {c.DisplayName} [{c.Tier.Label()}]");

    private string RenderAll() => RenderList(
        header: "Registry",
        emptyText: "(empty)",
        claims: _registry.All().ToList(),
        formatLine: c => $"  - {c.GetType().Name} [{c.Tier.Label()}]: {c.DisplayName}");

    private static string RenderList(
        string header,
        string emptyText,
        IReadOnlyList<Claim> claims,
        Func<Claim, string> formatLine)
    {
        if (claims.Count == 0)
            return $"{header}: {emptyText}\n";

        var sb = new StringBuilder();
        sb.AppendLine($"{header}: {claims.Count} Claim(s)");
        foreach (var c in claims)
            sb.AppendLine(formatLine(c));
        return sb.ToString();
    }
}
