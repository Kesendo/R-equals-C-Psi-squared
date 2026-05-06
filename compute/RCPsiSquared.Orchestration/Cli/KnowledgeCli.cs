using System.Text;
using KnowledgeTier = RCPsiSquared.Core.Knowledge.Tier;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Orchestration.Cli;

/// <summary>Layer 3 consumer 1 (CLI inspector). Renders deterministic ASCII output for the
/// supported <see cref="KnowledgeQuery"/> shapes. Output is keyed and memoised; the registry
/// is immutable, so the rendered string is cacheable as long as both registry and query are
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

    private string RenderTier(KnowledgeTier level)
    {
        var claims = _registry.AllOfTier(level);
        if (claims.Count == 0)
            return $"Tier {level.Label()}: (empty)\n";
        var sb = new StringBuilder();
        sb.AppendLine($"Tier {level.Label()}: {claims.Count} Claim(s)");
        foreach (var c in claims)
            sb.AppendLine($"  - {c.GetType().Name}: {c.DisplayName}");
        return sb.ToString();
    }

    private string RenderAncestors(Type childType)
    {
        var ancestors = _registry.AncestorsOf(childType);
        if (ancestors.Count == 0)
            return $"Ancestors of {childType.Name}: (no ancestors)\n";
        var sb = new StringBuilder();
        sb.AppendLine($"Ancestors of {childType.Name}:");
        foreach (var c in ancestors)
            sb.AppendLine($"  - {c.GetType().Name}: {c.DisplayName} [{c.Tier.Label()}]");
        return sb.ToString();
    }

    private string RenderDescendants(Type parentType)
    {
        var descendants = _registry.DescendantsOf(parentType);
        if (descendants.Count == 0)
            return $"Descendants of {parentType.Name}: (no descendants)\n";
        var sb = new StringBuilder();
        sb.AppendLine($"Descendants of {parentType.Name}:");
        foreach (var c in descendants)
            sb.AppendLine($"  - {c.GetType().Name}: {c.DisplayName} [{c.Tier.Label()}]");
        return sb.ToString();
    }

    private string RenderAll()
    {
        var sb = new StringBuilder();
        sb.AppendLine($"Registry: {_registry.All().Count()} Claim(s)");
        foreach (var c in _registry.All())
            sb.AppendLine($"  - {c.GetType().Name} [{c.Tier.Label()}]: {c.DisplayName}");
        return sb.ToString();
    }
}
