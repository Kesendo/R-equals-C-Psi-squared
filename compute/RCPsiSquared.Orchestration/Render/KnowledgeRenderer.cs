using System.Text;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Orchestration.Render;

/// <summary>Layer-3 consumer that renders a <see cref="ClaimRegistry"/> as Markdown,
/// grouped by <see cref="Tier"/>. Closes Pain C from the original Schicht 0+1 brainstorm:
/// auto-generation of human-readable summaries of the typed-knowledge graph from the live
/// registry, instead of hand-maintained markdown that drifts.
///
/// <para>Output structure:</para>
/// <list type="bullet">
///   <item>One H1 header for the registry as a whole, with totals.</item>
///   <item>One H2 section per non-empty Tier (in strength order, strongest first).</item>
///   <item>One bullet per claim: <c>**TypeName**: DisplayName</c> with the anchor in
///         italics on the next line.</item>
///   <item>An H2 "Edges" section summarising the inheritance edge count, plus the first
///         few edges as bullets.</item>
/// </list>
///
/// <para>The output is deterministic and idempotent. Rendering the same registry twice
/// produces byte-identical Markdown. Cache: keyed by registry reference; the registry is
/// immutable.</para></summary>
public sealed class KnowledgeRenderer
{
    private readonly ClaimRegistry _registry;
    private string? _cached;

    public KnowledgeRenderer(ClaimRegistry registry)
    {
        _registry = registry;
    }

    public string Render()
    {
        if (_cached is not null) return _cached;

        var sb = new StringBuilder();
        var allClaims = _registry.All().ToList();
        var allEdges = _registry.AllEdges().ToList();

        sb.AppendLine($"# Knowledge Registry");
        sb.AppendLine();
        sb.AppendLine($"{allClaims.Count} claim(s), {allEdges.Count} edge(s).");
        sb.AppendLine();

        foreach (var tier in TierStrength.AllByStrength)
        {
            var claimsAtTier = _registry.AllOfTier(tier);
            if (claimsAtTier.Count == 0) continue;

            sb.AppendLine($"## {tier.Label()}");
            sb.AppendLine();
            sb.AppendLine($"{claimsAtTier.Count} claim(s).");
            sb.AppendLine();

            foreach (var claim in claimsAtTier)
            {
                sb.AppendLine($"- **{claim.GetType().Name}**: {claim.DisplayName}");
                if (!string.IsNullOrWhiteSpace(claim.Anchor))
                    sb.AppendLine($"  - Anchor: {claim.Anchor}");
            }
            sb.AppendLine();
        }

        if (allEdges.Count > 0)
        {
            sb.AppendLine($"## Edges");
            sb.AppendLine();
            sb.AppendLine($"{allEdges.Count} edge(s) total.");
            sb.AppendLine();
            foreach (var edge in allEdges)
                sb.AppendLine($"- {edge.Parent.Name} → {edge.Child.Name}");
            sb.AppendLine();
        }

        _cached = sb.ToString();
        return _cached;
    }
}
