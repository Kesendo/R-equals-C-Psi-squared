using KnowledgeTier = RCPsiSquared.Core.Knowledge.Tier;

namespace RCPsiSquared.Orchestration.Cli;

/// <summary>Queries the <see cref="KnowledgeCli"/> can answer. Each subtype is one
/// concrete query shape so the CLI can switch over them exhaustively.</summary>
public abstract record KnowledgeQuery
{
    public sealed record Tier(KnowledgeTier Level) : KnowledgeQuery;
    public sealed record Ancestors(Type Of) : KnowledgeQuery;
    public sealed record Descendants(Type Of) : KnowledgeQuery;
    public sealed record All() : KnowledgeQuery;
}
