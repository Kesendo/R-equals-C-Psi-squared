using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Diagnostics.Knowledge;

namespace RCPsiSquared.Diagnostics.Tests.Knowledge;

public class ClaimRegistryInspectableNodeTests
{
    [Fact]
    public void Build_OneGroupPerNonEmptyTier_ChildrenAreTheClaims()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var node = ClaimRegistryInspectableNode.Build(registry);

        var groups = node.Children.ToList();
        var nonEmptyTiers = Enum.GetValues<Tier>()
            .Count(t => registry.AllOfTier(t).Count > 0);
        Assert.Equal(nonEmptyTiers, groups.Count);

        // Every claim in the registry surfaces under exactly one tier group.
        int leafCount = groups.Sum(g => g.Children.Count());
        Assert.Equal(registry.Count, leafCount);
    }

    [Fact]
    public void Build_EachGroupHoldsOnlyItsTiersClaims()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var node = ClaimRegistryInspectableNode.Build(registry);

        foreach (var tier in Enum.GetValues<Tier>())
        {
            var expected = registry.AllOfTier(tier);
            if (expected.Count == 0) continue;
            var group = node.Children.Single(g => g.DisplayName.StartsWith(tier.Label()));
            var claims = group.Children.Cast<Claim>().ToList();
            Assert.All(claims, c => Assert.Equal(tier, c.Tier));
            Assert.Equal(expected.Count, claims.Count);
        }
    }

    [Fact]
    public void Build_RootSummaryReportsTotalClaimCount()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var node = ClaimRegistryInspectableNode.Build(registry);
        Assert.Contains(registry.Count.ToString(), node.Summary);
    }
}
