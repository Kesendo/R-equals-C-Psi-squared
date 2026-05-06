using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Orchestration.Cli;
using RCPsiSquared.Runtime.ObjectManager;
using Knowledge = RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Orchestration.Tests.Cli;

public class KnowledgeCliTests
{
    private sealed class Foo : Claim
    {
        public Foo() : base("Foo", Tier.Tier1Derived,
            "compute/RCPsiSquared.Runtime.Tests/TestSupport/TestClaims.cs") { }
        public override string DisplayName => "Foo";
        public override string Summary => "synthetic Foo";
    }

    private sealed class Bar : Claim
    {
        public Bar(Foo parent) : base("Bar", Tier.Tier1Derived,
            "compute/RCPsiSquared.Runtime.Tests/TestSupport/TestClaims.cs") { _ = parent; }
        public override string DisplayName => "Bar";
        public override string Summary => "synthetic Bar (child of Foo)";
    }

    [Fact]
    public void Cli_TierQuery_EmptyRegistry_RendersEmpty()
    {
        var registry = new ClaimRegistryBuilder().Build();
        var cli = new KnowledgeCli(registry);
        var output = cli.Render(new KnowledgeQuery.Tier(Knowledge.Tier.Tier1Derived));
        Assert.Contains("(empty)", output);
    }

    [Fact]
    public void Cli_TierQuery_WithClaim_RendersTable()
    {
        var registry = new ClaimRegistryBuilder()
            .Register<Foo>(_ => new Foo())
            .Build();
        var cli = new KnowledgeCli(registry);
        var output = cli.Render(new KnowledgeQuery.Tier(Knowledge.Tier.Tier1Derived));
        Assert.Contains("Foo", output);
    }

    [Fact]
    public void Cli_AncestorsQuery_NoEdges_RendersEmpty()
    {
        var registry = new ClaimRegistryBuilder()
            .Register<Foo>(_ => new Foo())
            .Build();
        var cli = new KnowledgeCli(registry);
        var output = cli.Render(new KnowledgeQuery.Ancestors(typeof(Foo)));
        Assert.Contains("(no ancestors)", output);
    }

    [Fact]
    public void Cli_DescendantsQuery_OneEdge_RendersChild()
    {
        var registry = new ClaimRegistryBuilder()
            .Register<Foo>(_ => new Foo())
            .Register<Bar>(b => new Bar(b.Get<Foo>()))
            .Build();
        var cli = new KnowledgeCli(registry);
        var output = cli.Render(new KnowledgeQuery.Descendants(typeof(Foo)));
        Assert.Contains("Bar", output);
    }
}
