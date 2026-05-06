using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.Tests.TestSupport;

namespace RCPsiSquared.Runtime.Tests.ObjectManager;

public class ClaimRegistryBuilderTests
{
    private sealed class Child : Claim
    {
        public FooT1D Parent { get; }
        public Child(FooT1D parent)
            : base("synthetic Child", Tier.Tier1Derived,
                   "compute/RCPsiSquared.Runtime.Tests/TestSupport/TestClaims.cs")
        { Parent = parent; }
        public override string DisplayName => "Child";
        public override string Summary => "synthetic child of Foo";
    }

    [Fact]
    public void Build_SingleClaim_NoEdges_Succeeds()
    {
        var registry = new ClaimRegistryBuilder()
            .Register<FooT1D>(_ => new FooT1D())
            .Build();
        Assert.Single(registry.All());
    }

    [Fact]
    public void Build_TwoClaims_OneEdge_TopologicalOrderParentFirst()
    {
        var registry = new ClaimRegistryBuilder()
            .Register<Child>(b => new Child(b.Get<FooT1D>()))
            .Register<FooT1D>(_ => new FooT1D())
            .Build();

        Assert.Equal(2, registry.All().Count());
        Assert.Equal(typeof(FooT1D), registry.TopologicalOrder[0]);
        Assert.Equal(typeof(Child), registry.TopologicalOrder[1]);

        var edges = registry.AllEdges().ToList();
        Assert.Single(edges);
        Assert.Equal(typeof(FooT1D), edges[0].Parent);
        Assert.Equal(typeof(Child), edges[0].Child);
    }

    [Fact]
    public void Build_EdgeReason_DefaultIsClaimNames()
    {
        var registry = new ClaimRegistryBuilder()
            .Register<Child>(b => new Child(b.Get<FooT1D>()))
            .Register<FooT1D>(_ => new FooT1D())
            .Build();

        var edge = registry.AllEdges().Single();
        Assert.Contains("FooT1D", edge.Reason);
        Assert.Contains("Child", edge.Reason);
    }
}
