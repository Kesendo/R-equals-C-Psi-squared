using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.Tests.TestSupport;

namespace RCPsiSquared.Runtime.Tests.ObjectManager;

public class ClaimRegistryTests
{
    [Fact]
    public void Build_EmptyRegistry_ReturnsZeroClaims()
    {
        var registry = new ClaimRegistryBuilder().Build();
        Assert.Empty(registry.All());
    }

    [Fact]
    public void Build_SingleClaim_NoEdges_Succeeds()
    {
        var registry = new ClaimRegistryBuilder()
            .Register<FooT1D>(_ => new FooT1D())
            .Build();

        Assert.True(registry.Contains<FooT1D>());
        Assert.IsType<FooT1D>(registry.Get<FooT1D>());
        Assert.Single(registry.All());
        Assert.Empty(registry.AllEdges());
    }

    [Fact]
    public void Get_NotRegistered_Throws()
    {
        var registry = new ClaimRegistryBuilder().Build();
        Assert.Throws<KeyNotFoundException>((Action)(() => registry.Get<FooT1D>()));
    }

    [Fact]
    public void TryGet_NotRegistered_ReturnsFalse()
    {
        var registry = new ClaimRegistryBuilder().Build();
        Assert.False(registry.TryGet<FooT1D>(out var _));
    }

    [Fact]
    public void AllOfTier_FiltersByTier()
    {
        var registry = new ClaimRegistryBuilder()
            .Register<FooT1D>(_ => new FooT1D())
            .Register<BarT2E>(_ => new BarT2E())
            .Build();

        Assert.Single(registry.AllOfTier(Tier.Tier1Derived));
        Assert.Single(registry.AllOfTier(Tier.Tier2Empirical));
        Assert.Empty(registry.AllOfTier(Tier.Retracted));
    }

    private sealed class GrandChild : Claim
    {
        public GrandChild(Child child) : base("GrandChild", Tier.Tier1Derived,
            "compute/RCPsiSquared.Runtime.Tests/TestSupport/TestClaims.cs")
        { _ = child; }
        public override string DisplayName => "GrandChild";
        public override string Summary => "synthetic grandchild";
    }

    private sealed class Child : Claim
    {
        public Child(FooT1D parent) : base("Child", Tier.Tier1Derived,
            "compute/RCPsiSquared.Runtime.Tests/TestSupport/TestClaims.cs")
        { _ = parent; }
        public override string DisplayName => "Child";
        public override string Summary => "synthetic child";
    }

    [Fact]
    public void AncestorsOf_TransitiveClosure()
    {
        var registry = new ClaimRegistryBuilder()
            .Register<FooT1D>(_ => new FooT1D())
            .Register<Child>(b => new Child(b.Get<FooT1D>()))
            .Register<GrandChild>(b => new GrandChild(b.Get<Child>()))
            .Build();

        var ancestors = registry.AncestorsOf<GrandChild>().Select(c => c.GetType()).ToHashSet();
        Assert.Contains(typeof(Child), ancestors);
        Assert.Contains(typeof(FooT1D), ancestors);
        Assert.Equal(2, ancestors.Count);
    }

    [Fact]
    public void DescendantsOf_TransitiveClosure()
    {
        var registry = new ClaimRegistryBuilder()
            .Register<FooT1D>(_ => new FooT1D())
            .Register<Child>(b => new Child(b.Get<FooT1D>()))
            .Register<GrandChild>(b => new GrandChild(b.Get<Child>()))
            .Build();

        var descendants = registry.DescendantsOf<FooT1D>().Select(c => c.GetType()).ToHashSet();
        Assert.Contains(typeof(Child), descendants);
        Assert.Contains(typeof(GrandChild), descendants);
        Assert.Equal(2, descendants.Count);
    }

    [Fact]
    public void EdgesFrom_DirectChildren()
    {
        var registry = new ClaimRegistryBuilder()
            .Register<FooT1D>(_ => new FooT1D())
            .Register<Child>(b => new Child(b.Get<FooT1D>()))
            .Build();

        var edgesFromFoo = registry.EdgesFrom<FooT1D>();
        Assert.Single(edgesFromFoo);
        Assert.Equal(typeof(Child), edgesFromFoo[0].Child);
    }

    [Fact]
    public void EdgesInto_DirectParents()
    {
        var registry = new ClaimRegistryBuilder()
            .Register<FooT1D>(_ => new FooT1D())
            .Register<Child>(b => new Child(b.Get<FooT1D>()))
            .Build();

        var edgesIntoChild = registry.EdgesInto<Child>();
        Assert.Single(edgesIntoChild);
        Assert.Equal(typeof(FooT1D), edgesIntoChild[0].Parent);
    }
}
