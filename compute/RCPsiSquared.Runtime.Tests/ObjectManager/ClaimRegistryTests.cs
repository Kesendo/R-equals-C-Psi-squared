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
}
