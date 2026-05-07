using RCPsiSquared.Core.F86.Item1Derivation;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Runtime.F86Main;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.Tests.F86Main;

public class F86Item1RegistrationTests
{
    [Fact]
    public void RegisterF86Item1Light_BuildsTwoClaims()
    {
        var registry = new ClaimRegistryBuilder()
            .RegisterF86Item1Light(N: 5, n: 1, gammaZero: 0.05)
            .Build();

        Assert.Equal(2, registry.All().Count());
        Assert.True(registry.Contains<C2BlockShape>());
        Assert.True(registry.Contains<C2ChannelUniformAnalytical>());
    }

    [Fact]
    public void RegisterF86Item1Light_ChannelUniformDescendsFromBlockShape()
    {
        var registry = new ClaimRegistryBuilder()
            .RegisterF86Item1Light(N: 5, n: 1, gammaZero: 0.05)
            .Build();

        var ancestors = registry.AncestorsOf<C2ChannelUniformAnalytical>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(C2BlockShape), ancestors);
    }

    [Fact]
    public void RegisterF86Item1Light_BothTier1Derived()
    {
        var registry = new ClaimRegistryBuilder()
            .RegisterF86Item1Light(N: 5, n: 1, gammaZero: 0.05)
            .Build();

        Assert.Equal(Tier.Tier1Derived, registry.Get<C2BlockShape>().Tier);
        Assert.Equal(Tier.Tier1Derived, registry.Get<C2ChannelUniformAnalytical>().Tier);
    }

    [Fact]
    public void RegisterF86Item1Light_NonC2Block_Throws()
    {
        // c=3 block (N=5, n=2) is rejected; Item 1 is c=2 only.
        Assert.Throws<ArgumentException>(() =>
            new ClaimRegistryBuilder()
                .RegisterF86Item1Light(N: 5, n: 2, gammaZero: 0.05)
                .Build());
    }
}
