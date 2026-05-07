using RCPsiSquared.Core.F71;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Runtime.F71Family;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.Tests.F71Family;

public class F71FamilyRegistrationTests
{
    [Fact]
    public void RegisterF71Family_BuildsFourClaims()
    {
        var registry = new ClaimRegistryBuilder()
            .RegisterF71Family(N: 5)
            .Build();

        Assert.Equal(4, registry.All().Count());
        Assert.True(registry.Contains<C1MirrorIdentity>());
        Assert.True(registry.Contains<F71MirrorOperator>());
        Assert.True(registry.Contains<F71BondOrbitDecomposition>());
        Assert.True(registry.Contains<F86MirrorGeneralisationLink>());
    }

    [Fact]
    public void RegisterF71Family_TopologicalOrder_IdentityFirst()
    {
        var registry = new ClaimRegistryBuilder()
            .RegisterF71Family(N: 5)
            .Build();

        var order = registry.TopologicalOrder.ToList();
        Assert.Equal(typeof(C1MirrorIdentity), order[0]);
    }

    [Fact]
    public void RegisterF71Family_AllTier1Derived()
    {
        var registry = new ClaimRegistryBuilder()
            .RegisterF71Family(N: 5)
            .Build();

        Assert.All(registry.All(), c => Assert.Equal(Tier.Tier1Derived, c.Tier));
    }

    [Fact]
    public void RegisterF71Family_DescendantsOfIdentity_ReturnsThree()
    {
        var registry = new ClaimRegistryBuilder()
            .RegisterF71Family(N: 5)
            .Build();

        var descendants = registry.DescendantsOf<C1MirrorIdentity>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Equal(3, descendants.Count);
        Assert.Contains(typeof(F71MirrorOperator), descendants);
        Assert.Contains(typeof(F71BondOrbitDecomposition), descendants);
        Assert.Contains(typeof(F86MirrorGeneralisationLink), descendants);
    }
}
