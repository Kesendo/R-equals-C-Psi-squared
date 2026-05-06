using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.F1;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Runtime.F1Family;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.Tests.F1Family;

public class F1FamilyRegistrationTests
{
    private static ChainSystem DefaultChain(int N = 5) =>
        new ChainSystem(N: N, J: 1.0, GammaZero: 0.05,
            HType: HamiltonianType.XY, Topology: TopologyKind.Chain);

    [Fact]
    public void RegisterF1Family_BuildsThreeClaims()
    {
        var registry = new ClaimRegistryBuilder()
            .RegisterF1Family(DefaultChain())
            .Build();

        Assert.Equal(3, registry.All().Count());
        Assert.True(registry.Contains<ChainSystemPrimitive>());
        Assert.True(registry.Contains<F1PalindromeIdentity>());
        Assert.True(registry.Contains<PalindromeResidualScalingClaim>());
    }

    [Fact]
    public void RegisterF1Family_TopologicalOrder_PrimitiveFirst()
    {
        var registry = new ClaimRegistryBuilder()
            .RegisterF1Family(DefaultChain())
            .Build();

        var firstIndex = registry.TopologicalOrder.ToList().IndexOf(typeof(ChainSystemPrimitive));
        var f1Index = registry.TopologicalOrder.ToList().IndexOf(typeof(F1PalindromeIdentity));
        var f73Index = registry.TopologicalOrder.ToList().IndexOf(typeof(PalindromeResidualScalingClaim));

        Assert.True(firstIndex < f1Index);
        Assert.True(f1Index < f73Index);
    }

    [Fact]
    public void RegisterF1Family_Cli_AncestorsOfF73_ContainsF1AndChain()
    {
        var registry = new ClaimRegistryBuilder()
            .RegisterF1Family(DefaultChain())
            .Build();

        var ancestors = registry.AncestorsOf<PalindromeResidualScalingClaim>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(F1PalindromeIdentity), ancestors);
        Assert.Contains(typeof(ChainSystemPrimitive), ancestors);
    }
}
