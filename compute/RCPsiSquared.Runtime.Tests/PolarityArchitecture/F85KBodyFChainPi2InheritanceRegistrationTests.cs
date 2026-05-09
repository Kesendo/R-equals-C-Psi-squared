using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.F1Family;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class F85KBodyFChainPi2InheritanceRegistrationTests
{
    private static ChainSystem DefaultChain() =>
        new ChainSystem(N: 5, J: 1.0, GammaZero: 0.05,
            HType: HamiltonianType.Heisenberg, Topology: TopologyKind.Chain);

    private static ClaimRegistryBuilder BuildBaseRegistry() =>
        new ClaimRegistryBuilder()
            .RegisterF1Family(DefaultChain())
            .RegisterPi2Family()
            .RegisterPi2DyadicLadder()
            .RegisterF88PopcountCoherence()
            .RegisterF88StaticDyadicAnchor()
            .RegisterPi2OperatorSpaceMirror()
            .RegisterF49Pi2Inheritance()
            .RegisterF83AntiFractionPi2Inheritance();

    [Fact]
    public void RegisterF85_AddsClaim()
    {
        var registry = BuildBaseRegistry()
            .RegisterF85KBodyFChainPi2Inheritance()
            .Build();

        Assert.True(registry.Contains<F85KBodyFChainPi2Inheritance>());
    }

    [Fact]
    public void RegisterF85_TierIsTier1Derived()
    {
        var registry = BuildBaseRegistry()
            .RegisterF85KBodyFChainPi2Inheritance()
            .Build();

        Assert.Equal(Tier.Tier1Derived,
            registry.Get<F85KBodyFChainPi2Inheritance>().Tier);
    }

    [Fact]
    public void RegisterF85_AncestorsContainF49AndF83()
    {
        var registry = BuildBaseRegistry()
            .RegisterF85KBodyFChainPi2Inheritance()
            .Build();

        var ancestors = registry.AncestorsOf<F85KBodyFChainPi2Inheritance>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(F49Pi2Inheritance), ancestors);
        Assert.Contains(typeof(F83AntiFractionPi2Inheritance), ancestors);
    }

    [Theory]
    [InlineData(2, 4)]
    [InlineData(3, 14)]
    [InlineData(4, 40)]
    public void RegisterF85_Pi2OddCountTableAcrossRegistry(int k, long expected)
    {
        var registry = BuildBaseRegistry()
            .RegisterF85KBodyFChainPi2Inheritance()
            .Build();

        Assert.Equal(expected, registry.Get<F85KBodyFChainPi2Inheritance>().Pi2OddCount(k));
    }

    [Fact]
    public void RegisterF85_MatchesF83CoefficientsAcrossRegistry()
    {
        var registry = BuildBaseRegistry()
            .RegisterF85KBodyFChainPi2Inheritance()
            .Build();

        Assert.True(registry.Get<F85KBodyFChainPi2Inheritance>().MatchesF83Coefficients());
    }
}
