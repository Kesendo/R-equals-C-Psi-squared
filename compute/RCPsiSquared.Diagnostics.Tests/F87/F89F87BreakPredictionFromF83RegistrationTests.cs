using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Diagnostics.F87;
using RCPsiSquared.Runtime.F1Family;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Diagnostics.Tests.F87;

public class F89F87BreakPredictionFromF83RegistrationTests
{
    private static ChainSystem DefaultChain() =>
        new ChainSystem(N: 5, J: 1.0, GammaZero: 0.05,
            HType: HamiltonianType.Heisenberg, Topology: TopologyKind.Chain);

    private static ClaimRegistryBuilder BuildBaseRegistry() =>
        new ClaimRegistryBuilder()
            .RegisterF1Family(DefaultChain())
            .RegisterPi2Family()
            .RegisterF87Family()
            .RegisterF89F87TrulyInheritance();

    [Fact]
    public void RegisterF89F87BreakPredictionFromF83_AddsClaim()
    {
        var registry = BuildBaseRegistry()
            .RegisterF89F87BreakPredictionFromF83()
            .Build();

        Assert.True(registry.Contains<F89F87BreakPredictionFromF83>());
    }

    [Fact]
    public void RegisterF89F87BreakPredictionFromF83_AncestorsContainTrulyInheritanceAndTrichotomy()
    {
        var registry = BuildBaseRegistry()
            .RegisterF89F87BreakPredictionFromF83()
            .Build();

        var ancestors = registry.AncestorsOf<F89F87BreakPredictionFromF83>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(F89F87TrulyInheritance), ancestors);
        Assert.Contains(typeof(F87TrichotomyClassification), ancestors);
    }

    [Fact]
    public void RegisterF89F87BreakPredictionFromF83_TierIsTier1Derived()
    {
        var registry = BuildBaseRegistry()
            .RegisterF89F87BreakPredictionFromF83()
            .Build();

        Assert.Equal(Tier.Tier1Derived, registry.Get<F89F87BreakPredictionFromF83>().Tier);
    }

    [Fact]
    public void RegisterF89F87BreakPredictionFromF83_WithoutTrulyInheritance_Throws()
    {
        Assert.Throws<InvariantViolationException>(() =>
            new ClaimRegistryBuilder()
                .RegisterF1Family(DefaultChain())
                .RegisterPi2Family()
                .RegisterF87Family()
                // Missing: RegisterF89F87TrulyInheritance
                .RegisterF89F87BreakPredictionFromF83()
                .Build());
    }
}
