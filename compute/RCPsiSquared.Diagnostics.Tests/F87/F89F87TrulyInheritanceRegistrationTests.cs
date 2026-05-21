using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Diagnostics.F87;
using RCPsiSquared.Runtime.F1Family;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Diagnostics.Tests.F87;

public class F89F87TrulyInheritanceRegistrationTests
{
    private static ChainSystem DefaultChain() =>
        new ChainSystem(N: 5, J: 1.0, GammaZero: 0.05,
            HType: HamiltonianType.Heisenberg, Topology: TopologyKind.Chain);

    private static ClaimRegistryBuilder BuildBaseRegistry() =>
        new ClaimRegistryBuilder()
            .RegisterF1Family(DefaultChain())
            .RegisterPi2Family()
            .RegisterF87Family();

    [Fact]
    public void RegisterF89F87TrulyInheritance_AddsClaim()
    {
        var registry = BuildBaseRegistry()
            .RegisterF89F87TrulyInheritance()
            .Build();

        Assert.True(registry.Contains<F89F87TrulyInheritance>());
    }

    [Fact]
    public void RegisterF89F87TrulyInheritance_AncestorsContainTrichotomyClassification()
    {
        var registry = BuildBaseRegistry()
            .RegisterF89F87TrulyInheritance()
            .Build();

        var ancestors = registry.AncestorsOf<F89F87TrulyInheritance>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(F87TrichotomyClassification), ancestors);
    }

    [Fact]
    public void RegisterF89F87TrulyInheritance_TierIsTier1Derived()
    {
        var registry = BuildBaseRegistry()
            .RegisterF89F87TrulyInheritance()
            .Build();

        Assert.Equal(Tier.Tier1Derived, registry.Get<F89F87TrulyInheritance>().Tier);
    }

    [Fact]
    public void RegisterF89F87TrulyInheritance_WithoutF87Family_Throws()
    {
        Assert.Throws<InvariantViolationException>(() =>
            new ClaimRegistryBuilder()
                .RegisterF1Family(DefaultChain())
                .RegisterPi2Family()
                // Missing: RegisterF87Family
                .RegisterF89F87TrulyInheritance()
                .Build());
    }
}
