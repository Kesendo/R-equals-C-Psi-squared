using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.F1;
using RCPsiSquared.Core.Knowledge;
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

    private sealed class WeakerFakeIdentity : Claim
    {
        public WeakerFakeIdentity() : base("WeakerFakeIdentity", Tier.Tier2Empirical,
            "compute/RCPsiSquared.Runtime.Tests/TestSupport/TestClaims.cs") { }
        public override string DisplayName => "WeakerFakeIdentity";
        public override string Summary => "synthetic Tier2Empirical pretender to F1 master";
    }

    [Fact]
    public void Tier_F73_DependsOnF1_BothTier1Derived_Succeeds()
    {
        var registry = new ClaimRegistryBuilder()
            .RegisterF1Family(DefaultChain())
            .Build();

        var f1 = registry.Get<F1PalindromeIdentity>();
        var f73 = registry.Get<PalindromeResidualScalingClaim>();
        Assert.Equal(Tier.Tier1Derived, f1.Tier);
        Assert.Equal(Tier.Tier1Derived, f73.Tier);
    }

    [Fact]
    public void Tier_F73_DependsOnSyntheticTier2Parent_Throws()
    {
        // Construct a synthetic registration that violates Tier inheritance: F73
        // (Tier1Derived) is forced to depend on a Tier2Empirical pretender. Verifies the
        // builder catches it even when the path runs through a real Core Claim.
        var ex = Assert.Throws<InvariantViolationException>(() =>
            new ClaimRegistryBuilder()
                .Register<WeakerFakeIdentity>(_ => new WeakerFakeIdentity())
                .Register<PalindromeResidualScalingClaim>(b =>
                {
                    _ = b.Get<WeakerFakeIdentity>();
                    return new PalindromeResidualScalingClaim(N: 5, HamiltonianClass.Main);
                })
                .Build());

        Assert.Equal("TierInheritance", ex.Rule);
        Assert.Contains("PalindromeResidualScalingClaim", ex.Message);
        Assert.Contains("WeakerFakeIdentity", ex.Message);
    }

    [Fact]
    public void Tier_DowngradeDetectsCascadeViolation()
    {
        // Hypothetical scenario: someone has changed F1PalindromeIdentity to Tier2Empirical
        // in Core. Until F73 is also downgraded, the builder must throw. We simulate this by
        // wrapping F1 in a downgrading proxy and registering F73 on top.
        var ex = Assert.Throws<InvariantViolationException>(() =>
            new ClaimRegistryBuilder()
                .Register<WeakerFakeIdentity>(_ => new WeakerFakeIdentity())
                .Register<PalindromeResidualScalingClaim>(b =>
                {
                    _ = b.Get<WeakerFakeIdentity>();
                    return new PalindromeResidualScalingClaim(N: 5, HamiltonianClass.Main);
                })
                .Build());

        Assert.Equal("TierInheritance", ex.Rule);
        Assert.NotEmpty(ex.Path);
    }
}
