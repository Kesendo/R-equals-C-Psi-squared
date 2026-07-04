using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

/// <summary>Schicht-1 wiring tests for <see cref="SeedExistenceCountingClaim"/>: the
/// seed-existence counting theorem r(0⁺) − r(∞) = N − 1 (odd N) on the (1,2) block, with
/// typed parent edges to <see cref="AbsorptionTheoremClaim"/> (the −2/−6 rung diagonal) and
/// <see cref="ChiralKClaim"/> (the chiral pairing λ_{N+1−k} = −λ_k behind the 3-to-1
/// triple↔resonance bijection).</summary>
public class SeedExistenceCountingClaimRegistrationTests
{
    private static ClaimRegistryBuilder BuildBaseRegistry() =>
        new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterPi2DyadicLadder()
            .RegisterAbsorptionTheoremClaim()
            .RegisterChiralK();

    [Fact]
    public void RegisterSeedExistenceCounting_AddsClaim_Tier1Derived()
    {
        var registry = BuildBaseRegistry()
            .RegisterSeedExistenceCountingClaim()
            .Build();
        Assert.True(registry.Contains<SeedExistenceCountingClaim>());
        Assert.Equal(Tier.Tier1Derived, registry.Get<SeedExistenceCountingClaim>().Tier);
    }

    [Fact]
    public void RegisterSeedExistenceCounting_AncestorsContainBothParents()
    {
        var registry = BuildBaseRegistry()
            .RegisterSeedExistenceCountingClaim()
            .Build();
        var ancestors = registry.AncestorsOf<SeedExistenceCountingClaim>()
            .Select(c => c.GetType()).ToHashSet();
        Assert.Contains(typeof(AbsorptionTheoremClaim), ancestors);
        Assert.Contains(typeof(ChiralKClaim), ancestors);
    }
}
