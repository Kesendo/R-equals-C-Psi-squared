using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

/// <summary>Schicht-1 wiring tests for <see cref="BetaExoticExcludedAtN5Claim"/>: the β-exotic
/// excluded exactly at N = 5, both R-parities, by the certified disc-multiplicity reading. One typed
/// parent edge, to <see cref="SeedExistenceCountingClaim"/> — the theorem that FORCES the count-drops
/// whose character this claim pins, and whose own scope note names the β-exotic as its open item.
///
/// <para>The ancestor set must also reach the counting theorem's own parents, since a per-N
/// certificate about the count-drops rests on the pencil those parents define.</para></summary>
public class BetaExoticExcludedAtN5ClaimRegistrationTests
{
    private static ClaimRegistryBuilder BuildBaseRegistry() =>
        new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterPi2DyadicLadder()
            .RegisterAbsorptionTheoremClaim()
            .RegisterChiralK()
            .RegisterSeedExistenceCountingClaim();

    [Fact]
    public void RegisterBetaExoticExcludedAtN5_AddsClaim_Tier1Derived()
    {
        var registry = BuildBaseRegistry()
            .RegisterBetaExoticExcludedAtN5Claim()
            .Build();
        Assert.True(registry.Contains<BetaExoticExcludedAtN5Claim>());
        Assert.Equal(Tier.Tier1Derived, registry.Get<BetaExoticExcludedAtN5Claim>().Tier);
    }

    [Fact]
    public void RegisterBetaExoticExcludedAtN5_AncestorsReachTheCountingTheoremAndItsParents()
    {
        var registry = BuildBaseRegistry()
            .RegisterBetaExoticExcludedAtN5Claim()
            .Build();
        var ancestors = registry.AncestorsOf<BetaExoticExcludedAtN5Claim>()
            .Select(c => c.GetType()).ToHashSet();
        Assert.Contains(typeof(SeedExistenceCountingClaim), ancestors);
        Assert.Contains(typeof(AbsorptionTheoremClaim), ancestors);
        Assert.Contains(typeof(ChiralKClaim), ancestors);
    }
}
