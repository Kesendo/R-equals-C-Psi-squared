using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

/// <summary>Schicht-1 wiring tests for <see cref="BetaExoticPerNExclusionClaim"/>: the β-exotic
/// excluded exactly at N = 5 and N = 7, both R-parities, by the certified disc-multiplicity reading. One typed
/// parent edge, to <see cref="SeedExistenceCountingClaim"/>: the theorem that FORCES the count-drops
/// whose character this claim pins, and whose own scope note names the β-exotic as its open item.
///
/// <para>The ancestor set must also reach the counting theorem's own parents, since a per-N
/// certificate about the count-drops rests on the pencil those parents define.</para></summary>
public class BetaExoticPerNExclusionClaimRegistrationTests
{
    private static ClaimRegistryBuilder BuildBaseRegistry() =>
        new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterPi2DyadicLadder()
            .RegisterAbsorptionTheoremClaim()
            .RegisterChiralK()
            .RegisterSeedExistenceCountingClaim();

    [Fact]
    public void RegisterBetaExoticPerNExclusion_AddsClaim_Tier1Derived()
    {
        var registry = BuildBaseRegistry()
            .RegisterBetaExoticPerNExclusionClaim()
            .Build();
        Assert.True(registry.Contains<BetaExoticPerNExclusionClaim>());
        Assert.Equal(Tier.Tier1Derived, registry.Get<BetaExoticPerNExclusionClaim>().Tier);
    }

    [Fact]
    public void RegisterBetaExoticPerNExclusion_AncestorsReachTheCountingTheoremAndItsParents()
    {
        var registry = BuildBaseRegistry()
            .RegisterBetaExoticPerNExclusionClaim()
            .Build();
        var ancestors = registry.AncestorsOf<BetaExoticPerNExclusionClaim>()
            .Select(c => c.GetType()).ToHashSet();
        Assert.Contains(typeof(SeedExistenceCountingClaim), ancestors);
        Assert.Contains(typeof(AbsorptionTheoremClaim), ancestors);
        Assert.Contains(typeof(ChiralKClaim), ancestors);
    }
}
