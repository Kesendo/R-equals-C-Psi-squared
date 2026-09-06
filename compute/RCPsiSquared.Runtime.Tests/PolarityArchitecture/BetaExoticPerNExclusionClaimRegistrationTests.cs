using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

/// <summary>Schicht-1 wiring tests for <see cref="BetaExoticPerNExclusionClaim"/>: the β-exotic
/// excluded exactly at N = 5, N = 7, and N = 9, both R-parities, by certified disc-multiplicity readings. One typed
/// parent edge, to <see cref="SeedExistenceCountingClaim"/>: the endpoint-surplus theorem. Literal
/// drops and their character are additional per-N premises, not consequences of that parent.
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
