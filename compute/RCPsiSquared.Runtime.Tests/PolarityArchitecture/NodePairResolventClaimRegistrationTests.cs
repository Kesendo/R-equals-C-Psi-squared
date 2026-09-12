using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Knowledge;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class NodePairResolventClaimRegistrationTests
{
    private static readonly Lazy<RCPsiSquared.Runtime.ObjectManager.ClaimRegistry> Registry =
        new(KnowledgeRegistryFactory.BuildDefault);

    [Fact]
    public void DefaultRegistry_ContainsTheTier1DerivedClaim()
    {
        var registry = Registry.Value;

        Assert.True(registry.Contains<NodePairResolventClaim>());
        Assert.Equal(Tier.Tier1Derived, registry.Get<NodePairResolventClaim>().Tier);
    }

    [Fact]
    public void EdgeIntoClaim_IsExactlyTheBlindSeatPremise()
    {
        var parents = Registry.Value.EdgesInto<NodePairResolventClaim>()
            .Select(edge => edge.Parent)
            .ToArray();

        Assert.Equal(new[] { typeof(SeatCutBlindnessClaim) }, parents);
        Assert.Same(
            Registry.Value.Get<SeatCutBlindnessClaim>(),
            Registry.Value.Get<NodePairResolventClaim>().SeatBlindness);
    }
}
