using System.Linq;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Knowledge;
using Xunit;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class ThreeLadderHingeClaimRegistrationTests
{
    [Fact]
    public void BuildDefault_ContainsThreeLadderHingeClaim_WithBothParents()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        Assert.True(registry.Contains<ThreeLadderHingeClaim>());

        var claim = registry.Get<ThreeLadderHingeClaim>();
        Assert.Equal(Tier.Tier1Derived, claim.Tier);
        Assert.Equal(claim.Cases.Count, claim.PassCount);

        // the weld: both parents must be the SHARED registry instances
        Assert.Same(registry.Get<AbsorptionTheoremClaim>(), claim.Absorption);
        Assert.Same(registry.Get<MomentTowerPumpChannelClaim>(), claim.Moment);
    }

    [Fact]
    public void Ancestors_IncludeBothLadderParents()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var ancestors = registry.AncestorsOf<ThreeLadderHingeClaim>()
            .Select(c => c.GetType()).ToHashSet();
        Assert.Contains(typeof(AbsorptionTheoremClaim), ancestors);
        Assert.Contains(typeof(MomentTowerPumpChannelClaim), ancestors);
    }
}
