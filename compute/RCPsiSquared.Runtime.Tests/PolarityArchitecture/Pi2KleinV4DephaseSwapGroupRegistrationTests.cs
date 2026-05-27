using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Knowledge;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class Pi2KleinV4DephaseSwapGroupRegistrationTests
{
    [Fact]
    public void RegisterPi2KleinV4DephaseSwapGroup_AddsClaim()
    {
        var registry = new ClaimRegistryBuilder()
            .RegisterPi2KleinV4DephaseSwapGroup()
            .Build();

        Assert.True(registry.Contains<Pi2KleinV4DephaseSwapGroup>());
    }

    [Fact]
    public void RegisterPi2KleinV4DephaseSwapGroup_TierIsTier1Derived()
    {
        var registry = new ClaimRegistryBuilder()
            .RegisterPi2KleinV4DephaseSwapGroup()
            .Build();

        Assert.Equal(Tier.Tier1Derived,
            registry.Get<Pi2KleinV4DephaseSwapGroup>().Tier);
    }

    [Fact]
    public void RegisterPi2KleinV4DephaseSwapGroup_NoCtorParents()
    {
        // Standalone primitive following the Pi2DyadicLadderRegistration pattern; no
        // b.Get<X>() edges. Building with only this registration must succeed.
        var registry = new ClaimRegistryBuilder()
            .RegisterPi2KleinV4DephaseSwapGroup()
            .Build();

        var ancestors = registry.AncestorsOf<Pi2KleinV4DephaseSwapGroup>().ToList();
        Assert.Empty(ancestors);
    }

    [Fact]
    public void BuildDefault_ContainsPi2KleinV4DephaseSwapGroup()
    {
        // The Claim is wired into the production registry via KnowledgeRegistryFactory.
        var registry = KnowledgeRegistryFactory.BuildDefault();
        Assert.True(registry.Contains<Pi2KleinV4DephaseSwapGroup>());
    }
}
