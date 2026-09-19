using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class TwoReadingsClaimRegistrationTests
{
    [Fact]
    public void RegisterTwoReadingsClaim_AddsClaim()
    {
        var registry = new ClaimRegistryBuilder()
            .RegisterTwoReadingsClaim()
            .Build();

        Assert.True(registry.Contains<TwoReadingsClaim>());
    }

    [Fact]
    public void RegisterTwoReadingsClaim_TierIsOpenQuestion()
    {
        var registry = new ClaimRegistryBuilder()
            .RegisterTwoReadingsClaim()
            .Build();

        Assert.Equal(Tier.OpenQuestion, registry.Get<TwoReadingsClaim>().Tier);
    }

    [Fact]
    public void RegisterTwoReadingsClaim_HasNoTypedAncestors()
    {
        var registry = new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterTwoReadingsClaim()
            .Build();

        Assert.Empty(registry.AncestorsOf<TwoReadingsClaim>());
    }

    [Fact]
    public void RegisterTwoReadingsClaim_DoesNotRequirePi2Family()
    {
        var registry = new ClaimRegistryBuilder()
            .RegisterTwoReadingsClaim()
            .Build();

        Assert.True(registry.Contains<TwoReadingsClaim>());
    }
}
