using RCPsiSquared.Core.F86;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Runtime.F86Main;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.Tests.F86Main;

public class F86ExtendedRegistrationTests
{
    [Fact]
    public void RegisterF86Extended_BuildsFourClaims()
    {
        var registry = new ClaimRegistryBuilder()
            .RegisterF86Extended(gammaZero: 0.05)
            .Build();

        Assert.Equal(4, registry.All().Count());
        Assert.True(registry.Contains<PerBlockQPeakClaim>());
        Assert.True(registry.Contains<UniversalShapePrediction>());
        Assert.True(registry.Contains<SigmaZeroChromaticityScaling>());
        Assert.True(registry.Contains<RetractedClaim>());
    }

    [Fact]
    public void RegisterF86Extended_TierBandsCovered()
    {
        var registry = new ClaimRegistryBuilder()
            .RegisterF86Extended(gammaZero: 0.05)
            .Build();

        Assert.Equal(Tier.Tier1Candidate, registry.Get<PerBlockQPeakClaim>().Tier);
        Assert.Equal(Tier.Tier1Candidate, registry.Get<UniversalShapePrediction>().Tier);
        Assert.Equal(Tier.Tier2Empirical, registry.Get<SigmaZeroChromaticityScaling>().Tier);
        Assert.Equal(Tier.Retracted, registry.Get<RetractedClaim>().Tier);
    }

    [Fact]
    public void RegisterF86Extended_AllOfTier_FiltersCorrectly()
    {
        var registry = new ClaimRegistryBuilder()
            .RegisterF86Extended(gammaZero: 0.05)
            .Build();

        Assert.Equal(2, registry.AllOfTier(Tier.Tier1Candidate).Count);
        Assert.Single(registry.AllOfTier(Tier.Tier2Empirical));
        Assert.Single(registry.AllOfTier(Tier.Retracted));
        Assert.Empty(registry.AllOfTier(Tier.Tier1Derived));
    }
}
