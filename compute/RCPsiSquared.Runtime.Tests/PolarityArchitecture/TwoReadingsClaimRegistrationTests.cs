using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class TwoReadingsClaimRegistrationTests
{
    private static ClaimRegistryBuilder BuildBaseRegistry() =>
        new ClaimRegistryBuilder()
            .RegisterPi2Family();

    [Fact]
    public void RegisterTwoReadingsClaim_AddsClaim()
    {
        var registry = BuildBaseRegistry()
            .RegisterTwoReadingsClaim()
            .Build();

        Assert.True(registry.Contains<TwoReadingsClaim>());
    }

    [Fact]
    public void RegisterTwoReadingsClaim_TierIsTier1Derived()
    {
        var registry = BuildBaseRegistry()
            .RegisterTwoReadingsClaim()
            .Build();

        Assert.Equal(Tier.Tier1Derived, registry.Get<TwoReadingsClaim>().Tier);
    }

    [Fact]
    public void RegisterTwoReadingsClaim_AncestorsContainPolynomialFoundation()
    {
        var registry = BuildBaseRegistry()
            .RegisterTwoReadingsClaim()
            .Build();

        var ancestors = registry.AncestorsOf<TwoReadingsClaim>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(PolynomialFoundationClaim), ancestors);
    }

    [Fact]
    public void RegisterTwoReadingsClaim_WithoutPi2Family_Throws()
    {
        Assert.Throws<InvariantViolationException>(() =>
            new ClaimRegistryBuilder()
                // Missing: RegisterPi2Family
                .RegisterTwoReadingsClaim()
                .Build());
    }
}
