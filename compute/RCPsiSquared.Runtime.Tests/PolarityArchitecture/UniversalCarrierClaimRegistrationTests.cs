using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class UniversalCarrierClaimRegistrationTests
{
    private static ClaimRegistryBuilder BuildBaseRegistry() =>
        new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterPi2DyadicLadder()
            .RegisterAbsorptionTheoremClaim()
            .RegisterPolynomialDiscriminantAnchor();

    [Fact]
    public void RegisterUniversalCarrierClaim_AddsClaim()
    {
        var registry = BuildBaseRegistry()
            .RegisterUniversalCarrierClaim()
            .Build();

        Assert.True(registry.Contains<UniversalCarrierClaim>());
    }

    [Fact]
    public void RegisterUniversalCarrierClaim_TierIsTier1Derived()
    {
        var registry = BuildBaseRegistry()
            .RegisterUniversalCarrierClaim()
            .Build();

        Assert.Equal(Tier.Tier1Derived, registry.Get<UniversalCarrierClaim>().Tier);
    }

    [Fact]
    public void RegisterUniversalCarrierClaim_AncestorsContainAllThreeParents()
    {
        var registry = BuildBaseRegistry()
            .RegisterUniversalCarrierClaim()
            .Build();

        var ancestors = registry.AncestorsOf<UniversalCarrierClaim>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(AbsorptionTheoremClaim), ancestors);
        Assert.Contains(typeof(Pi2DyadicLadderClaim), ancestors);
        Assert.Contains(typeof(PolynomialDiscriminantAnchorClaim), ancestors);
    }

    [Fact]
    public void RegisterUniversalCarrierClaim_WithoutPolynomialDiscriminantAnchor_Throws()
    {
        Assert.Throws<InvariantViolationException>(() =>
            new ClaimRegistryBuilder()
                .RegisterPi2Family()
                .RegisterPi2DyadicLadder()
                .RegisterAbsorptionTheoremClaim()
                // Missing: RegisterPolynomialDiscriminantAnchor
                .RegisterUniversalCarrierClaim()
                .Build());
    }
}
