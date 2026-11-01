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
    public void RegisterUniversalCarrierClaim_AncestorsContainBothParentsAndNotTheDiscriminant()
    {
        var registry = BuildBaseRegistry()
            .RegisterUniversalCarrierClaim()
            .Build();

        var ancestors = registry.AncestorsOf<UniversalCarrierClaim>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(AbsorptionTheoremClaim), ancestors);
        Assert.Contains(typeof(Pi2DyadicLadderClaim), ancestors);
        // t_peak's four is the Absorption rung-2 rate, not the discriminant's operator-space four.
        Assert.DoesNotContain(typeof(PolynomialDiscriminantAnchorClaim), ancestors);
    }

    [Fact]
    public void RegisterUniversalCarrierClaim_WithoutAbsorptionTheorem_Throws()
    {
        Assert.Throws<InvariantViolationException>(() =>
            new ClaimRegistryBuilder()
                .RegisterPi2Family()
                .RegisterPi2DyadicLadder()
                // Missing: RegisterAbsorptionTheoremClaim
                .RegisterUniversalCarrierClaim()
                .Build());
    }
}
