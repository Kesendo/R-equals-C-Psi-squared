using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class PolynomialDiscriminantAnchorRegistrationTests
{
    private static ClaimRegistryBuilder BuildBaseRegistry() =>
        new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterPi2DyadicLadder();

    [Fact]
    public void RegisterPolynomialDiscriminantAnchor_AddsClaim()
    {
        var registry = BuildBaseRegistry()
            .RegisterPolynomialDiscriminantAnchor()
            .Build();
        Assert.True(registry.Contains<PolynomialDiscriminantAnchorClaim>());
    }

    [Fact]
    public void RegisterPolynomialDiscriminantAnchor_TierIsTier1Derived()
    {
        var registry = BuildBaseRegistry()
            .RegisterPolynomialDiscriminantAnchor()
            .Build();
        Assert.Equal(Tier.Tier1Derived, registry.Get<PolynomialDiscriminantAnchorClaim>().Tier);
    }

    [Fact]
    public void RegisterPolynomialDiscriminantAnchor_AncestorsContainAllThreeParents()
    {
        var registry = BuildBaseRegistry()
            .RegisterPolynomialDiscriminantAnchor()
            .Build();
        var ancestors = registry.AncestorsOf<PolynomialDiscriminantAnchorClaim>()
            .Select(c => c.GetType()).ToHashSet();
        Assert.Contains(typeof(PolynomialFoundationClaim), ancestors);
        Assert.Contains(typeof(QubitDimensionalAnchorClaim), ancestors);
        Assert.Contains(typeof(Pi2DyadicLadderClaim), ancestors);
    }

    [Fact]
    public void RegisterPolynomialDiscriminantAnchor_DiscriminantIsFour()
    {
        var registry = BuildBaseRegistry()
            .RegisterPolynomialDiscriminantAnchor()
            .Build();
        Assert.Equal(4.0, registry.Get<PolynomialDiscriminantAnchorClaim>().DiscriminantViaCoefficients, precision: 14);
        Assert.Equal(4.0, registry.Get<PolynomialDiscriminantAnchorClaim>().DiscriminantViaLadder, precision: 14);
    }

    [Fact]
    public void RegisterPolynomialDiscriminantAnchor_AllReadingsAgreeAcrossRegistry()
    {
        var registry = BuildBaseRegistry()
            .RegisterPolynomialDiscriminantAnchor()
            .Build();
        Assert.True(registry.Get<PolynomialDiscriminantAnchorClaim>().AllReadingsAgree());
    }

    [Fact]
    public void RegisterPolynomialDiscriminantAnchor_MirrorClosureHoldsAcrossRegistry()
    {
        var registry = BuildBaseRegistry()
            .RegisterPolynomialDiscriminantAnchor()
            .Build();
        Assert.True(registry.Get<PolynomialDiscriminantAnchorClaim>().MirrorClosureHolds());
    }

    [Fact]
    public void RegisterPolynomialDiscriminantAnchor_WithoutPi2DyadicLadder_Throws()
    {
        Assert.Throws<InvariantViolationException>(() =>
            new ClaimRegistryBuilder()
                .RegisterPi2Family()
                .RegisterPolynomialDiscriminantAnchor()
                .Build());
    }
}
