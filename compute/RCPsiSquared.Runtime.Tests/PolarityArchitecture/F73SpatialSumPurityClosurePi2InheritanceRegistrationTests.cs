using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class F73SpatialSumPurityClosurePi2InheritanceRegistrationTests
{
    private static ClaimRegistryBuilder BuildBaseRegistry() =>
        new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterPi2DyadicLadder()
            .RegisterF70DeltaNSelectionRulePi2Inheritance()
            .RegisterF72BlockDiagonalPurityPi2Inheritance();

    [Fact]
    public void RegisterF73_AddsClaim()
    {
        var registry = BuildBaseRegistry()
            .RegisterF73SpatialSumPurityClosurePi2Inheritance()
            .Build();

        Assert.True(registry.Contains<F73SpatialSumPurityClosurePi2Inheritance>());
    }

    [Fact]
    public void RegisterF73_TierIsTier1Derived()
    {
        var registry = BuildBaseRegistry()
            .RegisterF73SpatialSumPurityClosurePi2Inheritance()
            .Build();

        Assert.Equal(Tier.Tier1Derived,
            registry.Get<F73SpatialSumPurityClosurePi2Inheritance>().Tier);
    }

    [Fact]
    public void RegisterF73_AncestorsContainAllThreeParents()
    {
        var registry = BuildBaseRegistry()
            .RegisterF73SpatialSumPurityClosurePi2Inheritance()
            .Build();

        var ancestors = registry.AncestorsOf<F73SpatialSumPurityClosurePi2Inheritance>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(Pi2DyadicLadderClaim), ancestors);
        Assert.Contains(typeof(F70DeltaNSelectionRulePi2Inheritance), ancestors);
        Assert.Contains(typeof(F72BlockDiagonalPurityPi2Inheritance), ancestors);
    }

    [Fact]
    public void RegisterF73_AncestorsTransitivelyReachPolynomialFoundation()
    {
        var registry = BuildBaseRegistry()
            .RegisterF73SpatialSumPurityClosurePi2Inheritance()
            .Build();

        var ancestors = registry.AncestorsOf<F73SpatialSumPurityClosurePi2Inheritance>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(PolynomialFoundationClaim), ancestors);
    }

    [Fact]
    public void RegisterF73_VerifiedValueMatchesAcrossRegistry()
    {
        var registry = BuildBaseRegistry()
            .RegisterF73SpatialSumPurityClosurePi2Inheritance()
            .Build();

        Assert.Equal(9.157819e-3, registry.Get<F73SpatialSumPurityClosurePi2Inheritance>().VerifiedValueAtN5Gamma0p05T20(), precision: 6);
    }

    [Fact]
    public void RegisterF73_ClosureAtTZeroIsHalfAcrossRegistry()
    {
        var registry = BuildBaseRegistry()
            .RegisterF73SpatialSumPurityClosurePi2Inheritance()
            .Build();

        Assert.True(registry.Get<F73SpatialSumPurityClosurePi2Inheritance>().ClosureAtTZeroIsHalf());
    }
}
