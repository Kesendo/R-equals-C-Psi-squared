using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class F66PoleModesPi2InheritanceRegistrationTests
{
    private static ClaimRegistryBuilder BuildBaseRegistry() =>
        new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterPi2DyadicLadder();

    [Fact]
    public void RegisterF66_AddsClaim()
    {
        var registry = BuildBaseRegistry()
            .RegisterF66PoleModesPi2Inheritance()
            .Build();

        Assert.True(registry.Contains<F66PoleModesPi2Inheritance>());
    }

    [Fact]
    public void RegisterF66_TierIsTier1Derived()
    {
        var registry = BuildBaseRegistry()
            .RegisterF66PoleModesPi2Inheritance()
            .Build();

        Assert.Equal(Tier.Tier1Derived,
            registry.Get<F66PoleModesPi2Inheritance>().Tier);
    }

    [Fact]
    public void RegisterF66_AncestorsContainBothPi2Anchors()
    {
        var registry = BuildBaseRegistry()
            .RegisterF66PoleModesPi2Inheritance()
            .Build();

        var ancestors = registry.AncestorsOf<F66PoleModesPi2Inheritance>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(Pi2DyadicLadderClaim), ancestors);
        Assert.Contains(typeof(QubitDimensionalAnchorClaim), ancestors);
    }

    [Fact]
    public void RegisterF66_AncestorsTransitivelyReachPolynomialFoundation()
    {
        var registry = BuildBaseRegistry()
            .RegisterF66PoleModesPi2Inheritance()
            .Build();

        var ancestors = registry.AncestorsOf<F66PoleModesPi2Inheritance>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(PolynomialFoundationClaim), ancestors);
    }

    [Fact]
    public void RegisterF66_UpperPoleCoefficientIsTwo()
    {
        var registry = BuildBaseRegistry()
            .RegisterF66PoleModesPi2Inheritance()
            .Build();

        Assert.Equal(2.0, registry.Get<F66PoleModesPi2Inheritance>().UpperPoleCoefficient, precision: 14);
    }

    [Fact]
    public void RegisterF66_UpperPoleCoefficientMatchesQubitDimensionAcrossRegistry()
    {
        var registry = BuildBaseRegistry()
            .RegisterF66PoleModesPi2Inheritance()
            .Build();

        Assert.True(registry.Get<F66PoleModesPi2Inheritance>().UpperPoleCoefficientMatchesQubitDimension());
    }

    [Theory]
    [InlineData(0.05, 0.1)]
    [InlineData(0.5, 1.0)]
    [InlineData(2.5, 5.0)]
    public void RegisterF66_UpperPoleAlphaAcrossRegistry(double gammaZero, double expected)
    {
        var registry = BuildBaseRegistry()
            .RegisterF66PoleModesPi2Inheritance()
            .Build();

        Assert.Equal(expected, registry.Get<F66PoleModesPi2Inheritance>().UpperPoleAlpha(gammaZero), precision: 12);
    }

    [Theory]
    [InlineData(3, 4)]
    [InlineData(5, 6)]
    [InlineData(7, 8)]
    public void RegisterF66_EndpointMultiplicityIsNPlus1(int N, int expected)
    {
        var registry = BuildBaseRegistry()
            .RegisterF66PoleModesPi2Inheritance()
            .Build();

        Assert.Equal(expected, registry.Get<F66PoleModesPi2Inheritance>().EndpointMultiplicity(N));
    }

    [Theory]
    [InlineData(3)]
    [InlineData(5)]
    [InlineData(7)]
    public void RegisterF66_PalindromicWeightPairIsZeroAndN(int N)
    {
        var registry = BuildBaseRegistry()
            .RegisterF66PoleModesPi2Inheritance()
            .Build();

        var (l, u) = registry.Get<F66PoleModesPi2Inheritance>().PalindromicWeightPairOfPoles(N);
        Assert.Equal(0, l);
        Assert.Equal(N, u);
    }
}
