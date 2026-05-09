using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class F72BlockDiagonalPurityPi2InheritanceRegistrationTests
{
    private static ClaimRegistryBuilder BuildBaseRegistry() =>
        new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterPi2DyadicLadder()
            .RegisterF70DeltaNSelectionRulePi2Inheritance();

    [Fact]
    public void RegisterF72_AddsClaim()
    {
        var registry = BuildBaseRegistry()
            .RegisterF72BlockDiagonalPurityPi2Inheritance()
            .Build();

        Assert.True(registry.Contains<F72BlockDiagonalPurityPi2Inheritance>());
    }

    [Fact]
    public void RegisterF72_TierIsTier1Derived()
    {
        var registry = BuildBaseRegistry()
            .RegisterF72BlockDiagonalPurityPi2Inheritance()
            .Build();

        Assert.Equal(Tier.Tier1Derived,
            registry.Get<F72BlockDiagonalPurityPi2Inheritance>().Tier);
    }

    [Fact]
    public void RegisterF72_AncestorsContainBothLadderAndF70()
    {
        var registry = BuildBaseRegistry()
            .RegisterF72BlockDiagonalPurityPi2Inheritance()
            .Build();

        var ancestors = registry.AncestorsOf<F72BlockDiagonalPurityPi2Inheritance>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(Pi2DyadicLadderClaim), ancestors);
        Assert.Contains(typeof(F70DeltaNSelectionRulePi2Inheritance), ancestors);
        Assert.Contains(typeof(HalfAsStructuralFixedPointClaim), ancestors);
    }

    [Fact]
    public void RegisterF72_AncestorsTransitivelyReachPolynomialFoundation()
    {
        var registry = BuildBaseRegistry()
            .RegisterF72BlockDiagonalPurityPi2Inheritance()
            .Build();

        var ancestors = registry.AncestorsOf<F72BlockDiagonalPurityPi2Inheritance>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(PolynomialFoundationClaim), ancestors);
    }

    [Fact]
    public void RegisterF72_BaselineIsExactlyOneHalfAcrossRegistry()
    {
        var registry = BuildBaseRegistry()
            .RegisterF72BlockDiagonalPurityPi2Inheritance()
            .Build();

        Assert.Equal(0.5, registry.Get<F72BlockDiagonalPurityPi2Inheritance>().BaselineTraceSquared, precision: 14);
    }

    [Fact]
    public void RegisterF72_SingleSiteBlockCountIsTwoAcrossRegistry()
    {
        var registry = BuildBaseRegistry()
            .RegisterF72BlockDiagonalPurityPi2Inheritance()
            .Build();

        Assert.Equal(2.0, registry.Get<F72BlockDiagonalPurityPi2Inheritance>().SingleSiteBlockCount, precision: 14);
    }

    [Theory]
    [InlineData(1, 2)]
    [InlineData(2, 3)]
    [InlineData(3, 4)]
    public void RegisterF72_SubBlockCountAcrossRegistry(int kLocal, int expected)
    {
        var registry = BuildBaseRegistry()
            .RegisterF72BlockDiagonalPurityPi2Inheritance()
            .Build();

        Assert.Equal(expected, registry.Get<F72BlockDiagonalPurityPi2Inheritance>().SubBlockCountForKLocal(kLocal));
    }

    [Fact]
    public void RegisterF72_SingleSiteBlockCountMatchesF70PlusOne()
    {
        var registry = BuildBaseRegistry()
            .RegisterF72BlockDiagonalPurityPi2Inheritance()
            .Build();

        Assert.True(registry.Get<F72BlockDiagonalPurityPi2Inheritance>().SingleSiteBlockCountMatchesF70PlusOne());
    }
}
