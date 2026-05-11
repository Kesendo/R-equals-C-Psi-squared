using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class F75MirrorPairMiPi2InheritanceRegistrationTests
{
    private static ClaimRegistryBuilder BuildBaseRegistry() =>
        new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterPi2DyadicLadder()
            .RegisterAbsorptionTheoremClaim()
            .RegisterHalfIntegerMirror(N: 5)
            .RegisterF71MirrorSymmetryPi2Inheritance()
            .RegisterF66PoleModesPi2Inheritance()
            .RegisterF65XxChainSpectrumPi2Inheritance();

    [Fact]
    public void RegisterF75_AddsClaim()
    {
        var registry = BuildBaseRegistry()
            .RegisterF75MirrorPairMiPi2Inheritance()
            .Build();

        Assert.True(registry.Contains<F75MirrorPairMiPi2Inheritance>());
    }

    [Fact]
    public void RegisterF75_TierIsTier1Derived()
    {
        var registry = BuildBaseRegistry()
            .RegisterF75MirrorPairMiPi2Inheritance()
            .Build();

        Assert.Equal(Tier.Tier1Derived,
            registry.Get<F75MirrorPairMiPi2Inheritance>().Tier);
    }

    [Fact]
    public void RegisterF75_AncestorsContainAllThreeAnchors()
    {
        var registry = BuildBaseRegistry()
            .RegisterF75MirrorPairMiPi2Inheritance()
            .Build();

        var ancestors = registry.AncestorsOf<F75MirrorPairMiPi2Inheritance>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(Pi2DyadicLadderClaim), ancestors);
        Assert.Contains(typeof(BilinearApexClaim), ancestors);
        Assert.Contains(typeof(F71MirrorSymmetryPi2Inheritance), ancestors);
    }

    [Fact]
    public void RegisterF75_AncestorsTransitivelyReachHalfIntegerMirrorAndPolynomialFoundation()
    {
        // F75 → F71 → HalfIntegerMirror → ... → PolynomialFoundation
        var registry = BuildBaseRegistry()
            .RegisterF75MirrorPairMiPi2Inheritance()
            .Build();

        var ancestors = registry.AncestorsOf<F75MirrorPairMiPi2Inheritance>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(HalfIntegerMirrorClaim), ancestors);
        Assert.Contains(typeof(PolynomialFoundationClaim), ancestors);
    }

    [Fact]
    public void RegisterF75_MaxMIPerPairIsTwoAcrossRegistry()
    {
        var registry = BuildBaseRegistry()
            .RegisterF75MirrorPairMiPi2Inheritance()
            .Build();

        Assert.Equal(2.0, registry.Get<F75MirrorPairMiPi2Inheritance>().MaxMIPerPair, precision: 14);
    }

    [Fact]
    public void RegisterF75_MIAtBilinearApexEqualsTwo()
    {
        var registry = BuildBaseRegistry()
            .RegisterF75MirrorPairMiPi2Inheritance()
            .Build();

        Assert.True(registry.Get<F75MirrorPairMiPi2Inheritance>().MIAtBilinearApexEqualsMaxMIPerPair());
    }

    [Theory]
    [InlineData(5, 2, 1.245)]
    [InlineData(7, 4, 1.245)]
    public void RegisterF75_BondingModeMMTableAcrossRegistry(int N, int k, double expected)
    {
        var registry = BuildBaseRegistry()
            .RegisterF75MirrorPairMiPi2Inheritance()
            .Build();

        Assert.Equal(expected, registry.Get<F75MirrorPairMiPi2Inheritance>().BondingModeMMAtZero(N, k), precision: 2);
    }
}
