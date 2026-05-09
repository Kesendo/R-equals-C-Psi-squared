using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class F49cShadowCrossingPi2InheritanceRegistrationTests
{
    private static ClaimRegistryBuilder BuildBaseRegistry() =>
        new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterPi2DyadicLadder()
            .RegisterF88PopcountCoherence()
            .RegisterF88StaticDyadicAnchor()
            .RegisterPi2OperatorSpaceMirror();

    [Fact]
    public void RegisterF49c_AddsClaim()
    {
        var registry = BuildBaseRegistry()
            .RegisterF49cShadowCrossingPi2Inheritance()
            .Build();

        Assert.True(registry.Contains<F49cShadowCrossingPi2Inheritance>());
    }

    [Fact]
    public void RegisterF49c_TierIsTier1Derived()
    {
        var registry = BuildBaseRegistry()
            .RegisterF49cShadowCrossingPi2Inheritance()
            .Build();

        Assert.Equal(Tier.Tier1Derived,
            registry.Get<F49cShadowCrossingPi2Inheritance>().Tier);
    }

    [Fact]
    public void RegisterF49c_AncestorsContainBothPi2Anchors()
    {
        var registry = BuildBaseRegistry()
            .RegisterF49cShadowCrossingPi2Inheritance()
            .Build();

        var ancestors = registry.AncestorsOf<F49cShadowCrossingPi2Inheritance>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(Pi2DyadicLadderClaim), ancestors);
        Assert.Contains(typeof(Pi2OperatorSpaceMirrorClaim), ancestors);
    }

    [Theory]
    [InlineData(2, 4.0)]
    [InlineData(3, 16.0)]
    [InlineData(4, 64.0)]
    [InlineData(6, 1024.0)]
    public void RegisterF49c_FourPowerNMinus1FactorAcrossRegistry(int N, double expected)
    {
        var registry = BuildBaseRegistry()
            .RegisterF49cShadowCrossingPi2Inheritance()
            .Build();

        var f = registry.Get<F49cShadowCrossingPi2Inheritance>();
        Assert.Equal(expected, f.FourPowerNMinus1Factor(N), precision: 12);
        Assert.Equal(expected, f.MirrorPinnedFourPowerNMinus1(N), precision: 12);
    }

    [Theory]
    [InlineData(2, 1.0 / 8.0)]
    [InlineData(3, 2.0 / 48.0)]
    [InlineData(4, 3.0 / 256.0)]
    public void RegisterF49c_RSquaredAcrossRegistry(int N, double expected)
    {
        var registry = BuildBaseRegistry()
            .RegisterF49cShadowCrossingPi2Inheritance()
            .Build();

        Assert.Equal(expected, registry.Get<F49cShadowCrossingPi2Inheritance>().RSquared(N), precision: 12);
    }

    [Theory]
    [InlineData(2)]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    public void RegisterF49c_GapMatchesShadowCrossingMinusBalanced(int N)
    {
        var registry = BuildBaseRegistry()
            .RegisterF49cShadowCrossingPi2Inheritance()
            .Build();

        var f = registry.Get<F49cShadowCrossingPi2Inheritance>();
        Assert.Equal(f.RSquared(N) - f.F49ShadowBalancedRSquared(N), f.ShadowCrossingMinusBalancedGap(N), precision: 14);
    }

    [Fact]
    public void RegisterF49c_WithoutOperatorSpaceMirror_Throws()
    {
        Assert.Throws<InvariantViolationException>(() =>
            new ClaimRegistryBuilder()
                .RegisterPi2Family()
                .RegisterPi2DyadicLadder()
                // Missing: F88* + Pi2OperatorSpaceMirror
                .RegisterF49cShadowCrossingPi2Inheritance()
                .Build());
    }
}
