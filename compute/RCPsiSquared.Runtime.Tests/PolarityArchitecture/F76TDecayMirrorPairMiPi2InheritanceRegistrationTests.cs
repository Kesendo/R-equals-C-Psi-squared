using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class F76TDecayMirrorPairMiPi2InheritanceRegistrationTests
{
    private static ClaimRegistryBuilder BuildBaseRegistry() =>
        new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterPi2DyadicLadder()
            .RegisterAbsorptionTheoremClaim()
            .RegisterHalfIntegerMirror(N: 5)
            .RegisterF71MirrorSymmetryPi2Inheritance()
            .RegisterF66PoleModesPi2Inheritance()
            .RegisterF65XxChainSpectrumPi2Inheritance()
            .RegisterF75MirrorPairMiPi2Inheritance();

    [Fact]
    public void RegisterF76_AddsClaim()
    {
        var registry = BuildBaseRegistry()
            .RegisterF76TDecayMirrorPairMiPi2Inheritance()
            .Build();

        Assert.True(registry.Contains<F76TDecayMirrorPairMiPi2Inheritance>());
    }

    [Fact]
    public void RegisterF76_TierIsTier1Derived()
    {
        var registry = BuildBaseRegistry()
            .RegisterF76TDecayMirrorPairMiPi2Inheritance()
            .Build();

        Assert.Equal(Tier.Tier1Derived,
            registry.Get<F76TDecayMirrorPairMiPi2Inheritance>().Tier);
    }

    [Fact]
    public void RegisterF76_AncestorsContainBothParents()
    {
        var registry = BuildBaseRegistry()
            .RegisterF76TDecayMirrorPairMiPi2Inheritance()
            .Build();

        var ancestors = registry.AncestorsOf<F76TDecayMirrorPairMiPi2Inheritance>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(Pi2DyadicLadderClaim), ancestors);
        Assert.Contains(typeof(F75MirrorPairMiPi2Inheritance), ancestors);
    }

    [Fact]
    public void RegisterF76_AncestorsTransitivelyReachF71_F75_BilinearApex()
    {
        // F76 → F75 → {F71, BilinearApex} → ... → PolynomialFoundation
        var registry = BuildBaseRegistry()
            .RegisterF76TDecayMirrorPairMiPi2Inheritance()
            .Build();

        var ancestors = registry.AncestorsOf<F76TDecayMirrorPairMiPi2Inheritance>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(F71MirrorSymmetryPi2Inheritance), ancestors);
        Assert.Contains(typeof(BilinearApexClaim), ancestors);
        Assert.Contains(typeof(HalfIntegerMirrorClaim), ancestors);
        Assert.Contains(typeof(PolynomialFoundationClaim), ancestors);
    }

    [Fact]
    public void RegisterF76_DecayRateCoefficientIsFourAcrossRegistry()
    {
        var registry = BuildBaseRegistry()
            .RegisterF76TDecayMirrorPairMiPi2Inheritance()
            .Build();

        Assert.Equal(4.0, registry.Get<F76TDecayMirrorPairMiPi2Inheritance>().DecayRateCoefficient, precision: 14);
    }

    [Theory]
    [InlineData(5, 2, 0.05, 0.1, 0.936)]
    [InlineData(7, 2, 0.05, 0.1, 0.932)]
    public void RegisterF76_EnvelopeRatioAcrossRegistry(int N, int k, double gammaZero, double t, double expected)
    {
        var registry = BuildBaseRegistry()
            .RegisterF76TDecayMirrorPairMiPi2Inheritance()
            .Build();

        Assert.Equal(expected, registry.Get<F76TDecayMirrorPairMiPi2Inheritance>().EnvelopeRatioForBondingMode(N, k, gammaZero, t), precision: 2);
    }

    [Fact]
    public void RegisterF76_RecoversF75AtZero_AcrossRegistry()
    {
        var registry = BuildBaseRegistry()
            .RegisterF76TDecayMirrorPairMiPi2Inheritance()
            .Build();

        var f76 = registry.Get<F76TDecayMirrorPairMiPi2Inheritance>();
        Assert.True(f76.RecoversF75AtZero(0.25));
        Assert.True(f76.RecoversF75AtZero(0.4));
        Assert.True(f76.RecoversF75AtZero(0.5));
    }
}
