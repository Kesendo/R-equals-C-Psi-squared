using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class F77MmSaturationPi2InheritanceRegistrationTests
{
    private static ClaimRegistryBuilder BuildBaseRegistry() =>
        new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterPi2DyadicLadder()
            .RegisterHalfIntegerMirror(N: 5)
            .RegisterF71MirrorSymmetryPi2Inheritance()
            .RegisterF75MirrorPairMiPi2Inheritance();

    [Fact]
    public void RegisterF77_AddsClaim()
    {
        var registry = BuildBaseRegistry()
            .RegisterF77MmSaturationPi2Inheritance()
            .Build();

        Assert.True(registry.Contains<F77MmSaturationPi2Inheritance>());
    }

    [Fact]
    public void RegisterF77_TierIsTier1Derived()
    {
        var registry = BuildBaseRegistry()
            .RegisterF77MmSaturationPi2Inheritance()
            .Build();

        Assert.Equal(Tier.Tier1Derived,
            registry.Get<F77MmSaturationPi2Inheritance>().Tier);
    }

    [Fact]
    public void RegisterF77_AncestorsContainLadderAndHalfFixedPoint()
    {
        var registry = BuildBaseRegistry()
            .RegisterF77MmSaturationPi2Inheritance()
            .Build();

        var ancestors = registry.AncestorsOf<F77MmSaturationPi2Inheritance>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(Pi2DyadicLadderClaim), ancestors);
        Assert.Contains(typeof(HalfAsStructuralFixedPointClaim), ancestors);
    }

    [Fact]
    public void RegisterF77_AncestorsTransitivelyReachPolynomialFoundation()
    {
        var registry = BuildBaseRegistry()
            .RegisterF77MmSaturationPi2Inheritance()
            .Build();

        var ancestors = registry.AncestorsOf<F77MmSaturationPi2Inheritance>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(PolynomialFoundationClaim), ancestors);
    }

    [Fact]
    public void RegisterF77_SaturationBitIsExactlyOne()
    {
        var registry = BuildBaseRegistry()
            .RegisterF77MmSaturationPi2Inheritance()
            .Build();

        Assert.Equal(1.0, registry.Get<F77MmSaturationPi2Inheritance>().SaturationBit, precision: 14);
    }

    [Fact]
    public void RegisterF77_LandsOnSelfMirrorPivot()
    {
        var registry = BuildBaseRegistry()
            .RegisterF77MmSaturationPi2Inheritance()
            .Build();

        Assert.True(registry.Get<F77MmSaturationPi2Inheritance>().LandsOnSelfMirrorPivot);
    }

    [Fact]
    public void RegisterF77_InversionIdentityHoldsAcrossRegistry()
    {
        // Cross-registry verification: a_0 · a_2 = 1 = a_1 = SaturationBit.
        var registry = BuildBaseRegistry()
            .RegisterF77MmSaturationPi2Inheritance()
            .Build();

        Assert.True(registry.Get<F77MmSaturationPi2Inheritance>().InversionIdentityHolds);
    }

    [Theory]
    [InlineData(101)]
    [InlineData(1001)]
    [InlineData(10001)]
    public void RegisterF77_RescaledDeviationConvergesAcrossRegistry(int N)
    {
        var registry = BuildBaseRegistry()
            .RegisterF77MmSaturationPi2Inheritance()
            .Build();

        var f = registry.Get<F77MmSaturationPi2Inheritance>();
        Assert.Equal(f.AsymptoticCorrectionCoefficient, f.RescaledDeviation(N), precision: 10);
    }
}
