using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class F97CardioidHalfFixedPointPi2InheritanceRegistrationTests
{
    private static ClaimRegistryBuilder BuildBaseRegistry() =>
        new ClaimRegistryBuilder()
            .RegisterPi2Family();

    [Fact]
    public void RegisterF97_AddsClaim()
    {
        var registry = BuildBaseRegistry()
            .RegisterF97CardioidHalfFixedPointPi2Inheritance()
            .Build();

        Assert.True(registry.Contains<F97CardioidHalfFixedPointPi2Inheritance>());
    }

    [Fact]
    public void RegisterF97_TierIsTier1Derived()
    {
        var registry = BuildBaseRegistry()
            .RegisterF97CardioidHalfFixedPointPi2Inheritance()
            .Build();

        Assert.Equal(Tier.Tier1Derived,
            registry.Get<F97CardioidHalfFixedPointPi2Inheritance>().Tier);
    }

    [Fact]
    public void RegisterF97_AncestorsContainHalfAndQuarter()
    {
        var registry = BuildBaseRegistry()
            .RegisterF97CardioidHalfFixedPointPi2Inheritance()
            .Build();

        var ancestors = registry.AncestorsOf<F97CardioidHalfFixedPointPi2Inheritance>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(HalfAsStructuralFixedPointClaim), ancestors);
        Assert.Contains(typeof(QuarterAsBilinearMaxvalClaim), ancestors);
        Assert.Contains(typeof(NinetyDegreeMirrorMemoryClaim), ancestors);
        Assert.Contains(typeof(PolynomialFoundationClaim), ancestors);
    }

    [Fact]
    public void RegisterF97_FixedPointMagnitudeIsHalf_AtCanonicalAngles()
    {
        var registry = BuildBaseRegistry()
            .RegisterF97CardioidHalfFixedPointPi2Inheritance()
            .Build();

        var f = registry.Get<F97CardioidHalfFixedPointPi2Inheritance>();
        Assert.Equal(0.5, f.FixedPointMagnitude(0.0), precision: 14);
        Assert.Equal(0.5, f.FixedPointMagnitude(Math.PI / 3), precision: 14);
        Assert.Equal(0.5, f.FixedPointMagnitude(Math.PI), precision: 14);
    }

    [Fact]
    public void RegisterF97_AllDriftChecksPass()
    {
        var registry = BuildBaseRegistry()
            .RegisterF97CardioidHalfFixedPointPi2Inheritance()
            .Build();

        var f = registry.Get<F97CardioidHalfFixedPointPi2Inheritance>();
        Assert.True(f.MagnitudeInvariantAroundCardioid());
        Assert.True(f.AlgebraicIdentityHolds(Math.PI / 3));
        Assert.True(f.CuspAgreesWithF95Threshold());
        Assert.True(f.TailAtMinusThreeQuarters());
    }
}
