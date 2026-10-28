using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class F97CardioidHalfFixedPointRegistrationTests
{
    private static ClaimRegistry Build() => new ClaimRegistryBuilder()
        .RegisterPi2Family()
        .RegisterF97CardioidHalfFixedPoint()
        .Build();

    [Fact]
    public void RegistrationWiresTheQuarterParentAndNoPolarityAxis()
    {
        var registry = Build();
        var claim = registry.Get<F97CardioidHalfFixedPoint>();
        Assert.Same(registry.Get<QuarterAsBilinearMaxvalClaim>(), claim.Quarter);
        Assert.Contains(registry.AncestorsOf<F97CardioidHalfFixedPoint>(),
            ancestor => ancestor is QuarterAsBilinearMaxvalClaim);
        Assert.True(claim.CuspTakesQuartersValue());
        Assert.Equal(Math.Sqrt(5.0),
            F97CardioidHalfFixedPoint.MultiplierMagnitude(claim.OtherFixedPoint(Math.PI / 2)), precision: 14);
        Assert.False(claim is IZ2AxisClaim);
    }

    [Fact]
    public void RegistrationWithoutPi2Family_Throws()
    {
        Assert.ThrowsAny<Exception>(() => new ClaimRegistryBuilder()
            .RegisterF97CardioidHalfFixedPoint()
            .Build());
    }
}
