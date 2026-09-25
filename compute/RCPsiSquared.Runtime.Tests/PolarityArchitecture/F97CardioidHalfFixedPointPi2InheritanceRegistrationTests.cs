using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class F97CardioidHalfFixedPointPi2InheritanceRegistrationTests
{
    private static ClaimRegistry Build() => new ClaimRegistryBuilder()
        .RegisterPi2Family()
        .RegisterF97CardioidHalfFixedPointPi2Inheritance()
        .Build();

    [Fact]
    public void RegistrationWiresTheQuarterParentAndNoPolarityAxis()
    {
        var registry = Build();
        var claim = registry.Get<F97CardioidHalfFixedPointPi2Inheritance>();
        Assert.Same(registry.Get<QuarterAsBilinearMaxvalClaim>(), claim.Quarter);
        Assert.Contains(registry.AncestorsOf<F97CardioidHalfFixedPointPi2Inheritance>(),
            ancestor => ancestor is QuarterAsBilinearMaxvalClaim);
        Assert.True(claim.CuspTakesQuartersValue());
        Assert.Equal(Math.Sqrt(5.0),
            F97CardioidHalfFixedPointPi2Inheritance.MultiplierMagnitude(claim.OtherFixedPoint(Math.PI / 2)), precision: 14);
        Assert.False(claim is IZ2AxisClaim);
    }

    [Fact]
    public void RegistrationWithoutPi2Family_Throws()
    {
        Assert.ThrowsAny<Exception>(() => new ClaimRegistryBuilder()
            .RegisterF97CardioidHalfFixedPointPi2Inheritance()
            .Build());
    }
}
