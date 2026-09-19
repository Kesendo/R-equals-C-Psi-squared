using System.Numerics;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class F97CardioidHalfFixedPointPi2InheritanceRegistrationTests
{
    private static ClaimRegistry Build() => new ClaimRegistryBuilder()
        .RegisterF97CardioidHalfFixedPointPi2Inheritance()
        .Build();

    [Fact]
    public void RegistrationAddsExactParentlessClaim()
    {
        var registry = Build();
        var claim = registry.Get<F97CardioidHalfFixedPointPi2Inheritance>();
        Assert.True(claim.SelectedMultiplierIsMarginal(Math.PI / 2));
        Assert.Equal(Math.Sqrt(5.0), F97CardioidHalfFixedPointPi2Inheritance.MultiplierMagnitude(claim.OtherFixedPoint(Math.PI / 2)), precision: 14);
        Assert.True(claim.OriginRootsAreDistinct());
        Assert.Empty(registry.AncestorsOf<F97CardioidHalfFixedPointPi2Inheritance>());
        Assert.False(claim is IZ2AxisClaim);
    }
}
