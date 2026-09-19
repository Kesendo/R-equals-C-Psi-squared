using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class F94BornDeviationFourThirdsPi2InheritanceRegistrationTests
{
    private static ClaimRegistry Build() => new ClaimRegistryBuilder()
        .RegisterF94BornDeviationFourThirdsPi2Inheritance()
        .Build();

    [Fact]
    public void RegistrationAddsExactParentlessClaim()
    {
        var registry = Build();
        var claim = registry.Get<F94BornDeviationFourThirdsPi2Inheritance>();
        Assert.Equal(8.0 / 6.0, claim.Coefficient);
        Assert.Empty(registry.AncestorsOf<F94BornDeviationFourThirdsPi2Inheritance>());
        Assert.False(claim is IZ2AxisClaim);
    }
}
