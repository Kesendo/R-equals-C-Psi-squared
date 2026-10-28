using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class F94BornDeviationFourThirdsRegistrationTests
{
    private static ClaimRegistry Build() => new ClaimRegistryBuilder()
        .RegisterF94BornDeviationFourThirds()
        .Build();

    [Fact]
    public void RegistrationAddsExactParentlessClaim()
    {
        var registry = Build();
        var claim = registry.Get<F94BornDeviationFourThirds>();
        Assert.Equal(8.0 / 6.0, claim.Coefficient);
        Assert.Empty(registry.AncestorsOf<F94BornDeviationFourThirds>());
        Assert.False(claim is IZ2AxisClaim);
    }
}
