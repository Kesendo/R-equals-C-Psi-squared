using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class F96BornSubdominantSlopesRegistrationTests
{
    private static ClaimRegistry Build() => new ClaimRegistryBuilder()
        .RegisterF96BornSubdominantSlopes()
        .Build();

    [Fact]
    public void RegistrationAddsExactParentlessClaim()
    {
        var registry = Build();
        var claim = registry.Get<F96BornSubdominantSlopes>();
        Assert.Equal(Tier.Tier1Derived, claim.Tier);
        Assert.Equal(-16.0 / 9.0, claim.SlopeSingleFlipped, precision: 15);
        Assert.Equal(-8.0 / 3.0, claim.SlopeDoubleFlipped, precision: 15);
        Assert.Empty(registry.AncestorsOf<F96BornSubdominantSlopes>());
        Assert.False(claim is IZ2AxisClaim);
    }
}
