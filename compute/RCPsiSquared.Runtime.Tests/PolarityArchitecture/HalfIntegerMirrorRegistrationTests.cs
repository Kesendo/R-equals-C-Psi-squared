using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class HalfIntegerMirrorRegistrationTests
{
    [Fact]
    public void RegisterHalfIntegerMirror_AddsClaim()
    {
        var registry = new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterHalfIntegerMirror(N: 5)
            .Build();

        Assert.True(registry.Contains<HalfIntegerMirrorClaim>());
        Assert.Equal(9, registry.All().Count()); // 8 Pi2 + 1
    }

    [Fact]
    public void RegisterHalfIntegerMirror_OddN_HalfIntegerRegime()
    {
        var registry = new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterHalfIntegerMirror(N: 5)
            .Build();

        var claim = registry.Get<HalfIntegerMirrorClaim>();
        Assert.True(claim.IsHalfIntegerRegime);
        Assert.Equal(2.5, claim.WXY);
    }

    [Fact]
    public void RegisterHalfIntegerMirror_EvenN_IntegerRegime()
    {
        var registry = new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterHalfIntegerMirror(N: 4)
            .Build();

        var claim = registry.Get<HalfIntegerMirrorClaim>();
        Assert.False(claim.IsHalfIntegerRegime);
        Assert.Equal(2.0, claim.WXY);
    }

    [Fact]
    public void RegisterHalfIntegerMirror_AncestorContainsQubitDimensional()
    {
        var registry = new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterHalfIntegerMirror(N: 5)
            .Build();

        var ancestors = registry.AncestorsOf<HalfIntegerMirrorClaim>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(QubitDimensionalAnchorClaim), ancestors);
        Assert.Contains(typeof(PolynomialFoundationClaim), ancestors);
    }

    [Fact]
    public void RegisterHalfIntegerMirror_Tier1Derived()
    {
        var registry = new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterHalfIntegerMirror(N: 5)
            .Build();

        Assert.Equal(Tier.Tier1Derived, registry.Get<HalfIntegerMirrorClaim>().Tier);
    }
}
