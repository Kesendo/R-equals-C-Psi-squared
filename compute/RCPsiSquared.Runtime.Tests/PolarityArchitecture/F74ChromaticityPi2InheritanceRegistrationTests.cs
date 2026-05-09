using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class F74ChromaticityPi2InheritanceRegistrationTests
{
    private static ClaimRegistryBuilder BuildBaseRegistry() =>
        new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterPi2DyadicLadder();

    [Fact]
    public void RegisterF74_AddsClaim()
    {
        var registry = BuildBaseRegistry()
            .RegisterF74ChromaticityPi2Inheritance()
            .Build();
        Assert.True(registry.Contains<F74ChromaticityPi2Inheritance>());
    }

    [Fact]
    public void RegisterF74_TierIsTier1Derived()
    {
        var registry = BuildBaseRegistry()
            .RegisterF74ChromaticityPi2Inheritance()
            .Build();
        Assert.Equal(Tier.Tier1Derived, registry.Get<F74ChromaticityPi2Inheritance>().Tier);
    }

    [Fact]
    public void RegisterF74_AncestorsContainPi2DyadicLadder()
    {
        var registry = BuildBaseRegistry()
            .RegisterF74ChromaticityPi2Inheritance()
            .Build();
        var ancestors = registry.AncestorsOf<F74ChromaticityPi2Inheritance>()
            .Select(c => c.GetType()).ToHashSet();
        Assert.Contains(typeof(Pi2DyadicLadderClaim), ancestors);
    }

    [Theory]
    [InlineData(2, 5, 3)]
    [InlineData(0, 5, 1)]
    public void RegisterF74_ChromaticityAcrossRegistry(int n, int N, int expected)
    {
        var registry = BuildBaseRegistry()
            .RegisterF74ChromaticityPi2Inheritance()
            .Build();
        Assert.Equal(expected, registry.Get<F74ChromaticityPi2Inheritance>().Chromaticity(n, N));
    }

    [Fact]
    public void RegisterF74_MirrorPalindromicityHoldsAcrossRegistry()
    {
        var registry = BuildBaseRegistry()
            .RegisterF74ChromaticityPi2Inheritance()
            .Build();
        var f74 = registry.Get<F74ChromaticityPi2Inheritance>();
        for (int N = 2; N <= 10; N++)
            Assert.True(f74.MirrorPalindromicityHolds(N));
    }

    [Fact]
    public void RegisterF74_WithoutPi2DyadicLadder_Throws()
    {
        Assert.Throws<InvariantViolationException>(() =>
            new ClaimRegistryBuilder()
                .RegisterPi2Family()
                .RegisterF74ChromaticityPi2Inheritance()
                .Build());
    }
}
