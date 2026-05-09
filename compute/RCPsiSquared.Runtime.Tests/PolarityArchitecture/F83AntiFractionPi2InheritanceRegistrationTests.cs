using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class F83AntiFractionPi2InheritanceRegistrationTests
{
    private static ClaimRegistryBuilder BuildBaseRegistry() =>
        new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterPi2DyadicLadder();

    [Fact]
    public void RegisterF83_AddsClaim()
    {
        var registry = BuildBaseRegistry()
            .RegisterF83AntiFractionPi2Inheritance()
            .Build();

        Assert.True(registry.Contains<F83AntiFractionPi2Inheritance>());
    }

    [Fact]
    public void RegisterF83_TierIsTier1Derived()
    {
        var registry = BuildBaseRegistry()
            .RegisterF83AntiFractionPi2Inheritance()
            .Build();

        Assert.Equal(Tier.Tier1Derived,
            registry.Get<F83AntiFractionPi2Inheritance>().Tier);
    }

    [Fact]
    public void RegisterF83_AncestorsContainBilinearApex_FirstDirectEdge()
    {
        // Before F83, BilinearApexClaim had 0 descendants. F83 is the first
        // F-formula on the argmax-side of the bilinear-apex pair (per Tom
        // 2026-05-09 mirror-map check).
        var registry = BuildBaseRegistry()
            .RegisterF83AntiFractionPi2Inheritance()
            .Build();

        var ancestors = registry.AncestorsOf<F83AntiFractionPi2Inheritance>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(BilinearApexClaim), ancestors);
        Assert.Contains(typeof(QuarterAsBilinearMaxvalClaim), ancestors);
        Assert.Contains(typeof(Pi2DyadicLadderClaim), ancestors);
    }

    [Fact]
    public void RegisterF83_BilinearApexNowHasDescendants()
    {
        // The mirror-map detector showed BilinearApex with 0 descendants;
        // F83 fills the gap. Cross-registry verification.
        var registry = BuildBaseRegistry()
            .RegisterF83AntiFractionPi2Inheritance()
            .Build();

        var bilinearApexDescendants = registry.DescendantsOf<BilinearApexClaim>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(F83AntiFractionPi2Inheritance), bilinearApexDescendants);
    }

    [Fact]
    public void RegisterF83_AncestorsTransitivelyReachPolynomialFoundation()
    {
        var registry = BuildBaseRegistry()
            .RegisterF83AntiFractionPi2Inheritance()
            .Build();

        var ancestors = registry.AncestorsOf<F83AntiFractionPi2Inheritance>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(PolynomialFoundationClaim), ancestors);
    }

    [Fact]
    public void RegisterF83_MaximumAntiFractionIsOneHalfAcrossRegistry()
    {
        var registry = BuildBaseRegistry()
            .RegisterF83AntiFractionPi2Inheritance()
            .Build();

        Assert.Equal(0.5, registry.Get<F83AntiFractionPi2Inheritance>().MaximumAntiFraction, precision: 14);
    }

    [Fact]
    public void RegisterF83_QuarterCrossoverHoldsAcrossRegistry()
    {
        // Cross-registry verification of the BilinearApex ↔ QuarterAsBilinearMaxval
        // bridge: at r=1/2, anti-fraction = 1/4 = a_3.
        var registry = BuildBaseRegistry()
            .RegisterF83AntiFractionPi2Inheritance()
            .Build();

        Assert.True(registry.Get<F83AntiFractionPi2Inheritance>().QuarterCrossoverHolds());
    }

    [Theory]
    [InlineData(0.0, 0.5)]
    [InlineData(0.5, 0.25)]
    [InlineData(1.0, 1.0 / 6.0)]
    public void RegisterF83_AntiFractionAcrossRegistry(double r, double expected)
    {
        var registry = BuildBaseRegistry()
            .RegisterF83AntiFractionPi2Inheritance()
            .Build();

        Assert.Equal(expected, registry.Get<F83AntiFractionPi2Inheritance>().AntiFraction(r), precision: 12);
    }

    [Theory]
    [InlineData(32.0, 0.0, 3, 1024.0)]
    [InlineData(0.0, 32.0, 3, 2048.0)]
    [InlineData(16.0, 16.0, 3, 1536.0)]
    public void RegisterF83_MNormSquaredVerifiedTableAcrossRegistry(double hOddSq, double hEvenSq, int N, double expected)
    {
        var registry = BuildBaseRegistry()
            .RegisterF83AntiFractionPi2Inheritance()
            .Build();

        Assert.Equal(expected, registry.Get<F83AntiFractionPi2Inheritance>().MNormSquared(hOddSq, hEvenSq, N), precision: 8);
    }
}
