using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class F60GhzBornBelowFoldPi2InheritanceRegistrationTests
{
    private static ClaimRegistryBuilder BuildBaseRegistry() =>
        new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterPi2DyadicLadder();

    [Fact]
    public void RegisterF60_AddsClaim()
    {
        var registry = BuildBaseRegistry()
            .RegisterF60GhzBornBelowFoldPi2Inheritance()
            .Build();

        Assert.True(registry.Contains<F60GhzBornBelowFoldPi2Inheritance>());
    }

    [Fact]
    public void RegisterF60_TierIsTier1Derived()
    {
        var registry = BuildBaseRegistry()
            .RegisterF60GhzBornBelowFoldPi2Inheritance()
            .Build();

        Assert.Equal(Tier.Tier1Derived,
            registry.Get<F60GhzBornBelowFoldPi2Inheritance>().Tier);
    }

    [Fact]
    public void RegisterF60_AncestorsContainPolarityLayerOrigin()
    {
        // F60 sits directly on the 0.5-shift axis; PolarityLayerOriginClaim must be
        // an ancestor (the "1/2 off-diagonal" IS the ±0.5 polarity pair).
        var registry = BuildBaseRegistry()
            .RegisterF60GhzBornBelowFoldPi2Inheritance()
            .Build();

        var ancestors = registry.AncestorsOf<F60GhzBornBelowFoldPi2Inheritance>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(PolarityLayerOriginClaim), ancestors);
        Assert.Contains(typeof(Pi2DyadicLadderClaim), ancestors);
        Assert.Contains(typeof(QuarterAsBilinearMaxvalClaim), ancestors);
    }

    [Fact]
    public void RegisterF60_AncestorsTransitivelyReachPolynomialFoundation()
    {
        var registry = BuildBaseRegistry()
            .RegisterF60GhzBornBelowFoldPi2Inheritance()
            .Build();

        var ancestors = registry.AncestorsOf<F60GhzBornBelowFoldPi2Inheritance>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(PolynomialFoundationClaim), ancestors);
    }

    [Fact]
    public void RegisterF60_OffDiagonalElementIsHalfAcrossRegistry()
    {
        var registry = BuildBaseRegistry()
            .RegisterF60GhzBornBelowFoldPi2Inheritance()
            .Build();

        Assert.Equal(0.5, registry.Get<F60GhzBornBelowFoldPi2Inheritance>().OffDiagonalElement, precision: 14);
    }

    [Fact]
    public void RegisterF60_FoldPositionIsQuarterAcrossRegistry()
    {
        var registry = BuildBaseRegistry()
            .RegisterF60GhzBornBelowFoldPi2Inheritance()
            .Build();

        Assert.Equal(0.25, registry.Get<F60GhzBornBelowFoldPi2Inheritance>().FoldPosition, precision: 14);
    }

    [Theory]
    [InlineData(2, 1.0 / 3.0, false)]   // above fold (Bell+)
    [InlineData(3, 1.0 / 7.0, true)]    // below fold
    [InlineData(5, 1.0 / 31.0, true)]
    public void RegisterF60_FullVerifiedTableAcrossRegistry(int N, double expectedCpsi, bool expectedBelowFold)
    {
        var registry = BuildBaseRegistry()
            .RegisterF60GhzBornBelowFoldPi2Inheritance()
            .Build();

        var f = registry.Get<F60GhzBornBelowFoldPi2Inheritance>();
        Assert.Equal(expectedCpsi, f.CPsiAtZeroForGhz(N), precision: 12);
        Assert.Equal(expectedBelowFold, f.IsBornBelowFold(N));
    }

    [Fact]
    public void RegisterF60_BellPlusAboveFoldAcrossRegistry()
    {
        var registry = BuildBaseRegistry()
            .RegisterF60GhzBornBelowFoldPi2Inheritance()
            .Build();

        Assert.True(registry.Get<F60GhzBornBelowFoldPi2Inheritance>().BellPlusAboveFold());
    }
}
