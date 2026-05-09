using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class F62WStateBornBelowFoldPi2InheritanceRegistrationTests
{
    private static ClaimRegistryBuilder BuildBaseRegistry() =>
        new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterPi2DyadicLadder()
            .RegisterF88PopcountCoherence()
            .RegisterF88StaticDyadicAnchor()
            .RegisterPi2OperatorSpaceMirror()
            .RegisterPi2I4MemoryLoop()
            .RegisterF38Pi2InvolutionPi2Inheritance()
            .RegisterF63LCommutesPi2Pi2Inheritance()
            .RegisterF61BitAParityPi2Inheritance();

    [Fact]
    public void RegisterF62_AddsClaim()
    {
        var registry = BuildBaseRegistry()
            .RegisterF62WStateBornBelowFoldPi2Inheritance()
            .Build();

        Assert.True(registry.Contains<F62WStateBornBelowFoldPi2Inheritance>());
    }

    [Fact]
    public void RegisterF62_TierIsTier1Derived()
    {
        var registry = BuildBaseRegistry()
            .RegisterF62WStateBornBelowFoldPi2Inheritance()
            .Build();

        Assert.Equal(Tier.Tier1Derived,
            registry.Get<F62WStateBornBelowFoldPi2Inheritance>().Tier);
    }

    [Fact]
    public void RegisterF62_AncestorsContainQuarterAndF61()
    {
        var registry = BuildBaseRegistry()
            .RegisterF62WStateBornBelowFoldPi2Inheritance()
            .Build();

        var ancestors = registry.AncestorsOf<F62WStateBornBelowFoldPi2Inheritance>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(Pi2DyadicLadderClaim), ancestors);
        Assert.Contains(typeof(QuarterAsBilinearMaxvalClaim), ancestors);
        Assert.Contains(typeof(F61BitAParityPi2Inheritance), ancestors);
    }

    [Fact]
    public void RegisterF62_AncestorsTransitivelyReachF63AndPolynomialFoundation()
    {
        // F62 → F61 → F63 → F38 → ... → PolynomialFoundation
        var registry = BuildBaseRegistry()
            .RegisterF62WStateBornBelowFoldPi2Inheritance()
            .Build();

        var ancestors = registry.AncestorsOf<F62WStateBornBelowFoldPi2Inheritance>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(F63LCommutesPi2Pi2Inheritance), ancestors);
        Assert.Contains(typeof(F38Pi2InvolutionPi2Inheritance), ancestors);
        Assert.Contains(typeof(PolynomialFoundationClaim), ancestors);
    }

    [Theory]
    [InlineData(2, 1.0 / 3.0, false)]   // Bell+ regime
    [InlineData(3, 10.0 / 81.0, true)]
    [InlineData(5, 26.0 / 375.0, true)]
    public void RegisterF62_FullVerifiedTableAcrossRegistry(int N, double expectedCpsi, bool expectedBelowFold)
    {
        var registry = BuildBaseRegistry()
            .RegisterF62WStateBornBelowFoldPi2Inheritance()
            .Build();

        var f = registry.Get<F62WStateBornBelowFoldPi2Inheritance>();
        Assert.Equal(expectedCpsi, f.CPsiAtZeroForWState(N), precision: 12);
        Assert.Equal(expectedBelowFold, f.IsBornBelowFold(N));
    }

    [Fact]
    public void RegisterF62_FoldPositionIsQuarterAcrossRegistry()
    {
        var registry = BuildBaseRegistry()
            .RegisterF62WStateBornBelowFoldPi2Inheritance()
            .Build();

        Assert.Equal(0.25, registry.Get<F62WStateBornBelowFoldPi2Inheritance>().FoldPosition, precision: 14);
    }

    [Fact]
    public void RegisterF62_F60SiblingPairConsistentAcrossRegistry()
    {
        var registry = BuildBaseRegistry()
            .RegisterF60GhzBornBelowFoldPi2Inheritance()
            .RegisterF62WStateBornBelowFoldPi2Inheritance()
            .Build();

        var f60 = registry.Get<F60GhzBornBelowFoldPi2Inheritance>();
        var f62 = registry.Get<F62WStateBornBelowFoldPi2Inheritance>();

        Assert.Equal(f60.SmallestNBelowFold, f62.SmallestNBelowFold);
        Assert.Equal(f60.FoldPosition, f62.FoldPosition, precision: 14);
        Assert.True(f60.BellPlusAboveFold());
        Assert.True(f62.BellPlusAboveFold());
    }
}
