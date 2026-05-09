using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class F63LCommutesPi2Pi2InheritanceRegistrationTests
{
    private static ClaimRegistryBuilder BuildBaseRegistry() =>
        new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterPi2DyadicLadder()
            .RegisterF88PopcountCoherence()
            .RegisterF88StaticDyadicAnchor()
            .RegisterPi2OperatorSpaceMirror()
            .RegisterPi2I4MemoryLoop()
            .RegisterF38Pi2InvolutionPi2Inheritance();

    [Fact]
    public void RegisterF63_AddsClaim()
    {
        var registry = BuildBaseRegistry()
            .RegisterF63LCommutesPi2Pi2Inheritance()
            .Build();

        Assert.True(registry.Contains<F63LCommutesPi2Pi2Inheritance>());
    }

    [Fact]
    public void RegisterF63_TierIsTier1Derived()
    {
        var registry = BuildBaseRegistry()
            .RegisterF63LCommutesPi2Pi2Inheritance()
            .Build();

        Assert.Equal(Tier.Tier1Derived,
            registry.Get<F63LCommutesPi2Pi2Inheritance>().Tier);
    }

    [Fact]
    public void RegisterF63_AncestorsContainF38AndLadder()
    {
        var registry = BuildBaseRegistry()
            .RegisterF63LCommutesPi2Pi2Inheritance()
            .Build();

        var ancestors = registry.AncestorsOf<F63LCommutesPi2Pi2Inheritance>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(F38Pi2InvolutionPi2Inheritance), ancestors);
        Assert.Contains(typeof(Pi2DyadicLadderClaim), ancestors);
    }

    [Fact]
    public void RegisterF63_AncestorsTransitivelyReachPolynomialFoundation()
    {
        var registry = BuildBaseRegistry()
            .RegisterF63LCommutesPi2Pi2Inheritance()
            .Build();

        var ancestors = registry.AncestorsOf<F63LCommutesPi2Pi2Inheritance>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(PolynomialFoundationClaim), ancestors);
    }

    [Fact]
    public void RegisterF63_BlockCountIsFour()
    {
        var registry = BuildBaseRegistry()
            .RegisterF63LCommutesPi2Pi2Inheritance()
            .Build();

        Assert.Equal(4.0, registry.Get<F63LCommutesPi2Pi2Inheritance>().BlockCount, precision: 14);
    }

    [Theory]
    [InlineData(2, 4.0)]
    [InlineData(3, 16.0)]
    [InlineData(5, 256.0)]
    public void RegisterF63_PerBlockDimensionAcrossRegistry(int N, double expected)
    {
        var registry = BuildBaseRegistry()
            .RegisterF63LCommutesPi2Pi2Inheritance()
            .Build();

        Assert.Equal(expected, registry.Get<F63LCommutesPi2Pi2Inheritance>().PerBlockDimension(N), precision: 12);
    }

    [Theory]
    [InlineData(2, 4L, 6L)]
    [InlineData(3, 28L, 28L)]
    [InlineData(4, 122L, 124L)]
    [InlineData(5, 506L, 506L)]
    public void RegisterF63_MirrorTableAcrossRegistry(int N, long expectedEven, long expectedOdd)
    {
        var registry = BuildBaseRegistry()
            .RegisterF63LCommutesPi2Pi2Inheritance()
            .Build();

        var f = registry.Get<F63LCommutesPi2Pi2Inheritance>();
        Assert.Equal(expectedEven, f.MirrorEvenSector(N));
        Assert.Equal(expectedOdd, f.MirrorOddSector(N));
    }

    [Theory]
    [InlineData(2)]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    public void RegisterF63_MatchesF66EndpointMultiplicityAcrossN(int N)
    {
        var registry = BuildBaseRegistry()
            .RegisterF63LCommutesPi2Pi2Inheritance()
            .Build();

        Assert.True(registry.Get<F63LCommutesPi2Pi2Inheritance>().MatchesF66EndpointMultiplicity(N));
    }

    [Theory]
    [InlineData(1)]
    [InlineData(2)]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    public void RegisterF63_PerBlockDimensionAgreesWithF38(int N)
    {
        var registry = BuildBaseRegistry()
            .RegisterF63LCommutesPi2Pi2Inheritance()
            .Build();

        Assert.True(registry.Get<F63LCommutesPi2Pi2Inheritance>().PerBlockDimensionAgreesWithF38(N));
    }

    [Fact]
    public void RegisterF63_WithoutF38_Throws()
    {
        Assert.Throws<InvariantViolationException>(() =>
            new ClaimRegistryBuilder()
                .RegisterPi2Family()
                .RegisterPi2DyadicLadder()
                // Missing: F88* + Pi2OperatorSpaceMirror + Pi2I4MemoryLoop + F38
                .RegisterF63LCommutesPi2Pi2Inheritance()
                .Build());
    }
}
