using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class F61BitAParityPi2InheritanceRegistrationTests
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
            .RegisterF63LCommutesPi2Pi2Inheritance();

    [Fact]
    public void RegisterF61_AddsClaim()
    {
        var registry = BuildBaseRegistry()
            .RegisterF61BitAParityPi2Inheritance()
            .Build();

        Assert.True(registry.Contains<F61BitAParityPi2Inheritance>());
    }

    [Fact]
    public void RegisterF61_TierIsTier1Derived()
    {
        var registry = BuildBaseRegistry()
            .RegisterF61BitAParityPi2Inheritance()
            .Build();

        Assert.Equal(Tier.Tier1Derived,
            registry.Get<F61BitAParityPi2Inheritance>().Tier);
    }

    [Fact]
    public void RegisterF61_AncestorsContainF63AndLadder()
    {
        var registry = BuildBaseRegistry()
            .RegisterF61BitAParityPi2Inheritance()
            .Build();

        var ancestors = registry.AncestorsOf<F61BitAParityPi2Inheritance>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(F63LCommutesPi2Pi2Inheritance), ancestors);
        Assert.Contains(typeof(Pi2DyadicLadderClaim), ancestors);
    }

    [Fact]
    public void RegisterF61_AncestorsTransitivelyReachF38AndPolynomialFoundation()
    {
        // F61 → F63 → F38 → ... → PolynomialFoundation
        var registry = BuildBaseRegistry()
            .RegisterF61BitAParityPi2Inheritance()
            .Build();

        var ancestors = registry.AncestorsOf<F61BitAParityPi2Inheritance>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(F38Pi2InvolutionPi2Inheritance), ancestors);
        Assert.Contains(typeof(PolynomialFoundationClaim), ancestors);
    }

    [Fact]
    public void RegisterF61_BlockCountAgreesWithF63AcrossRegistry()
    {
        var registry = BuildBaseRegistry()
            .RegisterF61BitAParityPi2Inheritance()
            .Build();

        Assert.True(registry.Get<F61BitAParityPi2Inheritance>().BlockCountAgreesWithF63());
    }

    [Theory]
    [InlineData(2)]
    [InlineData(3)]
    [InlineData(5)]
    public void RegisterF61_PerBlockDimensionAgreesWithF63AcrossN(int N)
    {
        var registry = BuildBaseRegistry()
            .RegisterF61BitAParityPi2Inheritance()
            .Build();

        Assert.True(registry.Get<F61BitAParityPi2Inheritance>().PerBlockDimensionAgreesWithF63(N));
    }

    [Theory]
    [InlineData(2, 4.0)]
    [InlineData(3, 16.0)]
    [InlineData(5, 256.0)]
    public void RegisterF61_PerBlockDimensionAcrossRegistry(int N, double expected)
    {
        var registry = BuildBaseRegistry()
            .RegisterF61BitAParityPi2Inheritance()
            .Build();

        Assert.Equal(expected, registry.Get<F61BitAParityPi2Inheritance>().PerBlockDimension(N), precision: 12);
    }

    [Fact]
    public void RegisterF61_Z2AxisIsBitA()
    {
        var registry = BuildBaseRegistry()
            .RegisterF61BitAParityPi2Inheritance()
            .Build();

        Assert.Equal("bit_a (n_XY)", registry.Get<F61BitAParityPi2Inheritance>().Z2Axis);
    }

    [Fact]
    public void RegisterF61_WithoutF63_Throws()
    {
        Assert.Throws<InvariantViolationException>(() =>
            new ClaimRegistryBuilder()
                .RegisterPi2Family()
                .RegisterPi2DyadicLadder()
                // Missing: F88* + Pi2OperatorSpaceMirror + Pi2I4MemoryLoop + F38 + F63
                .RegisterF61BitAParityPi2Inheritance()
                .Build());
    }
}
