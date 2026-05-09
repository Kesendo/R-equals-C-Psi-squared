using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class F38Pi2InvolutionPi2InheritanceRegistrationTests
{
    private static ClaimRegistryBuilder BuildBaseRegistry() =>
        new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterPi2DyadicLadder()
            .RegisterF88PopcountCoherence()
            .RegisterF88StaticDyadicAnchor()
            .RegisterPi2OperatorSpaceMirror()
            .RegisterPi2I4MemoryLoop();

    [Fact]
    public void RegisterF38_AddsClaim()
    {
        var registry = BuildBaseRegistry()
            .RegisterF38Pi2InvolutionPi2Inheritance()
            .Build();

        Assert.True(registry.Contains<F38Pi2InvolutionPi2Inheritance>());
    }

    [Fact]
    public void RegisterF38_TierIsTier1Derived()
    {
        var registry = BuildBaseRegistry()
            .RegisterF38Pi2InvolutionPi2Inheritance()
            .Build();

        Assert.Equal(Tier.Tier1Derived,
            registry.Get<F38Pi2InvolutionPi2Inheritance>().Tier);
    }

    [Fact]
    public void RegisterF38_AncestorsContainAllThreePi2Anchors()
    {
        var registry = BuildBaseRegistry()
            .RegisterF38Pi2InvolutionPi2Inheritance()
            .Build();

        var ancestors = registry.AncestorsOf<F38Pi2InvolutionPi2Inheritance>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(Pi2DyadicLadderClaim), ancestors);
        Assert.Contains(typeof(Pi2OperatorSpaceMirrorClaim), ancestors);
        Assert.Contains(typeof(Pi2I4MemoryLoopClaim), ancestors);
    }

    [Fact]
    public void RegisterF38_AncestorsTransitivelyReachHalfAsStructuralFixedPoint()
    {
        // The "1/2 half-half balance" is anchored transitively via the Pi2DyadicLadder
        // (a_2 = 1/2) which has HalfAsStructuralFixedPoint as one of its known anchors.
        var registry = BuildBaseRegistry()
            .RegisterF38Pi2InvolutionPi2Inheritance()
            .Build();

        var ancestors = registry.AncestorsOf<F38Pi2InvolutionPi2Inheritance>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(HalfAsStructuralFixedPointClaim), ancestors);
    }

    [Theory]
    [InlineData(1, 4.0, 2L)]
    [InlineData(2, 16.0, 8L)]
    [InlineData(3, 64.0, 32L)]
    [InlineData(6, 4096.0, 2048L)]
    public void RegisterF38_DimensionsAcrossRegistry(int N, double fullDim, long eigenspaceDim)
    {
        var registry = BuildBaseRegistry()
            .RegisterF38Pi2InvolutionPi2Inheritance()
            .Build();

        var f = registry.Get<F38Pi2InvolutionPi2Inheritance>();
        Assert.Equal((long)fullDim, f.FullOperatorSpaceDimension(N));
        Assert.Equal(eigenspaceDim, f.EigenspaceDimension(N));
        Assert.Equal(fullDim, f.MirrorPinnedFullDimension(N), precision: 12);
        Assert.Equal((double)eigenspaceDim, f.EigenspaceDimensionViaLadder(N), precision: 10);
    }

    [Theory]
    [InlineData(1)]
    [InlineData(2)]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    public void RegisterF38_LadderFactorisationHoldsAcrossN(int N)
    {
        var registry = BuildBaseRegistry()
            .RegisterF38Pi2InvolutionPi2Inheritance()
            .Build();

        Assert.True(registry.Get<F38Pi2InvolutionPi2Inheritance>().LadderFactorisationHolds(N));
    }

    [Fact]
    public void RegisterF38_MemoryLoopClosureIsLive()
    {
        var registry = BuildBaseRegistry()
            .RegisterF38Pi2InvolutionPi2Inheritance()
            .Build();

        var f = registry.Get<F38Pi2InvolutionPi2Inheritance>();
        Assert.Equal(4, f.CyclicOrder);
        Assert.True(f.MemoryLoopClosesAtFour());
        Assert.Equal(new[] { +1, -1, +1, -1 }, f.Pi2EigenvaluesFromMemoryLoop());
    }

    [Fact]
    public void RegisterF38_WithoutMemoryLoop_Throws()
    {
        Assert.Throws<InvariantViolationException>(() =>
            new ClaimRegistryBuilder()
                .RegisterPi2Family()
                .RegisterPi2DyadicLadder()
                .RegisterF88PopcountCoherence()
                .RegisterF88StaticDyadicAnchor()
                .RegisterPi2OperatorSpaceMirror()
                // Missing: Pi2I4MemoryLoop
                .RegisterF38Pi2InvolutionPi2Inheritance()
                .Build());
    }

    [Fact]
    public void RegisterF38_WithoutOperatorSpaceMirror_Throws()
    {
        Assert.Throws<InvariantViolationException>(() =>
            new ClaimRegistryBuilder()
                .RegisterPi2Family()
                .RegisterPi2DyadicLadder()
                .RegisterPi2I4MemoryLoop()
                // Missing: F88* + Pi2OperatorSpaceMirror
                .RegisterF38Pi2InvolutionPi2Inheritance()
                .Build());
    }
}
