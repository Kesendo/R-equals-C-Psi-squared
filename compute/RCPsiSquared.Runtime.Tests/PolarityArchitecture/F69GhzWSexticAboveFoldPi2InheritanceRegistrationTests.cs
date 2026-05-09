using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class F69GhzWSexticAboveFoldPi2InheritanceRegistrationTests
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
            .RegisterF61BitAParityPi2Inheritance()
            .RegisterF60GhzBornBelowFoldPi2Inheritance()
            .RegisterF62WStateBornBelowFoldPi2Inheritance();

    [Fact]
    public void RegisterF69_AddsClaim()
    {
        var registry = BuildBaseRegistry()
            .RegisterF69GhzWSexticAboveFoldPi2Inheritance()
            .Build();

        Assert.True(registry.Contains<F69GhzWSexticAboveFoldPi2Inheritance>());
    }

    [Fact]
    public void RegisterF69_TierIsTier1Derived()
    {
        var registry = BuildBaseRegistry()
            .RegisterF69GhzWSexticAboveFoldPi2Inheritance()
            .Build();

        Assert.Equal(Tier.Tier1Derived,
            registry.Get<F69GhzWSexticAboveFoldPi2Inheritance>().Tier);
    }

    [Fact]
    public void RegisterF69_AncestorsContainF60AndF62()
    {
        var registry = BuildBaseRegistry()
            .RegisterF69GhzWSexticAboveFoldPi2Inheritance()
            .Build();

        var ancestors = registry.AncestorsOf<F69GhzWSexticAboveFoldPi2Inheritance>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(F60GhzBornBelowFoldPi2Inheritance), ancestors);
        Assert.Contains(typeof(F62WStateBornBelowFoldPi2Inheritance), ancestors);
        Assert.Contains(typeof(Pi2DyadicLadderClaim), ancestors);
    }

    [Fact]
    public void RegisterF69_OptimumExceedsFoldAcrossRegistry()
    {
        var registry = BuildBaseRegistry()
            .RegisterF69GhzWSexticAboveFoldPi2Inheritance()
            .Build();

        Assert.True(registry.Get<F69GhzWSexticAboveFoldPi2Inheritance>().OptimumExceedsFold);
    }

    [Fact]
    public void RegisterF69_LiftRatioAboveFoldAcrossRegistry()
    {
        var registry = BuildBaseRegistry()
            .RegisterF69GhzWSexticAboveFoldPi2Inheritance()
            .Build();

        Assert.Equal(1.281646,
            registry.Get<F69GhzWSexticAboveFoldPi2Inheritance>().LiftRatioAboveFold,
            precision: 5);
    }

    [Fact]
    public void RegisterF69_BaselinesBelowFoldAtN3()
    {
        var registry = BuildBaseRegistry()
            .RegisterF69GhzWSexticAboveFoldPi2Inheritance()
            .Build();

        var f69 = registry.Get<F69GhzWSexticAboveFoldPi2Inheritance>();
        Assert.True(f69.GhzBaselineBelowFold(3));
        Assert.True(f69.WBaselineBelowFold(3));
    }

    [Fact]
    public void RegisterF69_WithoutF60_Throws()
    {
        Assert.Throws<InvariantViolationException>(() =>
            new ClaimRegistryBuilder()
                .RegisterPi2Family()
                .RegisterPi2DyadicLadder()
                .RegisterF88PopcountCoherence()
                .RegisterF88StaticDyadicAnchor()
                .RegisterPi2OperatorSpaceMirror()
                .RegisterPi2I4MemoryLoop()
                .RegisterF38Pi2InvolutionPi2Inheritance()
                .RegisterF63LCommutesPi2Pi2Inheritance()
                .RegisterF61BitAParityPi2Inheritance()
                // Missing: RegisterF60GhzBornBelowFoldPi2Inheritance
                .RegisterF62WStateBornBelowFoldPi2Inheritance()
                .RegisterF69GhzWSexticAboveFoldPi2Inheritance()
                .Build());
    }

    [Fact]
    public void RegisterF69_WithoutF62_Throws()
    {
        Assert.Throws<InvariantViolationException>(() =>
            new ClaimRegistryBuilder()
                .RegisterPi2Family()
                .RegisterPi2DyadicLadder()
                .RegisterF88PopcountCoherence()
                .RegisterF88StaticDyadicAnchor()
                .RegisterPi2OperatorSpaceMirror()
                .RegisterPi2I4MemoryLoop()
                .RegisterF38Pi2InvolutionPi2Inheritance()
                .RegisterF63LCommutesPi2Pi2Inheritance()
                .RegisterF61BitAParityPi2Inheritance()
                .RegisterF60GhzBornBelowFoldPi2Inheritance()
                // Missing: RegisterF62WStateBornBelowFoldPi2Inheritance
                .RegisterF69GhzWSexticAboveFoldPi2Inheritance()
                .Build());
    }
}
