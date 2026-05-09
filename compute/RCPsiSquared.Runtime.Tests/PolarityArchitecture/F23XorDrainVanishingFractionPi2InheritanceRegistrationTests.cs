using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class F23XorDrainVanishingFractionPi2InheritanceRegistrationTests
{
    private static ClaimRegistryBuilder BuildBaseRegistry() =>
        new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterPi2DyadicLadder()
            .RegisterF88PopcountCoherence()
            .RegisterF88StaticDyadicAnchor()
            .RegisterPi2OperatorSpaceMirror();

    [Fact]
    public void RegisterF23_AddsClaim()
    {
        var registry = BuildBaseRegistry()
            .RegisterF23XorDrainVanishingFractionPi2Inheritance()
            .Build();

        Assert.True(registry.Contains<F23XorDrainVanishingFractionPi2Inheritance>());
    }

    [Fact]
    public void RegisterF23_TierIsTier1Derived()
    {
        var registry = BuildBaseRegistry()
            .RegisterF23XorDrainVanishingFractionPi2Inheritance()
            .Build();

        Assert.Equal(Tier.Tier1Derived,
            registry.Get<F23XorDrainVanishingFractionPi2Inheritance>().Tier);
    }

    [Fact]
    public void RegisterF23_AncestorsContainBothParents()
    {
        var registry = BuildBaseRegistry()
            .RegisterF23XorDrainVanishingFractionPi2Inheritance()
            .Build();

        var ancestors = registry.AncestorsOf<F23XorDrainVanishingFractionPi2Inheritance>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(Pi2DyadicLadderClaim), ancestors);
        Assert.Contains(typeof(Pi2OperatorSpaceMirrorClaim), ancestors);
    }

    [Theory]
    [InlineData(3, 0.0625)]
    [InlineData(5, 6.0 / 1024.0)]
    public void RegisterF23_XorDrainFractionAcrossRegistry(int N, double expected)
    {
        var registry = BuildBaseRegistry()
            .RegisterF23XorDrainVanishingFractionPi2Inheritance()
            .Build();

        Assert.Equal(expected, registry.Get<F23XorDrainVanishingFractionPi2Inheritance>().XorDrainFraction(N), precision: 12);
    }

    [Theory]
    [InlineData(1)]
    [InlineData(3)]
    [InlineData(6)]
    public void RegisterF23_MatchesMirrorTableAcrossRegistry(int N)
    {
        var registry = BuildBaseRegistry()
            .RegisterF23XorDrainVanishingFractionPi2Inheritance()
            .Build();

        Assert.True(registry.Get<F23XorDrainVanishingFractionPi2Inheritance>().MatchesMirrorTable(N));
    }

    [Fact]
    public void RegisterF23_MacroscopicallyNegligibleAtLargeN()
    {
        var registry = BuildBaseRegistry()
            .RegisterF23XorDrainVanishingFractionPi2Inheritance()
            .Build();

        var f23 = registry.Get<F23XorDrainVanishingFractionPi2Inheritance>();
        Assert.False(f23.IsMacroscopicallyNegligible(N: 5));
        Assert.True(f23.IsMacroscopicallyNegligible(N: 15));
    }

    [Fact]
    public void RegisterF23_WithoutMirror_Throws()
    {
        Assert.Throws<InvariantViolationException>(() =>
            new ClaimRegistryBuilder()
                .RegisterPi2Family()
                .RegisterPi2DyadicLadder()
                .RegisterF88PopcountCoherence()
                .RegisterF88StaticDyadicAnchor()
                // Missing: RegisterPi2OperatorSpaceMirror
                .RegisterF23XorDrainVanishingFractionPi2Inheritance()
                .Build());
    }
}
