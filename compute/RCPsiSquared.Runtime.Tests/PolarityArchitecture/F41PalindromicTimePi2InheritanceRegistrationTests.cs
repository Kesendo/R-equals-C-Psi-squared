using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.F1Family;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class F41PalindromicTimePi2InheritanceRegistrationTests
{
    private static ChainSystem DefaultChain() =>
        new ChainSystem(N: 5, J: 1.0, GammaZero: 0.05,
            HType: HamiltonianType.Heisenberg, Topology: TopologyKind.Chain);

    private static ClaimRegistryBuilder BuildBaseRegistry() =>
        new ClaimRegistryBuilder()
            .RegisterF1Family(DefaultChain())
            .RegisterPi2Family()
            .RegisterPi2DyadicLadder()
            .RegisterPi2I4MemoryLoop()
            .RegisterF1Pi2Inheritance();

    [Fact]
    public void RegisterF41_AddsClaim()
    {
        var registry = BuildBaseRegistry()
            .RegisterF41PalindromicTimePi2Inheritance()
            .Build();

        Assert.True(registry.Contains<F41PalindromicTimePi2Inheritance>());
    }

    [Fact]
    public void RegisterF41_TierIsTier1Derived()
    {
        var registry = BuildBaseRegistry()
            .RegisterF41PalindromicTimePi2Inheritance()
            .Build();

        Assert.Equal(Tier.Tier1Derived,
            registry.Get<F41PalindromicTimePi2Inheritance>().Tier);
    }

    [Fact]
    public void RegisterF41_AncestorsContainF1AndPi2DyadicLadder()
    {
        var registry = BuildBaseRegistry()
            .RegisterF41PalindromicTimePi2Inheritance()
            .Build();

        var ancestors = registry.AncestorsOf<F41PalindromicTimePi2Inheritance>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(F1Pi2Inheritance), ancestors);
        Assert.Contains(typeof(Pi2DyadicLadderClaim), ancestors);
        Assert.Contains(typeof(Pi2I4MemoryLoopClaim), ancestors);
    }

    [Theory]
    [InlineData(2, 1.0, Math.PI)]
    [InlineData(3, 1.0, 2.0 * Math.PI)]
    public void RegisterF41_PalindromicTimeAcrossRegistry(int N, double J, double expected)
    {
        var registry = BuildBaseRegistry()
            .RegisterF41PalindromicTimePi2Inheritance()
            .Build();

        Assert.Equal(expected, registry.Get<F41PalindromicTimePi2Inheritance>().PalindromicTime(N, J), precision: 12);
    }

    [Fact]
    public void RegisterF41_PeriodFrequencyProductIsTwoPiAcrossRegistry()
    {
        var registry = BuildBaseRegistry()
            .RegisterF41PalindromicTimePi2Inheritance()
            .Build();

        var f41 = registry.Get<F41PalindromicTimePi2Inheritance>();
        Assert.True(f41.PeriodFrequencyProductIsTwoPi(N: 3, J: 1.0));
        Assert.True(f41.PeriodFrequencyProductIsTwoPi(N: 5, J: 0.5));
    }

    [Fact]
    public void RegisterF41_WithoutF1_Throws()
    {
        Assert.Throws<InvariantViolationException>(() =>
            new ClaimRegistryBuilder()
                .RegisterF1Family(DefaultChain())
                .RegisterPi2Family()
                .RegisterPi2DyadicLadder()
                .RegisterPi2I4MemoryLoop()
                // Missing: RegisterF1Pi2Inheritance
                .RegisterF41PalindromicTimePi2Inheritance()
                .Build());
    }
}
