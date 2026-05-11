using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.F1Family;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class F68PalindromicPartnerPi2InheritanceRegistrationTests
{
    private static ChainSystem DefaultChain() =>
        new ChainSystem(N: 5, J: 1.0, GammaZero: 0.05,
            HType: HamiltonianType.Heisenberg, Topology: TopologyKind.Chain);

    private static ClaimRegistryBuilder BuildBaseRegistry() =>
        new ClaimRegistryBuilder()
            .RegisterF1Family(DefaultChain())
            .RegisterPi2Family()
            .RegisterPi2DyadicLadder()
            .RegisterAbsorptionTheoremClaim()
            .RegisterPi2I4MemoryLoop()
            .RegisterF1Pi2Inheritance()
            .RegisterF66PoleModesPi2Inheritance()
            .RegisterF65XxChainSpectrumPi2Inheritance()
            .RegisterF67BondingBellPairPi2Inheritance();

    [Fact]
    public void RegisterF68_AddsClaim()
    {
        var registry = BuildBaseRegistry()
            .RegisterF68PalindromicPartnerPi2Inheritance()
            .Build();

        Assert.True(registry.Contains<F68PalindromicPartnerPi2Inheritance>());
    }

    [Fact]
    public void RegisterF68_TierIsTier1Derived()
    {
        var registry = BuildBaseRegistry()
            .RegisterF68PalindromicPartnerPi2Inheritance()
            .Build();

        Assert.Equal(Tier.Tier1Derived,
            registry.Get<F68PalindromicPartnerPi2Inheritance>().Tier);
    }

    [Fact]
    public void RegisterF68_AncestorsContainF67AndF1()
    {
        var registry = BuildBaseRegistry()
            .RegisterF68PalindromicPartnerPi2Inheritance()
            .Build();

        var ancestors = registry.AncestorsOf<F68PalindromicPartnerPi2Inheritance>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(F67BondingBellPairPi2Inheritance), ancestors);
        Assert.Contains(typeof(F1Pi2Inheritance), ancestors);
        Assert.Contains(typeof(Pi2DyadicLadderClaim), ancestors);
        // Transitive: F65 (via F67), F66 (via F65), Pi2I4MemoryLoop (via F1).
        Assert.Contains(typeof(F65XxChainSpectrumPi2Inheritance), ancestors);
        Assert.Contains(typeof(F66PoleModesPi2Inheritance), ancestors);
        Assert.Contains(typeof(Pi2I4MemoryLoopClaim), ancestors);
    }

    [Theory]
    [InlineData(3, 0.05, 0.075)]
    [InlineData(5, 0.05, 0.1 - 0.05 / 6.0)]
    public void RegisterF68_PartnerRateAcrossRegistry(int N, double gammaZero, double expected)
    {
        var registry = BuildBaseRegistry()
            .RegisterF68PalindromicPartnerPi2Inheritance()
            .Build();

        Assert.Equal(expected, registry.Get<F68PalindromicPartnerPi2Inheritance>().PartnerRate(N, gammaZero), precision: 12);
    }

    [Fact]
    public void RegisterF68_PalindromicSumHoldsAcrossRegistry()
    {
        var registry = BuildBaseRegistry()
            .RegisterF68PalindromicPartnerPi2Inheritance()
            .Build();

        var f68 = registry.Get<F68PalindromicPartnerPi2Inheritance>();
        Assert.True(f68.PalindromicSumHolds(3, 0.05));
        Assert.True(f68.PalindromicSumHolds(5, 0.05));
        Assert.True(f68.PalindromicSumHolds(7, 0.1));
    }

    [Fact]
    public void RegisterF68_WithoutF67_Throws()
    {
        Assert.Throws<InvariantViolationException>(() =>
            new ClaimRegistryBuilder()
                .RegisterF1Family(DefaultChain())
                .RegisterPi2Family()
                .RegisterPi2DyadicLadder()
                .RegisterAbsorptionTheoremClaim()
                .RegisterPi2I4MemoryLoop()
                .RegisterF1Pi2Inheritance()
                .RegisterF66PoleModesPi2Inheritance()
                .RegisterF65XxChainSpectrumPi2Inheritance()
                // Missing: RegisterF67BondingBellPairPi2Inheritance
                .RegisterF68PalindromicPartnerPi2Inheritance()
                .Build());
    }
}
