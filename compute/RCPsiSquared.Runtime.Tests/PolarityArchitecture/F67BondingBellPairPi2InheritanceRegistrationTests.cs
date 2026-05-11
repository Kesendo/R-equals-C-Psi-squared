using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class F67BondingBellPairPi2InheritanceRegistrationTests
{
    private static ClaimRegistryBuilder BuildBaseRegistry() =>
        new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterPi2DyadicLadder()
            .RegisterAbsorptionTheoremClaim()
            .RegisterF66PoleModesPi2Inheritance()
            .RegisterF65XxChainSpectrumPi2Inheritance();

    [Fact]
    public void RegisterF67_AddsClaim()
    {
        var registry = BuildBaseRegistry()
            .RegisterF67BondingBellPairPi2Inheritance()
            .Build();

        Assert.True(registry.Contains<F67BondingBellPairPi2Inheritance>());
    }

    [Fact]
    public void RegisterF67_TierIsTier1Derived()
    {
        var registry = BuildBaseRegistry()
            .RegisterF67BondingBellPairPi2Inheritance()
            .Build();

        Assert.Equal(Tier.Tier1Derived,
            registry.Get<F67BondingBellPairPi2Inheritance>().Tier);
    }

    [Fact]
    public void RegisterF67_AncestorsContainF65AndPi2DyadicLadder()
    {
        var registry = BuildBaseRegistry()
            .RegisterF67BondingBellPairPi2Inheritance()
            .Build();

        var ancestors = registry.AncestorsOf<F67BondingBellPairPi2Inheritance>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(F65XxChainSpectrumPi2Inheritance), ancestors);
        Assert.Contains(typeof(Pi2DyadicLadderClaim), ancestors);
        // Transitive ancestors via F65: F66 + QubitDimensionalAnchor.
        Assert.Contains(typeof(F66PoleModesPi2Inheritance), ancestors);
    }

    [Theory]
    [InlineData(3, 1.0, 0.5)]
    [InlineData(5, 1.0, 1.0 / 6.0)]
    public void RegisterF67_BondingModeDecayRateAcrossRegistry(int N, double gammaZero, double expected)
    {
        var registry = BuildBaseRegistry()
            .RegisterF67BondingBellPairPi2Inheritance()
            .Build();

        Assert.Equal(expected, registry.Get<F67BondingBellPairPi2Inheritance>().BondingModeDecayRate(N, gammaZero), precision: 8);
    }

    [Fact]
    public void RegisterF67_PalindromicEquivalenceHolds()
    {
        var registry = BuildBaseRegistry()
            .RegisterF67BondingBellPairPi2Inheritance()
            .Build();

        Assert.True(registry.Get<F67BondingBellPairPi2Inheritance>().PalindromicAEqualsCAtBondingMode(7));
    }

    [Fact]
    public void RegisterF67_WithoutF65_Throws()
    {
        Assert.Throws<InvariantViolationException>(() =>
            new ClaimRegistryBuilder()
                .RegisterPi2Family()
                .RegisterPi2DyadicLadder()
                .RegisterAbsorptionTheoremClaim()
                .RegisterF66PoleModesPi2Inheritance()
                // Missing: RegisterF65XxChainSpectrumPi2Inheritance
                .RegisterF67BondingBellPairPi2Inheritance()
                .Build());
    }
}
