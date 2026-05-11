using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class AbsorptionTheoremClaimRegistrationTests
{
    private static ClaimRegistryBuilder BuildBaseRegistry() =>
        new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterPi2DyadicLadder();

    [Fact]
    public void RegisterAbsorptionTheoremClaim_AddsClaim()
    {
        var registry = BuildBaseRegistry()
            .RegisterAbsorptionTheoremClaim()
            .Build();

        Assert.True(registry.Contains<AbsorptionTheoremClaim>());
    }

    [Fact]
    public void RegisterAbsorptionTheoremClaim_TierIsTier1Derived()
    {
        var registry = BuildBaseRegistry()
            .RegisterAbsorptionTheoremClaim()
            .Build();

        Assert.Equal(Tier.Tier1Derived,
            registry.Get<AbsorptionTheoremClaim>().Tier);
    }

    [Fact]
    public void RegisterAbsorptionTheoremClaim_AncestorsContainPi2DyadicLadder()
    {
        var registry = BuildBaseRegistry()
            .RegisterAbsorptionTheoremClaim()
            .Build();

        var ancestors = registry.AncestorsOf<AbsorptionTheoremClaim>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(Pi2DyadicLadderClaim), ancestors);
    }

    [Fact]
    public void RegisterAbsorptionTheoremClaim_AbsorptionQuantumCoefficientIsTwo()
    {
        var registry = BuildBaseRegistry()
            .RegisterAbsorptionTheoremClaim()
            .Build();

        Assert.Equal(2.0,
            registry.Get<AbsorptionTheoremClaim>().AbsorptionQuantumCoefficient,
            precision: 14);
    }

    [Theory]
    [InlineData(0.05, 0.1)]
    [InlineData(1.0, 2.0)]
    public void RegisterAbsorptionTheoremClaim_AbsorptionQuantumAcrossRegistry(double gammaZero, double expected)
    {
        var registry = BuildBaseRegistry()
            .RegisterAbsorptionTheoremClaim()
            .Build();

        Assert.Equal(expected,
            registry.Get<AbsorptionTheoremClaim>().AbsorptionQuantum(gammaZero),
            precision: 12);
    }

    [Fact]
    public void RegisterAbsorptionTheoremClaim_F66BecomesDescendantWhenWired()
    {
        var registry = BuildBaseRegistry()
            .RegisterAbsorptionTheoremClaim()
            .RegisterF66PoleModesPi2Inheritance()
            .Build();

        var descendants = registry.DescendantsOf<AbsorptionTheoremClaim>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(F66PoleModesPi2Inheritance), descendants);
    }

    [Fact]
    public void RegisterAbsorptionTheoremClaim_F66AncestorsContainAbsorptionTheorem()
    {
        var registry = BuildBaseRegistry()
            .RegisterAbsorptionTheoremClaim()
            .RegisterF66PoleModesPi2Inheritance()
            .Build();

        var ancestors = registry.AncestorsOf<F66PoleModesPi2Inheritance>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(AbsorptionTheoremClaim), ancestors);
    }
}
