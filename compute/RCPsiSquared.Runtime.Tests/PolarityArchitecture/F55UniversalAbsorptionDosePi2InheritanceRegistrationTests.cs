using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class F55UniversalAbsorptionDosePi2InheritanceRegistrationTests
{
    private static ClaimRegistryBuilder BuildBaseRegistry() =>
        new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterPi2DyadicLadder()
            .RegisterAbsorptionTheoremClaim()
            .RegisterF50WeightOneDegeneracyPi2Inheritance();

    [Fact]
    public void RegisterF55_AddsClaim()
    {
        var registry = BuildBaseRegistry()
            .RegisterF55UniversalAbsorptionDosePi2Inheritance()
            .Build();

        Assert.True(registry.Contains<F55UniversalAbsorptionDosePi2Inheritance>());
    }

    [Fact]
    public void RegisterF55_TierIsTier1Derived()
    {
        var registry = BuildBaseRegistry()
            .RegisterF55UniversalAbsorptionDosePi2Inheritance()
            .Build();

        Assert.Equal(Tier.Tier1Derived,
            registry.Get<F55UniversalAbsorptionDosePi2Inheritance>().Tier);
    }

    [Fact]
    public void RegisterF55_AncestorsContainF50AndPi2DyadicLadder()
    {
        var registry = BuildBaseRegistry()
            .RegisterF55UniversalAbsorptionDosePi2Inheritance()
            .Build();

        var ancestors = registry.AncestorsOf<F55UniversalAbsorptionDosePi2Inheritance>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(F50WeightOneDegeneracyPi2Inheritance), ancestors);
        Assert.Contains(typeof(Pi2DyadicLadderClaim), ancestors);
    }

    [Theory]
    [InlineData(1, 2)]
    [InlineData(7, 8)]
    public void RegisterF55_ImmortalModeCountAcrossRegistry(int N, int expected)
    {
        var registry = BuildBaseRegistry()
            .RegisterF55UniversalAbsorptionDosePi2Inheritance()
            .Build();

        Assert.Equal(expected, registry.Get<F55UniversalAbsorptionDosePi2Inheritance>().ImmortalModeCount(N));
    }

    [Fact]
    public void RegisterF55_DerivationConsistencyHoldsAcrossRegistry()
    {
        var registry = BuildBaseRegistry()
            .RegisterF55UniversalAbsorptionDosePi2Inheritance()
            .Build();

        var f55 = registry.Get<F55UniversalAbsorptionDosePi2Inheritance>();
        Assert.True(f55.DerivationConsistencyHolds(0.05));
        Assert.True(f55.DerivationConsistencyHolds(1.0));
    }

    [Fact]
    public void RegisterF55_RateMinMatchesF50AcrossRegistry()
    {
        var registry = BuildBaseRegistry()
            .RegisterF55UniversalAbsorptionDosePi2Inheritance()
            .Build();

        var f55 = registry.Get<F55UniversalAbsorptionDosePi2Inheritance>();
        Assert.True(f55.RateMinMatchesF50(0.05));
        Assert.True(f55.RateMinMatchesF50(0.5));
    }

    [Fact]
    public void RegisterF55_WithoutF50_Throws()
    {
        Assert.Throws<InvariantViolationException>(() =>
            new ClaimRegistryBuilder()
                .RegisterPi2Family()
                .RegisterPi2DyadicLadder()
                .RegisterAbsorptionTheoremClaim()
                // Missing: RegisterF50WeightOneDegeneracyPi2Inheritance
                .RegisterF55UniversalAbsorptionDosePi2Inheritance()
                .Build());
    }
}
