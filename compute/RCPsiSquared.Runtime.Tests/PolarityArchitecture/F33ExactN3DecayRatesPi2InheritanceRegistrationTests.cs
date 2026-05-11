using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class F33ExactN3DecayRatesPi2InheritanceRegistrationTests
{
    private static ClaimRegistryBuilder BuildBaseRegistry() =>
        new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterPi2DyadicLadder()
            .RegisterAbsorptionTheoremClaim()
            .RegisterF50WeightOneDegeneracyPi2Inheritance();

    [Fact]
    public void RegisterF33_AddsClaim()
    {
        var registry = BuildBaseRegistry()
            .RegisterF33ExactN3DecayRatesPi2Inheritance()
            .Build();

        Assert.True(registry.Contains<F33ExactN3DecayRatesPi2Inheritance>());
    }

    [Fact]
    public void RegisterF33_TierIsTier1Derived()
    {
        var registry = BuildBaseRegistry()
            .RegisterF33ExactN3DecayRatesPi2Inheritance()
            .Build();

        Assert.Equal(Tier.Tier1Derived,
            registry.Get<F33ExactN3DecayRatesPi2Inheritance>().Tier);
    }

    [Fact]
    public void RegisterF33_AncestorsContainF50AndPi2DyadicLadder()
    {
        var registry = BuildBaseRegistry()
            .RegisterF33ExactN3DecayRatesPi2Inheritance()
            .Build();

        var ancestors = registry.AncestorsOf<F33ExactN3DecayRatesPi2Inheritance>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(F50WeightOneDegeneracyPi2Inheritance), ancestors);
        Assert.Contains(typeof(Pi2DyadicLadderClaim), ancestors);
    }

    [Theory]
    [InlineData(0.05, 0.1)]
    [InlineData(1.0, 2.0)]
    public void RegisterF33_Rate1AcrossRegistry(double gammaZero, double expected)
    {
        var registry = BuildBaseRegistry()
            .RegisterF33ExactN3DecayRatesPi2Inheritance()
            .Build();

        Assert.Equal(expected, registry.Get<F33ExactN3DecayRatesPi2Inheritance>().Rate1(gammaZero), precision: 14);
    }

    [Fact]
    public void RegisterF33_RatesAreRationalRatiosAcrossRegistry()
    {
        var registry = BuildBaseRegistry()
            .RegisterF33ExactN3DecayRatesPi2Inheritance()
            .Build();

        Assert.True(registry.Get<F33ExactN3DecayRatesPi2Inheritance>().RatesAreRationalRatios(0.05));
    }

    [Fact]
    public void RegisterF33_AbsorptionTheoremTableAcrossRegistry()
    {
        var registry = BuildBaseRegistry()
            .RegisterF33ExactN3DecayRatesPi2Inheritance()
            .Build();

        var table = registry.Get<F33ExactN3DecayRatesPi2Inheritance>().AbsorptionTheoremTable;
        Assert.Equal(3, table.Count);
    }

    [Fact]
    public void RegisterF33_WithoutF50_Throws()
    {
        Assert.Throws<InvariantViolationException>(() =>
            new ClaimRegistryBuilder()
                .RegisterPi2Family()
                .RegisterPi2DyadicLadder()
                .RegisterAbsorptionTheoremClaim()
                // Missing: RegisterF50WeightOneDegeneracyPi2Inheritance
                .RegisterF33ExactN3DecayRatesPi2Inheritance()
                .Build());
    }
}
