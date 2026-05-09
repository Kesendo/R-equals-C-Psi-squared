using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class F50WeightOneDegeneracyPi2InheritanceRegistrationTests
{
    private static ClaimRegistryBuilder BuildBaseRegistry() =>
        new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterPi2DyadicLadder();

    [Fact]
    public void RegisterF50_AddsClaim()
    {
        var registry = BuildBaseRegistry()
            .RegisterF50WeightOneDegeneracyPi2Inheritance()
            .Build();

        Assert.True(registry.Contains<F50WeightOneDegeneracyPi2Inheritance>());
    }

    [Fact]
    public void RegisterF50_TierIsTier1Derived()
    {
        var registry = BuildBaseRegistry()
            .RegisterF50WeightOneDegeneracyPi2Inheritance()
            .Build();

        Assert.Equal(Tier.Tier1Derived,
            registry.Get<F50WeightOneDegeneracyPi2Inheritance>().Tier);
    }

    [Fact]
    public void RegisterF50_AncestorsContainPi2DyadicLadder()
    {
        var registry = BuildBaseRegistry()
            .RegisterF50WeightOneDegeneracyPi2Inheritance()
            .Build();

        var ancestors = registry.AncestorsOf<F50WeightOneDegeneracyPi2Inheritance>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(Pi2DyadicLadderClaim), ancestors);
    }

    [Theory]
    [InlineData(2, 4)]
    [InlineData(5, 10)]
    [InlineData(7, 14)]
    public void RegisterF50_TotalDegeneracyAcrossRegistry(int N, int expected)
    {
        var registry = BuildBaseRegistry()
            .RegisterF50WeightOneDegeneracyPi2Inheritance()
            .Build();

        Assert.Equal(expected, registry.Get<F50WeightOneDegeneracyPi2Inheritance>().TotalDegeneracy(N));
    }

    [Fact]
    public void RegisterF50_WithoutPi2DyadicLadder_Throws()
    {
        Assert.Throws<InvariantViolationException>(() =>
            new ClaimRegistryBuilder()
                .RegisterPi2Family()
                // Missing: RegisterPi2DyadicLadder
                .RegisterF50WeightOneDegeneracyPi2Inheritance()
                .Build());
    }
}
