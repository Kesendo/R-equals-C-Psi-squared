using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class F96BornSubdominantSlopesPi2InheritanceRegistrationTests
{
    private static ClaimRegistryBuilder BuildBaseRegistry() =>
        new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterPi2DyadicLadder()
            .RegisterF94BornDeviationFourThirdsPi2Inheritance();

    [Fact]
    public void RegisterF96_AddsClaim()
    {
        var registry = BuildBaseRegistry()
            .RegisterF96BornSubdominantSlopesPi2Inheritance()
            .Build();

        Assert.True(registry.Contains<F96BornSubdominantSlopesPi2Inheritance>());
    }

    [Fact]
    public void RegisterF96_TierIsTier1Derived()
    {
        var registry = BuildBaseRegistry()
            .RegisterF96BornSubdominantSlopesPi2Inheritance()
            .Build();

        Assert.Equal(Tier.Tier1Derived,
            registry.Get<F96BornSubdominantSlopesPi2Inheritance>().Tier);
    }

    [Fact]
    public void RegisterF96_AncestorsContainF94()
    {
        // The 4/3 unit of F94 is the typed parent that F96 elaborates.
        var registry = BuildBaseRegistry()
            .RegisterF96BornSubdominantSlopesPi2Inheritance()
            .Build();

        var ancestors = registry.AncestorsOf<F96BornSubdominantSlopesPi2Inheritance>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(F94BornDeviationFourThirdsPi2Inheritance), ancestors);
    }

    [Fact]
    public void RegisterF96_SinglyFlippedSlopeIsMinus16Over9()
    {
        var registry = BuildBaseRegistry()
            .RegisterF96BornSubdominantSlopesPi2Inheritance()
            .Build();

        var f = registry.Get<F96BornSubdominantSlopesPi2Inheritance>();
        Assert.Equal(-16.0 / 9.0, f.SlopeSingleFlipped, precision: 15);
    }

    [Fact]
    public void RegisterF96_DoublyFlippedSlopeIsMinus8Over3()
    {
        var registry = BuildBaseRegistry()
            .RegisterF96BornSubdominantSlopesPi2Inheritance()
            .Build();

        var f = registry.Get<F96BornSubdominantSlopesPi2Inheritance>();
        Assert.Equal(-8.0 / 3.0, f.SlopeDoubleFlipped, precision: 15);
    }
}
