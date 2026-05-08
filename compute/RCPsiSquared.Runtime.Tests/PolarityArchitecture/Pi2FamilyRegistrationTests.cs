using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class Pi2FamilyRegistrationTests
{
    [Fact]
    public void RegisterPi2Family_BuildsEightClaims()
    {
        var registry = new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .Build();

        Assert.Equal(8, registry.All().Count());
        Assert.True(registry.Contains<PolynomialFoundationClaim>());
        Assert.True(registry.Contains<QubitDimensionalAnchorClaim>());
        Assert.True(registry.Contains<NinetyDegreeMirrorMemoryClaim>());
        Assert.True(registry.Contains<PolarityLayerOriginClaim>());
        Assert.True(registry.Contains<BilinearApexClaim>());
        Assert.True(registry.Contains<HalfAsStructuralFixedPointClaim>());
        Assert.True(registry.Contains<QuarterAsBilinearMaxvalClaim>());
        Assert.True(registry.Contains<KleinFourCellClaim>());
    }

    [Fact]
    public void RegisterPi2Family_TopologicalOrder_TrunkFirst()
    {
        var registry = new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .Build();

        var order = registry.TopologicalOrder.ToList();
        Assert.Equal(typeof(PolynomialFoundationClaim), order[0]);

        // QubitDim and 90 degrees both directly depend on Polynomial; either may be index 1 or 2.
        var qubitIdx = order.IndexOf(typeof(QubitDimensionalAnchorClaim));
        var ninetyIdx = order.IndexOf(typeof(NinetyDegreeMirrorMemoryClaim));
        Assert.True(qubitIdx >= 1 && qubitIdx <= 2);
        Assert.True(ninetyIdx >= 1 && ninetyIdx <= 2);

        // KleinFourCell is the deepest; must be after PolarityLayerOrigin.
        var kleinIdx = order.IndexOf(typeof(KleinFourCellClaim));
        var polarityIdx = order.IndexOf(typeof(PolarityLayerOriginClaim));
        Assert.True(polarityIdx < kleinIdx);
    }

    [Fact]
    public void RegisterPi2Family_AllTier1Derived()
    {
        var registry = new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .Build();

        Assert.All(registry.All(), c => Assert.Equal(Tier.Tier1Derived, c.Tier));
    }

    [Fact]
    public void RegisterPi2Family_PolarityLayerOrigin_DescendsFromQubitDimensional()
    {
        var registry = new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .Build();

        var ancestors = registry.AncestorsOf<PolarityLayerOriginClaim>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(QubitDimensionalAnchorClaim), ancestors);
        Assert.Contains(typeof(PolynomialFoundationClaim), ancestors);
    }

    [Fact]
    public void RegisterPi2Family_KleinFourCell_DescendsFromPolarityLayerOrigin()
    {
        var registry = new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .Build();

        var ancestors = registry.AncestorsOf<KleinFourCellClaim>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(PolarityLayerOriginClaim), ancestors);
        Assert.Contains(typeof(QubitDimensionalAnchorClaim), ancestors);
        Assert.Contains(typeof(PolynomialFoundationClaim), ancestors);
    }

    [Fact]
    public void RegisterPi2Family_NinetyDegreeMemory_DescendsFromPolynomialOnly()
    {
        var registry = new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .Build();

        var ancestors = registry.AncestorsOf<NinetyDegreeMirrorMemoryClaim>()
            .Select(c => c.GetType()).ToHashSet();

        // 90 degrees is on the second branch off the trunk: only Polynomial as ancestor.
        Assert.Single(ancestors);
        Assert.Contains(typeof(PolynomialFoundationClaim), ancestors);
    }
}
