using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class F88PopcountCoherenceRegistrationTests
{
    [Fact]
    public void RegisterF88_AddsDualParentClaim()
    {
        var registry = new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterF88PopcountCoherence()
            .Build();

        Assert.Equal(8, registry.All().Count()); // 7 Pi2 + 1 F88
        Assert.True(registry.Contains<PopcountCoherenceClaim>());
    }

    [Fact]
    public void RegisterF88_AncestorsContainBothKleinFourCellAndPolarityLayerOrigin()
    {
        var registry = new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterF88PopcountCoherence()
            .Build();

        var ancestors = registry.AncestorsOf<PopcountCoherenceClaim>()
            .Select(c => c.GetType()).ToHashSet();

        // Dual-parent: both KleinFourCell AND PolarityLayerOrigin must appear.
        Assert.Contains(typeof(KleinFourCellClaim), ancestors);
        Assert.Contains(typeof(PolarityLayerOriginClaim), ancestors);

        // Both transitively reach PolynomialFoundation.
        Assert.Contains(typeof(QubitDimensionalAnchorClaim), ancestors);
        Assert.Contains(typeof(PolynomialFoundationClaim), ancestors);
    }

    [Fact]
    public void RegisterF88_TopologicalOrder_AfterBothParents()
    {
        var registry = new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterF88PopcountCoherence()
            .Build();

        var order = registry.TopologicalOrder.ToList();
        var f88Idx = order.IndexOf(typeof(PopcountCoherenceClaim));
        var kleinIdx = order.IndexOf(typeof(KleinFourCellClaim));
        var polarityIdx = order.IndexOf(typeof(PolarityLayerOriginClaim));

        Assert.True(kleinIdx < f88Idx, "KleinFourCell must resolve before PopcountCoherence");
        Assert.True(polarityIdx < f88Idx, "PolarityLayerOrigin must resolve before PopcountCoherence");
    }
}
