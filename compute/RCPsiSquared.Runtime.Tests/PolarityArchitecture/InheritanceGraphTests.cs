using RCPsiSquared.Core.F86;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class InheritanceGraphTests
{
    private static ClaimRegistry BuildFullPolarityRegistry() =>
        new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterF86PolarityLink()
            .RegisterF88PopcountCoherence()
            .Build();

    [Fact]
    public void DescendantsOf_PolynomialFoundation_ReturnsAllOthers()
    {
        var registry = BuildFullPolarityRegistry();

        var descendants = registry.DescendantsOf<PolynomialFoundationClaim>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Equal(9, descendants.Count);
        Assert.Contains(typeof(QubitDimensionalAnchorClaim), descendants);
        Assert.Contains(typeof(NinetyDegreeMirrorMemoryClaim), descendants);
        Assert.Contains(typeof(PolarityLayerOriginClaim), descendants);
        Assert.Contains(typeof(BilinearApexClaim), descendants);
        Assert.Contains(typeof(HalfAsStructuralFixedPointClaim), descendants);
        Assert.Contains(typeof(QuarterAsBilinearMaxvalClaim), descendants);
        Assert.Contains(typeof(KleinFourCellClaim), descendants);
        Assert.Contains(typeof(PolarityInheritanceLink), descendants);
        Assert.Contains(typeof(PopcountCoherenceClaim), descendants);
    }

    [Fact]
    public void AncestorsOf_PopcountCoherence_TraversesDualPath()
    {
        var registry = BuildFullPolarityRegistry();

        var ancestors = registry.AncestorsOf<PopcountCoherenceClaim>()
            .Select(c => c.GetType()).ToHashSet();

        // Dual-parent: both direct parents
        Assert.Contains(typeof(KleinFourCellClaim), ancestors);
        Assert.Contains(typeof(PolarityLayerOriginClaim), ancestors);

        // Operator-path closure: KleinFourCell -> PolarityLayerOrigin -> QubitDim -> Polynomial
        // State-path closure: PolarityLayerOrigin -> QubitDim -> Polynomial
        Assert.Contains(typeof(QubitDimensionalAnchorClaim), ancestors);
        Assert.Contains(typeof(PolynomialFoundationClaim), ancestors);

        // 90 degree memory is on a different branch; should NOT appear.
        Assert.DoesNotContain(typeof(NinetyDegreeMirrorMemoryClaim), ancestors);
    }

    private sealed class WeakerFakePolarityLayerOrigin : Claim
    {
        public WeakerFakePolarityLayerOrigin()
            : base("WeakerFakePolarityLayerOrigin (synthetic Tier2Empirical)",
                   Tier.Tier2Empirical,
                   "compute/RCPsiSquared.Runtime.Tests/TestSupport/TestClaims.cs") { }
        public override string DisplayName => "WeakerFakePolarityLayerOrigin";
        public override string Summary => "synthetic downgrade pretender for cascade test";
    }

    [Fact]
    public void TierCascade_DowngradedPolarityLayerOrigin_FailsBuild()
    {
        // Synthetic test: a Tier2Empirical "fake" PolarityLayerOrigin replaces the real one;
        // PopcountCoherenceClaim (Tier1Derived) must throw TierInheritance because its
        // synthetic parent is too weak.
        var ex = Assert.Throws<InvariantViolationException>(() =>
            new ClaimRegistryBuilder()
                .Register<WeakerFakePolarityLayerOrigin>(_ => new WeakerFakePolarityLayerOrigin())
                .Register<KleinFourCellClaim>(b =>
                {
                    _ = b.Get<WeakerFakePolarityLayerOrigin>();
                    return new KleinFourCellClaim();
                })
                .Register<PopcountCoherenceClaim>(b =>
                {
                    _ = b.Get<KleinFourCellClaim>();
                    _ = b.Get<WeakerFakePolarityLayerOrigin>();
                    return new PopcountCoherenceClaim();
                })
                .Build());

        Assert.Equal("TierInheritance", ex.Rule);
        Assert.Contains("WeakerFakePolarityLayerOrigin", ex.Message);
        // Edge iteration order is not stable; either child name may surface first,
        // but the message must name at least one of them.
        Assert.True(
            ex.Message.Contains("KleinFourCellClaim") || ex.Message.Contains("PopcountCoherenceClaim"),
            $"expected at least one child name in message, got: {ex.Message}");
    }
}
