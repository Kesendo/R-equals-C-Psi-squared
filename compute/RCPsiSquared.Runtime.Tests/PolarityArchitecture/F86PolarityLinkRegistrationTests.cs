using RCPsiSquared.Core.F86;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class F86PolarityLinkRegistrationTests
{
    [Fact]
    public void RegisterF86PolarityLink_AddsCrossKBClaim()
    {
        var registry = new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterF86PolarityLink()
            .Build();

        Assert.Equal(10, registry.All().Count()); // 9 Pi2 + 1 F86
        Assert.True(registry.Contains<PolarityInheritanceLink>());
    }

    [Fact]
    public void RegisterF86PolarityLink_AncestorsContainsPolarityTrio()
    {
        var registry = new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterF86PolarityLink()
            .Build();

        var ancestors = registry.AncestorsOf<PolarityInheritanceLink>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(PolynomialFoundationClaim), ancestors);
        Assert.Contains(typeof(PolarityLayerOriginClaim), ancestors);
        Assert.Contains(typeof(HalfAsStructuralFixedPointClaim), ancestors);
        // QubitDim is also an ancestor (transitively).
        Assert.Contains(typeof(QubitDimensionalAnchorClaim), ancestors);
    }

    [Fact]
    public void RegisterF86PolarityLink_TierInheritance_VerifiedSurvives()
    {
        // PolarityInheritanceLink is Tier2Verified; its parents are Tier1Derived.
        // TierStrength check (parent at least as strong as child) passes.
        var registry = new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterF86PolarityLink()
            .Build();

        var link = registry.Get<PolarityInheritanceLink>();
        Assert.Equal(Tier.Tier2Verified, link.Tier);
    }
}
