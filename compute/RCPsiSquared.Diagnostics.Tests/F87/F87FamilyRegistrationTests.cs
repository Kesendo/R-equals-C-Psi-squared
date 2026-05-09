using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.F87;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Diagnostics.Tests.F87;

public class F87FamilyRegistrationTests
{
    [Fact]
    public void RegisterF87Family_BuildsThreeClaims()
    {
        var registry = new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterF87Family()
            .Build();

        Assert.Equal(12, registry.All().Count()); // 9 Pi2 + 3 F87
        Assert.True(registry.Contains<F87TrichotomyClassification>());
        Assert.True(registry.Contains<DissipatorResonanceLaw>());
        Assert.True(registry.Contains<DissipatorAxisSelectsPolarityClaim>());
    }

    [Fact]
    public void RegisterF87Family_AllTier1Derived()
    {
        var registry = new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterF87Family()
            .Build();

        Assert.Equal(Tier.Tier1Derived, registry.Get<F87TrichotomyClassification>().Tier);
        Assert.Equal(Tier.Tier1Derived, registry.Get<DissipatorResonanceLaw>().Tier);
        Assert.Equal(Tier.Tier1Derived, registry.Get<DissipatorAxisSelectsPolarityClaim>().Tier);
    }

    [Fact]
    public void RegisterF87Family_DissipatorAxis_AncestorContainsPolarityLayerOrigin()
    {
        // The cross-KB edge from DissipatorAxisSelectsPolarityClaim (F87) to
        // PolarityLayerOriginClaim (Pi2) is the architecture's first
        // Diagnostics-to-Runtime cross-KB inheritance.
        var registry = new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterF87Family()
            .Build();

        var ancestors = registry.AncestorsOf<DissipatorAxisSelectsPolarityClaim>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(PolarityLayerOriginClaim), ancestors);
        Assert.Contains(typeof(QubitDimensionalAnchorClaim), ancestors);
        Assert.Contains(typeof(PolynomialFoundationClaim), ancestors);
    }

    [Fact]
    public void RegisterF87Family_DissipatorResonance_AncestorContainsF87Trichotomy()
    {
        // DissipatorResonanceLaw is the empirical statement about F87-hardness, so
        // F87TrichotomyClassification is its typed parent (the F87 closed form is
        // the source of "hardness" being a thing).
        var registry = new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterF87Family()
            .Build();

        var ancestors = registry.AncestorsOf<DissipatorResonanceLaw>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(F87TrichotomyClassification), ancestors);
    }

    [Fact]
    public void RegisterF87Family_DissipatorAxis_AncestorContainsDissipatorResonanceAndF87Trichotomy()
    {
        // DissipatorAxisSelectsPolarityClaim's docstring states "witnesses live upstream"
        // → DissipatorResonanceLaw is its second typed parent (alongside PolarityLayerOrigin).
        // Transitively, F87TrichotomyClassification is also an ancestor through
        // DissipatorResonanceLaw.
        var registry = new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterF87Family()
            .Build();

        var ancestors = registry.AncestorsOf<DissipatorAxisSelectsPolarityClaim>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(DissipatorResonanceLaw), ancestors);
        Assert.Contains(typeof(F87TrichotomyClassification), ancestors);
    }
}
