using System.Linq;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Diagnostics.F87;
using RCPsiSquared.Diagnostics.Knowledge;

namespace RCPsiSquared.Diagnostics.Tests.F87;

/// <summary>Wiring audit for <see cref="WindowedConverseThresholdClaim"/>, the Tier1Derived
/// two-reflection spine of the F87 windowed-converse monomial theorem. Typed parent
/// <see cref="F87DiagonalCellBipartiteWitnessSet"/> (Tier1Derived ≥ Tier1Derived).</summary>
public class WindowedConverseThresholdClaimRegistrationTests
{
    [Fact]
    public void BuildDefault_ContainsWindowedConverseThresholdClaim()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        Assert.True(registry.Contains<WindowedConverseThresholdClaim>());
    }

    [Fact]
    public void TypedParent_IsPresentInRegistry()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        Assert.True(registry.Contains<F87DiagonalCellBipartiteWitnessSet>(),
            "typed parent F87DiagonalCellBipartiteWitnessSet must be registered");
    }

    [Fact]
    public void Ancestors_ContainBipartiteWitnessSet()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var ancestors = registry.AncestorsOf<WindowedConverseThresholdClaim>()
            .Select(c => c.GetType())
            .ToList();
        Assert.Contains(typeof(F87DiagonalCellBipartiteWitnessSet), ancestors);
    }

    [Fact]
    public void TierIsTier1Derived()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        Assert.Equal(Tier.Tier1Derived, registry.Get<WindowedConverseThresholdClaim>().Tier);
    }

    [Fact]
    public void SelfCheck_AllCasesPass()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var claim = registry.Get<WindowedConverseThresholdClaim>();
        Assert.Equal(claim.Cases.Count, claim.PassCount);
        Assert.NotEmpty(claim.Cases);
    }
}
