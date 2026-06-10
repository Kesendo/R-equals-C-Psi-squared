using System.Linq;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Diagnostics.F87;
using RCPsiSquared.Diagnostics.Knowledge;

namespace RCPsiSquared.Diagnostics.Tests.F87;

/// <summary>Wiring audit for <see cref="WindowedConverseAllGammaClaim"/>, the F87 windowed-converse
/// all-γ theorem. After Phase B it was "proven modulo R-deg + R-sign"; 2026-06-10 retired R-deg
/// (girth dichotomy) and resolved R-sign (Pascal-Gram positivity) the same day. It is now
/// Tier1Derived with NO residual, the Tier1Derived <see cref="WindowedConverseThresholdClaim"/>
/// as a typed parent.</summary>
public class WindowedConverseAllGammaClaimRegistrationTests
{
    [Fact]
    public void BuildDefault_ContainsWindowedConverseAllGammaClaim()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        Assert.True(registry.Contains<WindowedConverseAllGammaClaim>());
    }

    [Fact]
    public void TypedParents_ArePresentInRegistry()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        Assert.True(registry.Contains<F87DiagonalCellBipartiteWitnessSet>(),
            "typed parent F87DiagonalCellBipartiteWitnessSet must be registered");
        Assert.True(registry.Contains<WindowedConverseThresholdClaim>(),
            "typed parent WindowedConverseThresholdClaim (the Phase B spine) must be registered");
    }

    [Fact]
    public void Ancestors_ContainBothParents()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var ancestors = registry.AncestorsOf<WindowedConverseAllGammaClaim>()
            .Select(c => c.GetType())
            .ToList();
        Assert.Contains(typeof(F87DiagonalCellBipartiteWitnessSet), ancestors);
        Assert.Contains(typeof(WindowedConverseThresholdClaim), ancestors);
    }

    [Fact]
    public void TierIsTier1Derived()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        Assert.Equal(Tier.Tier1Derived, registry.Get<WindowedConverseAllGammaClaim>().Tier);
    }

    [Fact]
    public void Summary_StatesNoResidual_RDegRetired_AndRSignResolved()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var claim = registry.Get<WindowedConverseAllGammaClaim>();
        Assert.Contains("NO residual", claim.Summary);
        Assert.Contains("R-deg retired", claim.Summary);
        Assert.Contains("R-sign resolved", claim.Summary);
        Assert.Contains("Pascal-Gram", claim.Summary);
        Assert.DoesNotContain("modulo", claim.Summary);
    }
}
