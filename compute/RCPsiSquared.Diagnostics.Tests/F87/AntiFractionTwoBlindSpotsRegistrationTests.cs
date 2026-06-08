using System.Linq;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Diagnostics.F87;
using RCPsiSquared.Diagnostics.Knowledge;

namespace RCPsiSquared.Diagnostics.Tests.F87;

/// <summary>Wiring audit for the two-blind-spots connector: registered beneath both parents
/// (F89F87BreakPredictionFromF83 and AntiFractionObstructionOrthogonalityClaim, both Diagnostics,
/// both Tier1Derived), Tier1Derived (5 >= 5).</summary>
public class AntiFractionTwoBlindSpotsRegistrationTests
{
    [Fact]
    public void BuildDefault_ContainsConnector() =>
        Assert.True(KnowledgeRegistryFactory.BuildDefault().Contains<AntiFractionTwoBlindSpotsClaim>());

    [Fact]
    public void Connector_BothTypedParents_ArePresent()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        Assert.True(registry.Contains<F89F87BreakPredictionFromF83>());
        Assert.True(registry.Contains<AntiFractionObstructionOrthogonalityClaim>());
    }

    [Fact]
    public void Connector_Ancestors_ContainBothParents()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var ancestors = registry.AncestorsOf<AntiFractionTwoBlindSpotsClaim>().Select(c => c.GetType()).ToHashSet();
        Assert.Contains(typeof(F89F87BreakPredictionFromF83), ancestors);
        Assert.Contains(typeof(AntiFractionObstructionOrthogonalityClaim), ancestors);
    }

    [Fact]
    public void Connector_TierIsTier1Derived_AndSelfCheckPasses()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var claim = registry.Get<AntiFractionTwoBlindSpotsClaim>();
        Assert.Equal(Tier.Tier1Derived, claim.Tier);
        foreach (var c in claim.Cases) Assert.True(c.Passes);
    }
}
