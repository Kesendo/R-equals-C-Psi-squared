using System.Linq;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.F87;
using RCPsiSquared.Diagnostics.Knowledge;

namespace RCPsiSquared.Diagnostics.Tests.F87;

/// <summary>Wiring audit for the F81/F115 connector: registered in the default registry beneath both
/// parents (F83AntiFractionPi2Inheritance in Core, WindowedHardnessClaim in Diagnostics), Tier1Derived,
/// parent >= child (5 >= 5).</summary>
public class AntiFractionObstructionOrthogonalityRegistrationTests
{
    [Fact]
    public void BuildDefault_ContainsConnector()
    {
        Assert.True(KnowledgeRegistryFactory.BuildDefault()
            .Contains<AntiFractionObstructionOrthogonalityClaim>());
    }

    [Fact]
    public void Connector_BothTypedParents_ArePresent()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        Assert.True(registry.Contains<F83AntiFractionPi2Inheritance>());
        Assert.True(registry.Contains<WindowedHardnessClaim>());
    }

    [Fact]
    public void Connector_Ancestors_ContainBothParents()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var ancestors = registry.AncestorsOf<AntiFractionObstructionOrthogonalityClaim>()
            .Select(c => c.GetType()).ToHashSet();
        Assert.Contains(typeof(F83AntiFractionPi2Inheritance), ancestors);
        Assert.Contains(typeof(WindowedHardnessClaim), ancestors);
    }

    [Fact]
    public void Connector_TierIsTier1Derived_AndSelfCheckPasses()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var claim = registry.Get<AntiFractionObstructionOrthogonalityClaim>();
        Assert.Equal(Tier.Tier1Derived, claim.Tier);
        foreach (var c in claim.Cases) Assert.True(c.Passes);
    }
}
