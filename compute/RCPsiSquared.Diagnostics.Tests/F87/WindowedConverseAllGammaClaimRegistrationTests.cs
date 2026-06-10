using System.Linq;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Diagnostics.F87;
using RCPsiSquared.Diagnostics.Knowledge;

namespace RCPsiSquared.Diagnostics.Tests.F87;

/// <summary>Wiring audit for <see cref="WindowedConverseAllGammaClaim"/>, the F87 windowed-converse
/// residual lemma. After Phase B it was "proven modulo R-deg + R-sign"; since 2026-06-10 the girth
/// dichotomy retired R-deg, so it is "proven modulo R-sign" (t_ℓ = 0 branch only). It stays
/// Tier1Candidate with the Tier1Derived <see cref="WindowedConverseThresholdClaim"/> as a typed
/// parent, and R-sign is its one named open content.</summary>
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
    public void TierIsTier1Candidate()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        Assert.Equal(Tier.Tier1Candidate, registry.Get<WindowedConverseAllGammaClaim>().Tier);
    }

    [Fact]
    public void Summary_StatesProvenModuloRSign_AndRDegRetired()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var claim = registry.Get<WindowedConverseAllGammaClaim>();
        Assert.Contains("modulo", claim.Summary);
        Assert.Contains("R-sign", claim.Summary);
        Assert.Contains("R-deg retired", claim.Summary);
        Assert.Contains("outright", claim.Summary.ToLowerInvariant());
    }
}
