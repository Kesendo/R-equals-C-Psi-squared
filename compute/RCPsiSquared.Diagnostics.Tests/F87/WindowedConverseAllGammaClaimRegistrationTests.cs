using System.Linq;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Diagnostics.F87;
using RCPsiSquared.Diagnostics.Knowledge;

namespace RCPsiSquared.Diagnostics.Tests.F87;

/// <summary>Wiring audit for <see cref="WindowedConverseAllGammaClaim"/>, the F87 windowed-converse
/// all-γ residual lemma (Phase A "type the seam"). Single typed parent
/// <see cref="F87DiagonalCellBipartiteWitnessSet"/> (Tier1Derived ≥ this claim's Tier1Candidate);
/// the Summary must keep stating the open status so it cannot silently drift to a proven claim.</summary>
public class WindowedConverseAllGammaClaimRegistrationTests
{
    [Fact]
    public void BuildDefault_ContainsWindowedConverseAllGammaClaim()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        Assert.True(registry.Contains<WindowedConverseAllGammaClaim>());
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
        var ancestors = registry.AncestorsOf<WindowedConverseAllGammaClaim>()
            .Select(c => c.GetType())
            .ToList();
        Assert.Contains(typeof(F87DiagonalCellBipartiteWitnessSet), ancestors);
    }

    [Fact]
    public void TierIsTier1Candidate()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        Assert.Equal(Tier.Tier1Candidate, registry.Get<WindowedConverseAllGammaClaim>().Tier);
    }

    [Fact]
    public void Summary_StatesOpenResidual()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var claim = registry.Get<WindowedConverseAllGammaClaim>();
        Assert.Contains("OPEN", claim.Summary);
        Assert.Contains("Phase B", claim.Summary);
    }
}
