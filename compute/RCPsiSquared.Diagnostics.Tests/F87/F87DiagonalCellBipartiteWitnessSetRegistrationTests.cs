using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.F87;
using RCPsiSquared.Diagnostics.Knowledge;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Diagnostics.Tests.F87;

/// <summary>Wiring audit for <see cref="F87DiagonalCellBipartiteWitnessSet"/>: it surfaces the
/// F103 §7 diagonal-K (bipartite-chirality) criterion as a registered, registry-queryable
/// Claim, mirroring how <see cref="F87StandardWitnessSet"/> surfaces the canonical witnesses.
/// Built from the full <see cref="KnowledgeRegistryFactory.BuildDefault"/> registry (N=5), so
/// the criterion is also exercised against the registry's production chain.</summary>
public class F87DiagonalCellBipartiteWitnessSetRegistrationTests
{
    [Fact]
    public void BuildDefault_ContainsBipartiteWitnessSet()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        Assert.True(registry.Contains<F87DiagonalCellBipartiteWitnessSet>());
    }

    [Fact]
    public void BipartiteWitnessSet_TypedParents_ArePresentInRegistry()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();

        Assert.True(registry.Contains<F87TrichotomyClassification>(),
            "typed parent F87TrichotomyClassification must be registered");
        Assert.True(registry.Contains<ChiralKClaim>(),
            "typed parent ChiralKClaim must be registered");
    }

    [Fact]
    public void BipartiteWitnessSet_Ancestors_ContainTrichotomyAndChiralK()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();

        var ancestors = registry.AncestorsOf<F87DiagonalCellBipartiteWitnessSet>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(F87TrichotomyClassification), ancestors);
        Assert.Contains(typeof(ChiralKClaim), ancestors);
    }

    [Fact]
    public void BipartiteWitnessSet_TierIsTier1Candidate()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        Assert.Equal(Tier.Tier1Candidate, registry.Get<F87DiagonalCellBipartiteWitnessSet>().Tier);
    }

    [Fact]
    public void BipartiteWitnessSet_HasFourWitnesses()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var set = registry.Get<F87DiagonalCellBipartiteWitnessSet>();
        Assert.Equal(4, set.Witnesses.Count);
    }

    [Fact]
    public void BipartiteWitnessSet_AllFourWitnessesPass()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var set = registry.Get<F87DiagonalCellBipartiteWitnessSet>();

        foreach (var w in set.Witnesses)
        {
            Assert.True(w.CriterionAgrees,
                $"witness '{w.WitnessName}': criterion predicted {w.PredictedClass}, actual {w.ActualClass}");
            Assert.Equal(w.ExpectedClass, w.ActualClass);
            Assert.True(w.Matches,
                $"witness '{w.WitnessName}': expected {w.ExpectedClass}, actual {w.ActualClass}, bipartite={w.IsBipartite}");
        }
    }
}
