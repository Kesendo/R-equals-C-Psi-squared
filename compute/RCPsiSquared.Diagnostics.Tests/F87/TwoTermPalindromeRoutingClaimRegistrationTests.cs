using System.Linq;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.F87;
using RCPsiSquared.Diagnostics.Knowledge;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.F87;

/// <summary>Wiring audit for <see cref="TwoTermPalindromeRoutingClaim"/>: it surfaces the
/// Liouvillian-free two-term Q-pair router (<see cref="TwoTermPalindromeRouting"/>) as a
/// registered, registry-queryable Claim, the non-diagonal Q-pair counterpart to the diagonal
/// <see cref="F87DiagonalCellBipartiteWitnessSet"/>. Built from the full
/// <see cref="KnowledgeRegistryFactory.BuildDefault"/> registry, so the router is also exercised
/// against the spectral authority at registry-build time.</summary>
public class TwoTermPalindromeRoutingClaimRegistrationTests
{
    [Fact]
    public void BuildDefault_ContainsTwoTermPalindromeRoutingClaim()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        Assert.True(registry.Contains<TwoTermPalindromeRoutingClaim>());
    }

    [Fact]
    public void RoutingClaim_TypedParents_ArePresentInRegistry()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();

        Assert.True(registry.Contains<F87TrichotomyClassification>(),
            "typed parent F87TrichotomyClassification must be registered");
        Assert.True(registry.Contains<F87DiagonalCellBipartiteWitnessSet>(),
            "typed parent F87DiagonalCellBipartiteWitnessSet must be registered");
        Assert.True(registry.Contains<CrossoverMirrorSqrtNinetyClaim>(),
            "typed parent CrossoverMirrorSqrtNinetyClaim must be registered");
    }

    [Fact]
    public void RoutingClaim_Ancestors_ContainAllThreeTypedParents()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();

        var ancestors = registry.AncestorsOf<TwoTermPalindromeRoutingClaim>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(F87TrichotomyClassification), ancestors);
        Assert.Contains(typeof(F87DiagonalCellBipartiteWitnessSet), ancestors);
        Assert.Contains(typeof(CrossoverMirrorSqrtNinetyClaim), ancestors);
    }

    [Fact]
    public void RoutingClaim_TierIsTier2Empirical()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        Assert.Equal(Tier.Tier2Empirical, registry.Get<TwoTermPalindromeRoutingClaim>().Tier);
    }

    [Fact]
    public void RoutingClaim_SelfCheck_AllPairsMatchAuthority()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var claim = registry.Get<TwoTermPalindromeRoutingClaim>();

        Assert.True(claim.AllPairsMatchAuthority,
            $"router disagreed with the spectral authority on {claim.PairCount - claim.MatchCount} of {claim.PairCount} pairs");
        Assert.Equal(claim.PairCount, claim.MatchCount);
    }

    [Fact]
    public void RoutingClaim_CoversTheFullPairSet_InclSelfPairs()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var claim = registry.Get<TwoTermPalindromeRoutingClaim>();

        // 9 single bilinears taken as unordered pairs incl. self-pairs = C(9,2) + 9 = 45.
        Assert.Equal(45, claim.PairCount);
    }

    [Fact]
    public void RoutingClaim_SurfacesTheSoftAndHardAnchors()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var claim = registry.Get<TwoTermPalindromeRoutingClaim>();

        // The canonical soft anchor XX+XZ routes to the Uniform Q-family, and the authority agrees.
        Assert.Equal(TrichotomyClass.Soft, claim.SoftExample.RoutedFate);
        Assert.Equal(QFamily.Uniform, claim.SoftExample.Family);
        Assert.Equal(TrichotomyClass.Soft, claim.SoftExample.AuthorityFate);

        // The hard anchor XY+XZ routes to None_ and carries its Q7 Klein-cell reason.
        Assert.Equal(TrichotomyClass.Hard, claim.HardExample.RoutedFate);
        Assert.Equal(QFamily.None_, claim.HardExample.Family);
        Assert.Equal(TrichotomyClass.Hard, claim.HardExample.AuthorityFate);
        Assert.False(string.IsNullOrWhiteSpace(claim.HardExample.Reason));
    }
}
