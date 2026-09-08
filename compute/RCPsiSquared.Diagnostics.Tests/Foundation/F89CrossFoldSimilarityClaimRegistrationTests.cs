using System.Linq;
using RCPsiSquared.Core.F1;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Knowledge;
using RCPsiSquared.Runtime.ObjectManager;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>Wiring audit for <see cref="F89CrossFoldSimilarityClaim"/> (F89d): the (SE,DE)↔(SE,w_{N−2}) cross-
/// fold is an exact antiunitary similarity at the matrix level, so independently certified Jordan character and
/// gap transport across it. The similarity does not certify sampled candidates. Two typed parents, both
/// Tier1Derived: <see cref="F1PalindromeIdentity"/> (the mirror the
/// fold realises) and <see cref="F89BranchLocusPalindromeClaim"/> (the spectrum-level fold this upgrades).</summary>
public class F89CrossFoldSimilarityClaimRegistrationTests
{
    [Fact]
    public void BuildDefault_ContainsClaim()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        Assert.True(registry.Contains<F89CrossFoldSimilarityClaim>());
    }

    [Fact]
    public void Claim_IsTier1Derived()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        Assert.Equal(Tier.Tier1Derived, registry.Get<F89CrossFoldSimilarityClaim>().Tier);
    }

    [Fact]
    public void Claim_Ancestors_ContainBothTypedParents()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var ancestors = registry.AncestorsOf<F89CrossFoldSimilarityClaim>()
            .Select(c => c.GetType()).ToHashSet();
        Assert.Contains(typeof(F1PalindromeIdentity), ancestors);
        Assert.Contains(typeof(F89BranchLocusPalindromeClaim), ancestors);
    }

    [Fact]
    public void FoldImage_OfN7SampledProposal_IsMinus9p058()
    {
        // The affine fold sends the sampled λ proposal to −λ−2N; this arithmetic does not certify character.
        Assert.Equal(-9.058, F89CrossFoldSimilarityClaim.FoldImageReal(-4.942, 7), 3);
        Assert.Equal(-7.0, F89CrossFoldSimilarityClaim.FoldCentre(1.0, 7), 12);
    }

    [Fact]
    public void Claim_SeparatesExactSimilarityFromUncertifiedPositiveDeltaCharacter()
    {
        var claim = KnowledgeRegistryFactory.BuildDefault().Get<F89CrossFoldSimilarityClaim>();
        var surface = $"{claim.Name} {claim.DisplayName} {claim.Summary}";

        Assert.Contains("N=4", surface);
        Assert.Contains("N=5/N=6", surface);
        Assert.Contains("positive Delta", surface, System.StringComparison.OrdinalIgnoreCase);
        Assert.Contains("Uncertified", surface);
        Assert.Contains("conditional", surface, System.StringComparison.OrdinalIgnoreCase);
        Assert.DoesNotContain("defect-or-lift", surface, System.StringComparison.OrdinalIgnoreCase);
        Assert.DoesNotContain("N=7 real-q diabolic", surface, System.StringComparison.OrdinalIgnoreCase);
    }
}
