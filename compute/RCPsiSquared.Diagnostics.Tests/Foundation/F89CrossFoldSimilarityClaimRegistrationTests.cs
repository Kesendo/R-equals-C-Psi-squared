using System.Linq;
using RCPsiSquared.Core.F1;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Knowledge;
using RCPsiSquared.Runtime.ObjectManager;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>Wiring audit for <see cref="F89CrossFoldSimilarityClaim"/> (F89d): the (SE,DE)↔(SE,w_{N−2}) cross-
/// fold is an EXACT antiunitary similarity at the matrix level, so the diabolics pair across it with character
/// and gap preserved. Two typed parents, both Tier1Derived: <see cref="F1PalindromeIdentity"/> (the mirror the
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
    public void FoldImage_OfN7RealQDiabolic_IsMinus9p058()
    {
        // λ = −4.942 in (SE,DE) folds to −λ−2N = 4.942 − 14 = −9.058 in (SE,w_{N−2}) at N=7.
        Assert.Equal(-9.058, F89CrossFoldSimilarityClaim.FoldImageReal(-4.942, 7), 3);
        Assert.Equal(-7.0, F89CrossFoldSimilarityClaim.FoldCentre(1.0, 7), 12);
    }
}
