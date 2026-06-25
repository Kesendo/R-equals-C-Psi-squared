using System.Linq;
using RCPsiSquared.Core.F1;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Knowledge;
using RCPsiSquared.Runtime.ObjectManager;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>Wiring audit for <see cref="F89BranchLocusPalindromeClaim"/>: the path-3 octic branch locus
/// is a palindrome (mirror about Re λ = −4, forced by the F1 palindrome carried antiunitarily on the
/// block). Two typed parents, both Tier1Derived: <see cref="F1PalindromeIdentity"/> (the palindrome the
/// locus inherits) and <see cref="F89Path3OcticEpClaim"/> (the octic + the diabolic on the line).</summary>
public class F89BranchLocusPalindromeClaimRegistrationTests
{
    [Fact]
    public void BuildDefault_ContainsClaim()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        Assert.True(registry.Contains<F89BranchLocusPalindromeClaim>());
    }

    [Fact]
    public void Claim_IsTier1Derived()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        Assert.Equal(Tier.Tier1Derived, registry.Get<F89BranchLocusPalindromeClaim>().Tier);
    }

    [Fact]
    public void Claim_Ancestors_ContainBothTypedParents()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var ancestors = registry.AncestorsOf<F89BranchLocusPalindromeClaim>()
            .Select(c => c.GetType()).ToHashSet();
        Assert.Contains(typeof(F1PalindromeIdentity), ancestors);
        Assert.Contains(typeof(F89Path3OcticEpClaim), ancestors);
    }
}
