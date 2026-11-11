using System.Linq;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Knowledge;
using RCPsiSquared.Runtime.ObjectManager;
using Xunit;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

/// <summary>Wiring of <see cref="SignedMaskBalanceClaim"/> (F173): it lands in the default registry as Tier1Derived
/// with its typed parent PalindromeComplementConnectionClaim (Theorem 2), and its anchor is the page.</summary>
public class SignedMaskBalanceClaimRegistrationTests
{
    private static ClaimRegistry Default() => KnowledgeRegistryFactory.BuildDefault();

    [Fact]
    public void BuildDefault_ContainsClaim() =>
        Assert.True(Default().Contains<SignedMaskBalanceClaim>());

    [Fact]
    public void TierIsTier1Derived() =>
        Assert.Equal(Tier.Tier1Derived, Default().Get<SignedMaskBalanceClaim>().Tier);

    [Fact]
    public void Ancestors_ContainTheComplementConnection()
    {
        var ancestors = Default().AncestorsOf<SignedMaskBalanceClaim>().Select(c => c.GetType()).ToHashSet();
        Assert.Contains(typeof(PalindromeComplementConnectionClaim), ancestors);
    }

    [Fact]
    public void DisplayName_IsF173_AndTheAnchorIsThePage()
    {
        var c = Default().Get<SignedMaskBalanceClaim>();
        Assert.StartsWith("F173: ", c.DisplayName);
        Assert.StartsWith("docs/proofs/PROOF_SIGNED_MASK_BALANCE.md", c.Anchor);
    }
}
