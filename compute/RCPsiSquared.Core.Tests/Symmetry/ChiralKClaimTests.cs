using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using Xunit;

namespace RCPsiSquared.Core.Tests.Symmetry;

/// <summary>Witnesses for <see cref="ChiralKClaim"/>: the typed-Claim wrapper around the
/// K sublattice/chiral symmetry. Anchors the claim's Tier1Derived status and confirms it
/// is a sibling root (no typed parent) to <see cref="ChiralK"/>'s operator helper and
/// to <c>PROOF_K_PARTNERSHIP.md</c>.</summary>
public class ChiralKClaimTests
{
    [Fact]
    public void Build_Tier_IsTier1Derived()
    {
        var claim = new ChiralKClaim();
        Assert.Equal(Tier.Tier1Derived, claim.Tier);
    }

    [Fact]
    public void Build_Anchor_CitesProofAndHelper()
    {
        var claim = new ChiralKClaim();
        Assert.Contains("PROOF_K_PARTNERSHIP", claim.Anchor);
        Assert.Contains("ChiralK.cs", claim.Anchor);
    }

    [Fact]
    public void Build_Summary_NamesTheAntisymmetryAndSpectrumInversion()
    {
        var claim = new ChiralKClaim();
        Assert.Contains("K·H·K = −H", claim.Summary);
        Assert.Contains("E_{N+1−k} = −E_k", claim.Summary);
    }
}
