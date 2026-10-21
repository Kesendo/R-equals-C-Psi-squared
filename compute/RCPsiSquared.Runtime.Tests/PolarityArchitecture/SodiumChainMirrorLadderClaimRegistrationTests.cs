using System.Linq;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Knowledge;
using RCPsiSquared.Runtime.ObjectManager;
using Xunit;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

/// <summary>Wiring of <see cref="SodiumChainMirrorLadderClaim"/> (F167, the chain): it lands in the
/// default registry as Tier1Derived with its typed parents Spin32QuadrupolarMirrorsClaim (the n = 1
/// case) and PalindromeTwoEndCountClaim (F158, the full coupling's break); its closed forms read the
/// page's rungs.</summary>
public class SodiumChainMirrorLadderClaimRegistrationTests
{
    private static ClaimRegistry Default() => KnowledgeRegistryFactory.BuildDefault();

    [Fact]
    public void BuildDefault_ContainsClaim() =>
        Assert.True(Default().Contains<SodiumChainMirrorLadderClaim>());

    [Fact]
    public void TierIsTier1Derived() =>
        Assert.Equal(Tier.Tier1Derived, Default().Get<SodiumChainMirrorLadderClaim>().Tier);

    [Fact]
    public void Ancestors_ContainTheNucleusAndTheTwoEndCount()
    {
        var ancestors = Default().AncestorsOf<SodiumChainMirrorLadderClaim>()
            .Select(c => c.GetType()).ToHashSet();
        Assert.Contains(typeof(Spin32QuadrupolarMirrorsClaim), ancestors);
        Assert.Contains(typeof(PalindromeTwoEndCountClaim), ancestors);
    }

    [Fact]
    public void DisplayName_IsF167_AndTheAnchorIsTheChainPage()
    {
        var c = Default().Get<SodiumChainMirrorLadderClaim>();
        Assert.StartsWith("F167: the sodium chain", c.DisplayName);
        Assert.StartsWith("experiments/THE_SODIUM_CHAIN.md", c.Anchor);
    }

    [Fact]
    public void Rungs_ReadThePagesForms_AndTheSingleNucleusAtNEqualsOne()
    {
        // the page's units, J = (6, 4, 2) at two nuclei: s₀ = 36, rungs 3J₀ = 18 apart
        Assert.Equal(72, SodiumChainMirrorLadderClaim.TwiceRung(2, 0, 6, 4, 2));
        Assert.Equal(108, SodiumChainMirrorLadderClaim.TwiceRung(2, 1, 6, 4, 2));
        Assert.Equal(144, SodiumChainMirrorLadderClaim.TwiceRung(2, 2, 6, 4, 2));
        // rates per nucleus: 2s_B = 6Σ(J₁ + J₂) + 6Σ_B J₀
        var rates = new[] { (6L, 4L, 2L), (10L, 2L, 4L) };
        Assert.Equal(72, SodiumChainMirrorLadderClaim.LemmaConstant(rates));
        Assert.Equal(72 + 36, SodiumChainMirrorLadderClaim.TwiceCentre(rates, 1));
        Assert.Equal(72 + 60, SodiumChainMirrorLadderClaim.TwiceCentre(rates, 2));
        Assert.Equal(72 + 96, SodiumChainMirrorLadderClaim.TwiceCentre(rates, 3));
        // n = 1: the two rungs are the parent's dark and bright centres
        var one = new[] { (5L, 1L, 2L) };
        Assert.Equal(2 * Spin32QuadrupolarMirrorsClaim.DarkCentre(1, 2), SodiumChainMirrorLadderClaim.TwiceCentre(one, 0));
        Assert.Equal(2 * Spin32QuadrupolarMirrorsClaim.BrightCentre(5, 1, 2), SodiumChainMirrorLadderClaim.TwiceCentre(one, 1));
    }
}
