using System.Linq;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Knowledge;
using RCPsiSquared.Runtime.ObjectManager;
using Xunit;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

/// <summary>Wiring of <see cref="Spin32QuadrupolarMirrorsClaim"/> (F167): it lands in the default
/// registry as Tier1Derived with its typed parent PalindromeTwoEndCountClaim (F158, the second route
/// to the global-centre theorem); its closed forms read the page's centres.</summary>
public class Spin32QuadrupolarMirrorsClaimRegistrationTests
{
    private static ClaimRegistry Default() => KnowledgeRegistryFactory.BuildDefault();

    [Fact]
    public void BuildDefault_ContainsClaim() =>
        Assert.True(Default().Contains<Spin32QuadrupolarMirrorsClaim>());

    [Fact]
    public void TierIsTier1Derived() =>
        Assert.Equal(Tier.Tier1Derived, Default().Get<Spin32QuadrupolarMirrorsClaim>().Tier);

    [Fact]
    public void Ancestors_ContainTheTwoEndCount()
    {
        var ancestors = Default().AncestorsOf<Spin32QuadrupolarMirrorsClaim>()
            .Select(c => c.GetType()).ToHashSet();
        Assert.Contains(typeof(PalindromeTwoEndCountClaim), ancestors);
    }

    [Fact]
    public void DisplayName_IsF167_AndTheAnchorIsThePage()
    {
        var c = Default().Get<Spin32QuadrupolarMirrorsClaim>();
        Assert.StartsWith("F167: ", c.DisplayName);
        Assert.StartsWith("experiments/THE_SODIUM_NUCLEUS_TWO_MIRRORS.md", c.Anchor);
    }

    [Fact]
    public void Centres_ReadThePagesForms()
    {
        // c₁ = 3(J₁ + J₂), c₂ = c₁ + 3J₀, 2s = 3J₀ + 6J₁ + 6J₂; at J₀ = 0 the two centres coincide with s
        Assert.Equal(9, Spin32QuadrupolarMirrorsClaim.DarkCentre(1, 2));
        Assert.Equal(21, Spin32QuadrupolarMirrorsClaim.BrightCentre(4, 1, 2));
        Assert.Equal(30, Spin32QuadrupolarMirrorsClaim.TwiceGlobalCentre(4, 1, 2));
        Assert.Equal(2 * Spin32QuadrupolarMirrorsClaim.DarkCentre(1, 2), Spin32QuadrupolarMirrorsClaim.TwiceGlobalCentre(0, 1, 2));
    }
}
