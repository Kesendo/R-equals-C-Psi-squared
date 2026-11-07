using System.Linq;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Knowledge;
using RCPsiSquared.Runtime.ObjectManager;
using Xunit;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

/// <summary>Wiring of <see cref="QuditClockMirrorClaim"/> (F172): it lands in the default registry as Tier1Derived
/// with its typed parent PalindromeTwoEndCountClaim (F158); its closed forms read the page's rates and centres.</summary>
public class QuditClockMirrorClaimRegistrationTests
{
    private static ClaimRegistry Default() => KnowledgeRegistryFactory.BuildDefault();

    [Fact]
    public void BuildDefault_ContainsClaim() =>
        Assert.True(Default().Contains<QuditClockMirrorClaim>());

    [Fact]
    public void TierIsTier1Derived() =>
        Assert.Equal(Tier.Tier1Derived, Default().Get<QuditClockMirrorClaim>().Tier);

    [Fact]
    public void Ancestors_ContainTheTwoEndCount()
    {
        var ancestors = Default().AncestorsOf<QuditClockMirrorClaim>().Select(c => c.GetType()).ToHashSet();
        Assert.Contains(typeof(PalindromeTwoEndCountClaim), ancestors);
    }

    [Fact]
    public void DisplayName_IsF172_AndTheAnchorIsThePage()
    {
        var c = Default().Get<QuditClockMirrorClaim>();
        Assert.StartsWith("F172: ", c.DisplayName);
        Assert.StartsWith("docs/proofs/PROOF_QUDIT_CLOCK_MIRROR.md", c.Anchor);
    }

    [Fact]
    public void ClosedForms_ReadThePage()
    {
        // the immune column and the opposite one, where the cosine is exactly 1 and -1 in floating point; the half-shift
        // exists at even d only
        Assert.True(QuditClockMirrorClaim.PairRateOverGamma(0, 4) == 0.0);
        Assert.True(QuditClockMirrorClaim.PairRateOverGamma(2, 4) == 4.0);
        Assert.True(QuditClockMirrorClaim.PairRateOverGamma(3, 6) == 4.0);
        Assert.True(QuditClockMirrorClaim.HasHalfShift(6));
        Assert.False(QuditClockMirrorClaim.HasHalfShift(5));
        Assert.Equal(16, QuditClockMirrorClaim.PairTwiceCentre(2, 2));
        Assert.Equal(8, QuditClockMirrorClaim.LoneTwiceCentre(2, 2));
    }
}
