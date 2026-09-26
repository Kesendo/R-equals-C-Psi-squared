using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using Xunit;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class QuditMirrorProtectionScalingClaimTests
{
    private static QuditMirrorProtectionScalingClaim MakeClaim() => QuditMirrorProtectionScalingClaim.Build();

    [Fact]
    public void Tier_IsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived, MakeClaim().Tier);
    }

    [Fact]
    public void Parent_IsTheTier1DerivedCap()
    {
        var claim = MakeClaim();
        Assert.NotNull(claim.Cap);
        Assert.Equal(Tier.Tier1Derived, claim.Cap.Tier);
    }

    [Fact]
    public void Anchor_ReferencesProofAndQutritVerifier()
    {
        var anchor = MakeClaim().Anchor;
        Assert.Contains("PROOF_QUDIT_PARTIAL_PALINDROME.md", anchor);
        Assert.Contains("qudit_g2_split.py", anchor);
    }

    [Fact]
    public void Cap_IsTwoDToTheN()
    {
        // (2d)^N: d=2,N=2 -> 16; d=3,N=2 -> 36; d=3,N=3 -> 216; d=4,N=2 -> 64
        Assert.Equal(16L, QuditMirrorProtectionScalingClaim.ProductCap(2, 2));
        Assert.Equal(36L, QuditMirrorProtectionScalingClaim.ProductCap(3, 2));
        Assert.Equal(216L, QuditMirrorProtectionScalingClaim.ProductCap(3, 3));
        Assert.Equal(64L, QuditMirrorProtectionScalingClaim.ProductCap(4, 2));
        // From d = 6 the (0, 2) pairing beats the swap: 180 > 144 at N = 2.
        Assert.Equal(180L, QuditMirrorProtectionScalingClaim.ProductCap(6, 2));
    }

    [Fact]
    public void Total_IsDToTheTwoN()
    {
        // d^{2N}: d=2,N=2 -> 16; d=3,N=2 -> 81; d=3,N=3 -> 729; d=4,N=2 -> 256
        Assert.Equal(16L, QuditMirrorProtectionScalingClaim.TotalCoherences(2, 2));
        Assert.Equal(81L, QuditMirrorProtectionScalingClaim.TotalCoherences(3, 2));
        Assert.Equal(729L, QuditMirrorProtectionScalingClaim.TotalCoherences(3, 3));
        Assert.Equal(256L, QuditMirrorProtectionScalingClaim.TotalCoherences(4, 2));
    }

    [Fact]
    public void ProtectedFraction_IsTwoOverD_ToTheN_ExactlyForDUpToFive()
    {
        // Integer identity P(d, N)·d^N == 2^N·d^{2N}: the fraction is (2/d)^N with no division.
        for (int d = 2; d <= 5; d++)
            for (int n = 1; n <= 6; n++)
                Assert.True(QuditMirrorProtectionScalingClaim.FractionIsTwoOverDToTheN(d, n), $"d={d}, N={n}");
        // From d = 6 at N ≥ 2 the cap is larger, so the identity fails; N = 1 keeps it.
        Assert.True(QuditMirrorProtectionScalingClaim.FractionIsTwoOverDToTheN(6, 1));
        Assert.False(QuditMirrorProtectionScalingClaim.FractionIsTwoOverDToTheN(6, 2));
        Assert.False(QuditMirrorProtectionScalingClaim.FractionIsTwoOverDToTheN(7, 3));
    }

    [Fact]
    public void FullMirror_HoldsExactlyAtTheQubit()
    {
        // protected fraction = 1 ⟺ d = 2; below 1 for every qudit (the trunk root d²−2d=0)
        Assert.True(QuditMirrorProtectionScalingClaim.IsFullMirror(2));
        Assert.False(QuditMirrorProtectionScalingClaim.IsFullMirror(3));
        Assert.False(QuditMirrorProtectionScalingClaim.IsFullMirror(4));

        foreach (var n in new[] { 1, 2, 3, 4 })
            Assert.True(QuditMirrorProtectionScalingClaim.ProductCap(2, n) == QuditMirrorProtectionScalingClaim.TotalCoherences(2, n));

        foreach (var (d, n) in new[] { (3, 2), (3, 3), (4, 2), (5, 2), (6, 2), (7, 3) })
            Assert.True(QuditMirrorProtectionScalingClaim.ProtectedFraction(d, n) < 1.0,
                $"qudit d={d}, N={n} should be below full protection");
    }

    [Fact]
    public void ProtectedFraction_DecaysExponentiallyInN_AtTheQutrit()
    {
        // (2/3)^N strictly decreasing in N
        double f2 = QuditMirrorProtectionScalingClaim.ProtectedFraction(3, 2);
        double f3 = QuditMirrorProtectionScalingClaim.ProtectedFraction(3, 3);
        double f4 = QuditMirrorProtectionScalingClaim.ProtectedFraction(3, 4);
        Assert.True(f3 < f2 && f4 < f3, "qutrit protection (2/3)^N must strictly decrease with N");
    }
}
