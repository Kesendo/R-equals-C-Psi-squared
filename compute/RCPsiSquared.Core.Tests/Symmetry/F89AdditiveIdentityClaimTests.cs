using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class F89AdditiveIdentityClaimTests
{
    private static F89AdditiveIdentityClaim BuildClaim() =>
        new F89AdditiveIdentityClaim(new F89TopologyOrbitClosure(new Pi2DyadicLadderClaim()));

    [Fact]
    public void Tier_IsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived, BuildClaim().Tier);
    }

    [Theory]
    [InlineData(7, 0.0, 6.0 / 49.0)]      // N=7, t=0: (N-1)/N² = 6/49
    [InlineData(5, 0.0, 4.0 / 25.0)]      // N=5, t=0: 4/25
    [InlineData(11, 0.0, 10.0 / 121.0)]   // N=11, t=0
    public void BarePerSite_AtTimeZero_MatchesClosedForm(int n, double t, double expected)
    {
        Assert.Equal(expected, F89AdditiveIdentityClaim.BarePerSite(n, 0.05, t), precision: 14);
    }

    [Fact]
    public void BarePerSite_DecaysAt4Gamma()
    {
        // S_bare(t) / S_bare(0) = exp(-4γt). At γ=0.05, t=10: exp(-2) = 0.13534
        double ratio = F89AdditiveIdentityClaim.BarePerSite(7, 0.05, 10)
                     / F89AdditiveIdentityClaim.BarePerSite(7, 0.05, 0);
        Assert.Equal(Math.Exp(-2.0), ratio, precision: 12);
    }

    [Theory]
    [InlineData(1, 7, 0)]    // m=1, N=7: 0·N = 0
    [InlineData(2, 7, 7)]    // m=2, N=7: 1·N = 7
    [InlineData(3, 7, 14)]   // m=3, N=7: 2·N = 14
    [InlineData(2, 5, 5)]    // m=2, N=5: 1·5 = 5
    public void OvercountingCoefficient_IsCorrect(int m, int n, int expected)
    {
        Assert.Equal(expected, F89AdditiveIdentityClaim.OvercountingCoefficient(m, n));
    }

    [Theory]
    [InlineData(new int[] { 1 }, 7, 5)]            // (1) at N=7: 1 block of 2 sites, 5 bare
    [InlineData(new int[] { 2 }, 7, 4)]            // (2) at N=7: 1 block of 3 sites, 4 bare
    [InlineData(new int[] { 1, 2 }, 7, 2)]         // (1,2) at N=7: 2+3=5 block sites, 2 bare
    [InlineData(new int[] { 1, 1, 1 }, 7, 1)]      // (1,1,1) at N=7: 6 block sites, 1 bare
    [InlineData(new int[] { 6 }, 7, 0)]            // (6) at N=7: 7 block sites, 0 bare
    public void BareSiteCount_IsCorrect(int[] kValues, int n, int expected)
    {
        Assert.Equal(expected, F89AdditiveIdentityClaim.BareSiteCount(kValues, n));
    }

    [Fact]
    public void BareSiteCount_TopologyExceedsN_Throws()
    {
        // (3, 3) at N=5 needs 4+4 = 8 sites > 5
        Assert.Throws<ArgumentException>(
            () => F89AdditiveIdentityClaim.BareSiteCount(new[] { 3, 3 }, 5));
    }

    [Fact]
    public void Combine_AllIsolatedAtN7_MatchesS1Directly()
    {
        // (1) topology: m=1, so subtraction coefficient = (m-1)·N = 0.
        // Combine should equal S_(1)(t) exactly.
        // We use a stub S_path_k that returns a known value for k=1.
        double sP1(int k, int N, double g, double t)
        {
            // S_(1) at N=7, t=0: (N-1)/N = 6/7
            Assert.Equal(1, k);
            return (double)(N - 1) / N * Math.Exp(-4.0 * g * t);
        }

        var result = F89AdditiveIdentityClaim.Combine(new[] { 1 }, 7, 0.05, 0.0, sP1);
        Assert.Equal(6.0 / 7.0, result, precision: 14);
    }

    [Fact]
    public void Combine_MixedTopologyAtN7_AppliesSubtraction()
    {
        // (1, 2) at N=7: m=2, subtraction = 1·7·(N-1)/N²·exp(-4γt) = 7·6/49 = 6/7 at t=0.
        // S_(1)(0) = 6/7, S_(2)(0) = 6/7 (both at N=7).
        // Combine = 6/7 + 6/7 - 6/7 = 6/7. ✓ (every S(0) equals (N-1)/N by S(0) closed form)
        double sPathK(int k, int N, double g, double t)
        {
            return (double)(N - 1) / N * Math.Exp(-4.0 * g * t);
        }

        var result = F89AdditiveIdentityClaim.Combine(new[] { 1, 2 }, 7, 0.05, 0.0, sPathK);
        Assert.Equal(6.0 / 7.0, result, precision: 12);
    }

    [Fact]
    public void BarePerSite_NegativeGamma_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(
            () => F89AdditiveIdentityClaim.BarePerSite(7, -0.05, 0));
    }

    [Fact]
    public void BarePerSite_NLessThan2_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(
            () => F89AdditiveIdentityClaim.BarePerSite(1, 0.05, 0));
    }
}
