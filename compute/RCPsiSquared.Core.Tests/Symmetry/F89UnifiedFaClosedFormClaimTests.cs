using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class F89UnifiedFaClosedFormClaimTests
{
    private static F89UnifiedFaClosedFormClaim BuildClaim()
    {
        var f89 = new F89TopologyOrbitClosure(new Pi2DyadicLadderClaim());
        var atLock = new F89PathKAtLockMechanismClaim(f89);
        return new F89UnifiedFaClosedFormClaim(f89, atLock);
    }

    [Fact]
    public void Tier_IsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived, BuildClaim().Tier);
    }

    [Theory]
    [InlineData(3, new double[] { 47.0, 14.0 }, 9)]
    [InlineData(4, new double[] { 25.0, 10.0 }, 4)]
    [InlineData(5, new double[] { 129.0, 82.0, 13.0 }, 25)]
    [InlineData(6, new double[] { 80.0, 72.0, 17.0 }, 18)]
    public void PathPolynomial_MatchesEmpiricalCoefficients(int k, double[] expectedCoefsLowToHigh, int expectedDenom)
    {
        var (coefs, denom) = F89UnifiedFaClosedFormClaim.PathPolynomial(k);
        Assert.Equal(expectedCoefsLowToHigh.Length, coefs.Length);
        for (int i = 0; i < expectedCoefsLowToHigh.Length; i++)
            Assert.Equal(expectedCoefsLowToHigh[i], coefs[i], precision: 14);
        Assert.Equal(expectedDenom, denom);
    }

    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    public void PolynomialDegree_IsFloorNBlockHalfMinusOne(int k)
    {
        // poly degree = F_a count - 1 = floor(N_block/2) - 1 where N_block = k+1
        var (coefs, _) = F89UnifiedFaClosedFormClaim.PathPolynomial(k);
        int polyDegree = coefs.Length - 1;
        int expected = (k + 1) / 2 - 1;
        Assert.Equal(expected, polyDegree);
    }

    [Theory]
    // path-3, N=11, n=2 (y = √5−1 ≈ 1.236): sigs = (33+14√5)/[9·11²·10]
    [InlineData(3, 2, 11, 0.005904951)]
    // path-3, N=11, n=4 (y = -(1+√5) ≈ -3.236): sigs = (33-14√5)/[9·11²·10]
    [InlineData(3, 4, 11, 0.0001556521)]
    // path-4, N=11, n=2 (y = 2): sigs = 45/[4·11²·10]
    [InlineData(4, 2, 11, 0.009297520)]
    // path-4, N=11, n=4 (y = -2): sigs = 5/[4·11²·10]
    [InlineData(4, 4, 11, 0.001033058)]
    // path-5, N=11, n=2 (y = 4cos(2π/7) ≈ 2.494): sigs = (13y²+82y+129)/[25·11²·10]
    [InlineData(5, 2, 11, 0.0136979332)]
    // path-6, N=11, n=4 (y = 0, zero-mode): sigs = 80/[18·11²·10] = 40/[9·11²·10]
    [InlineData(6, 4, 11, 0.0036730946)]
    public void Sigma_AtPathNAndN_MatchesEmpirical(int k, int n, int blochN, double expected)
    {
        double sigs = F89UnifiedFaClosedFormClaim.Sigma(k, n, blochN);
        Assert.Equal(expected, sigs, precision: 8);
    }

    [Theory]
    [InlineData(3, 11, 0.0060606)]       // sum = 22 / [3·N²(N-1)] = 22/[3·121·10]
    [InlineData(4, 11, 0.0103306)]       // sum = 25 / [2·121·10]
    [InlineData(5, 11, 0.0159669)]       // sum = 483 / [25·121·10]
    [InlineData(6, 11, 0.0235078)]       // sum = 256 / [9·121·10]
    public void SigmaSum_AtPathAndN_MatchesRationalClosedForm(int k, int blochN, double expected)
    {
        double sum = F89UnifiedFaClosedFormClaim.SigmaSum(k, blochN);
        Assert.Equal(expected, sum, precision: 7);
    }

    [Fact]
    public void PathPolynomial_UnsupportedPath_Throws()
    {
        // path-8 is the first unsupported path now that path-7 is tabulated
        Assert.Throws<ArgumentOutOfRangeException>(
            () => F89UnifiedFaClosedFormClaim.PathPolynomial(8));
    }

    [Fact]
    public void Sigma_NTooSmallForBlock_Throws()
    {
        // path-3 needs N ≥ N_block + 1 = 5 to have at least one bare site
        Assert.Throws<ArgumentOutOfRangeException>(
            () => F89UnifiedFaClosedFormClaim.Sigma(3, 2, 4));
    }
}
