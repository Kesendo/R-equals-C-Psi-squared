using System.Numerics;
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
    [InlineData(8, new double[] { 110.0, 68.0, 54.0, 13.0 }, 32)]
    [InlineData(9, new double[] { 1476.0, 440.0, 288.0, 190.0, 31.0 }, 324)]
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
        // path-47 is the first unsupported path: k=3..46 tabulated via Chebyshev pipeline.
        // k=47 would require long-typed Denominator (D_47 = 4,632,608,768 > int.MaxValue).
        Assert.Throws<ArgumentOutOfRangeException>(
            () => F89UnifiedFaClosedFormClaim.PathPolynomial(47));
    }

    [Theory]
    [InlineData(10, 200)]
    [InlineData(11, 968)]
    [InlineData(12, 288)]
    [InlineData(15, 7200)]
    [InlineData(16, 2048)]
    [InlineData(20, 12800)]
    [InlineData(24, 73728)]
    [InlineData(25, 640000)]
    [InlineData(28, 401408)]
    [InlineData(31, 7872512)]
    [InlineData(32, 2097152)]
    [InlineData(33, 17842176)]
    [InlineData(40, 52428800)]
    [InlineData(45, 2123366400)]
    [InlineData(46, 1109393408)]
    public void PathPolynomial_AtExtendedPath_DenominatorMatchesPredictDenominator(int k, int expectedD)
    {
        var (_, d) = F89UnifiedFaClosedFormClaim.PathPolynomial(k);
        Assert.Equal(expectedD, d);
        Assert.Equal(F89UnifiedFaClosedFormClaim.PredictDenominator(k), d);
    }

    [Fact]
    public void Sigma_NTooSmallForBlock_Throws()
    {
        // path-3 needs N ≥ N_block + 1 = 5 to have at least one bare site
        Assert.Throws<ArgumentOutOfRangeException>(
            () => F89UnifiedFaClosedFormClaim.Sigma(3, 2, 4));
    }

    [Theory]
    // path-10, N=12, n=2 (y = 4cos(π/6) = 2√3 ≈ 3.464): coefs (1120, 128, 32, 152, 37), D=200.
    // P(y) = 37y⁴ + 152y³ + 32y² + 128y + 1120 ; σ = P(y_2)/(200·12²·11).
    [InlineData(10, 2, 12)]
    [InlineData(11, 2, 13)]
    [InlineData(15, 4, 17)]
    [InlineData(20, 6, 22)]
    [InlineData(32, 4, 34)]
    [InlineData(46, 2, 48)]
    public void Sigma_BeyondK9_LiftedRestriction_ProducesValue(int k, int n, int blochN)
    {
        // Restriction k ∈ {3..9} lifted 2026-05-15; Sigma now works for any k ≥ 3.
        double sigma = F89UnifiedFaClosedFormClaim.Sigma(k, n, blochN);
        Assert.True(sigma > 0, $"sigma path-{k} n={n} N={blochN} should be positive; got {sigma}");
        Assert.True(double.IsFinite(sigma), "sigma should be finite");
    }

    [Theory]
    [InlineData(50, 2, 52)]
    [InlineData(60, 4, 62)]
    public void Sigma_BeyondTabulationViaPipeline_Works(int k, int n, int blochN)
    {
        double sigma = F89UnifiedFaClosedFormClaim.Sigma(k, n, blochN);
        Assert.True(sigma > 0);
        Assert.True(double.IsFinite(sigma));
    }

    [Theory]
    [InlineData(10, 12)]
    [InlineData(20, 22)]
    [InlineData(32, 34)]
    public void SigmaSum_BeyondK9_LiftedRestriction_ProducesValue(int k, int blochN)
    {
        double sum = F89UnifiedFaClosedFormClaim.SigmaSum(k, blochN);
        Assert.True(sum > 0);
        Assert.True(double.IsFinite(sum));
    }

    [Theory]
    [InlineData(47)]
    [InlineData(50)]
    [InlineData(100)]
    public void PredictDenominatorBig_MatchesPipelineDenominator(int k)
    {
        var pipelineD = F89UnifiedFaClosedFormClaim.ComputePathPolynomialBig(k).Denominator;
        var formulaD = F89UnifiedFaClosedFormClaim.PredictDenominatorBig(k);
        Assert.Equal(pipelineD, formulaD);
    }

    [Theory]
    [InlineData(3, 9)]
    [InlineData(10, 200)]
    [InlineData(46, 1109393408)]
    public void PredictDenominatorBig_MatchesIntPredictDenominator(int k, int expected)
    {
        Assert.Equal(new BigInteger(expected), F89UnifiedFaClosedFormClaim.PredictDenominatorBig(k));
    }

    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    [InlineData(8)]
    [InlineData(9)]
    public void PredictDenominator_MatchesTabulatedPathPolynomial(int k)
    {
        var (_, tabulatedD) = F89UnifiedFaClosedFormClaim.PathPolynomial(k);
        int predictedD = F89UnifiedFaClosedFormClaim.PredictDenominator(k);
        Assert.Equal(tabulatedD, predictedD);
    }

    [Theory]
    [InlineData(10, 200)]      // k=10: 2³·5² = 200
    [InlineData(11, 968)]      // k=11: 2³·11² = 968
    [InlineData(12, 288)]      // k=12: 2⁵·3² = 288
    [InlineData(13, 2704)]     // k=13: 2⁴·13² = 2704
    [InlineData(14, 1568)]     // k=14: 2⁵·7² = 1568
    [InlineData(15, 7200)]     // k=15: 2⁵·3²·5² = 7200
    [InlineData(16, 2048)]     // k=16: 2¹¹·1 = 2048
    [InlineData(17, 18496)]    // k=17: 2⁶·17² = 18496
    [InlineData(24, 73728)]    // k=24: 2¹³·3² = 73728
    public void PredictDenominator_BeyondTabulated_MatchesProbeExtraction(int k, int expected)
    {
        int predicted = F89UnifiedFaClosedFormClaim.PredictDenominator(k);
        Assert.Equal(expected, predicted);
    }
}
