using System.Linq;
using RCPsiSquared.Core.F1;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Knowledge;
using RCPsiSquared.Runtime.ObjectManager;
using Xunit;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

/// <summary>Wiring of <see cref="CorrelatedZDephasingPriceClaim"/>: it lands in the default registry as
/// Tier1Derived with its two typed parents, F1PalindromeIdentity (the identity it generalises) and
/// AbsorptionTheoremClaim (the cell rates the shared field reads); its closed forms read the page's
/// numbers.</summary>
public class CorrelatedZDephasingPriceClaimRegistrationTests
{
    private static ClaimRegistry Default() => KnowledgeRegistryFactory.BuildDefault();

    [Fact]
    public void BuildDefault_ContainsClaim() =>
        Assert.True(Default().Contains<CorrelatedZDephasingPriceClaim>());

    [Fact]
    public void TierIsTier1Derived() =>
        Assert.Equal(Tier.Tier1Derived, Default().Get<CorrelatedZDephasingPriceClaim>().Tier);

    [Fact]
    public void Ancestors_ContainBothParents()
    {
        var ancestors = Default().AncestorsOf<CorrelatedZDephasingPriceClaim>()
            .Select(c => c.GetType()).ToHashSet();
        Assert.Contains(typeof(F1PalindromeIdentity), ancestors);
        Assert.Contains(typeof(AbsorptionTheoremClaim), ancestors);
    }

    [Fact]
    public void BlockSigma_ReadsThePagesExamples_TheNegativeOfThePairCentre()
    {
        // N = 4, uniform rate gamma = 10 (sigma = 40), eta = 3/10 (kappa = 3): sigma(2,2) = 28 = 2.8 gamma,
        // sigma(0,4) = 76 = 7.6 gamma; local noise (kappa = 0) gives sigma itself on every block.
        Assert.Equal(28.0, CorrelatedZDephasingPriceClaim.BlockSigma(4, 40.0, 3.0, 2, 2));
        Assert.Equal(76.0, CorrelatedZDephasingPriceClaim.BlockSigma(4, 40.0, 3.0, 0, 4));
        Assert.Equal(40.0, CorrelatedZDephasingPriceClaim.BlockSigma(4, 40.0, 0.0, 1, 3));
    }

    [Fact]
    public void ResidualNorm_ReadsTheSharedFieldForm()
    {
        // uniform shared field: sum_{k<l} a_kl^2 = N(N-1)/2 (gamma eta)^2, so 8 4^N times it is
        // 4^(N+1) N (N-1) (gamma eta)^2; at N = 3, gamma eta = 3: 4^4 * 6 * 9 = 13824.
        Assert.Equal(13824.0, CorrelatedZDephasingPriceClaim.ResidualNormSquared(3, 3 * 9.0));
    }

    private static long[,] Covariances(int n, System.Func<int, int, long> c)
    {
        var m = new long[n, n];
        for (int k = 0; k < n; k++)
            for (int l = 0; l < n; l++)
                if (k != l) m[k, l] = c(System.Math.Min(k, l), System.Math.Max(k, l));
        return m;
    }

    [Fact]
    public void SingleValuePopcounts_ReadsThePagesPatternsAndTheFourShapes()
    {
        // the page's N = 4 patterns: two disjoint shared pairs (equal row sums), the sign pattern s = (1, 1, 1, -1)
        // (additive), a pattern decaying with distance (neither); one common covariance qualifies everything
        Assert.Equal(new[] { 0, 1, 3, 4 }, CorrelatedZDephasingPriceClaim.SingleValuePopcounts(
            Covariances(4, (k, l) => (k, l) is (0, 1) or (2, 3) ? 3 : 0)));
        var s = new long[] { 1, 1, 1, -1 };
        Assert.Equal(new[] { 0, 2, 4 }, CorrelatedZDephasingPriceClaim.SingleValuePopcounts(
            Covariances(4, (k, l) => 3 * s[k] * s[l])));
        Assert.Equal(new[] { 0, 4 }, CorrelatedZDephasingPriceClaim.SingleValuePopcounts(
            Covariances(4, (k, l) => (l - k) switch { 1 => 3, 2 => 1, _ => 0 })));
        Assert.Equal(new[] { 0, 1, 2, 3, 4 }, CorrelatedZDephasingPriceClaim.SingleValuePopcounts(
            Covariances(4, (k, l) => 3)));
        // N = 6: a ring keeps the edge rows, a star the middle, the two together neither
        Assert.Equal(new[] { 0, 1, 5, 6 }, CorrelatedZDephasingPriceClaim.SingleValuePopcounts(
            Covariances(6, (k, l) => l - k == 1 || (k, l) == (0, 5) ? 1 : 0)));
        Assert.Equal(new[] { 0, 3, 6 }, CorrelatedZDephasingPriceClaim.SingleValuePopcounts(
            Covariances(6, (k, l) => k == 0 ? 1 : 0)));
        Assert.Equal(new[] { 0, 6 }, CorrelatedZDephasingPriceClaim.SingleValuePopcounts(
            Covariances(6, (k, l) => (k == 0 ? 1 : 0) + (l - k == 1 || (k, l) == (0, 5) ? 1 : 0))));
        // N = 3: equal row sums already force one number; N = 5: the middle does not exist
        Assert.Equal(new[] { 0, 3 }, CorrelatedZDephasingPriceClaim.SingleValuePopcounts(
            Covariances(3, (k, l) => k + l)));
        Assert.Equal(new[] { 0, 5 }, CorrelatedZDephasingPriceClaim.SingleValuePopcounts(
            Covariances(5, (k, l) => k == 0 ? 2 : 0)));
    }

    [Fact]
    public void ThirdMomentNoiseTerm_IsTheTriangleForm_PositiveOnAValidField_NegativeOnAnIndefiniteOne()
    {
        // 4 tr(C (2a + Lambda) C) = 12 sum_{k<l} a_kl^2 (a_kk + a_ll) + 48 sum_triangles, here at N = 4 by hand
        var a = new long[,] { { 3, 1, -2, 0 }, { 1, 5, 4, -1 }, { -2, 4, 2, 3 }, { 0, -1, 3, 7 } };
        long edges = 0, triangles = 0;
        for (int k = 0; k < 4; k++)
            for (int l = k + 1; l < 4; l++)
            {
                edges += a[k, l] * a[k, l] * (a[k, k] + a[l, l]);
                for (int m = l + 1; m < 4; m++) triangles += a[k, l] * a[l, m] * a[k, m];
            }
        Assert.Equal(new System.Numerics.BigInteger(12 * edges + 48 * triangles), CorrelatedZDephasingPriceClaim.ThirdMomentNoiseTerm(a));
        // a positive definite a with covariances (eigenvalues about 0.64, 3.31, 7.04): positive
        var psd = new long[,] { { 6, 2, -1 }, { 2, 3, 1 }, { -1, 1, 2 } };
        Assert.True(CorrelatedZDephasingPriceClaim.ThirdMomentNoiseTerm(psd) > 0);
        // the page's indefinite example: zero diagonal, every covariance -1 at N = 3
        var ind = new long[,] { { 0, -1, -1 }, { -1, 0, -1 }, { -1, -1, 0 } };
        Assert.Equal(new System.Numerics.BigInteger(-48), CorrelatedZDephasingPriceClaim.ThirdMomentNoiseTerm(ind));
    }

    [Fact]
    public void ThirdMoment_ReadsThePagesNumbers()
    {
        // the verifier's end-site case: N = 3, rates 10, covariance -3 between the end sites, the open XY chain at J = 5,
        // Tr(H^2 Z_k Z_l) = 2^3 (-2 J^2) = -400 on the bonds and 0 on the end pair: Tr((2L + 2 sigma)^3) = -1105920
        var aEnd = new long[,] { { 10, 0, -3 }, { 0, 10, 0 }, { -3, 0, 10 } };
        var t2 = new long[,] { { 0, -400, 0 }, { -400, 0, -400 }, { 0, -400, 0 } };
        Assert.Equal(new System.Numerics.BigInteger(-1105920 / 8),
            CorrelatedZDephasingPriceClaim.ThirdMoment(aEnd, t2, new long[3]));
        // the shared field on the open XY chain at N = 4, gamma = 10, kappa = -2, J = 3 against the page's closed form
        // -4^N 4 (N - 1) kappa [2N(N - 2) kappa^2 + 3N gamma kappa + 3J^2] = -915456
        var shared = new long[4, 4];
        var t2c = new long[4, 4];
        for (int k = 0; k < 4; k++)
            for (int l = 0; l < 4; l++)
            {
                shared[k, l] = k == l ? 10 : -2;
                if (System.Math.Abs(k - l) == 1) t2c[k, l] = 16 * (-2 * 9);
            }
        Assert.Equal(new System.Numerics.BigInteger(-915456), CorrelatedZDephasingPriceClaim.ThirdMoment(shared, t2c, new long[4]));
        // local noise with single-site fields: F117's face 6 4^N sum gamma_k c_k^2 with Tr(H Z_k) = 2^N c_k
        var local = new long[,] { { 2, 0, 0 }, { 0, 5, 0 }, { 0, 0, 1 } };
        var c = new long[] { 1, -2, 3 };
        var t1 = new long[] { 8 * c[0], 8 * c[1], 8 * c[2] };
        Assert.Equal(new System.Numerics.BigInteger(6 * 64 * (2 * 1 + 5 * 4 + 1 * 9)),
            CorrelatedZDephasingPriceClaim.ThirdMoment(local, t2, t1));
    }

    [Fact]
    public void FarEndBreakPopcounts_ReadThePagesCases()
    {
        // one common covariance: the corners for kappa > 0 (the GHZ coherence), the middle for kappa < 0
        Assert.Equal(new[] { 0, 4 }, CorrelatedZDephasingPriceClaim.FarEndBreakPopcounts(Covariances(4, (k, l) => 3)));
        Assert.Equal(new[] { 2 }, CorrelatedZDephasingPriceClaim.FarEndBreakPopcounts(Covariances(4, (k, l) => -2)));
        Assert.Equal(new[] { 2, 3 }, CorrelatedZDephasingPriceClaim.FarEndBreakPopcounts(Covariances(5, (k, l) => -2)));
        Assert.Equal(new[] { 0, 1, 4, 5 }, CorrelatedZDephasingPriceClaim.FarEndBreakPopcounts(Covariances(5, (k, l) => 3)));
        // the verifier's additive pattern y = (-1, -1, 0, 0) at N = 4: summed covariance -6, the middle breaks
        var y = new long[] { -1, -1, 0, 0 };
        Assert.Equal(new[] { 2 }, CorrelatedZDephasingPriceClaim.FarEndBreakPopcounts(Covariances(4, (k, l) => y[k] + y[l])));
        // the end-site covariance: only 0 and N qualify, and there K lies below sigma, so nothing breaks here
        Assert.Empty(CorrelatedZDephasingPriceClaim.FarEndBreakPopcounts(Covariances(3, (k, l) => (k, l) == (0, 2) ? -3 : 0)));
        Assert.Empty(CorrelatedZDephasingPriceClaim.FarEndBreakPopcounts(Covariances(5, (k, l) => (k, l) == (0, 4) ? -3 : 0)));
        // no covariance: nothing breaks
        Assert.Empty(CorrelatedZDephasingPriceClaim.FarEndBreakPopcounts(new long[4, 4]));
    }

    [Fact]
    public void ThirdMoment_RejectsMismatchedTablesAndAnAsymmetricMatrix()
    {
        var a = new long[,] { { 1, 0 }, { 0, 1 } };
        Assert.Throws<System.ArgumentException>(() => CorrelatedZDephasingPriceClaim.ThirdMoment(a, new long[3, 3], new long[2]));
        Assert.Throws<System.ArgumentException>(() => CorrelatedZDephasingPriceClaim.ThirdMoment(a, new long[2, 2], new long[3]));
        var asym = new long[,] { { 1, 2 }, { 0, 1 } };
        Assert.Throws<System.ArgumentException>(() => CorrelatedZDephasingPriceClaim.ThirdMomentNoiseTerm(asym));
    }

    [Fact]
    public void SingleValuePopcounts_RejectsAnAsymmetricOrNonSquareOrOversizedMatrix()
    {
        var m = new long[3, 3];
        m[0, 1] = 1;
        Assert.Throws<System.ArgumentException>(() => CorrelatedZDephasingPriceClaim.SingleValuePopcounts(m));
        Assert.Throws<System.ArgumentException>(() => CorrelatedZDephasingPriceClaim.SingleValuePopcounts(new long[3, 4]));
        var big = new long[3, 3];
        big[0, 1] = big[1, 0] = long.MinValue;
        Assert.Throws<System.ArgumentOutOfRangeException>(() => CorrelatedZDephasingPriceClaim.SingleValuePopcounts(big));
        var diag = new long[3, 3];
        diag[0, 0] = 5_000_000_000L;                          // the diagonal is ignored, so a large rate is no error
        Assert.Equal(new[] { 0, 1, 2, 3 }, CorrelatedZDephasingPriceClaim.SingleValuePopcounts(diag));
    }
}
