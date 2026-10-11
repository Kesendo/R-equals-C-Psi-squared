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
