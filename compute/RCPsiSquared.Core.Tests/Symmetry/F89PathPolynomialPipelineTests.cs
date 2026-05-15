using System.Numerics;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Symmetry;

/// <summary>Native C# F89 path polynomial pipeline (Chebyshev expansion + orbit
/// polynomial reduction) bit-exact match against the cached
/// <see cref="F89UnifiedFaClosedFormClaim.PathPolynomial"/> tabulation k=3..46,
/// plus computation past the int.MaxValue boundary k ≥ 47.</summary>
public class F89PathPolynomialPipelineTests
{
    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    [InlineData(8)]
    [InlineData(9)]
    [InlineData(10)]
    [InlineData(15)]
    [InlineData(20)]
    [InlineData(24)]
    [InlineData(25)]
    [InlineData(30)]
    [InlineData(32)]
    [InlineData(33)]
    [InlineData(40)]
    [InlineData(45)]
    [InlineData(46)]
    public void Compute_MatchesTabulatedPathPolynomial_BitExact(int k)
    {
        var (tabulatedCoefs, tabulatedD) = F89UnifiedFaClosedFormClaim.PathPolynomial(k);
        var (pipelineCoefs, pipelineD) = F89PathPolynomialPipeline.Compute(k);

        Assert.Equal(new BigInteger(tabulatedD), pipelineD);
        Assert.Equal(tabulatedCoefs.Length, pipelineCoefs.Length);
        for (int i = 0; i < tabulatedCoefs.Length; i++)
        {
            var expected = new BigInteger((long)tabulatedCoefs[i]);
            Assert.Equal(expected, pipelineCoefs[i]);
        }
    }

    [Fact]
    public void Compute_AllTabulatedPaths_BitExactBatch()
    {
        for (int k = 3; k <= 46; k++)
        {
            var (tabulatedCoefs, tabulatedD) = F89UnifiedFaClosedFormClaim.PathPolynomial(k);
            var (pipelineCoefs, pipelineD) = F89PathPolynomialPipeline.Compute(k);
            Assert.Equal(new BigInteger(tabulatedD), pipelineD);
            Assert.Equal(tabulatedCoefs.Length, pipelineCoefs.Length);
            for (int i = 0; i < tabulatedCoefs.Length; i++)
                Assert.Equal(new BigInteger((long)tabulatedCoefs[i]), pipelineCoefs[i]);
        }
    }

    [Theory]
    [InlineData(47)]
    [InlineData(48)]
    [InlineData(50)]
    public void Compute_BeyondTabulation_ProducesValidResult(int k)
    {
        var (coefs, d) = F89PathPolynomialPipeline.Compute(k);
        int FA = (k + 1) / 2;
        Assert.Equal(FA, coefs.Length);
        Assert.True(d > BigInteger.Zero);
        Assert.Equal(BigInteger.One, BigInteger.GreatestCommonDivisor(
            BigInteger.GreatestCommonDivisor(BigInteger.Abs(coefs[0]), BigInteger.Abs(coefs[^1])), d));
    }

    [Fact]
    public void Compute_K47_DenominatorExceedsIntMaxValue()
    {
        var (_, d) = F89PathPolynomialPipeline.Compute(47);
        Assert.True(d > new BigInteger(int.MaxValue),
            $"D_47 should exceed int.MaxValue ({int.MaxValue}); got {d}.");
        Assert.Equal(new BigInteger(4_632_608_768L), d);
    }

    [Fact]
    public void Compute_K3_Path3AlgebraicAnchor()
    {
        var (coefs, d) = F89PathPolynomialPipeline.Compute(3);
        Assert.Equal(new BigInteger(9), d);
        Assert.Equal(2, coefs.Length);
        Assert.Equal(new BigInteger(47), coefs[0]);
        Assert.Equal(new BigInteger(14), coefs[1]);
    }

    [Fact]
    public void Compute_KLessThan3_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => F89PathPolynomialPipeline.Compute(2));
        Assert.Throws<ArgumentOutOfRangeException>(() => F89PathPolynomialPipeline.Compute(0));
    }

    [Theory]
    [InlineData(3)]
    [InlineData(7)]
    [InlineData(15)]
    [InlineData(32)]
    [InlineData(50)]
    public void ComputePathPolynomialBig_DelegatesToPipeline(int k)
    {
        var (pipelineCoefs, pipelineD) = F89PathPolynomialPipeline.Compute(k);
        var (claimCoefs, claimD) = F89UnifiedFaClosedFormClaim.ComputePathPolynomialBig(k);
        Assert.Equal(pipelineD, claimD);
        Assert.Equal(pipelineCoefs.Length, claimCoefs.Length);
        for (int i = 0; i < pipelineCoefs.Length; i++)
            Assert.Equal(pipelineCoefs[i], claimCoefs[i]);
    }
}
