using System.Numerics;
using RCPsiSquared.Core.Numerics;

namespace RCPsiSquared.Core.Tests.Numerics;

public class RationalPolynomialTests
{
    [Fact]
    public void Constructor_TrimsTrailingZeros()
    {
        var p = new RationalPolynomial(new BigRational(1), new BigRational(2), BigRational.Zero);
        Assert.Equal(1, p.Degree);
    }

    [Fact]
    public void Addition_HandlesUnequalLengths()
    {
        var a = new RationalPolynomial(new BigRational(1), new BigRational(2));
        var b = new RationalPolynomial(new BigRational(3), new BigRational(4), new BigRational(5));
        var sum = a + b;
        Assert.Equal(2, sum.Degree);
        Assert.Equal(new BigInteger(4), sum[0].Numerator);
        Assert.Equal(new BigInteger(6), sum[1].Numerator);
        Assert.Equal(new BigInteger(5), sum[2].Numerator);
    }

    [Fact]
    public void Multiplication_ComputesProductDegree()
    {
        var a = new RationalPolynomial(new BigRational(1), new BigRational(1));  // 1 + y
        var b = new RationalPolynomial(new BigRational(-1), new BigRational(1)); // -1 + y
        var product = a * b;
        Assert.Equal(2, product.Degree);
        Assert.Equal(new BigInteger(-1), product[0].Numerator);
        Assert.True(product[1].IsZero);
        Assert.Equal(new BigInteger(1), product[2].Numerator);
    }

    [Fact]
    public void DivMod_ExactDivision_ReturnsZeroRemainder()
    {
        var dividend = new RationalPolynomial(new BigRational(-1), BigRational.Zero, new BigRational(1));  // y² - 1
        var divisor = new RationalPolynomial(new BigRational(-1), new BigRational(1));  // y - 1
        var (quotient, remainder) = dividend.DivMod(divisor);
        Assert.True(remainder.IsZero);
        Assert.Equal(1, quotient.Degree);
        Assert.Equal(new BigInteger(1), quotient[0].Numerator);
        Assert.Equal(new BigInteger(1), quotient[1].Numerator);
    }

    [Fact]
    public void Mod_PolynomialReduction_ComputesRemainder()
    {
        // (y³ + 1) mod (y² + 1) = ?
        // y³ + 1 = y · (y² + 1) - y + 1, so remainder = -y + 1.
        var dividend = new RationalPolynomial(new BigRational(1), BigRational.Zero, BigRational.Zero, new BigRational(1));
        var divisor = new RationalPolynomial(new BigRational(1), BigRational.Zero, new BigRational(1));
        var remainder = dividend.Mod(divisor);
        Assert.Equal(1, remainder.Degree);
        Assert.Equal(new BigInteger(1), remainder[0].Numerator);
        Assert.Equal(new BigInteger(-1), remainder[1].Numerator);
    }

    [Fact]
    public void DenominatorLcm_RationalCoefs_ReturnsLcm()
    {
        var p = new RationalPolynomial(
            new BigRational(1, 6),
            new BigRational(1, 4),
            new BigRational(1, 9));
        Assert.Equal(new BigInteger(36), p.DenominatorLcm());
    }
}
