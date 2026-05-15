using System.Numerics;
using RCPsiSquared.Core.Numerics;

namespace RCPsiSquared.Core.Tests.Numerics;

public class BigRationalTests
{
    [Fact]
    public void Constructor_Reduces_ToLowestTerms()
    {
        var r = new BigRational(new BigInteger(6), new BigInteger(8));
        Assert.Equal(new BigInteger(3), r.Numerator);
        Assert.Equal(new BigInteger(4), r.Denominator);
    }

    [Fact]
    public void Constructor_NegativeDenominator_NormalizesSign()
    {
        var r = new BigRational(new BigInteger(3), new BigInteger(-4));
        Assert.Equal(new BigInteger(-3), r.Numerator);
        Assert.Equal(new BigInteger(4), r.Denominator);
    }

    [Fact]
    public void Constructor_ZeroNumerator_NormalizesToOneDenominator()
    {
        var r = new BigRational(BigInteger.Zero, new BigInteger(5));
        Assert.True(r.IsZero);
        Assert.Equal(BigInteger.One, r.Denominator);
    }

    [Fact]
    public void Constructor_ZeroDenominator_Throws()
    {
        Assert.Throws<DivideByZeroException>(() => new BigRational(BigInteger.One, BigInteger.Zero));
    }

    [Fact]
    public void Addition_HandlesDifferentDenominators()
    {
        var a = new BigRational(1, 2);
        var b = new BigRational(1, 3);
        var sum = a + b;
        Assert.Equal(new BigInteger(5), sum.Numerator);
        Assert.Equal(new BigInteger(6), sum.Denominator);
    }

    [Fact]
    public void Multiplication_ReducesAtConstruction()
    {
        var a = new BigRational(2, 3);
        var b = new BigRational(3, 4);
        var product = a * b;
        Assert.Equal(new BigInteger(1), product.Numerator);
        Assert.Equal(new BigInteger(2), product.Denominator);
    }

    [Fact]
    public void Division_HandlesReciprocal()
    {
        var a = new BigRational(2, 3);
        var b = new BigRational(4, 9);
        var quotient = a / b;
        Assert.Equal(new BigInteger(3), quotient.Numerator);
        Assert.Equal(new BigInteger(2), quotient.Denominator);
    }

    [Fact]
    public void Lcm_ComputesCorrectly()
    {
        Assert.Equal(new BigInteger(12), BigRational.Lcm(new BigInteger(4), new BigInteger(6)));
        Assert.Equal(new BigInteger(15), BigRational.Lcm(new BigInteger(3), new BigInteger(5)));
        Assert.Equal(new BigInteger(100), BigRational.Lcm(new BigInteger(20), new BigInteger(25)));
    }
}
