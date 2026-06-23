using RCPsiSquared.Core.Numerics;
using Xunit;

namespace RCPsiSquared.Core.Tests.Numerics;

/// <summary>Exact Z[i] polynomial division by a monic divisor — the operation that isolates F_d
/// (charpoly ÷ AT-factor) and underlies the validation triple. Coefficients lowest-first.</summary>
public class GaussianPolynomialTests
{
    [Fact]
    public void DivMod_ByMonicLinear_ExactDivision()
    {
        // (x² − 1) / (x − 1) = x + 1, remainder 0.
        var (q, r) = GaussianPolynomial.DivMod(
            new GaussianInteger[] { -1, 0, 1 }, new GaussianInteger[] { -1, 1 });
        Assert.Equal(new GaussianInteger[] { 1, 1 }, q);
        Assert.Empty(r);
    }

    [Fact]
    public void DivMod_WithNonzeroRemainder()
    {
        // (x² + 1) / (x − 1) = x + 1, remainder 2.
        var (q, r) = GaussianPolynomial.DivMod(
            new GaussianInteger[] { 1, 0, 1 }, new GaussianInteger[] { -1, 1 });
        Assert.Equal(new GaussianInteger[] { 1, 1 }, q);
        Assert.Equal(new GaussianInteger[] { 2 }, r);
    }

    [Fact]
    public void DivMod_GaussianDivisor_FactorsXSquaredPlusOne()
    {
        // x² + 1 = (x − i)(x + i): dividing by (x − i) gives x + i, remainder 0.
        var (q, r) = GaussianPolynomial.DivMod(
            new GaussianInteger[] { 1, 0, 1 }, new GaussianInteger[] { -GaussianInteger.I, 1 });
        Assert.Equal(new GaussianInteger[] { GaussianInteger.I, 1 }, q);
        Assert.Empty(r);
    }

    // ---- Multiply (builds the AT factor F_a·F_b) ----

    [Fact]
    public void Multiply_RealFactors()
    {
        // (x + 1)(x − 1) = x² − 1.
        Assert.Equal(new GaussianInteger[] { -1, 0, 1 },
            GaussianPolynomial.Multiply(new GaussianInteger[] { 1, 1 }, new GaussianInteger[] { -1, 1 }));
    }

    [Fact]
    public void Multiply_GaussianConjugateFactors()
    {
        // (x + i)(x − i) = x² + 1.
        Assert.Equal(new GaussianInteger[] { 1, 0, 1 },
            GaussianPolynomial.Multiply(new GaussianInteger[] { GaussianInteger.I, 1 },
                                        new GaussianInteger[] { -GaussianInteger.I, 1 }));
    }

    // ---- AreCoprime (the resultant ≠ 0 leg of the validation triple) ----

    [Fact]
    public void AreCoprime_DistinctLinearFactors_True()
    {
        // (x − 1) and (x − 2) share no root.
        Assert.True(GaussianPolynomial.AreCoprime(
            new GaussianInteger[] { -1, 1 }, new GaussianInteger[] { -2, 1 }));
    }

    [Fact]
    public void AreCoprime_SharedFactor_False()
    {
        // (x−1)(x−2) = x²−3x+2 and (x−1)(x−3) = x²−4x+3 share the factor (x − 1).
        Assert.False(GaussianPolynomial.AreCoprime(
            new GaussianInteger[] { 2, -3, 1 }, new GaussianInteger[] { 3, -4, 1 }));
    }

    [Fact]
    public void AreCoprime_GaussianConjugates_True()
    {
        // (x − i) and (x + i) are coprime over Q(i).
        Assert.True(GaussianPolynomial.AreCoprime(
            new GaussianInteger[] { -GaussianInteger.I, 1 }, new GaussianInteger[] { GaussianInteger.I, 1 }));
    }
}
