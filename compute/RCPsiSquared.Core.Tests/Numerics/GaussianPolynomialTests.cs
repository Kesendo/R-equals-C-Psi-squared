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

    // ---- Derivative / Evaluate ----

    [Fact]
    public void Derivative_Quadratic_And_Cubic()
    {
        // d/dx(x² − 1) = 2x; d/dx(x³) = 3x².
        Assert.Equal(new GaussianInteger[] { 0, 2 }, GaussianPolynomial.Derivative(new GaussianInteger[] { -1, 0, 1 }));
        Assert.Equal(new GaussianInteger[] { 0, 0, 3 }, GaussianPolynomial.Derivative(new GaussianInteger[] { 0, 0, 0, 1 }));
        Assert.Empty(GaussianPolynomial.Derivative(new GaussianInteger[] { 5 }));   // constant → 0
    }

    [Fact]
    public void Evaluate_Horner()
    {
        // (x² + 1) at x = 2 is 5; at x = i is 0.
        Assert.Equal(new GaussianInteger(5, 0), GaussianPolynomial.Evaluate(new GaussianInteger[] { 1, 0, 1 }, 2));
        Assert.Equal(GaussianInteger.Zero, GaussianPolynomial.Evaluate(new GaussianInteger[] { 1, 0, 1 }, GaussianInteger.I));
    }

    // ---- Resultant (the value behind AreCoprime) ----

    [Fact]
    public void Resultant_DistinctLinears_And_SharedFactor()
    {
        // Res(x − 1, x − 2) = det[[1,−1],[1,−2]] = −1 (nonzero ⟹ coprime).
        Assert.Equal(new GaussianInteger(-1, 0),
            GaussianPolynomial.Resultant(new GaussianInteger[] { -1, 1 }, new GaussianInteger[] { -2, 1 }));
        // a shared root forces the resultant to vanish.
        Assert.Equal(GaussianInteger.Zero,
            GaussianPolynomial.Resultant(new GaussianInteger[] { 2, -3, 1 }, new GaussianInteger[] { 3, -4, 1 }));
    }

    // ---- Discriminant (monic; its zeros in q are the branch loci) ----

    [Fact]
    public void Discriminant_Quadratic_MatchesBSquaredMinus4C()
    {
        // disc(x² + bx + c) = b² − 4c: disc(x² − 1) = 4, disc(x² + 1) = −4, disc((x−1)²) = 0.
        Assert.Equal(new GaussianInteger(4, 0), GaussianPolynomial.Discriminant(new GaussianInteger[] { -1, 0, 1 }));
        Assert.Equal(new GaussianInteger(-4, 0), GaussianPolynomial.Discriminant(new GaussianInteger[] { 1, 0, 1 }));
        Assert.Equal(GaussianInteger.Zero, GaussianPolynomial.Discriminant(new GaussianInteger[] { 1, -2, 1 }));
    }

    // ---- ComposeLinear: p(α·x + β), the fold F_corner(−λ − 2N) ----

    [Fact]
    public void ComposeLinear_ReflectShiftFold()
    {
        // (x² + 1)(−x) = x² + 1 (even); x²∘(x+1) = (x+1)² = x²+2x+1; x∘(−x−4) = −x−4.
        Assert.Equal(new GaussianInteger[] { 1, 0, 1 },
            GaussianPolynomial.ComposeLinear(new GaussianInteger[] { 1, 0, 1 }, -1, 0));
        Assert.Equal(new GaussianInteger[] { 1, 2, 1 },
            GaussianPolynomial.ComposeLinear(new GaussianInteger[] { 0, 0, 1 }, 1, 1));
        Assert.Equal(new GaussianInteger[] { -4, -1 },
            GaussianPolynomial.ComposeLinear(new GaussianInteger[] { 0, 1 }, -1, -4));
        // (x² + 1)(−x − 4) = x² + 8x + 17 (the −λ−2N fold with 2N = 4).
        Assert.Equal(new GaussianInteger[] { 17, 8, 1 },
            GaussianPolynomial.ComposeLinear(new GaussianInteger[] { 1, 0, 1 }, -1, -4));
    }

    [Fact]
    public void ComposeLinear_AgreesWithEvaluate()
    {
        // p(αx+β) evaluated at x must equal p evaluated at αx+β, for a random-ish p and fold.
        var p = new GaussianInteger[] { 3, -1, 0, 2 };        // 2x³ − x + 3
        GaussianInteger alpha = -1, beta = -4, x = 5;
        var composed = GaussianPolynomial.ComposeLinear(p, alpha, beta);
        Assert.Equal(GaussianPolynomial.Evaluate(p, alpha * x + beta),
                     GaussianPolynomial.Evaluate(composed, x));
    }
}
