using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class F89Path2CardanoClaimTests
{
    private static F89Path2CardanoClaim BuildClaim() =>
        new F89Path2CardanoClaim(new F89TopologyOrbitClosure(new Pi2DyadicLadderClaim()));

    [Fact]
    public void Tier_IsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived, BuildClaim().Tier);
    }

    [Theory]
    [InlineData(0.05, -0.10, -0.30)]    // γ=0.05: -2γ=-0.10, -6γ=-0.30
    [InlineData(0.025, -0.05, -0.15)]
    [InlineData(1.0, -2.0, -6.0)]
    public void LinearFactorRates_AreNegativeTwoAndSixGamma(double gamma, double expectedLow, double expectedHigh)
    {
        var (low, high) = F89Path2CardanoClaim.LinearFactorRates(gamma);
        Assert.Equal(expectedLow, low, precision: 14);
        Assert.Equal(expectedHigh, high, precision: 14);
    }

    [Theory]
    [InlineData(1.0, 10.0, 60.0, 120.0)]   // q=1: a=10, b=28+32=60, c=24·5=120
    [InlineData(1.5, 10.0, 100.0, 240.0)]  // q=1.5: a=10, b=28+72=100, c=24·10=240
    [InlineData(0.0, 10.0, 28.0, 24.0)]    // q=0: a=10, b=28, c=24
    public void CubicCoefficients_MatchSymbolicForm(double q, double expectedA, double expectedB, double expectedC)
    {
        var (a, b, c) = F89Path2CardanoClaim.CubicCoefficients(q);
        Assert.Equal(expectedA, a, precision: 14);
        Assert.Equal(expectedB, b, precision: 14);
        Assert.Equal(expectedC, c, precision: 14);
    }

    [Fact]
    public void CubicEigenvaluesNumerical_AtStandardQ_MatchSympyOutput()
    {
        // At q=1.5: roots are -3.0448 (real), -3.4776 ± 8.169i (complex pair)
        var roots = F89Path2CardanoClaim.CubicEigenvaluesNumerical(1.5);
        Assert.Equal(3, roots.Length);

        // Find the real root (smallest |Im|)
        var real_root = roots.OrderBy(r => Math.Abs(r.Imaginary)).First();
        Assert.True(Math.Abs(real_root.Imaginary) < 1e-9, "Real root should have ~zero imaginary part");
        Assert.Equal(-3.0448, real_root.Real, precision: 3);

        // Find the complex pair
        var complex_roots = roots.OrderByDescending(r => Math.Abs(r.Imaginary)).Take(2).ToList();
        Assert.Equal(2, complex_roots.Count);
        Assert.Equal(-3.4776, complex_roots[0].Real, precision: 3);
        Assert.Equal(-3.4776, complex_roots[1].Real, precision: 3);
        Assert.Equal(8.169, Math.Abs(complex_roots[0].Imaginary), precision: 2);
        Assert.Equal(8.169, Math.Abs(complex_roots[1].Imaginary), precision: 2);
        // Conjugate pair: imaginary parts opposite
        Assert.Equal(0.0, complex_roots[0].Imaginary + complex_roots[1].Imaginary, precision: 9);
    }

    [Fact]
    public void VerifyAtStandardQ_ReturnsTrue()
    {
        Assert.True(F89Path2CardanoClaim.VerifyAtStandardQ());
    }

    [Theory]
    [InlineData(0.5)]
    [InlineData(1.0)]
    [InlineData(1.5)]
    [InlineData(2.0)]
    [InlineData(3.0)]
    public void CubicEigenvalues_AreActualRoots(double q)
    {
        // Verify each root satisfies μ³ + a·μ² + b·μ + c = 0
        var (a, b, c) = F89Path2CardanoClaim.CubicCoefficients(q);
        var roots = F89Path2CardanoClaim.CubicEigenvaluesNumerical(q);
        foreach (var root in roots)
        {
            var residual = root * root * root + a * root * root + b * root + c;
            Assert.True(residual.Magnitude < 1e-9,
                $"Root {root} at q={q}: residual {residual} should be ~0");
        }
    }

    [Fact]
    public void CubicEigenvalues_NegativeQ_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(
            () => F89Path2CardanoClaim.CubicCoefficients(-1.0));
    }

    [Fact]
    public void LinearFactorRates_NegativeGamma_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(
            () => F89Path2CardanoClaim.LinearFactorRates(-0.05));
    }
}
