using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

public class ComplexCuspSpiralTests
{
    [Fact]
    public void CircleRadius_IsAQuarter() => Assert.Equal(0.25, ComplexCuspSpiral.CircleRadius);

    [Fact]
    public void Magnitude_AtTimeZero_IsOneThird()
    {
        // |CΨ_com|(0) = f(1+f²)/6 with f = 1 → 2/6 = 1/3, the Bell+ start.
        Assert.Equal(1.0 / 3.0, ComplexCuspSpiral.Magnitude(0.05, 0.0), 12);
    }

    [Fact]
    public void Magnitude_IsTheSameF25LawAsTheLine()
    {
        // The radial law is exactly InteriorHorizon.BellPlusCpsi: the drift does not touch it.
        for (double t = 0.0; t <= 5.0; t += 0.5)
            Assert.Equal(InteriorHorizon.BellPlusCpsi(0.05, t), ComplexCuspSpiral.Magnitude(0.05, t), 12);
    }

    [Fact]
    public void Argument_IsLinearInTime_AtRateMinusOmega()
    {
        Assert.Equal(0.3, ComplexCuspSpiral.Argument(0.4, 0.3, 0.0), 12);   // φ₀ at t=0
        Assert.Equal(0.3 - 0.4 * 2.0, ComplexCuspSpiral.Argument(0.4, 0.3, 2.0), 12);
    }

    [Fact]
    public void ReIm_AtOmegaZero_StayOnTheRealAxis()
    {
        // φ₀ = 0, Ω = 0: the spiral is the real-axis line, Im ≡ 0, Re = magnitude (the 1D interior axis).
        for (double t = 0.0; t <= 5.0; t += 1.0)
        {
            Assert.Equal(ComplexCuspSpiral.Magnitude(0.05, t), ComplexCuspSpiral.Re(0.05, 0.0, 0.0, t), 12);
            Assert.Equal(0.0, ComplexCuspSpiral.Im(0.05, 0.0, 0.0, t), 12);
        }
    }

    [Fact]
    public void WindingRate_IsOmegaOverFourGamma()
    {
        Assert.Equal(2.0, ComplexCuspSpiral.WindingRate(0.05, 0.4), 12);   // the script's Ω/4γ = 2
        Assert.Equal(7.5, ComplexCuspSpiral.WindingRate(0.05, 1.5), 12);
        Assert.Equal(0.0, ComplexCuspSpiral.WindingRate(0.05, 0.0), 12);
    }

    [Theory]
    [InlineData(0.05)]
    [InlineData(0.5)]
    [InlineData(1.3)]
    public void MagnitudeEFoldRate_StartsAtEightGamma_AndOnlyApproachesFourGamma(double gamma)
    {
        // The denominator Ω/(4γ) is NOT divided by. At t = 0, f = 1 and the magnitude sheds
        // e-folds twice as fast as f does; 4γ is the f → 0 limit, approached and never attained.
        Assert.Equal(8.0 * gamma, ComplexCuspSpiral.MagnitudeEFoldRate(gamma, 0.0));

        double previous = double.PositiveInfinity;
        foreach (double k in new[] { 0.0, 0.25, 1.0, 2.0, 4.0 })   // t = k/γ, so f = e^{−4k}
        {
            double rate = ComplexCuspSpiral.MagnitudeEFoldRate(gamma, k / gamma);
            Assert.True(rate < previous, $"rate should fall monotonically; {rate:R} at k={k}");
            Assert.True(rate > 4.0 * gamma, $"rate should stay above 4γ; {rate:R} at k={k}");
            previous = rate;
        }

        // The limit is a limit of the formula, not of the arithmetic: past f² ≈ 1e-16 the ratio
        // (1+3f²)/(1+f²) rounds to 1 and the reading sits exactly on 4γ. That is where a test
        // asserting strict inequality would fail for a reason that is not about the physics.
        Assert.Equal(4.0 * gamma, ComplexCuspSpiral.MagnitudeEFoldRate(gamma, 100.0 / gamma));
    }

    [Theory]
    [InlineData(0.05)]
    [InlineData(0.5)]
    [InlineData(1.3)]
    public void MagnitudeEFoldRate_AtTheCrossing_IsF25sOwnCrossingDerivative(double gamma)
    {
        // Not a resemblance: the rate is |dCΨ/dt|/|CΨ| and |CΨ| = ¼ on the circle, so the factor
        // over 4γ is exactly |dCΨ/dt|/γ = 2f*(1+3f*²)/3 whenever f*(1+f*²) = 3/2. Both routes are
        // float, so this is rounding only, gated at 8 eps against a value of order unity.
        double tc = ComplexCuspSpiral.CrossingTime(gamma);
        double factor = ComplexCuspSpiral.MagnitudeEFoldRate(gamma, tc) / (4.0 * gamma);
        const double fStar = 0.8612240997395737;
        double f25 = 2.0 * fStar * (1.0 + 3.0 * fStar * fStar) / 3.0;
        Assert.True(System.Math.Abs(factor - f25) < 8.0 * 2.220446049250313e-16,
            $"factor {factor:R} vs F25 crossing derivative {f25:R}");
    }

    [Theory]
    [InlineData(0.05, 0.4)]
    [InlineData(0.5, 0.4)]
    [InlineData(0.05, 1.5)]
    [InlineData(1.3, -0.7)]
    public void PhasePerMagnitudeEFold_IsWindingRateOverAConstant(double gamma, double omega)
    {
        // The two readings differ by (−ln f*)/ln(4/3) = 1.9255760, at every γ and Ω. That the
        // ratio is a constant is the content; that it is not 1 is what the old label got wrong.
        double ratio = ComplexCuspSpiral.WindingRate(gamma, omega)
                     / ComplexCuspSpiral.PhasePerMagnitudeEFold(gamma, omega);
        Assert.Equal(1.9255759835293094, ratio, 12);
        Assert.True(System.Math.Abs(ratio - 1.0) > 0.9);
    }

    [Theory]
    [InlineData(0.05)]
    [InlineData(0.5)]
    [InlineData(0.137)]
    [InlineData(1.0 / 3.0)]
    public void CrossingTime_PutsTheMagnitudeExactlyOnTheQuarterCircle(double gamma)
    {
        // The Newton solve returns the double whose round trip lands on 0.25 with no residual at
        // all, so this is compared exactly rather than gated: a nonzero residual here would be a
        // finding about the solve, not a tolerance to widen. Four γ including a non-dyadic one.
        double tc = ComplexCuspSpiral.CrossingTime(gamma);
        Assert.Equal(0.25, ComplexCuspSpiral.Magnitude(gamma, tc));
    }

    [Fact]
    public void CrossingTime_IsTheF25CrossingKOverGamma()
    {
        // t_cross = K/γ with K = −ln(f*)/4 = 0.037350132494447214, F25's γ-invariant. ≈0.747 at
        // γ = 0.05.
        Assert.Equal(0.037350132494447214 / 0.05, ComplexCuspSpiral.CrossingTime(0.05), 12);
    }

    [Theory]
    [InlineData(0.05, 0.10)]
    [InlineData(0.137, 0.274)]
    [InlineData(1.0 / 3.0, 2.0 / 3.0)]
    [InlineData(0.7, 0.21)]
    public void CrossingTime_ScalesAsOneOverGamma(double gammaA, double gammaB)
    {
        // t_cross = K/γ: the product γ·t_cross is the same K at every rate. Non-dyadic ratios are
        // included, since a doubling alone is exact in binary and could not break the claim.
        double ka = gammaA * ComplexCuspSpiral.CrossingTime(gammaA);
        double kb = gammaB * ComplexCuspSpiral.CrossingTime(gammaB);
        Assert.Equal(ka, kb, 14);
    }

    [Fact]
    public void CrossingTime_AtRadiusOneThird_IsZero_AndAboveIsNaN()
    {
        Assert.Equal(0.0, ComplexCuspSpiral.CrossingTime(0.05, 1.0 / 3.0), 9);
        Assert.True(double.IsNaN(ComplexCuspSpiral.CrossingTime(0.05, 0.4)));  // > 1/3: never reached
    }

    [Fact]
    public void CrossingArgument_IsZeroAtOmegaZero_TheHeadOnCrossing()
    {
        // Ω = 0 crosses ¼ on the real axis (arg 0): exactly the 1D interior crossing.
        Assert.Equal(0.0, ComplexCuspSpiral.CrossingArgument(0.05, 0.0, 0.0), 12);
    }

    [Fact]
    public void CrossingArgument_IsPhiZeroMinusOmegaTCross()
    {
        double tc = ComplexCuspSpiral.CrossingTime(0.05);
        Assert.Equal(0.0 - 0.4 * tc, ComplexCuspSpiral.CrossingArgument(0.05, 0.4, 0.0), 12);
    }

    [Fact]
    public void WindingNumber_IsOmegaTMaxOverTwoPi()
    {
        // |Ω|·tMax/(2π) full turns; Ω=0.5 over tMax=4π is exactly one turn; Ω=0 is none.
        Assert.Equal(1.0, ComplexCuspSpiral.WindingNumber(0.5, 4.0 * System.Math.PI), 12);
        Assert.Equal(0.0, ComplexCuspSpiral.WindingNumber(0.0, 10.0), 12);
    }

    [Fact]
    public void ReIm_AtNonzeroOmega_ArePolarComponentsOfTheSpiral()
    {
        // Off the real axis (Ω ≠ 0) Re/Im are |CΨ|·cos/sin(arg), and Re²+Im² = |CΨ|²: the spiral
        // sits on the magnitude its radial law dictates, only rotated.
        const double g = 0.05, w = 0.4, p = 0.3;
        for (double t = 0.5; t <= 5.0; t += 1.5)
        {
            double mag = ComplexCuspSpiral.Magnitude(g, t);
            double arg = ComplexCuspSpiral.Argument(w, p, t);
            double re = ComplexCuspSpiral.Re(g, w, p, t);
            double im = ComplexCuspSpiral.Im(g, w, p, t);
            Assert.Equal(mag * System.Math.Cos(arg), re, 12);
            Assert.Equal(mag * System.Math.Sin(arg), im, 12);
            Assert.Equal(mag, System.Math.Sqrt(re * re + im * im), 12);
        }
    }
}
