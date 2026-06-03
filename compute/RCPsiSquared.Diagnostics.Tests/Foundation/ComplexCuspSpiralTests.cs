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

    [Fact]
    public void CrossingTime_PutsTheMagnitudeOnTheQuarterCircle()
    {
        double tc = ComplexCuspSpiral.CrossingTime(0.05);
        Assert.Equal(0.25, ComplexCuspSpiral.Magnitude(0.05, tc), 9);   // by construction |CΨ| = ¼ there
        Assert.True(tc > 0.7 && tc < 0.8, $"≈0.747 at γ=0.05; got {tc}");
    }

    [Fact]
    public void CrossingTime_ScalesAsOneOverGamma()
    {
        // t_cross = K/γ: doubling γ halves it (the crossing K is dimensionless).
        double a = ComplexCuspSpiral.CrossingTime(0.05);
        double b = ComplexCuspSpiral.CrossingTime(0.10);
        Assert.Equal(a / 2.0, b, 9);
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
