using System;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

public class OddHarmonicApproachTests
{
    [Fact]
    public void Anchors_AreCanonical()
    {
        Assert.Equal(0.25, OddHarmonicApproach.Cusp);
        Assert.Equal(0.75, OddHarmonicApproach.CrossingThreshold);
        Assert.Equal(1.0, OddHarmonicApproach.BellPlusS);
    }

    [Fact]
    public void Rates_AreFourAndTwelveGamma()
    {
        Assert.Equal(0.2, OddHarmonicApproach.CarrierRate(0.05), 12);
        Assert.Equal(0.6, OddHarmonicApproach.HarmonicRate(0.05), 12);
    }

    [Fact]
    public void Weights_BellPlus_AreOneSixthEach()
    {
        var (w0, w1) = OddHarmonicApproach.Weights(1.0);
        Assert.Equal(1.0 / 6.0, w0, 12);
        Assert.Equal(1.0 / 6.0, w1, 12);
        var (z0, z1) = OddHarmonicApproach.Weights(0.0);
        Assert.Equal(0.0, z0, 12);
        Assert.Equal(0.0, z1, 12);
    }

    [Fact]
    public void InitialCpsi_IsSOverThree_AndEqualsWeightSum()
    {
        foreach (double s in new[] { 0.3, 0.5, 0.75, 0.9, 1.0 })
        {
            Assert.Equal(s / 3.0, OddHarmonicApproach.InitialCpsi(s), 12);
            var (w0, w1) = OddHarmonicApproach.Weights(s);
            Assert.Equal(s / 3.0, w0 + w1, 12);
            Assert.Equal(OddHarmonicApproach.InitialCpsi(s), OddHarmonicApproach.Cpsi(s, 0.05, 0.0), 12);
        }
    }

    [Fact]
    public void BellPlusMember_ReproducesF25()
    {
        for (double t = 0.0; t <= 10.0; t += 0.5)
            Assert.Equal(InteriorHorizon.BellPlusCpsi(0.05, t), OddHarmonicApproach.Cpsi(1.0, 0.05, t), 12);
    }

    /// <summary>CΨ read off the density matrix rather than off the two-exponential: for
    /// |ψ(α)⟩ = cosα|00⟩ + sinα|11⟩ under Z-dephasing the populations are frozen and the single
    /// coherence is ρ₀₃ = cosα·sinα·e^{−4γt}, so the whole 4×4 state is known. This computes
    /// CΨ = Tr(ρ²)·ℓ₁(ρ)/(d−1) with d = 4 from that matrix, entry by entry, and is a second route
    /// to the number: it shares the dephasing solution and nothing of the weight algebra, so a
    /// wrong w₀, a wrong w₁, a swapped rate or a wrong normalisation all break it.</summary>
    private static double CpsiFromDensityMatrix(double s, double gamma, double t)
    {
        // s = sin2α fixes cos²α, sin²α = (1 ± sqrt(1−s²))/2 up to the branch; both branches give the
        // same CΨ, since only s enters below. Take the upper one.
        double root = Math.Sqrt(Math.Max(0.0, 1.0 - s * s));
        double p00 = 0.5 * (1.0 + root);
        double p11 = 0.5 * (1.0 - root);
        double coherence = 0.5 * s * Math.Exp(-4.0 * gamma * t);   // cosα·sinα·f = (s/2)·f

        double purity = p00 * p00 + p11 * p11 + 2.0 * coherence * coherence;
        double l1OffDiagonal = 2.0 * coherence;                    // |ρ₀₃| + |ρ₃₀|
        return purity * l1OffDiagonal / 3.0;                       // d − 1 = 3
    }

    [Theory]
    [InlineData(1.0, 0.05, 0.0)]
    [InlineData(1.0, 0.05, 3.7)]
    [InlineData(0.8, 0.3, 1.1)]
    [InlineData(0.5, 1.0, 0.4)]
    [InlineData(0.13, 0.05, 9.0)]
    public void Cpsi_MatchesTheDensityMatrixItClaimsToDescribe(double s, double gamma, double t)
    {
        // Both routes are float; the residual measured over these five points is at most 3 eps
        // against a value of order 0.1, so this is gated at 8 eps and not at a convenient decimal.
        double closedForm = OddHarmonicApproach.Cpsi(s, gamma, t);
        double fromMatrix = CpsiFromDensityMatrix(s, gamma, t);
        Assert.True(Math.Abs(closedForm - fromMatrix) < 8.0 * 2.220446049250313e-16,
            $"closed form {closedForm:R} vs density matrix {fromMatrix:R}");
    }

    [Fact]
    public void Cpsi_WouldNotMatchIfTheHarmonicRateWereAnythingElse()
    {
        // The control: the second exponential is f³, i.e. 12γ. Swapping it for 8γ (an even
        // multiple, the plausible wrong answer) moves the number where the gate above can see it.
        const double s = 0.8, gamma = 0.3, t = 1.1;
        var (w0, w1) = OddHarmonicApproach.Weights(s);
        double wrong = w0 * Math.Exp(-4.0 * gamma * t) + w1 * Math.Exp(-8.0 * gamma * t);
        Assert.True(Math.Abs(wrong - CpsiFromDensityMatrix(s, gamma, t)) > 1e-3);
    }

    [Fact]
    public void CrossingTime_AtZeroGamma_IsNaNBecauseTheCurveIsConstant()
    {
        Assert.True(double.IsNaN(OddHarmonicApproach.CrossingTime(1.0, 0.0)));
    }

    [Fact]
    public void CrossingTime_RefusesNegativeGamma()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => OddHarmonicApproach.CrossingTime(1.0, -0.05));
    }

    [Theory]
    [InlineData(0.80, 0.05, true)]
    [InlineData(1.00, 0.05, true)]
    [InlineData(0.75, 0.05, false)]
    [InlineData(0.50, 0.05, false)]
    [InlineData(1.00, 0.00, false)]
    [InlineData(1.00, -0.05, false)]
    public void HasDownwardCrossing_IffGammaPositiveAndSAboveThreeQuarters(
        double s,
        double gamma,
        bool expected)
    {
        Assert.Equal(expected, OddHarmonicApproach.HasDownwardCrossing(s, gamma));
    }

    [Theory]
    [InlineData(0.9, 0.05)]
    [InlineData(0.6, 0.3)]
    [InlineData(0.3, 1.0)]
    public void LogRate_FallsFromFourGammaOnePlusSSquared_ToFourGamma(double s, double gamma)
    {
        // The slowing is the observable's: -d ln CΨ/dt = 4γ(w₀ + 3w₁)/(w₀ + w₁) = 4γ(1 + s²) at t = 0,
        // and 4γ late, while the coherence factor decays at 4γ throughout. At t = 0 it is read by a
        // central difference of ln Cpsi, whose truncation error is h² times the third derivative of
        // ln CΨ over 6: the error law is a ratio of 100 per decade of h. (At s = 1 that third
        // derivative vanishes, the two rates carrying equal weight, so s = 1 is left out.)
        double LogRate(double t, double h) =>
            -(Math.Log(OddHarmonicApproach.Cpsi(s, gamma, t + h)) - Math.Log(OddHarmonicApproach.Cpsi(s, gamma, t - h))) / (2.0 * h);
        double scale = 1.0 / gamma;
        double expected = 4.0 * gamma * (1.0 + s * s);
        double e1 = Math.Abs(LogRate(0.0, 1e-2 * scale) - expected);
        double e2 = Math.Abs(LogRate(0.0, 1e-3 * scale) - expected);
        Assert.True(e1 / e2 > 90.0 && e1 / e2 < 110.0, $"s={s}, γ={gamma}: error ratio {e1 / e2} is not the h² law");

        // Late, at K = γt = 3: the exact rate is 4γ plus the 12γ term's share 8γ·w₁e^(−24)/(w₀ + w₁e^(−24)).
        // There the difference's truncation is negligible (the third derivative of ln CΨ carries e^(−24)) and
        // its rounding is about ε·(1 + |ln CΨ|)/h, so that is the gate, with margin 8.
        var (w0, w1) = OddHarmonicApproach.Weights(s);
        double tLate = 3.0 * scale, hLate = 1e-3 * scale, e24 = Math.Exp(-24.0);
        double exactLate = 4.0 * gamma + 8.0 * gamma * w1 * e24 / (w0 + w1 * e24);
        double late = LogRate(tLate, hLate);
        double roundingLaw = 2.220446049250313e-16 * (1.0 + Math.Abs(Math.Log(OddHarmonicApproach.Cpsi(s, gamma, tLate)))) / hLate;
        Assert.True(Math.Abs(late - exactLate) <= 8.0 * roundingLaw,
            $"s={s}: late rate {late:R} vs exact {exactLate:R}, rounding law {roundingLaw:E2}");
    }

    [Fact]
    public void HarmonicFraction_IsSSquaredOverTwo()
    {
        Assert.Equal(0.5, OddHarmonicApproach.HarmonicFraction(1.0), 12);
        foreach (double s in new[] { 0.3, 0.6, 0.9 })
        {
            Assert.Equal(0.5 * s * s, OddHarmonicApproach.HarmonicFraction(s), 12);
            var (w0, w1) = OddHarmonicApproach.Weights(s);
            Assert.Equal(w1 / (w0 + w1), OddHarmonicApproach.HarmonicFraction(s), 12);
        }
    }

    [Fact]
    public void CrossingTime_PutsCpsiOnTheCusp_AndScalesAsOneOverGamma()
    {
        double tc = OddHarmonicApproach.CrossingTime(1.0, 0.05);
        Assert.Equal(0.25, OddHarmonicApproach.Cpsi(1.0, 0.05, tc), 7);
        double tc2 = OddHarmonicApproach.CrossingTime(1.0, 0.10);
        Assert.Equal(tc / 2.0, tc2, 6);
    }

    [Fact]
    public void CrossingTime_IsNaN_BelowThreshold()
    {
        Assert.True(double.IsNaN(OddHarmonicApproach.CrossingTime(0.75, 0.05)));
        Assert.True(double.IsNaN(OddHarmonicApproach.CrossingTime(0.5, 0.05)));
    }
}
