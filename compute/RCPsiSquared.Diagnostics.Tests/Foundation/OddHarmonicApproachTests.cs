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

    [Fact]
    public void Crosses_IffSAboveThreeQuarters()
    {
        Assert.True(OddHarmonicApproach.Crosses(0.8));
        Assert.True(OddHarmonicApproach.Crosses(1.0));
        Assert.False(OddHarmonicApproach.Crosses(0.75));
        Assert.False(OddHarmonicApproach.Crosses(0.5));
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
