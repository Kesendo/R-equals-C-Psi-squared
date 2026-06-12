using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

public class QuarterEnvelopeTests
{
    [Fact]
    public void Maxima_ParabolicApex_ReportsAboveSample()
    {
        // One interior peak at i=1; its true apex (parabola through 0.8,0.9,0.7) sits ABOVE 0.9.
        double[] cpsi = { 0.8, 0.9, 0.7 };
        double[] t = { 0.0, 1.0, 2.0 };
        var e = QuarterEnvelope.Of(cpsi, t);
        Assert.Single(e.Maxima);
        Assert.Equal(1.0, e.Maxima[0].Time, 12);
        Assert.Equal(0.904167, e.Maxima[0].ApexValue, 5);   // parabolic apex > the 0.9 sample
        Assert.True(e.IsNonIncreasing);                      // a single maximum is trivially non-increasing
    }

    [Fact]
    public void MonotoneDecreasing_LeftEndpointOnly_NonIncreasing()
    {
        double[] cpsi = { 0.9, 0.7, 0.5 };
        double[] t = { 0.0, 1.0, 2.0 };
        var e = QuarterEnvelope.Of(cpsi, t);
        Assert.Single(e.Maxima);                  // the descending-start t=0 endpoint
        Assert.Equal(0.9, e.Maxima[0].ApexValue, 12);   // endpoint keeps its raw value
        Assert.True(e.IsNonIncreasing);
        Assert.Equal(0, e.RiseCount);
    }

    [Fact]
    public void PredecessorRise_Detected_NotRunningMax()
    {
        // Three peaks 0.30, 0.28, 0.29 (symmetric neighbours so apex == sample). The third beats its
        // PREDECESSOR (0.28) but NOT the running max (0.30): predecessor semantics must flag it.
        double[] cpsi = { 0.30, 0.1, 0.28, 0.1, 0.29, 0.1 };
        double[] t = { 0, 1, 2, 3, 4, 5 };
        var e = QuarterEnvelope.Of(cpsi, t);
        Assert.Equal(1, e.RiseCount);
        Assert.False(e.IsNonIncreasing);
        Assert.Equal(4.0, e.FirstRiseTime!.Value, 12);
        Assert.Equal(0.01, e.MaxRiseMagnitude, 9);
    }

    [Fact]
    public void RiseTol_Gate_SuppressesSmallRise()
    {
        double[] cpsi = { 0.30, 0.1, 0.28, 0.1, 0.29, 0.1 };
        double[] t = { 0, 1, 2, 3, 4, 5 };
        var e = QuarterEnvelope.Of(cpsi, t, threshold: 0.25, riseTol: 0.02);
        Assert.Equal(0, e.RiseCount);   // 0.01 rise is below the 0.02 tolerance
        Assert.True(e.IsNonIncreasing);
    }

    [Fact]
    public void EnvelopeFold_IsLastAbsorbingDownCrossing()
    {
        // Crosses ¼ down (s=1), back up (0.3 at s=2), then down-for-good (s=3): the fold is the s=3 crossing.
        double[] cpsi = { 0.4, 0.2, 0.3, 0.1, 0.05 };
        double[] t = { 0, 1, 2, 3, 4 };
        var e = QuarterEnvelope.Of(cpsi, t);
        Assert.NotNull(e.EnvelopeFoldTime);
        Assert.Equal(2.25, e.EnvelopeFoldTime!.Value, 9);   // interp: 2 + (0.3-0.25)/(0.3-0.1)
    }

    [Fact]
    public void EnvelopeFold_NullWhenNeverSettlesBelow()
    {
        double[] cpsi = { 0.4, 0.2, 0.3 };   // ends ABOVE ¼
        double[] t = { 0, 1, 2 };
        var e = QuarterEnvelope.Of(cpsi, t);
        Assert.Null(e.EnvelopeFoldTime);
    }
}
