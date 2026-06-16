using System.Numerics;
using RCPsiSquared.Core.Numerics;
using Xunit;

namespace RCPsiSquared.Core.Tests.Numerics;

public class EigenvalueClockExtractionTests
{
    [Fact]
    public void SlowestRate_IsTheGap_OmegaIsTheImAtIt()
    {
        // slowest nonzero rate 0.1 (the ±2i pair); a faster real mode at 0.5; a steady mode at 0.
        var spec = new[]
        {
            new Complex(-0.1, 2.0), new Complex(-0.1, -2.0),
            new Complex(-0.5, 0.0), new Complex(0.0, 0.0),
        };
        var (gap, omega) = EigenvalueClockExtraction.ExtractClockFromSpectrum(spec);
        Assert.Equal(0.1, gap, 12);
        Assert.Equal(2.0, omega, 12);
    }

    [Fact]
    public void OverdampedGap_HasZeroOmega()
    {
        // the slowest mode is real (overdamped) at 0.05; the oscillating pair is faster (0.3) and ignored
        var spec = new[]
        {
            new Complex(-0.05, 0.0),
            new Complex(-0.3, 5.0), new Complex(-0.3, -5.0),
        };
        var (gap, omega) = EigenvalueClockExtraction.ExtractClockFromSpectrum(spec);
        Assert.Equal(0.05, gap, 12);
        Assert.Equal(0.0, omega, 12);
    }

    [Fact]
    public void NoDecayingMode_IsAStoppedClock()
    {
        var spec = new[] { new Complex(0.0, 1.0), new Complex(0.0, -1.0) };   // all steady (rate 0)
        var (gap, omega) = EigenvalueClockExtraction.ExtractClockFromSpectrum(spec);
        Assert.Equal(0.0, gap, 12);
        Assert.Equal(0.0, omega, 12);
    }
}
