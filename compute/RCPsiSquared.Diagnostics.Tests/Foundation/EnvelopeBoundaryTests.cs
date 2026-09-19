using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>The finite envelope_n4_rise atlas, gate-first. These tests preserve the named N=3/4/5,
/// Q, K-window and grid rows and reuse the witness's detector. They test finite observations and one
/// same-(N,Q,K) rescaling comparison; they do not infer an all-Q/all-N threshold, an absence theorem,
/// or a mechanism. See experiments/ENVELOPE_RISE_BOUNDARY.md.</summary>
public class EnvelopeBoundaryTests
{
    private const double Bar = EnvelopeTheoremWitness.RiseReportingBar;

    [Fact]
    public void SameNQKRescaling_AgreesToSixDecimals_OnTheNamedPair()
    {
        // The clock movement certifies the dimensionless lenses are pure (Q,K)-observables: under L→rL
        // (J,γ both scaled by r, the dose window K=γt held fixed) the global CΨ(K) curve is invariant. So
        // the sampled curves should agree at matched (N,Q,K). This pair at Q=500 is related by r=2:
        // (5,0.01,tMax=25) vs (10,0.02,tMax=12.5), both K_max=0.25. Agreement to six decimals is
        // the finite numerical check; it is not an every-Q identity inferred from samples.
        var a = EnvelopeTheoremWitness.GlobalReading(4, 5.0, 0.01, 25.0, 1600);
        var b = EnvelopeTheoremWitness.GlobalReading(4, 10.0, 0.02, 12.5, 1600);
        Assert.Equal(a.RiseCount, b.RiseCount);
        Assert.Equal(a.MaxRiseMagnitude, b.MaxRiseMagnitude, 6);
    }

    [Fact]
    public void NamedRows_N3AtQ2000ResolvesNoRise_N4AtQ500ResolvesARise()
    {
        // One N=3 null sample at Q=2000 and one N=4 positive sample at Q=500. They invite wider
        // all-Q/all-N classification; they do not prove an N floor or identify the internal mechanism.
        Assert.Equal(0, EnvelopeTheoremWitness.GlobalRiseCount(3, 20.0, 0.01, 25.0, 1600));   // Q=2000
        Assert.True(EnvelopeTheoremWitness.GlobalReading(4, 5.0, 0.01, 25.0, 1600).MaxRiseMagnitude
                    > Bar);                                                                    // Q=500: N=4 rises
    }

    [Fact]
    public void NamedN4Rows_Q13AndQ40_GiveDifferentFiniteClassifications()
    {
        // On this fixed K-window/grid the Q=40 row exceeds the reporting bar and Q=13 resolves no rise.
        Assert.True(EnvelopeTheoremWitness.GlobalReading(4, 0.40, 0.01, 25.0, 1600).MaxRiseMagnitude
                    > Bar);                                                                    // Q=40: rises
        Assert.Equal(0, EnvelopeTheoremWitness.GlobalRiseCount(4, 0.13, 0.01, 25.0, 1600));    // Q=13: holds
    }

    [Fact]
    public void NamedN4AndN5Rows_DifferAtQ40_AndN5RisesAtQ500()
    {
        // Three retained atlas rows: N=4/Q=40 exceeds the bar, N=5/Q=40 does not, and N=5/Q=500 does.
        // No monotone-in-N law or mechanism follows from these three samples.
        double n4_q40 = EnvelopeTheoremWitness.GlobalReading(4, 0.40, 0.01, 25.0, 1600).MaxRiseMagnitude;
        double n5_q40 = EnvelopeTheoremWitness.GlobalReading(5, 0.40, 0.01, 25.0, 1600).MaxRiseMagnitude;
        double n5_q500 = EnvelopeTheoremWitness.GlobalReading(5, 5.0, 0.01, 25.0, 1600).MaxRiseMagnitude;
        Assert.True(n4_q40 > Bar);    // N=4 rises at Q=40
        Assert.True(n5_q40 < Bar);    // N=5/Q=40 is below the reporting bar in this named window
        Assert.True(n5_q500 > Bar);   // but N=5 does rise at Q=500
    }
}
