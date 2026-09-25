using System.Linq;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>The finite envelope_n4_rise atlas, gate-first. These tests preserve the named N=3/4/5,
/// Q, K-window and grid rows on the XY chain and reuse the witness's detector. They test finite
/// observations and one same-(N,Q,K) rescaling comparison; they do not infer an all-Q/all-N threshold,
/// an absence theorem, or a mechanism. See experiments/ENVELOPE_RISE_BOUNDARY.md.</summary>
public class EnvelopeBoundaryTests
{
    private const double Bar = EnvelopeTheoremWitness.RiseReportingBar;

    [Fact]
    public void SameNQKRescaling_AgreesToEigensolverRounding_OnTheNamedPair()
    {
        // (J, γ, tMax) = (5, 0.01, 25) and (10, 0.02, 12.5) share (N, Q, K_max) = (4, 500, 0.25). Doubling
        // J and γ doubles the Liouvillian entry for entry and halving tMax halves every sample time, both
        // exact in binary64, so the K-grids coincide exactly and the two curves differ only through the
        // eigensolver's rounding on a doubled matrix. Error model: that rounding is O(κ(R)·ε) on the
        // curve, measured 1.4e-14 on CΨ and 1.7e-15 on the largest rise; the gate 1e-12 sits two decades
        // above it. The mutation J = 10.2 (Q = 510) moves the curve by more than 1e-6.
        var a = new Symphony(n: 4, j: 5.0, gamma: 0.01, initialState: InitialStateKind.BellPair,
            tMax: 25.0, tPoints: 1600);
        var b = new Symphony(n: 4, j: 10.0, gamma: 0.02, initialState: InitialStateKind.BellPair,
            tMax: 12.5, tPoints: 1600);
        for (int i = 0; i < a.TimeGrid.Count; i++)
            Assert.True(a.Gamma * a.TimeGrid[i] == b.Gamma * b.TimeGrid[i], $"K-grid sample {i}");

        var curveA = a.States.Select(Symphony.Cpsi).ToArray();
        var curveB = b.States.Select(Symphony.Cpsi).ToArray();
        double curveDiff = TempoCertificationMovement.MaxAbsDiff(curveA, curveB);
        Assert.True(curveDiff <= 1e-12, $"curve difference {curveDiff:R}");

        // The witness's detector on both curves: the same count, the largest rise within the model.
        var ra = QuarterEnvelope.Of(curveA, a.TimeGrid.ToArray(), riseTol: Bar);
        var rb = QuarterEnvelope.Of(curveB, b.TimeGrid.ToArray(), riseTol: Bar);
        Assert.Equal(ra.RiseCount, rb.RiseCount);
        Assert.True(ra.RiseCount > 0);
        Assert.True(System.Math.Abs(ra.MaxRiseMagnitude - rb.MaxRiseMagnitude) <= 1e-12,
            $"largest-rise difference {System.Math.Abs(ra.MaxRiseMagnitude - rb.MaxRiseMagnitude):R}");

        var mutated = new Symphony(n: 4, j: 10.2, gamma: 0.02, initialState: InitialStateKind.BellPair,
            tMax: 12.5, tPoints: 1600);
        Assert.True(TempoCertificationMovement.MaxAbsDiff(curveA,
            mutated.States.Select(Symphony.Cpsi).ToArray()) > 1e-6);
    }

    [Fact]
    public void NamedRows_N3AtQ2000MaximaFallStrictly_N4AtQ500ResolvesARise()
    {
        // One N=3 null sample at Q=2000 and one N=4 positive sample at Q=500. The N=3 row is raw: all 225
        // maxima fall, the closest pair by 4.0e-4, so the largest positive predecessor delta is exactly
        // zero. They invite wider all-Q/all-N classification; they do not prove an N floor or identify
        // the internal mechanism.
        var n3 = EnvelopeTheoremWitness.GlobalReading(3, 20.0, 0.01, 25.0, 1600);           // Q=2000
        Assert.Equal(0, n3.RiseCount);
        Assert.True(n3.IsNonIncreasing);
        Assert.True(n3.MaxRiseMagnitude == 0.0, $"raw max rise {n3.MaxRiseMagnitude:R}");
        Assert.True(EnvelopeTheoremWitness.GlobalReading(4, 5.0, 0.01, 25.0, 1600).MaxRiseMagnitude
                    > Bar);                                                                    // Q=500: N=4 rises
    }

    [Fact]
    public void NamedN4Rows_Q13AndQ40_GiveDifferentFiniteClassifications()
    {
        // On this fixed K-window/grid the Q=40 row exceeds the reporting bar, and the Q=13 row has a
        // single maximum, so no predecessor exists to rise above.
        Assert.True(EnvelopeTheoremWitness.GlobalReading(4, 0.40, 0.01, 25.0, 1600).MaxRiseMagnitude
                    > Bar);                                                                    // Q=40: rises
        var q13 = EnvelopeTheoremWitness.GlobalReading(4, 0.13, 0.01, 25.0, 1600);
        Assert.Single(q13.Maxima);                                                             // Q=13
        Assert.True(q13.MaxRiseMagnitude == 0.0);
    }

    [Fact]
    public void NamedN4AndN5Rows_DifferAtQ40_AndN5RisesAtQ500()
    {
        // Three retained atlas rows: N=4/Q=40 exceeds the bar, N=5/Q=40 does not, and N=5/Q=500 does.
        // No monotone-in-N law or mechanism follows from these three samples. The N=5/Q=40 row does
        // carry a raw 3.1e-4 rise that refinement keeps, so the bar is a reporting convention here.
        double n4_q40 = EnvelopeTheoremWitness.GlobalReading(4, 0.40, 0.01, 25.0, 1600).MaxRiseMagnitude;
        double n5_q40 = EnvelopeTheoremWitness.GlobalReading(5, 0.40, 0.01, 25.0, 1600).MaxRiseMagnitude;
        double n5_q500 = EnvelopeTheoremWitness.GlobalReading(5, 5.0, 0.01, 25.0, 1600).MaxRiseMagnitude;
        Assert.True(n4_q40 > Bar);    // N=4 rises at Q=40
        Assert.True(n5_q40 < Bar);    // N=5/Q=40 is below the reporting bar in this named window
        Assert.True(n5_q40 > 0.0);    // yet it carries a raw rise
        Assert.True(n5_q500 > Bar);   // and N=5 rises above the bar at Q=500
    }
}
