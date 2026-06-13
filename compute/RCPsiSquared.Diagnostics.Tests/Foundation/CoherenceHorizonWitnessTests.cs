using System;
using System.Collections.Generic;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

public class CoherenceHorizonWitnessTests
{
    [Theory]
    [InlineData(2, 1.0)]
    [InlineData(3, 1.41421)]   // √2
    [InlineData(4, 1.8785)]
    [InlineData(5, 2.37217)]
    public void Horizon_MatchesTheLiveComputedThreshold(int n, double expectedQStar)
    {
        var w = new CoherenceHorizonWitness();
        double q = w.Horizon(n);                  // computed live by bisecting Symphony.Clock.Omega
        Assert.Equal(expectedQStar, q, 2);        // the precise live Q*(N) at 2 decimals
    }

    /// <summary>The literal "quantum = carbon" check: the live-computed Horizon(n) equals the actual
    /// Frost-Hückel coherent↔incoherent doc value from docs/carbon/FROST_CIRCLE_AS_THE_CLOCK_FACE.md
    /// (√2 / 1.879 / 2.372 at N=3/4/5) at 2-decimal tolerance. N=2 is excluded: it is the EP base, the
    /// rung the polyene layer (N≥3) cannot reach, so it has no carbon doc value.</summary>
    [Theory]
    [InlineData(3, 1.41421356)]   // √2
    [InlineData(4, 1.879)]
    [InlineData(5, 2.372)]
    public void Horizon_EqualsCarbonDocValue(int n, double carbonDoc)
    {
        var w = new CoherenceHorizonWitness();
        Assert.Equal(carbonDoc, w.Horizon(n), 2);  // quantum coherence horizon = carbon Frost-Hückel threshold
    }

    [Fact]
    public void N2_IsTheExceptionalPointBase()
    {
        var w = new CoherenceHorizonWitness();
        Assert.Equal(1.0, w.Horizon(2), 3);  // the EP, the rung carbon (N≥3) cannot reach
    }

    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    public void EpModes_CoalescerIsThe02CoherenceWithRigidityToZero(int n)
    {
        var w = new CoherenceHorizonWitness();
        var ep = w.EpModes(n);
        // the coalescing gap mode is a genuine EP: rigidity collapses
        Assert.True(ep.Coalescer.Rigidity < 0.05,
            $"N={n}: coalescer rigidity {ep.Coalescer.Rigidity:F4} should be ~0");
        // and it is the {0,2}-coherence: weight split between n_diff 0 and 2, mean ≈ 1 (within bisection
        // / cluster noise; the {0,2} block has mean exactly 1, the live value drifts a few % off Q*).
        Assert.InRange(ep.CoalescerMeanNDiff, 0.8, 1.2);
        Assert.True(ep.CoalescerHist.GetValueOrDefault(0) > 0.3 && ep.CoalescerHist.GetValueOrDefault(2) > 0.3,
            $"N={n}: coalescer histogram should be ~{{0:1/2, 2:1/2}}, got 0={ep.CoalescerHist.GetValueOrDefault(0):F2} 2={ep.CoalescerHist.GetValueOrDefault(2):F2}");
    }

    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    public void EpModes_BandEdgeIsTheCoLocatedSurvivor(int n)
    {
        var w = new CoherenceHorizonWitness();
        var ep = w.EpModes(n);
        double bandIm = 2.0 * System.Math.Cos(System.Math.PI / (n + 1));
        // the band edge oscillates at 2cos(π/(N+1)) and does NOT coalesce: rigidity stays high
        Assert.Equal(bandIm, System.Math.Abs(ep.BandEdge.Lambda.Imaginary), 2);
        Assert.True(ep.BandEdgeR > 0.5,
            $"N={n}: band-edge rigidity {ep.BandEdgeR:F3} should stay near 1 (the survivor)");
    }

    [Fact]
    public void SqrtScalingRatio_IsConstant_AtN4_GenuineSecondOrderEp()
    {
        // Near a 2nd-order EP, Im² ∝ (Q−Q*): the ratio Im²/(Q−Q*) is the same at two offsets.
        var w = new CoherenceHorizonWitness();
        double r1 = w.SqrtScalingRatio(4, 0.03);
        double r2 = w.SqrtScalingRatio(4, 0.06);
        Assert.True(r1 > 0 && r2 > 0, $"ratios should be positive, got {r1}, {r2}");
        Assert.True(System.Math.Abs(r1 - r2) / r1 < 0.30,
            $"√-scaling ratio should be ~constant (2nd-order EP); got {r1:F3} vs {r2:F3}");
    }
}
