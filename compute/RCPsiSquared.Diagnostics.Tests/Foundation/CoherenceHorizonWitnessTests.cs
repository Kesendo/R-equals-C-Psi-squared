using System;
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
}
