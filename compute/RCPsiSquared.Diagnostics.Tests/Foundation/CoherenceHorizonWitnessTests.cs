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
    public void Horizon_MatchesTheCarbonCoherentIncoherentThreshold(int n, double carbon)
    {
        var w = new CoherenceHorizonWitness();
        double q = w.Horizon(n);             // computed live by bisecting Symphony.Clock.Omega
        Assert.Equal(carbon, q, 2);          // 2 decimals: our quantum Q*(N) = carbon's threshold
    }

    [Fact]
    public void N2_IsTheExceptionalPointBase()
    {
        var w = new CoherenceHorizonWitness();
        Assert.Equal(1.0, w.Horizon(2), 3);  // the EP, the rung carbon (N≥3) cannot reach
    }
}
