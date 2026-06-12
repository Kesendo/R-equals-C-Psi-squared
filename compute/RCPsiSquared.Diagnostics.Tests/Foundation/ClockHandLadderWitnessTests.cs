using System;
using System.Linq;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

public class ClockHandLadderWitnessTests
{
    private const double Tol = 1e-9;

    [Theory]
    [InlineData(3, 1.4142135623730951)] // sqrt(2)
    [InlineData(4, 1.6180339887498949)] // golden ratio phi
    [InlineData(5, 1.7320508075688772)] // sqrt(3)
    public void OmegaMem_AtNge3_IsTheF2bBandEdge_2JcosPiOverNPlus1(int n, double expected)
    {
        var w = new ClockHandLadderWitness();
        double omega = w.OmegaMem(n); // live from Symphony.Clock at J=1
        Assert.Equal(expected, omega, 6);
        double bandEdge = 2.0 * 1.0 * Math.Cos(Math.PI / (n + 1));
        Assert.True(Math.Abs(omega - bandEdge) < Tol,
            $"N={n}: live Omega {omega} must equal F2b band edge {bandEdge}");
    }
}
