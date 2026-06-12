using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Diagnostics.Foundation;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

public class ReadoutFisherTests
{
    [Fact]
    public void BondingState_IsNormalized_AndSineShaped()
    {
        var psi = ReadoutFisher.BondingState(4);
        Assert.Equal(16, psi.Count);
        Assert.Equal(1.0, psi.ConjugateDotProduct(psi).Real, 12);
        // amplitude on |1000⟩ (site 0 = MSB ⇒ index 8) ∝ sin(π·1/5)
        double a0 = psi[8].Real, a1 = psi[4].Real; // |0100⟩ ⇒ index 4
        Assert.Equal(Math.Sin(Math.PI * 1 / 5) / Math.Sin(Math.PI * 2 / 5), a0 / a1, 10);
    }

    [Fact]
    public void Trajectory_PreservesTrace_AtEveryPoint()
    {
        var ts = ReadoutFisher.KGrid(gamma: 0.05, kMax: 1.0, points: 20);
        var rhos = ReadoutFisher.Trajectory(n: 3, j: 1.0, gamma: 0.05,
                                            defectBond: null, deltaJ: 0.0, times: ts);
        Assert.Equal(20, rhos.Count);
        foreach (var rho in rhos)
            Assert.Equal(1.0, rho.Trace().Real, 9);
    }
}
