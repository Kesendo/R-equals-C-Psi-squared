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

    [Fact]
    public void Probs_Normalize_InEveryBasis()
    {
        var ts = ReadoutFisher.KGrid(0.05, 1.0, 5);
        var rho = ReadoutFisher.Trajectory(3, 1.0, 0.05, null, 0.0, ts)[2];
        foreach (var basis in new[] { ReadoutBasis.Z, ReadoutBasis.X, ReadoutBasis.Y })
        {
            var p = ReadoutFisher.Probs(rho, basis);
            Assert.Equal(1.0, p.Sum(), 9);
            Assert.All(p, x => Assert.True(x >= 0));
        }
    }

    // Reference values: two independent implementations agreed on these 2026-06-12
    // (recorded in the design spec). N=4, J=1, bonding state, δJ=0.02 forward
    // difference, strength at bond 0, K-window (0,1], 80 points.
    [Theory]
    [InlineData(ReadoutBasis.Z, 0.05, 1.594, 0.02)]   // Q=20
    [InlineData(ReadoutBasis.Z, 1.00, 0.0715, 0.005)] // Q=1 (the exceptional point)
    [InlineData(ReadoutBasis.X, 0.05, 0.669, 0.02)]   // Q=20
    [InlineData(ReadoutBasis.X, 1.00, 0.0004, 0.0002)]// Q=1
    public void StrengthFi_MatchesReferenceValues(ReadoutBasis basis, double gamma,
                                                  double expected, double tol)
    {
        double fi = ReadoutFisher.StrengthFiMax(n: 4, j: 1.0, gamma: gamma,
            defectBond: 0, deltaJ: 0.02, basis: basis, kMax: 1.0, points: 80);
        Assert.InRange(fi, expected - tol, expected + tol);
    }

    [Fact]
    public void XAndYBases_ReadIdentically_TheModelSymmetry()
    {
        double fx = ReadoutFisher.StrengthFiMax(4, 1.0, 0.2, 0, 0.02, ReadoutBasis.X, 1.0, 80);
        double fy = ReadoutFisher.StrengthFiMax(4, 1.0, 0.2, 0, 0.02, ReadoutBasis.Y, 1.0, 80);
        Assert.Equal(fx, fy, 8);
    }
}
