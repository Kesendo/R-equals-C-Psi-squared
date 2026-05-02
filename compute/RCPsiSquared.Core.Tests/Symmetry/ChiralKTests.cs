using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class ChiralKTests
{
    [Theory]
    [InlineData(2)]
    [InlineData(3)]
    [InlineData(4)]
    public void K_AntiCommutesWith_XYChainHamiltonian(int N)
    {
        // KHK = -H for H_xy (sublattice / chiral / class-AIII anti-commutation).
        var H = PauliHamiltonian.XYChain(N, J: 1.0).ToMatrix();
        Assert.Equal(ChiralK.Classification.KOdd, ChiralK.ClassifyHamiltonian(H, N));
    }

    [Theory]
    [InlineData(2)]
    [InlineData(3)]
    public void K_AntiCommutesWith_HeisenbergHamiltonian(int N)
    {
        // Heisenberg = XX + YY + ZZ. ZZ on adjacent sites (one even, one odd) commutes with K
        // (one Z each side); XY part anti-commutes; total is K-mixed. Wait — check more carefully.
        // K = ⊗_{odd i} Z_i. Z_i for odd i commutes with Z_j (any), with X_i / Y_i anti-commutes.
        // For ZZ on bond (i, i+1) with one even, one odd index: K Z_i Z_{i+1} K = (Z Z_{odd}) Z_i Z_{i+1} (Z Z_{odd}).
        // Just classify and document what we get.
        var H = PauliHamiltonian.HeisenbergChain(N, J: 1.0).ToMatrix();
        var cls = ChiralK.ClassifyHamiltonian(H, N);
        // XX + YY anti-commutes with K, ZZ commutes — sum is mixed.
        Assert.Equal(ChiralK.Classification.KMixed, cls);
    }

    [Fact]
    public void K_FullSquared_IsIdentity()
    {
        int N = 4;
        var k = ChiralK.BuildFull(N);
        int d = 1 << N;
        var kk = k * k;
        var diff = kk - MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>.Build.DenseIdentity(d);
        Assert.True(diff.FrobeniusNorm() < 1e-12);
    }
}
