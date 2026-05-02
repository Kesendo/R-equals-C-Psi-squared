using System.Numerics;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.States;
using RCPsiSquared.Diagnostics.DZero;

namespace RCPsiSquared.Diagnostics.Tests.DZero;

public class DZeroTests
{
    [Fact]
    public void StationaryModes_OnXYChain_HasNPlusOneKernelDimensions()
    {
        // Per F4: uniform XY chain + Z-dephasing has N+1 stationary modes (sector projectors).
        var chain = new ChainSystem(N: 3, J: 1.0, GammaZero: 0.05);
        var sm = StationaryModes.Compute(chain);
        Assert.Equal(chain.N + 1, sm.KernelDimension);
        // All kernel eigenvalues should be ~ 0
        foreach (var lambda in sm.Eigenvalues)
            Assert.True(lambda.Magnitude < 1e-9, $"kernel eigenvalue {lambda} not near zero");
    }

    [Fact]
    public void DZeroDecomposition_Sums_RhoD0_PlusRhoD2_EqualsRho()
    {
        var chain = new ChainSystem(N: 2, J: 1.0, GammaZero: 0.05);
        var psi = PolarityState.Build(N: 2, signs: new[] { +1, -1 });
        var rho = DensityMatrix.FromStateVector(psi);
        var dec = DZeroDecomposition.Decompose(rho, chain);
        var sum = dec.RhoD0 + dec.RhoD2;
        Assert.True((sum - rho).FrobeniusNorm() < 1e-10);
    }

    [Fact]
    public void SectorPopulations_PolarityState_HasBinomialDistribution()
    {
        // |+⟩^N in computational basis has equal amplitudes 1/√(2^N) on every basis state,
        // so p_n = (#k with popcount n) / 2^N = C(N, n) / 2^N (binomial).
        int N = 4;
        var psi = PolarityState.Uniform(N, sign: +1);
        var sec = SectorPopulations.FromStateVector(psi);
        Assert.Equal(N, sec.N);
        // Σ p_n = 1
        Assert.Equal(1.0, sec.P.Sum(), 10);
        // ⟨n⟩ = N/2 by symmetry of |+⟩^N
        Assert.Equal(N / 2.0, sec.MeanN, 10);
    }

    [Fact]
    public void SectorPopulations_HammingWeightVacuum_IsDelta()
    {
        // |0…0⟩ has popcount 0 → p_0 = 1, others 0.
        int N = 3;
        var psi = MathNet.Numerics.LinearAlgebra.Vector<Complex>.Build.Dense(1 << N);
        psi[0] = Complex.One;
        var sec = SectorPopulations.FromStateVector(psi);
        Assert.Equal(1.0, sec.P[0], 10);
        for (int n = 1; n <= N; n++) Assert.Equal(0.0, sec.P[n], 10);
        Assert.Equal(0.0, sec.MeanN, 10);
        Assert.Equal(0.0, sec.Entropy, 10);
    }
}
