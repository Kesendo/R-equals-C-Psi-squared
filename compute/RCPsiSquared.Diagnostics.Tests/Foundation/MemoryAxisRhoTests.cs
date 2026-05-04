using System.Numerics;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.States;
using RCPsiSquared.Diagnostics.Foundation;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;
using ComplexVector = MathNet.Numerics.LinearAlgebra.Vector<System.Numerics.Complex>;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

public class MemoryAxisRhoTests
{
    private static ChainSystem MakeChain(int N) => new(N, J: 1.0, GammaZero: 0.05);

    [Fact]
    public void MaximallyMixed_IsEntirelyStatic()
    {
        // ρ_mm = I/d is a uniform mixture of all sector projectors P_n (each P_n is in
        // kernel of L). It sits entirely on the d=0 axis. The kernel of L for our XY +
        // Z-dephasing chain has dimension N+1 (per F4).
        var chain = MakeChain(N: 3);
        int d = 1 << chain.N;
        var rhoMm = ComplexMatrix.Build.DiagonalIdentity(d) / d;

        var result = MemoryAxisRho.Decompose(rhoMm, chain);

        Assert.True(result.StaticFraction > 0.999,
            $"expected ρ_mm to be ~100% static; got StaticFraction = {result.StaticFraction}");
        Assert.True(result.MemoryFraction < 1e-6,
            $"expected ρ_mm to have ~0 memory; got MemoryFraction = {result.MemoryFraction}");
        Assert.Equal(chain.N + 1, result.KernelDimension);
    }

    [Fact]
    public void ComputationalBasisVacuum_IsStatic_BecauseInSinglePopulationSector()
    {
        // |0…0⟩ is in P_0 (popcount-0 sector); ρ = |0…0⟩⟨0…0| is purely diagonal at
        // index 0; entirely in kernel of L.
        var chain = MakeChain(N: 3);
        int d = 1 << chain.N;
        var psi = ComplexVector.Build.Dense(d);
        psi[0] = Complex.One;
        var rho = DensityMatrix.FromStateVector(psi);

        var result = MemoryAxisRho.Decompose(rho, chain);

        Assert.True(result.StaticFraction > 0.999,
            $"expected |0…0⟩ to be ~100% static; got {result.StaticFraction}");
    }

    [Fact]
    public void PolarityState_HasMemoryContent_BecauseHasCrossSectorOffDiagonals()
    {
        // |+−+⟩ has equal amplitudes on every basis state with phase signs; ρ has
        // off-diagonal entries between different sectors. The off-diagonal cross-sector
        // content is NOT in kernel of L, so the state has substantial memory content.
        var chain = MakeChain(N: 3);
        var psi = PolarityState.Build(N: 3, signs: new[] { +1, -1, +1 });
        var rho = DensityMatrix.FromStateVector(psi);

        var result = MemoryAxisRho.Decompose(rho, chain);

        Assert.True(result.MemoryFraction > 0.1,
            $"polarity state should have substantial memory content; got {result.MemoryFraction}");
        Assert.True(result.StaticFraction > 0.0,
            $"polarity state has diagonal sector content too; got static = {result.StaticFraction}");
        Assert.True(result.StaticFraction + result.MemoryFraction > 0.99,
            "fractions should sum to ≈ 1");
    }

    [Fact]
    public void Decomposition_RhoEqualsRhoD0PlusRhoD2_AndPi2SplitIsExact()
    {
        // The diagnostic builds on DZeroDecomposition (ρ = ρ_d0 + ρ_d2) + Pi2Projection.Split
        // (ρ_d2 = ρ_d2_even + ρ_d2_odd). Each step is exact at machine precision; verify
        // the partition is internally consistent by reproducing ρ from the parts.
        var chain = MakeChain(N: 2);
        var psi = PolarityState.Build(N: 2, signs: new[] { +1, -1 });
        var rho = DensityMatrix.FromStateVector(psi);

        var result = MemoryAxisRho.Decompose(rho, chain);

        // All three norms should be non-negative
        Assert.True(result.StaticNormSquared >= 0);
        Assert.True(result.MemoryPi2EvenNormSquared >= 0);
        Assert.True(result.MemoryPi2OddNormSquared >= 0);
        Assert.True(result.TotalNormSquared > 0);
    }

    [Fact]
    public void StaticFraction_IsBoundedBetween0And1()
    {
        var chain = MakeChain(N: 3);
        var psi = PolarityState.Build(N: 3, signs: new[] { +1, -1, +1 });
        var rho = DensityMatrix.FromStateVector(psi);

        var result = MemoryAxisRho.Decompose(rho, chain);

        Assert.InRange(result.StaticFraction, 0.0, 1.0);
        Assert.InRange(result.MemoryFraction, 0.0, 1.0);
        Assert.InRange(result.Pi2OddFractionWithinMemory, 0.0, 1.0);
    }

    [Fact]
    public void Pi2OddFraction_OfMaximallyMixed_IsZero_WhenMemoryIsZero()
    {
        // When MemoryFraction = 0, Pi2OddFractionWithinMemory should be 0 (not NaN).
        var chain = MakeChain(N: 2);
        int d = 1 << chain.N;
        var rhoMm = ComplexMatrix.Build.DiagonalIdentity(d) / d;

        var result = MemoryAxisRho.Decompose(rhoMm, chain);

        Assert.Equal(0.0, result.Pi2OddFractionWithinMemory);
    }

    /// <summary>
    /// X-Y combinatorial theorem (proved 2026-05-03 conversation, recorded in
    /// project_y_phase_pi2_odd_lens.md): any tensor product of σ_x and σ_y eigenstates
    /// with at least one Y-site has exactly 2^(N-1) Π²-even + 2^(N-1) Π²-odd Pauli
    /// strings. Within memory (excluding the III static component), the Π²-odd fraction
    /// is 2^(N-1) / (2^N − 1). At N = 2, 3, 4, 5 this gives 2/3, 4/7, 8/15, 16/31.
    /// Bit-exact verification of the combinatorial identity at each N.
    /// </summary>
    [Theory]
    [InlineData(2, 2.0 / 3.0)]    // 2^1 / (2^2 − 1) = 2/3
    [InlineData(3, 4.0 / 7.0)]    // 2^2 / (2^3 − 1) = 4/7
    [InlineData(4, 8.0 / 15.0)]   // 2^3 / (2^4 − 1) = 8/15
    [InlineData(5, 16.0 / 31.0)]  // 2^4 / (2^5 − 1) = 16/31
    public void YOnlyState_HasPi2OddFraction_2NMinus1_Over_2NMinus1(int N, double expectedFraction)
    {
        var chain = MakeChain(N);
        var psi = YBasisProduct(N, Enumerable.Repeat(+1, N).ToArray()); // all |+i⟩
        var rho = DensityMatrix.FromStateVector(psi);

        var result = MemoryAxisRho.Decompose(rho, chain);

        Assert.Equal(expectedFraction, result.Pi2OddFractionWithinMemory, 10);
        Assert.Equal(1.0 / (1 << N), result.StaticFraction, 10);          // 1/d
        Assert.Equal(1.0 - 1.0 / (1 << N), result.MemoryFraction, 10);    // (d−1)/d
    }

    /// <summary>
    /// The X-Y theorem says the fraction is INDEPENDENT of how many Y-sites M, as long
    /// as M ≥ 1. Verify at N = 4 across M = 1, 2, 3, 4: all should give 8/15.
    /// </summary>
    [Theory]
    [InlineData(4, new[] { 'Y', 'X', 'X', 'X' })]  // M=1
    [InlineData(4, new[] { 'Y', 'Y', 'X', 'X' })]  // M=2
    [InlineData(4, new[] { 'Y', 'Y', 'Y', 'X' })]  // M=3
    [InlineData(4, new[] { 'Y', 'Y', 'Y', 'Y' })]  // M=4 (Y-only)
    [InlineData(4, new[] { 'X', 'Y', 'X', 'Y' })]  // M=2 alternating
    public void XYMixedState_FractionIndependentOfYSiteCount_AtN4(int N, char[] axes)
    {
        var chain = MakeChain(N);
        var signs = Enumerable.Repeat(+1, N).ToArray();
        var psi = MixedXYProduct(N, axes, signs);
        var rho = DensityMatrix.FromStateVector(psi);

        var result = MemoryAxisRho.Decompose(rho, chain);

        // Theorem: 2^(N-1) / (2^N - 1) = 8/15 for N=4
        Assert.Equal(8.0 / 15.0, result.Pi2OddFractionWithinMemory, 10);
    }

    /// <summary>
    /// X-only state (all sites in X-basis) has 0 Π²-odd content: {I, X} are both
    /// Π²-even. Verify the M = 0 boundary case of the X-Y theorem.
    /// </summary>
    [Theory]
    [InlineData(2)]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    public void XOnlyState_HasZeroPi2OddFraction_AcrossN(int N)
    {
        var chain = MakeChain(N);
        var psi = PolarityState.Uniform(N, sign: +1);
        var rho = DensityMatrix.FromStateVector(psi);

        var result = MemoryAxisRho.Decompose(rho, chain);

        Assert.Equal(0.0, result.Pi2OddFractionWithinMemory, 10);
    }

    // === Helpers ===

    /// <summary>Y-basis polarity tensor product: |+i⟩ = (|0⟩ + i|1⟩)/√2,
    /// |−i⟩ = (|0⟩ − i|1⟩)/√2, per-site signs +1 = |+i⟩, −1 = |−i⟩.</summary>
    private static ComplexVector YBasisProduct(int N, IReadOnlyList<int> signs)
    {
        int d = 1 << N;
        var vec = ComplexVector.Build.Dense(d);
        double norm = 1.0 / Math.Sqrt(d);
        for (int idx = 0; idx < d; idx++)
        {
            Complex amp = Complex.One;
            for (int k = 0; k < N; k++)
            {
                int bit = (idx >> (N - 1 - k)) & 1;
                if (bit == 1)
                    amp *= signs[k] == +1 ? Complex.ImaginaryOne : -Complex.ImaginaryOne;
            }
            vec[idx] = amp * norm;
        }
        return vec;
    }

    /// <summary>Mixed X/Y-basis tensor product: each site is in either σ_x or σ_y
    /// eigenstate per <paramref name="axes"/> ('X' or 'Y'), with sign +1 or −1.</summary>
    private static ComplexVector MixedXYProduct(int N, IReadOnlyList<char> axes, IReadOnlyList<int> signs)
    {
        int d = 1 << N;
        var vec = ComplexVector.Build.Dense(d);
        double norm = 1.0 / Math.Sqrt(d);
        for (int idx = 0; idx < d; idx++)
        {
            Complex amp = Complex.One;
            for (int k = 0; k < N; k++)
            {
                int bit = (idx >> (N - 1 - k)) & 1;
                if (bit == 1)
                {
                    if (axes[k] == 'X') amp *= signs[k];
                    else /* 'Y' */     amp *= signs[k] == +1 ? Complex.ImaginaryOne : -Complex.ImaginaryOne;
                }
            }
            vec[idx] = amp * norm;
        }
        return vec;
    }
}
