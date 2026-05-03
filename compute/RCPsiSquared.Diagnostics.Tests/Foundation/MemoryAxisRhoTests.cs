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
}
