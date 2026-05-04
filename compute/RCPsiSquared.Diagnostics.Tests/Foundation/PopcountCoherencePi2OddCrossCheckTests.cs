using System.Numerics;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.States;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.DZero;
using RCPsiSquared.Diagnostics.Foundation;
using ComplexVector = MathNet.Numerics.LinearAlgebra.Vector<System.Numerics.Complex>;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>Validates <see cref="PopcountCoherencePi2Odd"/>'s closed-form
/// prediction against the numerical <see cref="MemoryAxisRho"/> reading on
/// popcount-coherence states |ψ⟩ = (|p⟩ + |q⟩)/√2. Python counterpart:
/// <c>simulations/_pi2_odd_landscape_sweep.py</c>.</summary>
public class PopcountCoherencePi2OddCrossCheckTests
{
    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    public void ClosedForm_MatchesMemoryAxisRho_AcrossAllPopcountCoherencePairs(int N)
    {
        var chain = new ChainSystem(N: N, J: 1.0, GammaZero: 0.05);
        var sm = StationaryModes.Compute(chain);
        int d = 1 << N;

        for (int np = 0; np < N; np++)
        {
            int nq = np + 1;
            int p = -1, q = -1;
            for (int x = 0; x < d && p < 0; x++)
            {
                if (System.Numerics.BitOperations.PopCount((uint)x) != np) continue;
                for (int y = 0; y < d; y++)
                {
                    if (System.Numerics.BitOperations.PopCount((uint)y) != nq) continue;
                    if (System.Numerics.BitOperations.PopCount((uint)(x ^ y)) != 1) continue;
                    p = x; q = y;
                    break;
                }
            }
            if (p < 0) continue;

            var psi = ComplexVector.Build.Dense(d);
            psi[p] = 1.0 / Math.Sqrt(2);
            psi[q] = 1.0 / Math.Sqrt(2);
            var rho = DensityMatrix.FromStateVector(psi);

            var reading = MemoryAxisRho.Decompose(rho, chain, sm);
            double predicted = PopcountCoherencePi2Odd.Pi2OddInMemory(N, np, nq);
            double predictedStatic = PopcountCoherencePi2Odd.StaticFraction(N, np, nq);

            Assert.True(Math.Abs(reading.StaticFraction - predictedStatic) < 1e-9,
                $"N={N} (n_p={np}, n_q={nq}): static numeric={reading.StaticFraction:G15}, predicted={predictedStatic:G15}");
            Assert.True(Math.Abs(reading.Pi2OddFractionWithinMemory - predicted) < 1e-9,
                $"N={N} (n_p={np}, n_q={nq}): Π²-odd/mem numeric={reading.Pi2OddFractionWithinMemory:G15}, predicted={predicted:G15}");
        }
    }

    [Fact]
    public void PopcountMirror_AtN5_GivesExact10Over19()
    {
        var chain = new ChainSystem(N: 5, J: 1.0, GammaZero: 0.05);
        var sm = StationaryModes.Compute(chain);
        var psi = ComplexVector.Build.Dense(32);
        psi[0b00011] = 1.0 / Math.Sqrt(2);
        psi[0b00111] = 1.0 / Math.Sqrt(2);
        var rho = DensityMatrix.FromStateVector(psi);

        var reading = MemoryAxisRho.Decompose(rho, chain, sm);
        Assert.Equal(10.0 / 19, reading.Pi2OddFractionWithinMemory, 10);
        Assert.Equal(10.0 / 19, PopcountCoherencePi2Odd.Pi2OddInMemory(5, 2, 3), 10);
    }
}
