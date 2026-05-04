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
/// popcount-coherence pair states |ψ⟩ = (|p⟩ + |q⟩)/√2 across all anchor
/// categories (HD = N Π²-classical, popcount-mirror, intra-mirror,
/// K-intermediate, generic). Python counterpart:
/// <c>simulations/_pi2_odd_general_closed_form.py</c>.</summary>
public class PopcountCoherencePi2OddCrossCheckTests
{
    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    public void ClosedForm_MatchesMemoryAxisRho_AllCompatiblePopcountCoherencePairs(int N)
    {
        var chain = new ChainSystem(N: N, J: 1.0, GammaZero: 0.05);
        var sm = StationaryModes.Compute(chain);
        int d = 1 << N;

        for (int np = 0; np <= N; np++)
        {
            for (int nq = np; nq <= N; nq++)
            {
                int hdMin = (np == nq) ? 2 : Math.Abs(nq - np);
                int hdMax = Math.Min(np + nq, N);
                int parity = (nq - np) % 2;
                for (int hd = hdMin; hd <= hdMax; hd++)
                {
                    if ((hd % 2) != parity) continue;

                    var (p, q) = FirstPair(N, np, nq, hd);
                    if (p < 0) continue;

                    var psi = ComplexVector.Build.Dense(d);
                    psi[p] = 1.0 / Math.Sqrt(2);
                    psi[q] = 1.0 / Math.Sqrt(2);
                    var rho = DensityMatrix.FromStateVector(psi);

                    var reading = MemoryAxisRho.Decompose(rho, chain, sm);
                    double predictedStatic = PopcountCoherencePi2Odd.StaticFraction(N, np, nq);
                    double predicted = PopcountCoherencePi2Odd.Pi2OddInMemory(N, np, nq, hd);

                    Assert.True(Math.Abs(reading.StaticFraction - predictedStatic) < 1e-9,
                        $"N={N} (n_p={np}, n_q={nq}, HD={hd}): static numeric={reading.StaticFraction:G15}, predicted={predictedStatic:G15}");
                    Assert.True(Math.Abs(reading.Pi2OddFractionWithinMemory - predicted) < 1e-9,
                        $"N={N} (n_p={np}, n_q={nq}, HD={hd}): Π²-odd/mem numeric={reading.Pi2OddFractionWithinMemory:G15}, predicted={predicted:G15}");
                }
            }
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
        Assert.Equal(10.0 / 19, PopcountCoherencePi2Odd.Pi2OddInMemory(5, 2, 3, hd: 1), 10);
    }

    [Theory]
    [InlineData(2, 0b00, 0b11)]              // GHZ_2 / Bell+ at popcount-(0, 2) HD = 2 = N
    [InlineData(3, 0b000, 0b111)]            // GHZ_3 popcount-(0, 3) HD = 3 = N
    [InlineData(4, 0b0000, 0b1111)]          // GHZ_4 popcount-(0, 4) HD = 4 = N
    [InlineData(2, 0b01, 0b10)]              // Singlet/Triplet popcount-(1, 1) HD = 2 = N
    public void HdComplementStates_ArePi2Classical(int N, int p, int q)
    {
        var chain = new ChainSystem(N: N, J: 1.0, GammaZero: 0.05);
        var sm = StationaryModes.Compute(chain);
        int d = 1 << N;
        var psi = ComplexVector.Build.Dense(d);
        psi[p] = 1.0 / Math.Sqrt(2);
        psi[q] = 1.0 / Math.Sqrt(2);
        var rho = DensityMatrix.FromStateVector(psi);

        var reading = MemoryAxisRho.Decompose(rho, chain, sm);
        Assert.Equal(0.0, reading.Pi2OddFractionWithinMemory, 10);

        int np = System.Numerics.BitOperations.PopCount((uint)p);
        int nq = System.Numerics.BitOperations.PopCount((uint)q);
        int hd = System.Numerics.BitOperations.PopCount((uint)(p ^ q));
        Assert.Equal(0.0, PopcountCoherencePi2Odd.Pi2OddInMemory(N, np, nq, hd), 10);
    }

    private static (int p, int q) FirstPair(int N, int np, int nq, int hd)
    {
        int d = 1 << N;
        for (int x = 0; x < d; x++)
        {
            if (System.Numerics.BitOperations.PopCount((uint)x) != np) continue;
            for (int y = 0; y < d; y++)
            {
                if (System.Numerics.BitOperations.PopCount((uint)y) != nq) continue;
                if (x == y) continue;
                if (System.Numerics.BitOperations.PopCount((uint)(x ^ y)) == hd)
                    return (x, y);
            }
        }
        return (-1, -1);
    }
}
