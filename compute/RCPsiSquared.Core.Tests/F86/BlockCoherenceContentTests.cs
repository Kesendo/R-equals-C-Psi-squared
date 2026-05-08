using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.F86;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;
using ComplexVector = MathNet.Numerics.LinearAlgebra.Vector<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Tests.F86;

public class BlockCoherenceContentTests
{
    private static ComplexMatrix OuterProduct(ComplexVector psi)
    {
        var col = psi.ToColumnMatrix();
        return col * col.ConjugateTranspose();
    }

    private static ComplexVector StateVector(params Complex[] amplitudes) =>
        ComplexVector.Build.DenseOfArray(amplitudes);

    [Fact]
    public void DickeMaxCoherent_AtN2_n0_Saturates_OneQuarter()
    {
        // (|D_0⟩ + |D_1⟩)/√2 = (|00⟩ + (|01⟩+|10⟩)/√2)/√2 saturates Theorem 2 at n=0 on 2 qubits.
        double s = 1.0 / Math.Sqrt(2);
        double t = 0.5; // 1/2 amplitude on each of |01⟩, |10⟩ inside |D_1⟩, then ×1/√2 outer
        var psi = StateVector(s, t, t, 0);
        var rho = OuterProduct(psi);
        Assert.Equal(0.25, BlockCoherenceContent.Compute(rho, n: 0), precision: 12);
    }

    [Fact]
    public void DickeMaxCoherent_AtN2_n1_Saturates_OneQuarter()
    {
        // (|D_1⟩ + |D_2⟩)/√2 = ((|01⟩+|10⟩)/√2 + |11⟩)/√2 saturates Theorem 2 at n=1 on 2 qubits.
        double s = 1.0 / Math.Sqrt(2);
        double t = 0.5;
        var psi = StateVector(0, t, t, s);
        var rho = OuterProduct(psi);
        Assert.Equal(0.25, BlockCoherenceContent.Compute(rho, n: 1), precision: 12);
    }

    [Fact]
    public void BellPlus_HasZero_CBlock_OnBothBlocks()
    {
        // |Φ+⟩ = (|00⟩ + |11⟩)/√2 sits in popcount-{0, 2}; no popcount-1 entries, so both
        // (popcount-0, popcount-1) and (popcount-1, popcount-2) blocks are empty.
        double s = 1.0 / Math.Sqrt(2);
        var psi = StateVector(s, 0, 0, s);
        var rho = OuterProduct(psi);
        Assert.Equal(0.0, BlockCoherenceContent.Compute(rho, n: 0), precision: 14);
        Assert.Equal(0.0, BlockCoherenceContent.Compute(rho, n: 1), precision: 14);
    }

    [Fact]
    public void Singlet_HasZero_CBlock_BothSectorsAreInPopcount1()
    {
        // |S⟩ = (|01⟩ − |10⟩)/√2 sits entirely in popcount-1; no cross-popcount entries.
        double s = 1.0 / Math.Sqrt(2);
        var psi = StateVector(0, s, -s, 0);
        var rho = OuterProduct(psi);
        Assert.Equal(0.0, BlockCoherenceContent.Compute(rho, n: 0), precision: 14);
        Assert.Equal(0.0, BlockCoherenceContent.Compute(rho, n: 1), precision: 14);
    }

    [Fact]
    public void MaximallyMixed_HasZero_CBlock_AllOffDiagonalsVanish()
    {
        var rho = ComplexMatrix.Build.DenseIdentity(4) * (1.0 / 4.0);
        Assert.Equal(0.0, BlockCoherenceContent.Compute(rho, n: 0), precision: 14);
        Assert.Equal(0.0, BlockCoherenceContent.Compute(rho, n: 1), precision: 14);
    }

    [Fact]
    public void Theorem2_CeilingHolds_OnRandomDensityMatrices()
    {
        // Sample 50 random pure states across 2-qubit and 3-qubit Hilbert spaces;
        // every C_block at every (n, n+1) block must satisfy C_block ≤ 1/4 + tol.
        var rng = new Random(42);
        const double tol = 1e-12;
        foreach (int N in new[] { 2, 3 })
        {
            int dim = 1 << N;
            for (int trial = 0; trial < 50; trial++)
            {
                var amps = new Complex[dim];
                for (int i = 0; i < dim; i++)
                    amps[i] = new Complex(rng.NextDouble() * 2 - 1, rng.NextDouble() * 2 - 1);
                var psi = ComplexVector.Build.DenseOfArray(amps);
                psi = psi / psi.L2Norm();
                var rho = OuterProduct(psi);
                for (int n = 0; n < N; n++)
                {
                    double cb = BlockCoherenceContent.Compute(rho, n);
                    Assert.True(cb <= BlockCoherenceContent.Quarter + tol,
                        $"Theorem 2 violation at N={N} n={n} trial={trial}: C_block={cb}");
                }
            }
        }
    }

    [Fact]
    public void Theorem2_CeilingHolds_OnRandomMixedStates()
    {
        // Mixed states: convex combinations of two random pure states. Same ceiling.
        var rng = new Random(7);
        const double tol = 1e-12;
        int N = 2;
        int dim = 1 << N;
        for (int trial = 0; trial < 50; trial++)
        {
            ComplexVector RandomPsi()
            {
                var a = new Complex[dim];
                for (int i = 0; i < dim; i++)
                    a[i] = new Complex(rng.NextDouble() * 2 - 1, rng.NextDouble() * 2 - 1);
                var v = ComplexVector.Build.DenseOfArray(a);
                return v / v.L2Norm();
            }
            double w = rng.NextDouble();
            var rho = w * OuterProduct(RandomPsi()) + (1 - w) * OuterProduct(RandomPsi());
            for (int n = 0; n < N; n++)
            {
                double cb = BlockCoherenceContent.Compute(rho, n);
                Assert.True(cb <= BlockCoherenceContent.Quarter + tol,
                    $"Theorem 2 violation on mixed state at N={N} n={n} trial={trial}: C_block={cb}");
            }
        }
    }

    [Fact]
    public void Compute_RejectsNonSquareMatrix()
    {
        var rho = ComplexMatrix.Build.Dense(4, 3);
        Assert.Throws<ArgumentException>(() => BlockCoherenceContent.Compute(rho, n: 0));
    }

    [Fact]
    public void Compute_RejectsNonPowerOfTwoDimension()
    {
        var rho = ComplexMatrix.Build.DenseIdentity(3);
        Assert.Throws<ArgumentException>(() => BlockCoherenceContent.Compute(rho, n: 0));
    }

    [Fact]
    public void Compute_RejectsOutOfRangeN()
    {
        var rho = ComplexMatrix.Build.DenseIdentity(4);
        Assert.Throws<ArgumentOutOfRangeException>(() => BlockCoherenceContent.Compute(rho, n: -1));
        Assert.Throws<ArgumentOutOfRangeException>(() => BlockCoherenceContent.Compute(rho, n: 2));
    }

    [Fact]
    public void Quarter_EqualsExactly_OneOverFour()
    {
        Assert.Equal(0.25, BlockCoherenceContent.Quarter, precision: 15);
    }
}
