using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.States;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;
using ComplexVector = MathNet.Numerics.LinearAlgebra.Vector<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Tests.States;

public class PartialTraceTests
{
    private static ComplexMatrix Build(int dim, Func<int, int, Complex> f)
    {
        var m = Matrix<Complex>.Build.Dense(dim, dim);
        for (int i = 0; i < dim; i++)
            for (int j = 0; j < dim; j++)
                m[i, j] = f(i, j);
        return m;
    }

    [Fact]
    public void KeepAll_IsIdentity()
    {
        // Tracing out nothing returns the original matrix.
        int N = 2;
        var rho = Build(4, (i, j) => new Complex(i + j, i - j));
        var reduced = PartialTrace.Of(rho, N, new[] { 0, 1 });
        Assert.Equal(4, reduced.RowCount);
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
                Assert.Equal(rho[i, j], reduced[i, j]);
    }

    [Fact]
    public void KeepNone_IsScalarTrace()
    {
        // Tracing out everything returns the (full trace) 1×1 matrix.
        int N = 2;
        var rho = Build(4, (i, j) => i == j ? new Complex(i + 1, 0) : Complex.Zero);
        var reduced = PartialTrace.Of(rho, N, Array.Empty<int>());
        Assert.Equal(1, reduced.RowCount);
        Assert.Equal(new Complex(1 + 2 + 3 + 4, 0), reduced[0, 0]);
    }

    [Fact]
    public void ProductState_PartialTrace_ReturnsFactor()
    {
        // For ρ = ρ_A ⊗ ρ_B (product), Tr_B[ρ] = Tr(ρ_B) · ρ_A.
        // Take ρ_A = diag(0.3, 0.7), ρ_B = diag(0.4, 0.6); both trace to 1.
        // ρ = diag(0.3·0.4, 0.3·0.6, 0.7·0.4, 0.7·0.6) = diag(0.12, 0.18, 0.28, 0.42).
        // Keep [0] → ρ_A = diag(0.3, 0.7).
        // Keep [1] → ρ_B = diag(0.4, 0.6).
        int N = 2;
        var rho = Matrix<Complex>.Build.Dense(4, 4);
        rho[0, 0] = 0.12; rho[1, 1] = 0.18; rho[2, 2] = 0.28; rho[3, 3] = 0.42;
        var rhoA = PartialTrace.Of(rho, N, new[] { 0 });
        var rhoB = PartialTrace.Of(rho, N, new[] { 1 });
        Assert.Equal(0.3, rhoA[0, 0].Real, 12);
        Assert.Equal(0.7, rhoA[1, 1].Real, 12);
        Assert.Equal(0.4, rhoB[0, 0].Real, 12);
        Assert.Equal(0.6, rhoB[1, 1].Real, 12);
    }

    [Fact]
    public void BellPlusOnTwoOfThree_DoesNotAffectThirdSite()
    {
        // |0⟩ ⊗ |Φ⁺⟩ = (|000⟩ + |011⟩)/√2 on 3 qubits, basis order 0=msb, 2=lsb.
        // Trace out qubit 0 → ρ_(1,2) = |Φ⁺⟩⟨Φ⁺| (the Bell+ on qubits 1, 2).
        int N = 3;
        var psi = ComplexVector.Build.Dense(8);
        double s = 1.0 / Math.Sqrt(2);
        psi[0b000] = s;
        psi[0b011] = s;
        var rho = DensityMatrix.FromStateVector(psi);

        var rhoBell = PartialTrace.Of(rho, N, new[] { 1, 2 });
        Assert.Equal(4, rhoBell.RowCount);
        // Expected Φ⁺ density matrix: 1/2 on (|00⟩+|11⟩)(⟨00|+⟨11|).
        Assert.Equal(0.5, rhoBell[0, 0].Real, 12);
        Assert.Equal(0.5, rhoBell[0, 3].Real, 12);
        Assert.Equal(0.5, rhoBell[3, 0].Real, 12);
        Assert.Equal(0.5, rhoBell[3, 3].Real, 12);
        Assert.Equal(0.0, rhoBell[1, 1].Real, 12);
        Assert.Equal(0.0, rhoBell[2, 2].Real, 12);
    }

    [Fact]
    public void PreservesTrace()
    {
        // Tr(Tr_B[ρ]) = Tr(ρ) for any ρ and any partition.
        int N = 4;
        int d = 1 << N;
        var rho = Matrix<Complex>.Build.Dense(d, d);
        for (int i = 0; i < d; i++) rho[i, i] = new Complex((i + 1) / 136.0, 0);  // trace = 1
        var traceFull = 0.0;
        for (int i = 0; i < d; i++) traceFull += rho[i, i].Real;
        Assert.Equal(1.0, traceFull, 12);

        foreach (int[] keep in new[] { new[] { 0 }, new[] { 1, 2 }, new[] { 0, 2 }, new[] { 1, 3 } })
        {
            var reduced = PartialTrace.Of(rho, N, keep);
            double traceReduced = 0.0;
            for (int i = 0; i < reduced.RowCount; i++) traceReduced += reduced[i, i].Real;
            Assert.Equal(1.0, traceReduced, 12);
        }
    }

    [Fact]
    public void RejectsNonPowerOfTwoMatrix()
    {
        var bad = Matrix<Complex>.Build.Dense(3, 3);
        Assert.Throws<ArgumentException>(() => PartialTrace.Of(bad, 2, new[] { 0 }));
    }

    [Fact]
    public void RejectsKeepIndexOutOfRange()
    {
        var rho = Matrix<Complex>.Build.Dense(4, 4);
        Assert.Throws<ArgumentOutOfRangeException>(() => PartialTrace.Of(rho, 2, new[] { 2 }));
        Assert.Throws<ArgumentOutOfRangeException>(() => PartialTrace.Of(rho, 2, new[] { -1 }));
    }

    [Fact]
    public void RejectsDuplicateKeepIndex()
    {
        var rho = Matrix<Complex>.Build.Dense(4, 4);
        Assert.Throws<ArgumentException>(() => PartialTrace.Of(rho, 2, new[] { 0, 0 }));
    }
}
