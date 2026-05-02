using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Core.Pauli;

namespace RCPsiSquared.Core.Tests.Lindblad;

public class LindbladianBuilderTests
{
    [Fact]
    public void Identity_DensityMatrix_IsStationary_ForAnyL()
    {
        // L · vec(I_d / d) = 0 for any Lindbladian, since L preserves trace and I/d is the
        // maximally mixed state (a fixed point of any unital channel).
        // For 1-qubit: L_z (I/2) = -i[H, I/2] + γ(Z·I/2·Z - I/2) = 0 + γ(I/2 - I/2) = 0.
        int N = 2;
        int d = 1 << N;
        var H = PauliHamiltonian.XYChain(N, J: 1.0).ToMatrix();
        var L = PauliDephasingDissipator.BuildZ(H, new[] { 0.05, 0.05 });

        var rhoI = Matrix<Complex>.Build.DenseIdentity(d) / d;
        // I/d is symmetric so vec_F and vec_C agree; the test is convention-agnostic. We pick
        // one flatten and check L·vec(I/d) = 0 (any unital channel fixes the maximally mixed).
        var vecRhoI = FlattenSymmetric(rhoI);
        var result = L * vecRhoI;
        for (int i = 0; i < result.Count; i++)
            Assert.True(result[i].Magnitude < 1e-10, $"L · vec(I/d) [{i}] = {result[i]}, expected 0");
    }

    [Fact]
    public void HamiltonianOnly_L_IsAntiHermitian()
    {
        // For pure unitary evolution dρ/dt = -i[H, ρ], the Liouvillian L = -i[H, ·] is
        // anti-Hermitian: L^† = -L (in operator sense on vec space).
        int N = 2;
        var H = PauliHamiltonian.XYChain(N, J: 1.0).ToMatrix();
        var L = LindbladianBuilder.Build(H, Array.Empty<Matrix<Complex>>());
        var sum = L + L.ConjugateTranspose();
        Assert.True(sum.FrobeniusNorm() < 1e-10, $"L + L† should be 0 for unitary L; got Frobenius {sum.FrobeniusNorm()}");
    }

    [Fact]
    public void Hermiticity_RejectsNonHermitianH()
    {
        var nonH = Matrix<Complex>.Build.DenseOfArray(new Complex[,]
        {
            { 0, 1 },
            { 0, 0 },
        });
        Assert.Throws<ArgumentException>(() => LindbladianBuilder.Build(nonH, Array.Empty<Matrix<Complex>>()));
    }

    /// <summary>Flatten a SYMMETRIC d×d matrix. vec_F and vec_C agree for symmetric inputs, so
    /// either convention gives the same vector. Caller responsibility to ensure symmetry — only
    /// safe for tests checking convention-agnostic identities (here: L · vec(I/d) = 0).</summary>
    private static MathNet.Numerics.LinearAlgebra.Vector<Complex> FlattenSymmetric(Matrix<Complex> m)
    {
        int d = m.RowCount;
        var v = MathNet.Numerics.LinearAlgebra.Vector<Complex>.Build.Dense(d * d);
        for (int i = 0; i < d; i++)
            for (int j = 0; j < d; j++)
                v[i * d + j] = m[i, j];
        return v;
    }
}
