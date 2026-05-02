using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.Pauli;

namespace RCPsiSquared.Core.Tests.Pauli;

public class PauliBasisTests
{
    [Theory]
    [InlineData(1)]
    [InlineData(2)]
    [InlineData(3)]
    public void DensityMatrix_RoundTrip_ThroughPauliBasis(int N)
    {
        // Build a random Hermitian density matrix, decompose to Pauli basis, reconstruct, compare.
        int d = 1 << N;
        var rho = Matrix<Complex>.Build.Random(d, d);
        rho = (rho + rho.ConjugateTranspose()) / 2.0;  // Hermitianise

        var pauliVec = PauliBasis.ToPauliVector(rho, N);
        var rhoBack = PauliBasis.FromPauliVector(pauliVec, N);

        for (int i = 0; i < d; i++)
            for (int j = 0; j < d; j++)
            {
                Assert.Equal(rho[i, j].Real, rhoBack[i, j].Real, 10);
                Assert.Equal(rho[i, j].Imaginary, rhoBack[i, j].Imaginary, 10);
            }
    }

    [Fact]
    public void Identity_DecomposesTo_Single_I_Component()
    {
        // For ρ = I_d (d×d identity): vec[I-string] = Tr(I_d) / d = 1, vec[else] = Tr(σ_k) / d = 0.
        int N = 2;
        int d = 1 << N;
        var rho = Matrix<Complex>.Build.DenseIdentity(d);

        var vec = PauliBasis.ToPauliVector(rho, N);
        Assert.Equal(1.0, vec[0].Real, 10);
        for (int k = 1; k < vec.Count; k++)
        {
            Assert.Equal(0.0, vec[k].Real, 10);
            Assert.Equal(0.0, vec[k].Imaginary, 10);
        }
    }

    [Fact]
    public void SingleQubit_BlochComponents_RecoverFromExpectations()
    {
        // 1-qubit Bloch parametrisation: ρ = (I + bX X + bY Y + bZ Z) / 2.
        // PauliBasis.ToPauliVector returns vec = (1/2, bX/2, bZ/2, bY/2) in I,X,Z,Y order.
        var bx = 0.3; var by = -0.4; var bz = 0.5;
        var rho = (Matrix<Complex>.Build.DenseIdentity(2)
            + bx * RCPsiSquared.Core.Pauli.PauliMatrix.Of(RCPsiSquared.Core.Pauli.PauliLetter.X)
            + by * RCPsiSquared.Core.Pauli.PauliMatrix.Of(RCPsiSquared.Core.Pauli.PauliLetter.Y)
            + bz * RCPsiSquared.Core.Pauli.PauliMatrix.Of(RCPsiSquared.Core.Pauli.PauliLetter.Z)) / 2.0;

        var vec = PauliBasis.ToPauliVector(rho, N: 1);
        Assert.Equal(0.5, vec[(int)RCPsiSquared.Core.Pauli.PauliLetter.I].Real, 10);
        Assert.Equal(bx / 2.0, vec[(int)RCPsiSquared.Core.Pauli.PauliLetter.X].Real, 10);
        Assert.Equal(bz / 2.0, vec[(int)RCPsiSquared.Core.Pauli.PauliLetter.Z].Real, 10);
        Assert.Equal(by / 2.0, vec[(int)RCPsiSquared.Core.Pauli.PauliLetter.Y].Real, 10);
    }
}
