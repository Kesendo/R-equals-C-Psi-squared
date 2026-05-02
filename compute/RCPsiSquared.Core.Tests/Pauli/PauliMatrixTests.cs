using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.Pauli;

namespace RCPsiSquared.Core.Tests.Pauli;

public class PauliMatrixTests
{
    [Theory]
    [InlineData(PauliLetter.X)]
    [InlineData(PauliLetter.Y)]
    [InlineData(PauliLetter.Z)]
    public void SigmaSquared_EqualsIdentity(PauliLetter letter)
    {
        var sigma = PauliMatrix.Of(letter);
        var sigmaSq = sigma * sigma;
        var i = Matrix<Complex>.Build.DenseIdentity(2);
        AssertMatrixApproxEqual(i, sigmaSq);
    }

    [Fact]
    public void SigmaX_Times_SigmaY_Equals_iSigmaZ()
    {
        var product = PauliMatrix.Of(PauliLetter.X) * PauliMatrix.Of(PauliLetter.Y);
        var iZ = new Complex(0, 1) * PauliMatrix.Of(PauliLetter.Z);
        AssertMatrixApproxEqual(iZ, product);
    }

    [Fact]
    public void SigmaY_Times_SigmaZ_Equals_iSigmaX()
    {
        var product = PauliMatrix.Of(PauliLetter.Y) * PauliMatrix.Of(PauliLetter.Z);
        var iX = new Complex(0, 1) * PauliMatrix.Of(PauliLetter.X);
        AssertMatrixApproxEqual(iX, product);
    }

    [Fact]
    public void SigmaZ_Times_SigmaX_Equals_iSigmaY()
    {
        var product = PauliMatrix.Of(PauliLetter.Z) * PauliMatrix.Of(PauliLetter.X);
        var iY = new Complex(0, 1) * PauliMatrix.Of(PauliLetter.Y);
        AssertMatrixApproxEqual(iY, product);
    }

    [Theory]
    [InlineData(PauliLetter.X)]
    [InlineData(PauliLetter.Y)]
    [InlineData(PauliLetter.Z)]
    public void SigmaTrace_IsZero(PauliLetter letter)
    {
        var sigma = PauliMatrix.Of(letter);
        Assert.Equal(0.0, sigma.Trace().Real, 12);
        Assert.Equal(0.0, sigma.Trace().Imaginary, 12);
    }

    [Theory]
    [InlineData(PauliLetter.X, PauliLetter.Y)]
    [InlineData(PauliLetter.X, PauliLetter.Z)]
    [InlineData(PauliLetter.Y, PauliLetter.Z)]
    public void SigmaSigma_TraceIsZero_ForDistinctNonIdentity(PauliLetter a, PauliLetter b)
    {
        var product = PauliMatrix.Of(a) * PauliMatrix.Of(b);
        Assert.Equal(0.0, product.Trace().Real, 12);
        Assert.Equal(0.0, product.Trace().Imaginary, 12);
    }

    [Theory]
    [InlineData(PauliLetter.X)]
    [InlineData(PauliLetter.Y)]
    [InlineData(PauliLetter.Z)]
    public void TrSigmaSquared_EqualsTwo(PauliLetter letter)
    {
        var sigma = PauliMatrix.Of(letter);
        var trace = (sigma * sigma).Trace();
        Assert.Equal(2.0, trace.Real, 12);
    }

    private static void AssertMatrixApproxEqual(Matrix<Complex> expected, Matrix<Complex> actual)
    {
        Assert.Equal(expected.RowCount, actual.RowCount);
        Assert.Equal(expected.ColumnCount, actual.ColumnCount);
        for (int i = 0; i < expected.RowCount; i++)
            for (int j = 0; j < expected.ColumnCount; j++)
            {
                Assert.Equal(expected[i, j].Real, actual[i, j].Real, 12);
                Assert.Equal(expected[i, j].Imaginary, actual[i, j].Imaginary, 12);
            }
    }
}
