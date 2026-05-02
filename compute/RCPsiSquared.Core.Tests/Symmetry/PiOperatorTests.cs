using System.Numerics;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class PiOperatorTests
{
    [Theory]
    [InlineData(PauliLetter.I, PauliLetter.X, 1.0, 0.0)]   // Z-deph: I → X, phase 1
    [InlineData(PauliLetter.X, PauliLetter.I, 1.0, 0.0)]   // X → I, phase 1
    [InlineData(PauliLetter.Y, PauliLetter.Z, 0.0, 1.0)]   // Y → Z, phase i
    [InlineData(PauliLetter.Z, PauliLetter.Y, 0.0, 1.0)]   // Z → Y, phase i
    public void ActOnLetter_ZDephasing_FollowsKleinFlip(PauliLetter input, PauliLetter expected, double phaseRe, double phaseIm)
    {
        var (newLetter, phase) = PiOperator.ActOnLetter(input, PauliLetter.Z);
        Assert.Equal(expected, newLetter);
        Assert.Equal(phaseRe, phase.Real, 12);
        Assert.Equal(phaseIm, phase.Imaginary, 12);
    }

    [Fact]
    public void Pi_IsOrderFour_ZDephasing()
    {
        // Π² is signed-permutation diagonal in the Pauli basis; Π⁴ = I.
        var pi = PiOperator.BuildFull(N: 2);
        var pi4 = pi * pi * pi * pi;
        long d2 = 1L << 4;
        var diff = pi4 - MathNet.Numerics.LinearAlgebra.Matrix<Complex>.Build.SparseIdentity((int)d2);
        Assert.True(diff.FrobeniusNorm() < 1e-10, $"Π⁴ should equal I; got Frobenius {diff.FrobeniusNorm()}");
    }

    [Theory]
    [InlineData("XX", PauliLetter.Z, +1)] // bit_b sum = 0+0 = 0 even → +1
    [InlineData("XY", PauliLetter.Z, -1)] // bit_b sum = 0+1 = 1 odd → -1
    [InlineData("YY", PauliLetter.Z, +1)] // bit_b sum = 1+1 = 2 even → +1
    [InlineData("ZZ", PauliLetter.Z, +1)] // bit_b sum = 1+1 = 2 even → +1
    [InlineData("YZ", PauliLetter.Z, +1)] // bit_b sum = 1+1 = 2 even → +1
    public void Pi2_Eigenvalue_MatchesBitBParity_ForZDephasing(string label, PauliLetter dephase, int expected)
    {
        var letters = PauliLabel.Parse(label);
        Assert.Equal(expected, PiOperator.SquaredEigenvalue(letters, dephase));
    }
}
