using RCPsiSquared.Core.Numerics;
using Xunit;

namespace RCPsiSquared.Core.Tests.Numerics;

/// <summary>Division-free Samuelson-Berkowitz characteristic polynomial over Z[i] — the engine
/// that turns the path-k block into its exact charpoly. Known-answer gates (coefficients
/// lowest-first, monic).</summary>
public class GaussianMatrixCharpolyTests
{
    [Fact]
    public void Characteristic_RealTwoByTwo_IsLambdaSquaredMinusTraceLambdaPlusDet()
    {
        // [[1,2],[3,4]]: det(λI − A) = λ² − 5λ − 2.
        var a = new GaussianInteger[,] { { 1, 2 }, { 3, 4 } };
        Assert.Equal(new GaussianInteger[] { -2, -5, 1 }, GaussianMatrixCharpoly.Characteristic(a));
    }

    [Fact]
    public void Characteristic_GaussianUpperTriangular_IsProductOfDiagonalRoots()
    {
        // [[i,1],[0,2]]: (λ − i)(λ − 2) = λ² − (2 + i)λ + 2i.  (the i in row 0 exercises Z[i] arithmetic)
        var a = new GaussianInteger[,] { { GaussianInteger.I, 1 }, { 0, 2 } };
        Assert.Equal(new GaussianInteger[] { new(0, 2), new(-2, -1), 1 },
                     GaussianMatrixCharpoly.Characteristic(a));
    }

    [Fact]
    public void Characteristic_TridiagonalThreeByThree_MatchesDirectDeterminant()
    {
        // [[2,1,0],[1,2,1],[0,1,2]]: (λ−2)(λ²−4λ+2) = λ³ − 6λ² + 10λ − 4.
        // Nonzero below-diagonal + above-diagonal entries exercise the R·Msub^k·S Toeplitz terms.
        var a = new GaussianInteger[,] { { 2, 1, 0 }, { 1, 2, 1 }, { 0, 1, 2 } };
        Assert.Equal(new GaussianInteger[] { -4, 10, -6, 1 }, GaussianMatrixCharpoly.Characteristic(a));
    }
}
