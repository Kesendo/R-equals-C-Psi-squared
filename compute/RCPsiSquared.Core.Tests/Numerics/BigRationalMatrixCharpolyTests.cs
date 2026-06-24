using RCPsiSquared.Core.Numerics;
using Xunit;

namespace RCPsiSquared.Core.Tests.Numerics;

/// <summary>Berkowitz characteristic polynomial over BigRational (lowest-first, monic).</summary>
public class BigRationalMatrixCharpolyTests
{
    [Fact]
    public void Characteristic_TridiagonalInteger_MatchesDirectDeterminant()
    {
        // [[2,1,0],[1,2,1],[0,1,2]]: det(λI − A) = λ³ − 6λ² + 10λ − 4.
        var a = new BigRational[,]
        {
            { 2, 1, 0 },
            { 1, 2, 1 },
            { 0, 1, 2 },
        };
        Assert.Equal(new BigRational[] { -4, 10, -6, 1 }, BigRationalMatrixCharpoly.Characteristic(a));
    }

    [Fact]
    public void Characteristic_DiagonalRational_IsProductOfRootFactors()
    {
        // diag(1/2, 1/3): (λ − 1/2)(λ − 1/3) = λ² − 5/6 λ + 1/6.
        var a = new BigRational[,]
        {
            { new BigRational(1, 2), BigRational.Zero },
            { BigRational.Zero, new BigRational(1, 3) },
        };
        Assert.Equal(
            new[] { new BigRational(1, 6), new BigRational(-5, 6), BigRational.One },
            BigRationalMatrixCharpoly.Characteristic(a));
    }
}
