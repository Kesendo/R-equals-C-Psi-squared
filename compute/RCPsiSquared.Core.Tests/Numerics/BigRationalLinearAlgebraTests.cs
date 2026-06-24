using System.Linq;
using RCPsiSquared.Core.Numerics;
using Xunit;

namespace RCPsiSquared.Core.Tests.Numerics;

/// <summary>Exact BigRational linear algebra — the foundation for the F89 full-D rate-confined
/// invariant-subspace AT reconstruction (nullspace + restriction-solve).</summary>
public class BigRationalLinearAlgebraTests
{
    private static BigRational[,] R(long[,] m)
    {
        var a = new BigRational[m.GetLength(0), m.GetLength(1)];
        for (int i = 0; i < m.GetLength(0); i++)
            for (int j = 0; j < m.GetLength(1); j++) a[i, j] = m[i, j];
        return a;
    }

    [Fact]
    public void Multiply_TwoByTwo()
    {
        var p = BigRationalLinearAlgebra.Multiply(R(new long[,] { { 1, 2 }, { 3, 4 } }),
                                                  R(new long[,] { { 5, 6 }, { 7, 8 } }));
        Assert.Equal((BigRational)19, p[0, 0]);
        Assert.Equal((BigRational)22, p[0, 1]);
        Assert.Equal((BigRational)43, p[1, 0]);
        Assert.Equal((BigRational)50, p[1, 1]);
    }

    [Fact]
    public void Transpose_TwoByThree()
    {
        var t = BigRationalLinearAlgebra.Transpose(R(new long[,] { { 1, 2, 3 }, { 4, 5, 6 } }));
        Assert.Equal(3, t.GetLength(0));
        Assert.Equal(2, t.GetLength(1));
        Assert.Equal((BigRational)2, t[1, 0]);
        Assert.Equal((BigRational)6, t[2, 1]);
    }

    [Fact]
    public void Nullspace_RankDeficient_HasOneVectorInTheKernel()
    {
        // [[1,2,3],[2,4,6]] has rank 1 ⟹ 2-dim nullspace; every basis vector v satisfies A·v = 0.
        var a = R(new long[,] { { 1, 2, 3 }, { 2, 4, 6 } });
        var ns = BigRationalLinearAlgebra.Nullspace(a);
        Assert.Equal(2, ns.Count);
        foreach (var v in ns)
        {
            Assert.False(v.All(x => x.IsZero));                 // nonzero
            for (int row = 0; row < 2; row++)
            {
                BigRational dot = BigRational.Zero;
                for (int j = 0; j < 3; j++) dot += a[row, j] * v[j];
                Assert.True(dot.IsZero);                        // A·v = 0
            }
        }
    }

    [Fact]
    public void Nullspace_FullColumnRank_IsEmpty()
    {
        Assert.Empty(BigRationalLinearAlgebra.Nullspace(R(new long[,] { { 1, 0 }, { 0, 1 } })));
    }

    [Fact]
    public void Solve_RecoversTheKnownSolution()
    {
        // [[2,1],[1,3]] · x = [[5],[10]] ⟹ x = [[1],[3]].
        var x = BigRationalLinearAlgebra.Solve(R(new long[,] { { 2, 1 }, { 1, 3 } }),
                                               R(new long[,] { { 5 }, { 10 } }));
        Assert.Equal((BigRational)1, x[0, 0]);
        Assert.Equal((BigRational)3, x[1, 0]);
    }
}
