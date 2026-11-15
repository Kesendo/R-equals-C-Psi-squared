using RCPsiSquared.Core.Numerics;
using Xunit;
using M = RCPsiSquared.Core.Numerics.GaussianRationalMatrix;

namespace RCPsiSquared.Core.Tests.Numerics;

/// <summary>The exact matrix helpers the quadrupolar generator is built from. The sandwich is the
/// one convention everything rests on: row-major vec, index d·i + j for |i⟩⟨j|.</summary>
public class GaussianRationalMatrixTests
{
    private static GaussianRational[,] G(long[,] re, long[,]? im = null)
    {
        int n = re.GetLength(0), p = re.GetLength(1);
        var m = M.Zero(n, p);
        for (int i = 0; i < n; i++)
            for (int j = 0; j < p; j++) m[i, j] = new GaussianRational(re[i, j], im is null ? 0 : im[i, j]);
        return m;
    }

    private static GaussianRational[] Vec(GaussianRational[,] x)
    {
        int d = x.GetLength(0);
        var v = new GaussianRational[d * d];
        for (int i = 0; i < d; i++)
            for (int j = 0; j < d; j++) v[d * i + j] = x[i, j];
        return v;
    }

    [Fact]
    public void Sandwich_ActsOnTheRowMajorVec_AsAThenRhoThenB()
    {
        var a = G(new long[,] { { 1, 2 }, { 3, 4 } }, new long[,] { { 0, 1 }, { -1, 0 } });
        var b = G(new long[,] { { 2, 0 }, { 1, -1 } }, new long[,] { { 1, 0 }, { 0, 2 } });
        var rho = G(new long[,] { { 5, -1 }, { 7, 3 } }, new long[,] { { 1, 2 }, { 0, -3 } });

        var direct = Vec(M.Multiply(M.Multiply(a, rho), b));
        var s = M.Sandwich(a, b);
        var viaSuper = new GaussianRational[4];
        var v = Vec(rho);
        for (int r = 0; r < 4; r++)
        {
            var acc = GaussianRational.Zero;
            for (int c = 0; c < 4; c++) acc += s[r, c] * v[c];
            viaSuper[r] = acc;
        }
        Assert.Equal(direct, viaSuper);
    }

    [Fact]
    public void Adjoint_IsTheConjugateTranspose_AndIdentityIsNeutral()
    {
        var a = G(new long[,] { { 1, 2 }, { 3, 4 } }, new long[,] { { 0, 1 }, { -1, 5 } });
        var adj = M.Adjoint(a);
        Assert.Equal(new GaussianRational(3, 1), adj[0, 1]);
        Assert.Equal(new GaussianRational(4, -5), adj[1, 1]);
        Assert.True(M.AreEqual(M.Multiply(a, M.Identity(2)), a));
        Assert.True(M.AreEqual(M.Subtract(M.Add(a, a), a), a));
        Assert.True(M.AreEqual(M.Transpose(M.Transpose(a)), a));
        Assert.False(M.AreEqual(a, M.Scale(a, 2)));
    }
}
