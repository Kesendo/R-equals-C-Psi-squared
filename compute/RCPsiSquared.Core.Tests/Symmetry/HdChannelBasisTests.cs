using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class HdChannelBasisTests
{
    [Theory]
    [InlineData(5, 1)]   // c=2
    [InlineData(7, 2)]   // c=3
    [InlineData(7, 3)]   // c=4
    [InlineData(8, 1)]   // c=2
    public void P_Columns_AreOrthonormal(int N, int n)
    {
        var block = new CoherenceBlock(N, n, gammaZero: 0.05);
        var basis = HdChannelBasis.Build(block);

        var gram = basis.P.ConjugateTranspose() * basis.P;
        var identity = Matrix<Complex>.Build.DenseIdentity(basis.C);
        for (int i = 0; i < basis.C; i++)
            for (int j = 0; j < basis.C; j++)
                Assert.Equal(identity[i, j].Real, gram[i, j].Real, 12);
    }

    [Theory]
    [InlineData(5, 1, new[] { 1, 3 })]
    [InlineData(7, 3, new[] { 1, 3, 5, 7 })]
    [InlineData(8, 4, new[] { 1, 3, 5, 7 })]
    public void HammingDistances_AreCorrect(int N, int n, int[] expected)
    {
        var block = new CoherenceBlock(N, n, gammaZero: 0.05);
        var basis = HdChannelBasis.Build(block);
        Assert.Equal(expected, basis.HammingDistances);
    }

    [Theory]
    [InlineData(5, 1)]
    [InlineData(6, 2)]
    [InlineData(7, 2)]
    [InlineData(7, 3)]
    [InlineData(8, 1)]
    public void M_H_Total_IsDiagonal_InChannelUniformBasis(int N, int n)
    {
        // Empirical structural fact: P^† · M_H_total · P is diagonal across all tested cases.
        var block = new CoherenceBlock(N, n, gammaZero: 0.05);
        var hd = HdChannelBasis.Build(block);

        var mhTotal = block.Decomposition.MhPerBond[0].Clone();
        for (int b = 1; b < block.Decomposition.NumBonds; b++)
            mhTotal = mhTotal + block.Decomposition.MhPerBond[b];

        var projected = hd.P.ConjugateTranspose() * mhTotal * hd.P;

        for (int i = 0; i < hd.C; i++)
        {
            for (int j = 0; j < hd.C; j++)
            {
                if (i == j) continue;
                double mag = projected[i, j].Magnitude;
                Assert.True(mag < 1e-10,
                    $"P†·M_H·P[{i},{j}] should be 0 in channel-uniform basis; got |z|={mag}");
            }
        }
    }
}
