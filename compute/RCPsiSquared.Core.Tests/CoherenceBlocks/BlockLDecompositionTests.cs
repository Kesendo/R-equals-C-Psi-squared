using System.Numerics;
using RCPsiSquared.Core.CoherenceBlocks;

namespace RCPsiSquared.Core.Tests.CoherenceBlocks;

public class BlockLDecompositionTests
{
    [Fact]
    public void D_IsDiagonal_WithMinusTwoGammaTimesHD()
    {
        // c=2 N=5: HD ∈ {1, 3}. With γ₀=0.05, diagonal entries are ∈ {-0.1, -0.3}.
        var block = new CoherenceBlock(N: 5, n: 1, gammaZero: 0.05);
        var decomp = block.Decomposition;

        for (int i = 0; i < block.Basis.MTotal; i++)
        {
            for (int j = 0; j < block.Basis.MTotal; j++)
            {
                if (i == j) continue;
                Assert.Equal(0.0, decomp.D[i, j].Real, 12);
                Assert.Equal(0.0, decomp.D[i, j].Imaginary, 12);
            }
        }

        for (int i = 0; i < block.Basis.MTotal; i++)
        {
            double v = decomp.D[i, i].Real;
            Assert.True(Math.Abs(v - (-0.1)) < 1e-12 || Math.Abs(v - (-0.3)) < 1e-12,
                $"D[{i},{i}].Real = {v}, expected -0.1 or -0.3");
            Assert.Equal(0.0, decomp.D[i, i].Imaginary, 12);
        }
    }

    [Fact]
    public void NumBonds_EqualsNMinusOne()
    {
        var block = new CoherenceBlock(N: 7, n: 2, gammaZero: 0.05);
        Assert.Equal(6, block.Decomposition.NumBonds);
        Assert.All(block.Decomposition.MhPerBond, mh =>
        {
            Assert.Equal(block.Basis.MTotal, mh.RowCount);
            Assert.Equal(block.Basis.MTotal, mh.ColumnCount);
        });
    }

    [Fact]
    public void MhPerBond_OffDiagonalsArePureImaginaryUnits()
    {
        // Each per-bond M_H matrix has entries ±i on bond-flip transitions, 0 elsewhere.
        var block = new CoherenceBlock(N: 5, n: 1, gammaZero: 0.05);
        var decomp = block.Decomposition;

        foreach (var mh in decomp.MhPerBond)
        {
            for (int i = 0; i < block.Basis.MTotal; i++)
            {
                for (int j = 0; j < block.Basis.MTotal; j++)
                {
                    var z = mh[i, j];
                    Assert.Equal(0.0, z.Real, 12);
                    double im = z.Imaginary;
                    Assert.True(
                        Math.Abs(im) < 1e-12 ||
                        Math.Abs(im - 1.0) < 1e-12 ||
                        Math.Abs(im - (-1.0)) < 1e-12,
                        $"M_H[{i},{j}] imaginary part = {im}, expected 0, ±1");
                }
            }
        }
    }

    [Fact]
    public void AssembleUniform_AtJZero_EqualsD()
    {
        var block = new CoherenceBlock(N: 5, n: 1, gammaZero: 0.05);
        var L = block.Decomposition.AssembleUniform(0.0);
        for (int i = 0; i < block.Basis.MTotal; i++)
        {
            for (int j = 0; j < block.Basis.MTotal; j++)
            {
                Assert.Equal(block.Decomposition.D[i, j].Real, L[i, j].Real, 12);
                Assert.Equal(block.Decomposition.D[i, j].Imaginary, L[i, j].Imaginary, 12);
            }
        }
    }

    [Fact]
    public void AssembleUniform_MatchesAssembleAtUniformBondList()
    {
        var block = new CoherenceBlock(N: 5, n: 1, gammaZero: 0.05);
        var uniform = block.Decomposition.AssembleUniform(0.7);
        var perBond = block.Decomposition.AssembleAt(Enumerable.Repeat(0.7, block.Decomposition.NumBonds).ToArray());
        for (int i = 0; i < block.Basis.MTotal; i++)
        {
            for (int j = 0; j < block.Basis.MTotal; j++)
            {
                Assert.Equal(uniform[i, j].Real, perBond[i, j].Real, 12);
                Assert.Equal(uniform[i, j].Imaginary, perBond[i, j].Imaginary, 12);
            }
        }
    }

    [Fact]
    public void AssembleAt_RejectsWrongBondCount()
    {
        var block = new CoherenceBlock(N: 5, n: 1, gammaZero: 0.05);
        Assert.Throws<ArgumentException>(() => block.Decomposition.AssembleAt(new[] { 1.0 }));
    }
}
