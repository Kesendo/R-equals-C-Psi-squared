using System;
using System.Numerics;
using RCPsiSquared.Core.F89PathK;
using Xunit;

namespace RCPsiSquared.Core.Tests.F89PathK;

/// <summary>
/// The CSR emitter must be entry-for-entry the same operator as the dense builders it
/// transcribes: WeightCoherenceBlock.Build (full block) and BuildReflectionSectorColumnMajor
/// (R-parity sector). Everything downstream (matvec, BiCGStab, witness) rests on this pin.
/// </summary>
public class WeightCoherenceSectorCsrTests
{
    private static Complex[,] DenseFromCsr(WeightCoherenceSectorCsr.Csr m)
    {
        var d = new Complex[m.Dim, m.Dim];
        for (int r = 0; r < m.Dim; r++)
            for (int k = m.RowPtr[r]; k < m.RowPtr[r + 1]; k++)
                d[r, m.ColIdx[k]] += m.Values[k];
        return d;
    }

    [Theory]
    [InlineData(5, 2, 3, 0.7, 0.0)]
    [InlineData(5, 1, 2, 0.620878, 0.0)]
    [InlineData(6, 2, 4, 1.3, 0.25)]   // complex q exercised too
    public void BuildFull_MatchesDenseBuilder(int n, int wKet, int wBra, double qRe, double qIm)
    {
        var q = new Complex(qRe, qIm);
        var dense = WeightCoherenceBlock.Build(n, wKet, wBra, q);
        var csr = WeightCoherenceSectorCsr.BuildFull(n, wKet, wBra, q);
        var back = DenseFromCsr(csr);
        Assert.Equal(dense.GetLength(0), csr.Dim);
        double worst = 0.0;
        for (int r = 0; r < csr.Dim; r++)
            for (int c = 0; c < csr.Dim; c++)
                worst = Math.Max(worst, (dense[r, c] - back[r, c]).Magnitude);
        Assert.True(worst <= 1e-15, $"full block ({wKet},{wBra}) N={n}: worst entry diff {worst}");
    }

    [Theory]
    [InlineData(5, 1, 2, true)]
    [InlineData(5, 1, 2, false)]
    [InlineData(5, 2, 3, true)]
    [InlineData(5, 2, 3, false)]
    public void BuildReflectionSector_MatchesDenseSectorBuilder(int n, int wKet, int wBra, bool odd)
    {
        var q = new Complex(0.620878, 0.0);
        var (aFlat, d) = WeightCoherenceBlock.BuildReflectionSectorColumnMajor(n, wKet, wBra, q, odd);
        var csr = WeightCoherenceSectorCsr.BuildReflectionSector(n, wKet, wBra, q, odd);
        Assert.Equal(d, csr.Dim);
        var back = DenseFromCsr(csr);
        double worst = 0.0;
        for (int c = 0; c < d; c++)
            for (int r = 0; r < d; r++)
                worst = Math.Max(worst, (aFlat[c * d + r] - back[r, c]).Magnitude);
        Assert.True(worst <= 1e-15, $"sector ({wKet},{wBra}) odd={odd}: worst entry diff {worst}");
    }

    [Fact]
    public void BuildFull_RowPtrIsWellFormedAndSorted()
    {
        var csr = WeightCoherenceSectorCsr.BuildFull(5, 2, 3, new Complex(0.7, 0.0));
        Assert.Equal(0, csr.RowPtr[0]);
        Assert.Equal(csr.Values.Length, csr.RowPtr[csr.Dim]);
        for (int r = 0; r < csr.Dim; r++)
            for (int k = csr.RowPtr[r] + 1; k < csr.RowPtr[r + 1]; k++)
                Assert.True(csr.ColIdx[k - 1] < csr.ColIdx[k], "column indices sorted, no duplicates");
    }
}
