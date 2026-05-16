using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.BlockSpectrum.JordanWigner;
using RCPsiSquared.Core.Knowledge;
using Xunit;
using Xunit.Abstractions;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Tests.BlockSpectrum.JordanWigner;

/// <summary>Witnesses for <see cref="JwSlaterPairSparseLBuilder"/>: builds the L = L_H + L_D
/// Liouvillian block directly in the JW Slater-pair basis as a sparse CSR matrix, without
/// ever materialising the dense computational-basis L or running the dense U^†·L·U product.
/// Validated against <see cref="JwSlaterPairLProjection"/> (the dense reference path verified
/// in Phase 2 step 1) element-by-element at small N.</summary>
public class JwSlaterPairSparseLBuilderTests
{
    private readonly ITestOutputHelper _out;
    public JwSlaterPairSparseLBuilderTests(ITestOutputHelper output) => _out = output;

    [Theory]
    [InlineData(4, 1, 2)]
    [InlineData(5, 1, 2)]
    [InlineData(5, 2, 2)]
    [InlineData(5, 2, 3)]
    [InlineData(6, 2, 3)]
    [InlineData(6, 3, 3)]
    [InlineData(7, 3, 4)]
    public void Build_UniformGamma_MatchesDenseProjection(int N, int pCol, int pRow)
    {
        // The sparse build must agree element-by-element with the dense U^†·L·U reference.
        // Tolerance 1e-10 matches the construction tolerance of JwSlaterPairBasis.
        var gamma = Enumerable.Repeat(0.1, N).ToArray();

        var sparse = JwSlaterPairSparseLBuilder.Build(N, pCol, pRow, gamma);
        var sparseDense = sparse.ToDense();

        var basis = JwSlaterPairBasis.Build(N, pCol, pRow);
        var refProj = JwSlaterPairLProjection.Build(basis, gamma);

        double diff = (refProj.LJw - sparseDense).FrobeniusNorm();
        double refNorm = refProj.LJw.FrobeniusNorm();
        _out.WriteLine($"N={N}, (p_c,p_r)=({pCol},{pRow}), dim={sparse.SectorDim}, " +
                       $"nnz={sparse.NnzTotal}, max-nnz/row={sparse.MaxNnzPerRow}, " +
                       $"‖sparse − dense‖_F / ‖dense‖_F = {diff / refNorm:G3}");
        Assert.True(diff / refNorm < 1e-10,
            $"Sparse/dense mismatch: ‖diff‖_F = {diff:G3}, ‖ref‖_F = {refNorm:G3}");
    }

    [Theory]
    [InlineData(4, 1, 2)]
    [InlineData(5, 2, 3)]
    [InlineData(6, 3, 3)]
    public void Build_NonUniformGamma_MatchesDenseProjection(int N, int pCol, int pRow)
    {
        // Non-uniform γ stresses the per-site dissipator weighting in the sparse build.
        var gamma = Enumerable.Range(0, N).Select(i => 0.05 + 0.03 * i).ToArray();

        var sparse = JwSlaterPairSparseLBuilder.Build(N, pCol, pRow, gamma);
        var sparseDense = sparse.ToDense();

        var basis = JwSlaterPairBasis.Build(N, pCol, pRow);
        var refProj = JwSlaterPairLProjection.Build(basis, gamma);

        double diff = (refProj.LJw - sparseDense).FrobeniusNorm();
        Assert.True(diff < 1e-10, $"Non-uniform γ sparse/dense mismatch: ‖diff‖_F = {diff:G3}");
    }

    [Theory]
    [InlineData(5, 2, 2)]
    [InlineData(6, 3, 3)]
    public void Build_NnzWithinTheoreticalBound(int N, int pCol, int pRow)
    {
        // Same Slater-swap-reach bound as JwSlaterPairLProjection.TheoreticalNnzBoundPerRow:
        // max-nnz/row (including diagonal) ≤ (1 + p_r(N−p_r)) · (1 + p_c(N−p_c)).
        var gamma = Enumerable.Repeat(0.1, N).ToArray();
        var sparse = JwSlaterPairSparseLBuilder.Build(N, pCol, pRow, gamma);

        int boundIncludingDiag = (1 + pRow * (N - pRow)) * (1 + pCol * (N - pCol));
        Assert.True(sparse.MaxNnzPerRow <= boundIncludingDiag,
            $"max-nnz/row {sparse.MaxNnzPerRow} exceeds bound {boundIncludingDiag}");
    }

    [Theory]
    [InlineData(4, 1, 2)]
    [InlineData(5, 2, 3)]
    public void Build_ZeroGamma_OnlyHamiltonianDiagonal(int N, int pCol, int pRow)
    {
        // At γ=0, L = L_H is purely diagonal in JW basis; the sparse build emits at most
        // one nnz per row (the Hamiltonian diagonal). Some rows can be exactly zero when the
        // dispersion symmetry ε_L = ε_K holds (e.g. N=5 has ε_3=0 so any L,K containing
        // mode 3 may produce a zero diagonal), so we assert MaxNnzPerRow ≤ 1, not ==.
        var zeroGamma = new double[N];
        var sparse = JwSlaterPairSparseLBuilder.Build(N, pCol, pRow, zeroGamma);

        Assert.True(sparse.MaxNnzPerRow <= 1,
            $"At γ=0 each row should have ≤ 1 nnz (diagonal only); got max={sparse.MaxNnzPerRow}.");
        // Every nnz must sit on the diagonal (column == row).
        for (int alpha = 0; alpha < sparse.SectorDim; alpha++)
            for (int e = sparse.RowPtr[alpha]; e < sparse.RowPtr[alpha + 1]; e++)
                Assert.Equal(alpha, sparse.ColIdx[e]);
    }

    [Fact]
    public void Build_Tier_IsTier1Derived()
    {
        var gamma = Enumerable.Repeat(0.1, 4).ToArray();
        var sparse = JwSlaterPairSparseLBuilder.Build(N: 4, pCol: 1, pRow: 2, gamma);
        Assert.Equal(Tier.Tier1Derived, sparse.Tier);
    }

    [Fact]
    public void Build_RejectsGammaWrongLength()
    {
        var wrongGamma = new double[3];
        Assert.Throws<ArgumentException>(() =>
            JwSlaterPairSparseLBuilder.Build(N: 4, pCol: 1, pRow: 2, wrongGamma));
    }

    [Fact]
    public void Build_AtN10HalfFilling_ConstructsWithinSparsityBudget()
    {
        // The N=10 push target: (p_c=5, p_r=5) max-block, sectorDim 63504, expected
        // nnz/row ≤ 676 (= 26·26 including diagonal). Verifies the sparse build
        // succeeds at the dim that exceeds JwSlaterPairBasis.MaxSectorDimForDenseWitness.
        var gamma = Enumerable.Repeat(0.05, 10).ToArray();
        var sparse = JwSlaterPairSparseLBuilder.Build(N: 10, pCol: 5, pRow: 5, gamma);

        Assert.Equal(63504, sparse.SectorDim);
        int boundIncludingDiag = (1 + 5 * 5) * (1 + 5 * 5);  // 676
        Assert.True(sparse.MaxNnzPerRow <= boundIncludingDiag,
            $"N=10 max-nnz/row {sparse.MaxNnzPerRow} exceeds bound {boundIncludingDiag}");

        _out.WriteLine($"N=10 (5,5): sectorDim={sparse.SectorDim}, nnz={sparse.NnzTotal:N0}, " +
                       $"max-nnz/row={sparse.MaxNnzPerRow}, mean-nnz/row={sparse.MeanNnzPerRow:F1}, " +
                       $"density={(double)sparse.MaxNnzPerRow / sparse.SectorDim:P2}");
    }
}
