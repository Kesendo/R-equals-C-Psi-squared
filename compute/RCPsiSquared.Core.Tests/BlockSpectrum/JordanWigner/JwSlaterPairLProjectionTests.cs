using RCPsiSquared.Core.BlockSpectrum.JordanWigner;
using RCPsiSquared.Core.Knowledge;
using Xunit;
using Xunit.Abstractions;

namespace RCPsiSquared.Core.Tests.BlockSpectrum.JordanWigner;

/// <summary>Witnesses for <see cref="JwSlaterPairLProjection"/>: the full L-block of the
/// chain XY + Z-dephasing Liouvillian transformed to the Slater-pair basis is sparse with
/// off-diagonal support bounded by (1 + p_r(N−p_r))·(1 + p_c(N−p_c)) − 1 — the Slater-pair
/// reach of a single Z_l ⊗ Z_l action. This is the structural foundation of the N=10
/// sparse-eig path (Phase 2 of the N=10 push).</summary>
public class JwSlaterPairLProjectionTests
{
    private readonly ITestOutputHelper _out;
    public JwSlaterPairLProjectionTests(ITestOutputHelper output) => _out = output;

    [Theory]
    [InlineData(4, 1, 2)]  // c=2 F86 block, dim 24, nnz-bound 19
    [InlineData(5, 1, 2)]  // c=2 N=5, dim 50, nnz-bound 27
    [InlineData(5, 2, 2)]  // dim 100, nnz-bound 48
    [InlineData(5, 2, 3)]  // dim 100, nnz-bound 48
    [InlineData(6, 2, 3)]  // dim 300, nnz-bound 90 (sym)
    [InlineData(6, 3, 3)]  // dim 400, nnz-bound 99
    public void Build_UniformGamma_OffDiagonalNnzWithinTheoreticalBound(int N, int pCol, int pRow)
    {
        var basis = JwSlaterPairBasis.Build(N, pCol, pRow);
        var gamma = Enumerable.Repeat(0.1, N).ToArray();
        var proj = JwSlaterPairLProjection.Build(basis, gamma);

        int bound = (1 + pRow * (N - pRow)) * (1 + pCol * (N - pCol)) - 1;
        _out.WriteLine($"N={N}, (p_c,p_r)=({pCol},{pRow}), dim={basis.U.RowCount}, " +
                       $"theoretical-nnz-bound/row={bound}, " +
                       $"empirical-max-nnz/row={proj.MaxOffDiagonalNonzeroPerRow}, " +
                       $"empirical-mean-nnz/row={proj.MeanOffDiagonalNonzeroPerRow:F2}");
        Assert.Equal(bound, proj.TheoreticalNnzBoundPerRow);
        Assert.True(proj.MaxOffDiagonalNonzeroPerRow <= proj.TheoreticalNnzBoundPerRow,
            $"Empirical max-nnz/row {proj.MaxOffDiagonalNonzeroPerRow} exceeds theoretical bound {bound}");
        Assert.True(proj.RespectsTheoreticalBound);
    }

    [Theory]
    [InlineData(4, 1, 2)]
    [InlineData(5, 2, 2)]
    [InlineData(6, 3, 3)]
    public void Build_FrobeniusSquared_PreservedUnderUnitaryTransform(int N, int pCol, int pRow)
    {
        // Unitary similarity preserves Frobenius² (basis change does not affect operator norm).
        var basis = JwSlaterPairBasis.Build(N, pCol, pRow);
        var gamma = Enumerable.Repeat(0.05, N).ToArray();
        var proj = JwSlaterPairLProjection.Build(basis, gamma);

        double lJwFrob = proj.LJw.FrobeniusNorm();
        double lBlockFrob = proj.LBlockComputational.FrobeniusNorm();
        Assert.Equal(lBlockFrob, lJwFrob, precision: 8);
    }

    [Theory]
    [InlineData(4, 1, 2)]
    [InlineData(5, 2, 3)]
    public void Build_ZeroGamma_LJwIsDiagonal(int N, int pCol, int pRow)
    {
        // At γ=0, L = L_H, which is diagonal in JW basis (already verified by basis witness).
        // L_JW off-diagonal must be all zero.
        var basis = JwSlaterPairBasis.Build(N, pCol, pRow);
        var zeroGamma = new double[N];
        var proj = JwSlaterPairLProjection.Build(basis, zeroGamma);

        Assert.Equal(0, proj.MaxOffDiagonalNonzeroPerRow);
        Assert.Equal(0.0, proj.MeanOffDiagonalNonzeroPerRow, precision: 12);
    }

    [Theory]
    [InlineData(5, 2, 2)]
    [InlineData(6, 2, 3)]
    public void Build_UniformGamma_DensityImprovesWithN(int N, int pCol, int pRow)
    {
        // Sparsity ratio (empirical-max-nnz / sectorDim) should be strictly < 1 — i.e. the
        // dissipator is genuinely sparse, not a full-rank perturbation.
        var basis = JwSlaterPairBasis.Build(N, pCol, pRow);
        var gamma = Enumerable.Repeat(0.1, N).ToArray();
        var proj = JwSlaterPairLProjection.Build(basis, gamma);

        int sectorDim = basis.U.RowCount;
        double density = (double)proj.MaxOffDiagonalNonzeroPerRow / sectorDim;
        Assert.True(density < 1.0,
            $"Density {density:P1} should be < 100% (sparse, not full)");
    }

    [Fact]
    public void Build_Tier_IsTier1Derived()
    {
        var basis = JwSlaterPairBasis.Build(N: 4, pCol: 1, pRow: 2);
        var gamma = Enumerable.Repeat(0.1, 4).ToArray();
        var proj = JwSlaterPairLProjection.Build(basis, gamma);
        Assert.Equal(Tier.Tier1Derived, proj.Tier);
    }

    [Fact]
    public void Build_RejectsGammaWrongLength()
    {
        var basis = JwSlaterPairBasis.Build(N: 4, pCol: 1, pRow: 2);
        var wrongGamma = new double[3];
        Assert.Throws<ArgumentException>(() => JwSlaterPairLProjection.Build(basis, wrongGamma));
    }

    [Fact]
    public void Reconnaissance_PrintsSparsityProgression()
    {
        // Reconnaissance: how does the sparsity ratio scale across N, p_c, p_r? Confirms the
        // claim that density improves with N (N=10 expected ~1%, N=6 ~25%).
        _out.WriteLine("  N | p_c | p_r | sectorDim | bound | empirical | density");
        _out.WriteLine("  --|-----|-----|-----------|-------|-----------|--------");
        foreach (var (N, pc, pr) in new[]
        {
            (4, 1, 2), (5, 1, 2), (5, 2, 2), (5, 2, 3),
            (6, 2, 3), (6, 3, 3), (7, 3, 4)
        })
        {
            var basis = JwSlaterPairBasis.Build(N, pc, pr);
            var gamma = Enumerable.Repeat(0.1, N).ToArray();
            var proj = JwSlaterPairLProjection.Build(basis, gamma);
            int dim = basis.U.RowCount;
            double density = (double)proj.MaxOffDiagonalNonzeroPerRow / dim;
            _out.WriteLine($"  {N} | {pc,3} | {pr,3} | {dim,9} | {proj.TheoreticalNnzBoundPerRow,5} | " +
                           $"{proj.MaxOffDiagonalNonzeroPerRow,9} | {density,7:P2}");
        }
    }
}
