using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.F89PathK;

/// <summary>Per-site (0,1) coherence reduction from the n_block-site block density matrix.
///
/// <para>Given a block-reduced ρ on n_block sites,
/// <c>⟨0_l|Tr_{block\\{l}}(ρ)|1_l⟩</c> isolates the (|0⟩,|1⟩) matrix element on site l after
/// tracing out the other (n_block − 1) block sites. Used by F89 eigenmode-projection workflows
/// to extract per-site amplitude trajectories from block dynamics.</para>
///
/// <para><see cref="PerSiteReductionMatrix"/> packs the reduction as a (n_block × d²) sparse-but-dense
/// matrix acting on the column-major vec(ρ) of the block, so a single matvec yields the per-site
/// amplitude vector for arbitrary block ρ.</para>
/// </summary>
public static class F89BlockSiteReduction
{
    public static Complex ReduceBlockToSite01(ComplexMatrix rhoBlock, int l, int nBlock)
    {
        if (nBlock < 1)
            throw new ArgumentOutOfRangeException(nameof(nBlock), nBlock, "nBlock must be ≥ 1.");
        if (l < 0 || l >= nBlock)
            throw new ArgumentOutOfRangeException(nameof(l), l, $"l must be in [0, {nBlock - 1}].");
        int d = 1 << nBlock;
        if (rhoBlock.RowCount != d || rhoBlock.ColumnCount != d)
            throw new ArgumentException($"rhoBlock must be {d}x{d}; got {rhoBlock.RowCount}x{rhoBlock.ColumnCount}.");

        var bitPos = F89BlockInitialRho.BlockBitPos(nBlock);
        var other = OtherSites(l, nBlock);
        Complex val = Complex.Zero;
        int sweep = 1 << (nBlock - 1);
        for (int c = 0; c < sweep; c++)
        {
            int idx0 = 0;
            for (int i = 0; i < nBlock - 1; i++)
            {
                int bit = (c >> (nBlock - 2 - i)) & 1;
                idx0 += bitPos[other[i]] * bit;
            }
            int idx1 = idx0 + bitPos[l];
            val += rhoBlock[idx0, idx1];
        }
        return val;
    }

    /// <summary>Build the (n_block × d²) reduction matrix w[l, idx] such that
    /// <c>Σ_idx w[l, idx] · vec_C(ρ)[idx] = ReduceBlockToSite01(ρ, l, n_block)</c>, where
    /// <c>vec_C(ρ)[b·d + a] = ρ[a, b]</c> is the column-major vectorisation matching Python
    /// <c>simulations/_f89_pathk_lib.py</c>.</summary>
    public static ComplexMatrix PerSiteReductionMatrix(int nBlock)
    {
        if (nBlock < 1)
            throw new ArgumentOutOfRangeException(nameof(nBlock), nBlock, "nBlock must be ≥ 1.");
        int d = 1 << nBlock;
        var bitPos = F89BlockInitialRho.BlockBitPos(nBlock);
        var w = Matrix<Complex>.Build.Dense(nBlock, d * d);
        int sweep = 1 << (nBlock - 1);
        for (int l = 0; l < nBlock; l++)
        {
            var other = OtherSites(l, nBlock);
            for (int c = 0; c < sweep; c++)
            {
                int idx0 = 0;
                for (int i = 0; i < nBlock - 1; i++)
                {
                    int bit = (c >> (nBlock - 2 - i)) & 1;
                    idx0 += bitPos[other[i]] * bit;
                }
                int idx1 = idx0 + bitPos[l];
                // Column-major vec: vec(M)[b·d + a] = M[a, b]; here a = idx0, b = idx1.
                w[l, idx1 * d + idx0] = Complex.One;
            }
        }
        return w;
    }

    private static int[] OtherSites(int l, int nBlock)
    {
        var other = new int[nBlock - 1];
        int k = 0;
        for (int s = 0; s < nBlock; s++) if (s != l) other[k++] = s;
        return other;
    }
}
