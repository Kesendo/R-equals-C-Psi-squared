namespace RCPsiSquared.Core.BlockSpectrum;

/// <summary>Builds the basis permutation that block-diagonalises the XY+Z-dephasing
/// Liouvillian by the joint label (popcount_col, popcount_row) ∈ {0..N}².
///
/// <para>Convention matches <see cref="Lindblad.LindbladianBuilder"/> and
/// <see cref="Lindblad.PauliDephasingDissipator"/>: the Liouville-space flat index is
/// <c>flat = row · d + col</c> (row-major vec(ρ)), where d = 2^N. With this convention
/// <c>L = −i·(H ⊗ I − I ⊗ H^T) + Σ_l γ_l·(P_l ⊗ P_l* − I ⊗ I)</c> correctly represents
/// <c>L·vec(ρ) = vec(−i[H,ρ] + dissipator)</c>. The popcount labels then partition the
/// 4^N flat indices into (N+1)² disjoint sectors.</para>
///
/// <para>Returned permutation is the index ordering that brings L into block-diagonal
/// form: <c>L_perm[i, j] = L[perm[i], perm[j]]</c> is block-diagonal, with block (p_c, p_r)
/// occupying contiguous rows/cols starting at <see cref="SectorRange.Offset"/> with size
/// <see cref="SectorRange.Size"/>.</para>
///
/// <para>Supports N ∈ [1, 12]: at N=12 Liouville space is 16M flat indices (~64 MB for an
/// int permutation), already past the typical block-decomposition use case.</para></summary>
public static class JointPopcountSectorBuilder
{
    public sealed record SectorRange(int PCol, int PRow, int Offset, int Size);

    public sealed class Decomposition
    {
        public int N { get; }
        public int D { get; }   // 2^N
        public int[] Permutation { get; }
        public IReadOnlyList<SectorRange> SectorRanges { get; }
        public Decomposition(int n, int d, int[] perm, IReadOnlyList<SectorRange> sectors)
        { N = n; D = d; Permutation = perm; SectorRanges = sectors; }
    }

    public static Decomposition Build(int N)
    {
        if (N < 1 || N > 12) throw new ArgumentOutOfRangeException(nameof(N), N, "Supported N range: 1..12.");
        int d = 1 << N;
        int liouvilleDim = d * d;

        // Label each flat index by (popcount_col, popcount_row). Convention follows
        // RCPsiSquared.Core.Lindblad: flat = row * d + col → row = flat / d, col = flat % d.
        var labels = new (int pCol, int pRow)[liouvilleDim];
        for (int flat = 0; flat < liouvilleDim; flat++)
        {
            int row = flat / d;
            int col = flat % d;
            labels[flat] = (PopCount(col), PopCount(row));
        }

        // Sort flat indices by (pCol, pRow); stable so within-sector order is original ascending
        var indices = Enumerable.Range(0, liouvilleDim).ToArray();
        Array.Sort(indices, (a, b) =>
        {
            int c = labels[a].pCol.CompareTo(labels[b].pCol);
            if (c != 0) return c;
            c = labels[a].pRow.CompareTo(labels[b].pRow);
            if (c != 0) return c;
            return a.CompareTo(b);
        });

        // Build sector ranges
        var sectors = new List<SectorRange>();
        int start = 0;
        while (start < liouvilleDim)
        {
            var (pc, pr) = labels[indices[start]];
            int end = start + 1;
            while (end < liouvilleDim && labels[indices[end]] == (pc, pr)) end++;
            sectors.Add(new SectorRange(pc, pr, start, end - start));
            start = end;
        }

        return new Decomposition(N, d, indices, sectors);
    }

    private static int PopCount(int x) => System.Numerics.BitOperations.PopCount((uint)x);
}
