using System;
using System.Collections.Generic;
using System.Numerics;

namespace RCPsiSquared.Core.F89PathK;

/// <summary>
/// CSR emission of the (wKet,wBra) weight-coherence block pencil L(q) = A + qC and of its
/// R-parity sectors, entry-identical to <see cref="WeightCoherenceBlock.Build"/> /
/// <see cref="WeightCoherenceBlock.BuildReflectionSectorColumnMajor"/> (the dense builders remain
/// the source of truth; the equivalence is pinned in WeightCoherenceSectorCsrTests). nnz per row is
/// O(N), so the N=11 deferred sectors (dim 54k..107k) fit in ~100 MB where dense storage would need
/// 47..182 GB. Basis convention: full block row = ketIndex * nBra + braIndex, exactly the dense
/// builder's ordering; sector basis = reflection fixed points (even sector only) then 2-cycle
/// representatives, exactly BuildReflectionSectorColumnMajor's ordering.
/// </summary>
public static class WeightCoherenceSectorCsr
{
    /// <summary>Compressed-sparse-row storage of a square complex operator: for row r the entries
    /// are Values[RowPtr[r]..RowPtr[r+1]) at columns ColIdx[.], columns ascending, no duplicates.</summary>
    public sealed record Csr(int Dim, int[] RowPtr, int[] ColIdx, Complex[] Values);

    /// <summary>The full (wKet,wBra) chain coherence block at complex q (γ = 1, Δ = 0), emitted as CSR,
    /// entry-identical to <see cref="WeightCoherenceBlock.Build"/>. Full-basis index = ketPos * nBra + braPos.</summary>
    public static Csr BuildFull(int n, int wKet, int wBra, Complex q)
    {
        var kets = WeightCoherenceBlock.Configs(n, wKet);
        var bras = WeightCoherenceBlock.Configs(n, wBra);
        var ketIndex = Index(kets);
        var braIndex = Index(bras);
        int nb = bras.Count, dim = kets.Count * nb;
        var cols = new List<(int Row, Complex Val)>[dim];
        var acc = new Dictionary<int, Complex>();
        for (int col = 0; col < dim; col++)
        {
            acc.Clear();
            ApplyColumn(n, col, kets, bras, ketIndex, braIndex, nb, q, Complex.One, acc);
            var list = new List<(int, Complex)>(acc.Count);
            foreach (var kv in acc) list.Add((kv.Key, kv.Value));
            cols[col] = list;
        }
        return CscToCsr(dim, cols);
    }

    /// <summary>One R-parity sector of the (wKet,wBra) chain block at complex q (γ = 1, Δ = 0), emitted as CSR,
    /// entry-identical to <see cref="WeightCoherenceBlock.BuildReflectionSectorColumnMajor"/>. Transcribes its
    /// sector-coordinate wrapping verbatim: basis vectors are the reflection fixed points e_t (even sector only)
    /// and the 2-cycle combinations (e_t ± e_{Rt})/√2 over orbit reps t &lt; Rt, in increasing t; the full-block
    /// column rule is applied to each sector basis vector and projected back onto the real orthonormal sector basis.</summary>
    public static Csr BuildReflectionSector(int n, int wKet, int wBra, Complex q, bool odd)
    {
        var kets = WeightCoherenceBlock.Configs(n, wKet);
        var bras = WeightCoherenceBlock.Configs(n, wBra);
        int nb = bras.Count;
        var ketIndex = Index(kets);
        var braIndex = Index(bras);
        var perm = WeightCoherenceBlock.ReflectionPermutation(n, wKet, wBra);

        // orbit reps: fixed points (even only), then 2-cycle reps t < perm[t]; orbitOf = sector row or −1
        var reps = new List<int>();
        var orbitOf = new int[perm.Length];
        Array.Fill(orbitOf, -1);
        for (int i = 0; i < perm.Length; i++)
        {
            if (perm[i] == i) { if (!odd) { orbitOf[i] = reps.Count; reps.Add(i); } }
            else if (i < perm[i]) { orbitOf[i] = orbitOf[perm[i]] = reps.Count; reps.Add(i); }
        }
        int d = reps.Count;
        double sSign = odd ? -1.0 : 1.0;
        double inv2 = 1.0 / Math.Sqrt(2.0);

        var cols = new List<(int Row, Complex Val)>[d];
        var acc = new Dictionary<int, Complex>();
        for (int col = 0; col < d; col++)
        {
            acc.Clear();
            int t = reps[col];
            bool fixedPt = perm[t] == t;
            double wNorm = fixedPt ? 1.0 : inv2;
            ApplyColumn(n, t, kets, bras, ketIndex, braIndex, nb, q, wNorm, acc);
            if (!fixedPt) ApplyColumn(n, perm[t], kets, bras, ketIndex, braIndex, nb, q, sSign * wNorm, acc);

            var sectorCol = new Dictionary<int, Complex>();
            foreach (var (fullRow, val) in acc)
            {
                int row = orbitOf[fullRow];
                if (row < 0) continue;                            // row absent from this sector (odd ∌ fixed)
                int rep = reps[row];
                // coefficient of e_fullRow inside the row's REAL sector basis vector (no conjugation)
                double coeff = perm[rep] == rep ? 1.0
                             : (fullRow == rep ? inv2 : sSign * inv2);
                AddTo(sectorCol, row, coeff * val);
            }
            var list = new List<(int, Complex)>(sectorCol.Count);
            foreach (var kv in sectorCol) list.Add((kv.Key, kv.Value));
            cols[col] = list;
        }
        return CscToCsr(d, cols);
    }

    private static Dictionary<int, int> Index(List<int> configs)
    {
        var index = new Dictionary<int, int>(configs.Count);
        for (int i = 0; i < configs.Count; i++) index[configs[i]] = i;
        return index;
    }

    /// <summary>The exact column action of L on full-basis index fullCol, scaled by weight (the same hop
    /// rule as <see cref="WeightCoherenceBlock.Build"/>: diagonal −2·n_diff, ket hops −2iq, bra hops +2iq,
    /// nearest-neighbour, Pauli-excluded). Accumulates into acc keyed by full-basis row index.</summary>
    private static void ApplyColumn(int n, int fullCol, List<int> kets, List<int> bras,
        Dictionary<int, int> ketIndex, Dictionary<int, int> braIndex, int nb, Complex q, Complex weight,
        Dictionary<int, Complex> acc)
    {
        int kc = kets[fullCol / nb], bc = bras[fullCol % nb];
        AddTo(acc, fullCol, weight * new Complex(-2.0 * BitOperations.PopCount((uint)(kc ^ bc)), 0));
        for (int site = 0; site < n; site++)
            if ((kc & (1 << site)) != 0)
                foreach (int s2 in new[] { site - 1, site + 1 })
                    if (s2 >= 0 && s2 < n && (kc & (1 << s2)) == 0)
                        AddTo(acc, ketIndex[(kc & ~(1 << site)) | (1 << s2)] * nb + braIndex[bc],
                            weight * Complex.ImaginaryOne * -2.0 * q);
        for (int site = 0; site < n; site++)
            if ((bc & (1 << site)) != 0)
                foreach (int s2 in new[] { site - 1, site + 1 })
                    if (s2 >= 0 && s2 < n && (bc & (1 << s2)) == 0)
                        AddTo(acc, ketIndex[kc] * nb + braIndex[(bc & ~(1 << site)) | (1 << s2)],
                            weight * Complex.ImaginaryOne * 2.0 * q);
    }

    private static void AddTo(Dictionary<int, Complex> acc, int row, Complex v)
        => acc[row] = acc.TryGetValue(row, out var cur) ? cur + v : v;

    /// <summary>Assemble CSR from per-column entry lists (a CSC view). Iterating columns ascending makes
    /// each row's column indices land in ascending order, so ColIdx per row is sorted with no duplicates.</summary>
    private static Csr CscToCsr(int dim, List<(int Row, Complex Val)>[] cols)
    {
        var rowPtr = new int[dim + 1];
        foreach (var col in cols) foreach (var (r, _) in col) rowPtr[r + 1]++;
        for (int i = 0; i < dim; i++) rowPtr[i + 1] += rowPtr[i];
        int nnz = rowPtr[dim];
        var colIdx = new int[nnz];
        var values = new Complex[nnz];
        var cursor = (int[])rowPtr.Clone();
        for (int c = 0; c < dim; c++)                 // ascending c => sorted colIdx per row
            foreach (var (r, v) in cols[c])
            { colIdx[cursor[r]] = c; values[cursor[r]] = v; cursor[r]++; }
        return new Csr(dim, rowPtr, colIdx, values);
    }
}
