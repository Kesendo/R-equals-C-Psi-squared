using System;
using System.Numerics;
using System.Threading.Tasks;
using RCPsiSquared.Core.F89PathK;

namespace RCPsiSquared.Core.Numerics;

/// <summary>
/// Numeric kernels on the <see cref="WeightCoherenceSectorCsr.Csr"/> sparse operator: the plain
/// matvec y = M·x, the shifted matvec y = (M − shift·I)·x, and the exact Hermitian (conjugate)
/// transpose. These are the primitives every downstream consumer of the weight-coherence sectors
/// needs (the BiCGStab inner solve, the shift-invert census). The matvec loop shape is the one used
/// by <see cref="BlockSpectrum.SparseShiftInvertArnoldi"/>'s ApplyShiftedMatvec (a Parallel.For over
/// rows, each row summing values[e]·x[colIdx[e]] and writing only y[row], so the parallelism is
/// race-free); here it is exposed as a standalone reusable kernel rather than a private helper.
/// </summary>
public static class CsrOps
{
    /// <summary>y = M·x. Parallel over rows; each row writes only y[row], so no synchronisation is
    /// needed. x and y must be length m.Dim; y is fully overwritten.</summary>
    public static void Multiply(WeightCoherenceSectorCsr.Csr m, Complex[] x, Complex[] y)
    {
        if (m is null) throw new ArgumentNullException(nameof(m));
        if (x is null) throw new ArgumentNullException(nameof(x));
        if (y is null) throw new ArgumentNullException(nameof(y));
        if (x.Length != m.Dim) throw new ArgumentException($"x length {x.Length} != dim {m.Dim}", nameof(x));
        if (y.Length != m.Dim) throw new ArgumentException($"y length {y.Length} != dim {m.Dim}", nameof(y));

        int[] rowPtr = m.RowPtr;
        int[] colIdx = m.ColIdx;
        Complex[] values = m.Values;
        Parallel.For(0, m.Dim, row =>
        {
            Complex sum = Complex.Zero;
            int start = rowPtr[row];
            int end = rowPtr[row + 1];
            for (int e = start; e < end; e++)
                sum += values[e] * x[colIdx[e]];
            y[row] = sum;
        });
    }

    /// <summary>y = (M − shift·I)·x. The diagonal shift is applied on top of the CSR matvec, exactly
    /// as <see cref="BlockSpectrum.SparseShiftInvertArnoldi"/>'s ApplyShiftedMatvec does. Parallel over
    /// rows; each row writes only y[row].</summary>
    public static void MultiplyShifted(WeightCoherenceSectorCsr.Csr m, Complex shift, Complex[] x, Complex[] y)
    {
        if (m is null) throw new ArgumentNullException(nameof(m));
        if (x is null) throw new ArgumentNullException(nameof(x));
        if (y is null) throw new ArgumentNullException(nameof(y));
        if (x.Length != m.Dim) throw new ArgumentException($"x length {x.Length} != dim {m.Dim}", nameof(x));
        if (y.Length != m.Dim) throw new ArgumentException($"y length {y.Length} != dim {m.Dim}", nameof(y));

        int[] rowPtr = m.RowPtr;
        int[] colIdx = m.ColIdx;
        Complex[] values = m.Values;
        Parallel.For(0, m.Dim, row =>
        {
            Complex sum = -shift * x[row];
            int start = rowPtr[row];
            int end = rowPtr[row + 1];
            for (int e = start; e < end; e++)
                sum += values[e] * x[colIdx[e]];
            y[row] = sum;
        });
    }

    /// <summary>The exact Hermitian (conjugate) transpose Mᴴ: (Mᴴ)[i,j] = conj(M[j,i]). Computed by a
    /// counting transpose in O(nnz): count the entries in each source column (= destination row), prefix-
    /// sum into the destination RowPtr, then scatter each source entry (r, c, v) to destination (c, r,
    /// conj(v)). Iterating the source rows in ascending order makes each destination row's column indices
    /// (the source row numbers) land already sorted with no duplicates, matching the CSR invariant of
    /// <see cref="WeightCoherenceSectorCsr.BuildFull"/>. No eigensolver, no tolerance: the conjugation is
    /// an exact sign flip on the imaginary part, so HermitianTranspose∘HermitianTranspose is the identity
    /// bit-for-bit.</summary>
    public static WeightCoherenceSectorCsr.Csr HermitianTranspose(WeightCoherenceSectorCsr.Csr m)
    {
        if (m is null) throw new ArgumentNullException(nameof(m));

        int dim = m.Dim;
        int[] rowPtr = m.RowPtr;
        int[] colIdx = m.ColIdx;
        Complex[] values = m.Values;
        int nnz = values.Length;

        // destination RowPtr: count entries per source column, then prefix-sum
        var tRowPtr = new int[dim + 1];
        for (int k = 0; k < nnz; k++) tRowPtr[colIdx[k] + 1]++;
        for (int i = 0; i < dim; i++) tRowPtr[i + 1] += tRowPtr[i];

        var tColIdx = new int[nnz];
        var tValues = new Complex[nnz];
        var cursor = (int[])tRowPtr.Clone();
        for (int r = 0; r < dim; r++)                     // ascending source rows => sorted dest columns
        {
            int start = rowPtr[r];
            int end = rowPtr[r + 1];
            for (int k = start; k < end; k++)
            {
                int c = colIdx[k];                        // source (r, c) -> destination (c, r)
                int dest = cursor[c]++;
                tColIdx[dest] = r;
                tValues[dest] = Complex.Conjugate(values[k]);
            }
        }
        return new WeightCoherenceSectorCsr.Csr(dim, tRowPtr, tColIdx, tValues);
    }
}
