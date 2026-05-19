using System.Numerics;
using System.Runtime.InteropServices;
using MathNet.Numerics.LinearAlgebra;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.BlockSpectrum;

/// <summary>Build a single Liouville-space block of the XY+Z-dephasing Liouvillian L
/// directly, without ever materialising the full (4^N) × (4^N) matrix.
///
/// <para>For a sector with Liouville-space flat indices <c>S = { f_0, f_1, … }</c>
/// (each <c>f = r·d + c</c>, d = 2^N), this routine returns the sector block
/// <c>B[i,j] = L[f_i, f_j]</c> by computing matrix elements on the fly from the
/// Hilbert-space Hamiltonian <c>H</c> (size 2^N × 2^N, cheap even at N=8) plus the
/// per-site Z-dephasing rates.</para>
///
/// <para>Memory profile: a single block of size <c>S</c> costs <c>S² × 16 bytes</c>.
/// For N=8 the largest joint-popcount sector has C(8,4)² = 4900 entries → ~384 MB per
/// block; with F71 refinement that halves to ~190 MB. The full L would cost
/// (4^8)² × 16 = 68.7 GB, explicitly avoided here.</para>
///
/// <para>Convention: matches <see cref="Lindblad.PauliDephasingDissipator.BuildZ"/>
/// row-major <c>flat = row · d + col</c> layout, with <c>L = -i (H⊗I − I⊗Hᵀ) +
/// Σ_l γ_l (Z_l ⊗ Z_l) − Σ_l γ_l · I</c>. The Z⊗Z term is diagonal in (row, col)
/// pair so contributes only along <c>flat = flat'</c>.</para></summary>
public static class PerBlockLiouvillianBuilder
{
    /// <summary>Build the block <c>B[i,j] = L[flatIndices[i], flatIndices[j]]</c>
    /// of the Z-dephasing Liouvillian for the given Hilbert-space Hamiltonian and
    /// per-site γ rates.</summary>
    /// <param name="H">Hilbert-space Hamiltonian, dense 2^N × 2^N.</param>
    /// <param name="gammaPerSite">Per-site Z-dephasing rates (length N).</param>
    /// <param name="flatIndices">Liouville-space flat indices defining the block.</param>
    /// <returns>Dense complex matrix of size flatIndices.Length × flatIndices.Length.</returns>
    public static ComplexMatrix BuildBlockZ(ComplexMatrix H, IReadOnlyList<double> gammaPerSite, IReadOnlyList<int> flatIndices)
    {
        if (H is null) throw new ArgumentNullException(nameof(H));
        if (gammaPerSite is null) throw new ArgumentNullException(nameof(gammaPerSite));
        if (flatIndices is null) throw new ArgumentNullException(nameof(flatIndices));

        int d = H.RowCount;
        int N = (int)Math.Round(Math.Log2(d));
        if ((1 << N) != d) throw new ArgumentException($"H dim {d} is not a power of 2", nameof(H));
        if (gammaPerSite.Count != N)
            throw new ArgumentException($"gamma list length {gammaPerSite.Count} != N={N}", nameof(gammaPerSite));

        int blockSize = flatIndices.Count;
        var B = Matrix<Complex>.Build.Dense(blockSize, blockSize);

        // Pre-decompose flat indices into (row, col).
        var rowOf = new int[blockSize];
        var colOf = new int[blockSize];
        for (int i = 0; i < blockSize; i++)
        {
            int f = flatIndices[i];
            rowOf[i] = f / d;
            colOf[i] = f % d;
        }

        // Diagonal Z-dephasing contribution: for flat = row*d + col,
        //   (Z_l ⊗ Z_l)[flat, flat] = (-1)^{bit_l(row) + bit_l(col)} = +1 if row,col match in bit l, else -1.
        //   Identity term contributes -1 per site with rate γ_l.
        // Net diagonal contribution at flat:
        //   diagDeph = Σ_l γ_l ( (-1)^{bit_l(row) ⊕ bit_l(col)} − 1 )
        //            = Σ_l (-2 γ_l) · [bit_l(row) != bit_l(col)]
        // (i.e. each disagreement contributes -2 γ_l).
        // This formula is bit-exact equivalent to the full PauliDephasingDissipator.BuildZ
        // diagonal sum over the N sites.
        for (int i = 0; i < blockSize; i++)
        {
            int r = rowOf[i];
            int c = colOf[i];
            int diff = r ^ c;
            double diag = 0;
            for (int l = 0; l < N; l++)
                if (((diff >> l) & 1) != 0) diag += -2.0 * gammaPerSite[l];
            B[i, i] = new Complex(diag, 0);
        }

        // Off-diagonal & diagonal Hamiltonian contribution:
        //   L[(r1,c1), (r2,c2)] = −i · ( H[r1, r2] · δ_{c1, c2} − δ_{r1, r2} · H[c2, c1] )
        // Iterate over (i, j) pairs and add the H contribution; the Z-deph diagonal
        // already added above is independent.
        // Optimisation: bucket block indices by row and by col so we only visit pairs
        // sharing one coordinate. For an XY chain this is O(blockSize · N) rather than
        // O(blockSize²) for the brute scan.
        var byRow = new Dictionary<int, List<int>>();
        var byCol = new Dictionary<int, List<int>>();
        for (int i = 0; i < blockSize; i++)
        {
            if (!byRow.TryGetValue(rowOf[i], out var listR)) { listR = new List<int>(); byRow[rowOf[i]] = listR; }
            listR.Add(i);
            if (!byCol.TryGetValue(colOf[i], out var listC)) { listC = new List<int>(); byCol[colOf[i]] = listC; }
            listC.Add(i);
        }

        // First Hamiltonian term: −i · H[r1, r2] · δ_{c1, c2}.
        // For each j in block, fix c2 = colOf[j], r2 = rowOf[j]; sum over i with c1 == c2:
        //   B[i, j] += −i · H[r1=rowOf[i], r2=rowOf[j]].
        // Iterate i over the bucket byCol[c2].
        var minusI = new Complex(0, -1);
        for (int j = 0; j < blockSize; j++)
        {
            int r2 = rowOf[j];
            int c2 = colOf[j];
            if (byCol.TryGetValue(c2, out var bucket))
            {
                foreach (int i in bucket)
                {
                    int r1 = rowOf[i];
                    Complex h = H[r1, r2];
                    if (h.Magnitude > 0)
                        B[i, j] += minusI * h;
                }
            }
        }

        // Second Hamiltonian term: −i · ( − δ_{r1, r2} · H[c2, c1] ) = +i · H[c2, c1] · δ_{r1, r2}.
        // For each j, fix r2 = rowOf[j], c2 = colOf[j]; sum over i with r1 == r2:
        //   B[i, j] += +i · H[c2, c1=colOf[i]].
        var plusI = new Complex(0, 1);
        for (int j = 0; j < blockSize; j++)
        {
            int r2 = rowOf[j];
            int c2 = colOf[j];
            if (byRow.TryGetValue(r2, out var bucket))
            {
                foreach (int i in bucket)
                {
                    int c1 = colOf[i];
                    Complex h = H[c2, c1];
                    if (h.Magnitude > 0)
                        B[i, j] += plusI * h;
                }
            }
        }

        return B;
    }

    /// <summary>Build the block <c>B[i,j] = L[flatIndices[i], flatIndices[j]]</c> into a freshly
    /// allocated native (unmanaged) Complex buffer in column-major layout
    /// (<c>buf[j*blockSize + i] = B[i,j]</c>), bypassing both MathNet's managed wrapper and the
    /// LP64 2 GB single-array marshalling ceiling. Memory math: <c>blockSize² × 16 B</c>; at the
    /// N=9 chain's largest joint-popcount sector (15 876² = ~4 GB) this is the only path that
    /// passes the LAPACK boundary.
    ///
    /// <para><b>Ownership transfer.</b> The returned <see cref="IntPtr"/> points at a buffer
    /// allocated via <see cref="NativeMemory.AllocZeroed(nuint)"/> and owned by the caller; the
    /// caller MUST release it with <see cref="NativeMemory.Free(void*)"/> in a <c>finally</c>
    /// block paired with the LAPACK call. Forgetting to free leaks the entire block (~4 GB at
    /// N=9 chain). The LAPACK convention is that <c>zgeev</c> destroys its input matrix; the
    /// buffer therefore has single-use semantics after handing it to MklDirect.</para>
    ///
    /// <para>The math is identical to <see cref="BuildBlockZ"/>; only the storage backend
    /// differs (native column-major <c>Complex*</c> instead of MathNet <c>Matrix&lt;Complex&gt;</c>
    /// in its <c>DenseColumnMajorMatrixStorage</c>). The two methods produce bit-exact agreement
    /// at block sizes where both are valid; this is the parity witness for the LP64 bridge.</para>
    ///
    /// <para>Anchors: <see cref="LiouvillianBlockSpectrum.ComputeSpectrumPerBlock"/> dispatches
    /// to this path for blocks above the LP64 ceiling (block size > 11 585); the parity test
    /// <c>PerBlockLiouvillianBuilderNativeMemoryParityTests</c> verifies bit-exactness vs the
    /// MathNet path at small N where both routes are exercisable.</para></summary>
    /// <param name="H">Hilbert-space Hamiltonian, dense 2^N × 2^N (cheap even at N=9: 512×512).</param>
    /// <param name="gammaPerSite">Per-site Z-dephasing rates (length N).</param>
    /// <param name="flatIndices">Liouville-space flat indices defining the block.</param>
    /// <returns>An <see cref="IntPtr"/> at a column-major <c>Complex[blockSize², row-major (i, j) =
    /// j*blockSize + i]</c> native buffer; caller frees with
    /// <see cref="NativeMemory.Free(void*)"/>.</returns>
    public static unsafe IntPtr BuildBlockZIntoNativeMemory(
        ComplexMatrix H, IReadOnlyList<double> gammaPerSite, IReadOnlyList<int> flatIndices)
    {
        if (H is null) throw new ArgumentNullException(nameof(H));
        if (gammaPerSite is null) throw new ArgumentNullException(nameof(gammaPerSite));
        if (flatIndices is null) throw new ArgumentNullException(nameof(flatIndices));

        int d = H.RowCount;
        int N = (int)Math.Round(Math.Log2(d));
        if ((1 << N) != d) throw new ArgumentException($"H dim {d} is not a power of 2", nameof(H));
        if (gammaPerSite.Count != N)
            throw new ArgumentException($"gamma list length {gammaPerSite.Count} != N={N}", nameof(gammaPerSite));

        int blockSize = flatIndices.Count;
        long elements = (long)blockSize * blockSize;
        long bytes = elements * sizeof(Complex);

        // Native allocation: zeroed so the unwritten upper-triangle off-diagonals start at 0+0i,
        // matching MathNet's Build.Dense initial state. AllocZeroed is critical: the partial-fill
        // loops below only touch the cells where a contribution exists.
        var ptr = NativeMemory.AllocZeroed((nuint)bytes);
        var data = (Complex*)ptr;

        try
        {
            // Pre-decompose flat indices into (row, col); identical to BuildBlockZ.
            var rowOf = new int[blockSize];
            var colOf = new int[blockSize];
            for (int i = 0; i < blockSize; i++)
            {
                int f = flatIndices[i];
                rowOf[i] = f / d;
                colOf[i] = f % d;
            }

            // Diagonal Z-dephasing: same formula as BuildBlockZ. Column-major position of (i,i)
            // is i * blockSize + i; per element accumulation as if into B[i, i].
            for (int i = 0; i < blockSize; i++)
            {
                int r = rowOf[i];
                int c = colOf[i];
                int diff = r ^ c;
                double diag = 0;
                for (int l = 0; l < N; l++)
                    if (((diff >> l) & 1) != 0) diag += -2.0 * gammaPerSite[l];
                data[(long)i * blockSize + i] = new Complex(diag, 0);
            }

            // Bucket block indices by row and by col; same bucket structure as BuildBlockZ.
            var byRow = new Dictionary<int, List<int>>();
            var byCol = new Dictionary<int, List<int>>();
            for (int i = 0; i < blockSize; i++)
            {
                if (!byRow.TryGetValue(rowOf[i], out var listR)) { listR = new List<int>(); byRow[rowOf[i]] = listR; }
                listR.Add(i);
                if (!byCol.TryGetValue(colOf[i], out var listC)) { listC = new List<int>(); byCol[colOf[i]] = listC; }
                listC.Add(i);
            }

            // First Hamiltonian term: −i · H[r1, r2] · δ_{c1, c2}; columns indexed by j.
            // Column-major position of (i, j) is j * blockSize + i.
            var minusI = new Complex(0, -1);
            for (int j = 0; j < blockSize; j++)
            {
                int r2 = rowOf[j];
                int c2 = colOf[j];
                long colBase = (long)j * blockSize;
                if (byCol.TryGetValue(c2, out var bucket))
                {
                    foreach (int i in bucket)
                    {
                        int r1 = rowOf[i];
                        Complex h = H[r1, r2];
                        if (h.Magnitude > 0)
                            data[colBase + i] += minusI * h;
                    }
                }
            }

            // Second Hamiltonian term: +i · H[c2, c1] · δ_{r1, r2}.
            var plusI = new Complex(0, 1);
            for (int j = 0; j < blockSize; j++)
            {
                int r2 = rowOf[j];
                int c2 = colOf[j];
                long colBase = (long)j * blockSize;
                if (byRow.TryGetValue(r2, out var bucket))
                {
                    foreach (int i in bucket)
                    {
                        int c1 = colOf[i];
                        Complex h = H[c2, c1];
                        if (h.Magnitude > 0)
                            data[colBase + i] += plusI * h;
                    }
                }
            }

            return (IntPtr)data;
        }
        catch
        {
            // If anything throws after AllocZeroed succeeds, free immediately so the caller's
            // try/finally is the sole owner from this point on. Without this, exceptions before
            // return leak the buffer.
            NativeMemory.Free(data);
            throw;
        }
    }
}
