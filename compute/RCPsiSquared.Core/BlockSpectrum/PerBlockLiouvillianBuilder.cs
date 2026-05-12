using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.BlockSpectrum;

/// <summary>Build a single Liouville-space block of the XY+Z-dephasing Liouvillian L
/// directly, without ever materialising the full (4^N) × (4^N) matrix.
///
/// <para>For a sector with Liouville-space flat indices <c>S = { f_0, f_1, … }</c>
/// (each <c>f = r·d + c</c>, d = 2^N), this routine returns the sector block
/// <c>B[i,j] = L[f_i, f_j]</c> by computing matrix elements on the fly from the
/// Hilbert-space Hamiltonian <c>H</c> (size 2^N × 2^N — cheap even at N=8) plus the
/// per-site Z-dephasing rates.</para>
///
/// <para>Memory profile: a single block of size <c>S</c> costs <c>S² × 16 bytes</c>.
/// For N=8 the largest joint-popcount sector has C(8,4)² = 4900 entries → ~384 MB per
/// block; with F71 refinement that halves to ~190 MB. The full L would cost
/// (4^8)² × 16 = 68.7 GB — explicitly avoided here.</para>
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
}
