using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.F86.JordanWigner;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.BlockSpectrum.JordanWigner;

/// <summary>Builds the L = L_H + L_D Liouvillian block of the chain XY + Z-dephasing
/// dynamics directly in the JW Slater-pair basis of joint-popcount sector (p_c, p_r) as a
/// CSR sparse complex matrix, bypassing the dense computational-basis L and the
/// <see cref="JwSlaterPairLProjection"/> dense U^†·L·U product. This is the N=10 push
/// enabler: at sector dimension 63504 (max block at N=10) the dense path is ~64 GB and
/// physically infeasible on commodity hardware; this primitive produces the same operator
/// in ~700 MB of sparse storage.
///
/// <para>Construction is via explicit matrix-element formulas derived from the JW expansion
/// of Z_l = 1 − 2·Σ_{k,k'} ψ_k(l)·ψ_{k'}(l)·η_k†·η_{k'}:</para>
/// <list type="bullet">
///   <item><b>Diagonal</b> at α = (L, K):
///       <c>L_JW[α, α] = −i·(ε_L − ε_K) + Σ_l γ_l·[(1 − 2n_l^L)(1 − 2n_l^K) − 1]</c>
///       with <c>ε_M = Σ_{k∈M} ε_k</c> and <c>n_l^M = Σ_{k∈M} ψ_k(l)²</c>.</item>
///   <item><b>Col-only swap</b> L' = L, K' = K \ {k'_0} ∪ {k_0}:
///       <c>L_JW[α', α] = −2·sign_K·Σ_l γ_l·(1 − 2n_l^L)·ψ_{k_0}(l)·ψ_{k'_0}(l)</c>.</item>
///   <item><b>Row-only swap</b> L' = L \ {j'_0} ∪ {j_0}, K' = K:
///       <c>L_JW[α', α] = −2·sign_L·Σ_l γ_l·(1 − 2n_l^K)·ψ_{j_0}(l)·ψ_{j'_0}(l)</c>.</item>
///   <item><b>Both swap</b>:
///       <c>L_JW[α', α] = 4·sign_L·sign_K·Σ_l γ_l·ψ_{j_0}(l)·ψ_{j'_0}(l)·ψ_{k_0}(l)·ψ_{k'_0}(l)</c>.</item>
/// </list>
///
/// <para>Signs follow the standard ordered Slater convention
/// <c>|L⟩ = η_{l_1}†·...·η_{l_p}†|0⟩</c> with <c>l_1 &lt; ... &lt; l_p</c>: for a swap
/// k'_0 → k_0, sign = (−1)^(m + q) where m is the 1-indexed position of k'_0 in L and q is
/// the 1-indexed position k_0 would take if inserted into sorted L \ {k'_0}.</para>
///
/// <para><b>Tier outcome: <see cref="Tier.Tier1Derived"/>.</b> Element-by-element
/// reconstruction of the same operator <see cref="JwSlaterPairLProjection"/> builds via
/// dense U^†·L·U; validated against the dense reference at small N (4–7) to relative
/// Frobenius residual &lt; 1e-10 in the test suite. The sparsity bound matches:
/// max-nnz/row ≤ (1 + p_r·(N − p_r))·(1 + p_c·(N − p_c)) including the diagonal.</para>
///
/// <para>Anchor: <see cref="JwSlaterPairBasis"/> (basis-transform reference) +
/// <see cref="JwSlaterPairLProjection"/> (dense validation) +
/// <see cref="XyJordanWignerModes"/> (sine-mode basis + dispersion);
/// textbook free-fermion algebra for η_k† η_{k'} acting on ordered Slater states.</para>
/// </summary>
public sealed class JwSlaterPairSparseLBuilder : Claim
{
    /// <summary>Magnitude threshold for inclusion in the CSR sparse pattern; entries with
    /// |v| ≤ <see cref="NonzeroThreshold"/> are dropped to suppress FP round-off noise from
    /// the per-l accumulation (typical noise floor ~1e-15·γ).</summary>
    public const double NonzeroThreshold = 1e-12;

    public int N { get; }
    public int PCol { get; }
    public int PRow { get; }
    public double J { get; }
    public IReadOnlyList<double> GammaPerSite { get; }
    public int SectorDim { get; }

    /// <summary>CSR row-pointer array of length <see cref="SectorDim"/> + 1.
    /// <c>RowPtr[α + 1] − RowPtr[α]</c> is the nnz count of row α.</summary>
    public int[] RowPtr { get; }

    /// <summary>CSR column-index array of length <see cref="NnzTotal"/>.
    /// Sorted ascending within each row's slice.</summary>
    public int[] ColIdx { get; }

    /// <summary>CSR value array of length <see cref="NnzTotal"/>, aligned with
    /// <see cref="ColIdx"/>.</summary>
    public Complex[] Values { get; }

    public long NnzTotal => Values.LongLength;
    public int MaxNnzPerRow { get; }
    public double MeanNnzPerRow => (double)NnzTotal / SectorDim;

    public static JwSlaterPairSparseLBuilder Build(int N, int pCol, int pRow,
        IReadOnlyList<double> gammaPerSite, double J = 1.0)
    {
        if (N < 1) throw new ArgumentOutOfRangeException(nameof(N), N, "N must be ≥ 1.");
        if (pCol < 0 || pCol > N) throw new ArgumentOutOfRangeException(nameof(pCol), pCol, $"pCol must be in [0, {N}].");
        if (pRow < 0 || pRow > N) throw new ArgumentOutOfRangeException(nameof(pRow), pRow, $"pRow must be in [0, {N}].");
        if (gammaPerSite is null) throw new ArgumentNullException(nameof(gammaPerSite));
        if (gammaPerSite.Count != N)
            throw new ArgumentException(
                $"gammaPerSite length {gammaPerSite.Count} != N {N}", nameof(gammaPerSite));

        var modes = XyJordanWignerModes.Build(N, J);
        var rowSubsets = JwSlaterPairBasis.EnumerateSortedSubsets(N, pRow);
        var colSubsets = JwSlaterPairBasis.EnumerateSortedSubsets(N, pCol);
        int mr = rowSubsets.Count, mc = colSubsets.Count;
        int sectorDim = mr * mc;

        // Precompute sine-mode matrix as flat double[N*N] for fast index access:
        // psi[k_index * N + site] where k_index is 0-indexed (k_index = k - 1).
        var psi = new double[N * N];
        for (int k = 1; k <= N; k++)
            for (int site = 0; site < N; site++)
                psi[(k - 1) * N + site] = modes.SineMode(k, site);

        var dispersion = modes.Dispersion;

        var (rowSubsetToIndex, rowEps, rowOccupation) = BuildSubsetTables(rowSubsets, N, psi, dispersion);
        var (colSubsetToIndex, colEps, colOccupation) = BuildSubsetTables(colSubsets, N, psi, dispersion);

        var rowSwaps = BuildSwapTable(rowSubsets, rowSubsetToIndex, N);
        var colSwaps = BuildSwapTable(colSubsets, colSubsetToIndex, N);

        var valuesList = new List<Complex>(capacity: sectorDim * 8);
        var colIdxList = new List<int>(capacity: sectorDim * 8);
        var rowPtr = new int[sectorDim + 1];
        int maxNnz = 0;

        var perRow = new SortedDictionary<int, Complex>();

        for (int alpha = 0; alpha < sectorDim; alpha++)
        {
            int Lidx = alpha / mc;
            int Kidx = alpha % mc;

            perRow.Clear();

            // (1) Diagonal element at (alpha, alpha). Site index l is MSB-convention
            // (matches XyJordanWignerModes.SineMode); PerBlockLiouvillianBuilder.BuildBlockZ
            // applies gammaPerSite[bit_l] with bit_l = LSB-bit-l = N − 1 − site_l, so the
            // sparse build must mirror that mapping: gammaPerSite[N - 1 - l].
            double imagH = -(rowEps[Lidx] - colEps[Kidx]);
            double diagD = 0.0;
            for (int l = 0; l < N; l++)
            {
                double nL = rowOccupation[Lidx * N + l];
                double nK = colOccupation[Kidx * N + l];
                diagD += gammaPerSite[N - 1 - l] * ((1.0 - 2.0 * nL) * (1.0 - 2.0 * nK) - 1.0);
            }
            var diag = new Complex(diagD, imagH);
            if (diag.Magnitude > NonzeroThreshold)
                perRow[alpha] = diag;

            // (2) Col-only swaps: L' = L, K' = K \ {kprime_0} ∪ {k_0}.
            foreach (var swap in colSwaps[Kidx])
            {
                double sum = 0.0;
                for (int l = 0; l < N; l++)
                {
                    double nL = rowOccupation[Lidx * N + l];
                    sum += gammaPerSite[N - 1 - l] * (1.0 - 2.0 * nL)
                         * psi[(swap.KEnter - 1) * N + l] * psi[(swap.KLeave - 1) * N + l];
                }
                double val = -2.0 * swap.Sign * sum;
                if (Math.Abs(val) > NonzeroThreshold)
                {
                    int beta = Lidx * mc + swap.TargetIndex;
                    Accumulate(perRow, beta, new Complex(val, 0.0));
                }
            }

            // (3) Row-only swaps: L' = L \ {jprime_0} ∪ {j_0}, K' = K.
            foreach (var swap in rowSwaps[Lidx])
            {
                double sum = 0.0;
                for (int l = 0; l < N; l++)
                {
                    double nK = colOccupation[Kidx * N + l];
                    sum += gammaPerSite[N - 1 - l] * (1.0 - 2.0 * nK)
                         * psi[(swap.KEnter - 1) * N + l] * psi[(swap.KLeave - 1) * N + l];
                }
                double val = -2.0 * swap.Sign * sum;
                if (Math.Abs(val) > NonzeroThreshold)
                {
                    int beta = swap.TargetIndex * mc + Kidx;
                    Accumulate(perRow, beta, new Complex(val, 0.0));
                }
            }

            // (4) Both-swap: L → L', K → K'.
            foreach (var rowSwap in rowSwaps[Lidx])
            {
                foreach (var colSwap in colSwaps[Kidx])
                {
                    double sum = 0.0;
                    for (int l = 0; l < N; l++)
                    {
                        sum += gammaPerSite[N - 1 - l]
                             * psi[(rowSwap.KEnter - 1) * N + l] * psi[(rowSwap.KLeave - 1) * N + l]
                             * psi[(colSwap.KEnter - 1) * N + l] * psi[(colSwap.KLeave - 1) * N + l];
                    }
                    double val = 4.0 * rowSwap.Sign * colSwap.Sign * sum;
                    if (Math.Abs(val) > NonzeroThreshold)
                    {
                        int beta = rowSwap.TargetIndex * mc + colSwap.TargetIndex;
                        Accumulate(perRow, beta, new Complex(val, 0.0));
                    }
                }
            }

            int rowNnz = perRow.Count;
            if (rowNnz > maxNnz) maxNnz = rowNnz;
            foreach (var kv in perRow)
            {
                colIdxList.Add(kv.Key);
                valuesList.Add(kv.Value);
            }
            rowPtr[alpha + 1] = colIdxList.Count;
        }

        return new JwSlaterPairSparseLBuilder(N, pCol, pRow, J, gammaPerSite.ToArray(),
            sectorDim, rowPtr, colIdxList.ToArray(), valuesList.ToArray(), maxNnz);
    }

    private static void Accumulate(SortedDictionary<int, Complex> map, int beta, Complex value)
    {
        if (map.TryGetValue(beta, out var existing))
            map[beta] = existing + value;
        else
            map[beta] = value;
    }

    private static (Dictionary<long, int> subsetToIndex, double[] eps, double[] occupation)
        BuildSubsetTables(IReadOnlyList<int[]> subsets, int N, double[] psi, IReadOnlyList<double> dispersion)
    {
        int m = subsets.Count;
        var subsetToIndex = new Dictionary<long, int>(capacity: m);
        var eps = new double[m];
        var occupation = new double[m * N];
        for (int idx = 0; idx < m; idx++)
        {
            var subset = subsets[idx];
            subsetToIndex[ToBitmask(subset)] = idx;
            double e = 0.0;
            foreach (var k in subset) e += dispersion[k - 1];
            eps[idx] = e;
            for (int l = 0; l < N; l++)
            {
                double n = 0.0;
                foreach (var k in subset)
                {
                    double v = psi[(k - 1) * N + l];
                    n += v * v;
                }
                occupation[idx * N + l] = n;
            }
        }
        return (subsetToIndex, eps, occupation);
    }

    private static List<SwapEntry>[] BuildSwapTable(IReadOnlyList<int[]> subsets,
        Dictionary<long, int> subsetToIndex, int N)
    {
        int m = subsets.Count;
        var table = new List<SwapEntry>[m];
        var member = new bool[N + 1];  // 1-indexed
        for (int idx = 0; idx < m; idx++)
        {
            var subset = subsets[idx];
            for (int i = 1; i <= N; i++) member[i] = false;
            foreach (var k in subset) member[k] = true;

            var entries = new List<SwapEntry>(capacity: subset.Length * (N - subset.Length));
            // For each (kLeave ∈ subset, kEnter ∉ subset), compute target subset bitmask + sign.
            for (int posLeave = 0; posLeave < subset.Length; posLeave++)
            {
                int kLeave = subset[posLeave];
                int m1Indexed = posLeave + 1;  // 1-indexed position of kLeave in subset
                for (int kEnter = 1; kEnter <= N; kEnter++)
                {
                    if (member[kEnter]) continue;
                    // q: 1-indexed insertion position of kEnter in (subset \ {kLeave}) sorted.
                    int q = 1;
                    for (int i = 0; i < subset.Length; i++)
                    {
                        if (subset[i] == kLeave) continue;
                        if (subset[i] < kEnter) q++;
                    }
                    int sign = (((m1Indexed + q) & 1) == 0) ? 1 : -1;
                    long targetMask = ToBitmask(subset) ^ (1L << (kLeave - 1)) ^ (1L << (kEnter - 1));
                    int targetIndex = subsetToIndex[targetMask];
                    entries.Add(new SwapEntry(targetIndex, kEnter, kLeave, sign));
                }
            }
            table[idx] = entries;
        }
        return table;
    }

    private static long ToBitmask(int[] sortedSubset)
    {
        long mask = 0L;
        foreach (var k in sortedSubset) mask |= 1L << (k - 1);
        return mask;
    }

    /// <summary>Convert CSR to a dense <see cref="ComplexMatrix"/>. Intended for small-N
    /// validation against <see cref="JwSlaterPairLProjection"/>; will OOM at large
    /// sectorDim (the entire purpose of the sparse build is to avoid the dense
    /// SectorDim²·16 byte storage cost).</summary>
    public ComplexMatrix ToDense()
    {
        var dense = Matrix<Complex>.Build.Dense(SectorDim, SectorDim);
        for (int alpha = 0; alpha < SectorDim; alpha++)
            for (int e = RowPtr[alpha]; e < RowPtr[alpha + 1]; e++)
                dense[alpha, ColIdx[e]] = Values[e];
        return dense;
    }

    private JwSlaterPairSparseLBuilder(
        int n, int pCol, int pRow, double j, double[] gamma,
        int sectorDim, int[] rowPtr, int[] colIdx, Complex[] values, int maxNnz)
        : base($"L_JW direct sparse build on Slater-pair sector (p_c={pCol}, p_r={pRow}, N={n}): " +
               $"dim {sectorDim}, nnz {values.LongLength}, max-nnz/row {maxNnz} ≤ theoretical " +
               $"{(1 + pRow * (n - pRow)) * (1 + pCol * (n - pCol))} = (1 + p_r·(N−p_r))·(1 + p_c·(N−p_c)).",
               Tier.Tier1Derived,
               "compute/RCPsiSquared.Core/BlockSpectrum/JordanWigner/JwSlaterPairBasis.cs (basis + subset enumeration) + " +
               "compute/RCPsiSquared.Core/BlockSpectrum/JordanWigner/JwSlaterPairLProjection.cs (dense validation reference) + " +
               "compute/RCPsiSquared.Core/F86/JordanWigner/XyJordanWignerModes.cs (sine modes + dispersion); " +
               "free-fermion algebra for η_k†·η_{k'} on ordered Slater states with sign (−1)^(m+q).")
    {
        N = n; PCol = pCol; PRow = pRow; J = j;
        GammaPerSite = gamma;
        SectorDim = sectorDim;
        RowPtr = rowPtr;
        ColIdx = colIdx;
        Values = values;
        MaxNnzPerRow = maxNnz;
    }

    public override string DisplayName =>
        $"L_JW sparse build (p_c={PCol}, p_r={PRow}, N={N}); dim {SectorDim}, nnz {NnzTotal:N0}";

    public override string Summary =>
        $"nnz={NnzTotal:N0}, max-nnz/row={MaxNnzPerRow}, mean={MeanNnzPerRow:F1} ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return InspectableNode.RealScalar("N", N);
            yield return InspectableNode.RealScalar("p_c", PCol);
            yield return InspectableNode.RealScalar("p_r", PRow);
            yield return InspectableNode.RealScalar("J", J);
            yield return InspectableNode.RealScalar("sector dimension", SectorDim);
            yield return InspectableNode.RealScalar("nnz total", NnzTotal);
            yield return InspectableNode.RealScalar("max nnz/row", MaxNnzPerRow);
            yield return InspectableNode.RealScalar("mean nnz/row", MeanNnzPerRow, "F2");
        }
    }

    /// <summary>One source → target Slater swap: <c>k_enter</c> replaces <c>k_leave</c>
    /// (both 1-indexed mode labels), producing the subset indexed by <c>TargetIndex</c>
    /// with fermion sign <c>Sign ∈ {+1, −1}</c>.</summary>
    private readonly record struct SwapEntry(int TargetIndex, int KEnter, int KLeave, int Sign);
}
