using System.Numerics;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.BlockSpectrum;

/// <summary>Sparse-CSR variant of <see cref="KleinFourGroupSelfPairedRefinement.BuildSubBlockL"/>:
/// for one Klein character χ on the (N/2, N/2) self-paired sector, builds the L_χ
/// sub-block as a CSR triplet (RowPtr, ColIdx, Values) without ever materialising the
/// dense matrix. At N=10 (5, 5) the dense path trips the .NET/MKL int32 marshaling cap
/// (4.2 GB > 2 GB) — sparse storage drops the memory ~600× to ~7 MB per sub-block
/// (mean nnz/row ≈ 10.7, ~170 k total nnz at dim 16132), making the sub-block streamable
/// through <see cref="SparseShiftInvertArnoldi.Run"/> for top-K slow-mode extraction.
///
/// <para>The sparsity pattern is structurally the same as in the dense
/// <see cref="KleinFourGroupSelfPairedRefinement.BuildSubBlockL"/>:
/// L_χ[α, β] = (1/√(|O_α|·|O_β|)) · Σ_{g_α, g_β} χ(g_α)·χ(g_β) · L[g_α·x_α, g_β·x_β].
/// The underlying L matrix-element formula (chain XY + uniform Z-dephasing) is
/// element-wise sparse with O(N) non-zeros per (i, j) row in the computational basis,
/// inflated by a constant Klein-orbit factor — ≤ 4·4·O(N) ≈ 64 per sub-block row at N=10.
/// Empirical mean nnz/row is ~10, reflecting partial cancellation in the character sum.</para>
///
/// <para><b>Tier outcome: <see cref="Tier.Tier1Derived"/>.</b> Element-wise reconstruction
/// of the same operator that <see cref="KleinFourGroupSelfPairedRefinement.BuildSubBlockL"/>
/// builds densely; validated against the dense reference at small N (4..8) to relative
/// Frobenius residual &lt; 1e-10 in the test suite.</para>
///
/// <para>Anchor: <see cref="KleinFourGroupSelfPairedRefinement"/> (orbit + character source)
/// + <see cref="JordanWigner.JwSlaterPairSparseLBuilder"/> (CSR convention template).</para>
/// </summary>
public sealed class KleinFourGroupSelfPairedSparseLBuilder : Claim
{
    /// <summary>Magnitude threshold for inclusion in the CSR sparse pattern; entries with
    /// |v| ≤ <see cref="NonzeroThreshold"/> are dropped to suppress FP round-off noise.</summary>
    public const double NonzeroThreshold = 1e-12;

    public int N { get; }
    public KleinCharacter Character { get; }
    public IReadOnlyList<double> GammaPerSite { get; }
    public int SectorDim { get; }

    public int[] RowPtr { get; }
    public int[] ColIdx { get; }
    public Complex[] Values { get; }

    public long NnzTotal => Values.LongLength;
    public int MaxNnzPerRow { get; }
    public double MeanNnzPerRow => (double)NnzTotal / SectorDim;

    /// <summary>Build the CSR sparse L sub-block with uniform J = 1. Convenience wrapper that
    /// forwards to the per-bond overload with <c>bondJ = [1, 1, ..., 1]</c> of length N − 1.
    /// Existing callers continue to compile unchanged; pass an explicit <c>bondJ</c> list to
    /// access non-uniform per-bond couplings.</summary>
    public static KleinFourGroupSelfPairedSparseLBuilder Build(int N,
        KleinCharacter character, IReadOnlyList<double> gammaPerSite)
    {
        // Match the per-bond overload's factory guard up-front so `new double[N - 1]` never
        // sees a negative size (would silently throw OverflowException). The bare form then
        // matches sibling KleinFourGroupSelfPairedRefinement.BuildSubBlockL for consistency.
        if (N < 2) throw new ArgumentOutOfRangeException(nameof(N), N, "N must be ≥ 2.");
        var bondJ = new double[N - 1];
        for (int b = 0; b < bondJ.Length; b++) bondJ[b] = 1.0;
        return Build(N, character, gammaPerSite, bondJ);
    }

    /// <summary>Build the CSR sparse L sub-block with per-bond J profile. The underlying
    /// chain-XY Hamiltonian is <c>H = Σ_b (J_b/2)·(X_b X_{b+1} + Y_b Y_{b+1})</c>;
    /// <paramref name="bondJ"/> must have length N − 1. Use this overload for F100-territory
    /// experiments; uniform J = 1 callers can use the scalar overload
    /// <see cref="Build(int, KleinCharacter, IReadOnlyList{double})"/>.</summary>
    public static KleinFourGroupSelfPairedSparseLBuilder Build(int N,
        KleinCharacter character, IReadOnlyList<double> gammaPerSite,
        IReadOnlyList<double> bondJ)
    {
        if (N < 2) throw new ArgumentOutOfRangeException(nameof(N), N, "N must be ≥ 2.");
        if ((N & 1) != 0) throw new ArgumentException(
            $"N must be even (self-paired sector requires N/2 ∈ ℤ); got N={N}.", nameof(N));
        if (gammaPerSite is null) throw new ArgumentNullException(nameof(gammaPerSite));
        if (gammaPerSite.Count != N)
            throw new ArgumentException(
                $"gammaPerSite length {gammaPerSite.Count} != N {N}", nameof(gammaPerSite));
        if (bondJ is null) throw new ArgumentNullException(nameof(bondJ));
        if (bondJ.Count != N - 1)
            throw new ArgumentException(
                $"bondJ length {bondJ.Count} != N - 1 = {N - 1}", nameof(bondJ));

        var refinement = KleinFourGroupSelfPairedRefinement.Build(N);
        var ordered = refinement.Orbits
            .Where(o => o.SurvivingCharacters.Contains(character))
            .ToArray();
        int dim = ordered.Length;

        int[] bondFlip = Enumerable.Range(0, N - 1).Select(b => (1 << b) | (1 << (b + 1))).ToArray();

        // Pre-index orbit members for fast (row, col) → orbit-index lookup.
        var memberToOrbitIdx = new Dictionary<(int, int), (int orbitIdx, int gIdx)>(dim * 4);
        for (int orbitIdx = 0; orbitIdx < dim; orbitIdx++)
        {
            var o = ordered[orbitIdx];
            for (int g = 0; g < o.Size; g++)
                memberToOrbitIdx[o.Members[g]] = (orbitIdx, g);
        }

        var rowPtr = new int[dim + 1];
        var colIdxList = new List<int>(dim * 16);
        var valuesList = new List<Complex>(dim * 16);
        int maxNnz = 0;

        var perRow = new SortedDictionary<int, Complex>();

        for (int alpha = 0; alpha < dim; alpha++)
        {
            var oRow = ordered[alpha];
            double invSqrtSizeRow = 1.0 / Math.Sqrt(oRow.Size);
            perRow.Clear();

            // Walk every L-coupled (rowβ, colβ) reachable from each (rowα, colα) = oRow.Members[gRow].
            // For chain XY + Z-deph: L[(rowα, colα), (rowβ, colβ)] is nonzero only when:
            //   (a) colα == colβ AND rowα = rowβ ⊕ bondFlip[b] (ket-side hop), value −i·J
            //   (b) rowα == rowβ AND colα = colβ ⊕ bondFlip[b] (bra-side hop), value +i·J
            //   (c) rowα == rowβ AND colα == colβ                 (diagonal),  value −2·Σ_l γ_l·(n_l^row ⊕ n_l^col)
            // For each non-zero L element, find which orbit β contains (rowβ, colβ) and
            // accumulate its contribution to L_χ[α, β] with character sign χ(g_α)·χ(g_β).
            for (int gRow = 0; gRow < oRow.Size; gRow++)
            {
                double chiRow = ChiSign(character, gRow, oRow.GroupElementIndices);
                var (rowA, colA) = oRow.Members[gRow];

                // (a) ket-side hops: rowβ = rowA ⊕ bondFlip[b], colβ = colA
                for (int b = 0; b < N - 1; b++)
                {
                    int bitB = (rowA >> b) & 1;
                    int bitBp1 = (rowA >> (b + 1)) & 1;
                    if (bitB == bitBp1) continue;
                    int rowBeta = rowA ^ bondFlip[b];
                    AccumulateCoupling(rowBeta, colA, new Complex(0, -bondJ[b]), chiRow, perRow,
                        memberToOrbitIdx, ordered, character, invSqrtSizeRow);
                }

                // (b) bra-side hops: rowβ = rowA, colβ = colA ⊕ bondFlip[b]
                for (int b = 0; b < N - 1; b++)
                {
                    int bitB = (colA >> b) & 1;
                    int bitBp1 = (colA >> (b + 1)) & 1;
                    if (bitB == bitBp1) continue;
                    int colBeta = colA ^ bondFlip[b];
                    AccumulateCoupling(rowA, colBeta, new Complex(0, +bondJ[b]), chiRow, perRow,
                        memberToOrbitIdx, ordered, character, invSqrtSizeRow);
                }

                // (c) diagonal dissipator: rowβ = rowA, colβ = colA
                int xorIJ = rowA ^ colA;
                double diss = 0.0;
                for (int l = 0; l < N; l++)
                    if (((xorIJ >> l) & 1) != 0) diss -= 2 * gammaPerSite[l];
                if (diss != 0.0)
                    AccumulateCoupling(rowA, colA, new Complex(diss, 0), chiRow, perRow,
                        memberToOrbitIdx, ordered, character, invSqrtSizeRow);
            }

            int rowNnz = 0;
            foreach (var kv in perRow)
            {
                if (kv.Value.Magnitude <= NonzeroThreshold) continue;
                colIdxList.Add(kv.Key);
                valuesList.Add(kv.Value);
                rowNnz++;
            }
            if (rowNnz > maxNnz) maxNnz = rowNnz;
            rowPtr[alpha + 1] = colIdxList.Count;
        }

        return new KleinFourGroupSelfPairedSparseLBuilder(N, character, gammaPerSite.ToArray(),
            dim, rowPtr, colIdxList.ToArray(), valuesList.ToArray(), maxNnz);
    }

    /// <summary>Accumulate a single L-coupling into <c>perRow[βIndex] += chiRow · chiCol · value · invSqrtSize</c>
    /// after mapping (rowβ, colβ) to its orbit β. If the orbit's surviving characters don't
    /// include χ, the contribution is zero.</summary>
    private static void AccumulateCoupling(int rowBeta, int colBeta, Complex value, double chiRow,
        SortedDictionary<int, Complex> perRow,
        Dictionary<(int, int), (int orbitIdx, int gIdx)> memberToOrbitIdx,
        KleinOrbit[] ordered,
        KleinCharacter character, double invSqrtSizeRow)
    {
        if (!memberToOrbitIdx.TryGetValue((rowBeta, colBeta), out var orbitInfo)) return;
        // Find orbit β within the ordered (χ-surviving) list; if its orbit isn't in this
        // character's sub-block, skip.
        int orbitGlobalIdx = orbitInfo.orbitIdx;
        if (orbitGlobalIdx < 0 || orbitGlobalIdx >= ordered.Length) return;
        var oCol = ordered[orbitGlobalIdx];
        // The orbit could be in the surviving list at a different position than orbitGlobalIdx
        // (memberToOrbitIdx indexes against the ORDERED array, so orbitInfo.orbitIdx IS the
        // position in ordered[]).
        double chiCol = ChiSign(character, orbitInfo.gIdx, oCol.GroupElementIndices);
        double invSqrtSizeCol = 1.0 / Math.Sqrt(oCol.Size);
        Complex contrib = value * (chiRow * chiCol * invSqrtSizeRow * invSqrtSizeCol);
        if (perRow.TryGetValue(orbitGlobalIdx, out var existing))
            perRow[orbitGlobalIdx] = existing + contrib;
        else
            perRow[orbitGlobalIdx] = contrib;
    }

    private static double ChiSign(KleinCharacter chi, int gIndex, IReadOnlyList<KleinGroupElement> groupElements)
    {
        var g = groupElements[gIndex];
        bool flipF71 = (g & KleinGroupElement.F71) == KleinGroupElement.F71;
        bool flipXN = (g & KleinGroupElement.XN) == KleinGroupElement.XN;
        int sign = 1;
        if (flipF71 && (chi == KleinCharacter.MinusPlus || chi == KleinCharacter.MinusMinus)) sign = -sign;
        if (flipXN && (chi == KleinCharacter.PlusMinus || chi == KleinCharacter.MinusMinus)) sign = -sign;
        return sign;
    }

    private KleinFourGroupSelfPairedSparseLBuilder(int n, KleinCharacter character,
        double[] gammaPerSite, int dim, int[] rowPtr, int[] colIdx, Complex[] values, int maxNnz)
        : base($"Klein sub-block CSR sparse L on self-paired (N/2, N/2) sector at N={n}, " +
               $"character {character}: dim {dim}, nnz {values.LongLength}, max-nnz/row {maxNnz}, " +
               $"mean-nnz/row {(double)values.LongLength / dim:F2}.",
               Tier.Tier1Derived,
               "compute/RCPsiSquared.Core/BlockSpectrum/KleinFourGroupSelfPairedRefinement.cs (orbit + character source; dense reference for validation) + " +
               "compute/RCPsiSquared.Core/BlockSpectrum/JordanWigner/JwSlaterPairSparseLBuilder.cs (CSR convention template); " +
               "chain XY + Z-dephasing Liouvillian L matrix element formula.")
    {
        N = n; Character = character;
        GammaPerSite = gammaPerSite;
        SectorDim = dim;
        RowPtr = rowPtr;
        ColIdx = colIdx;
        Values = values;
        MaxNnzPerRow = maxNnz;
    }

    public override string DisplayName =>
        $"Klein sparse L sub-block χ={Character} (N={N}, dim {SectorDim}, nnz {NnzTotal:N0})";

    public override string Summary =>
        $"χ={Character}, dim={SectorDim}, nnz={NnzTotal:N0}, max-nnz/row={MaxNnzPerRow}, " +
        $"mean={MeanNnzPerRow:F1} ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return InspectableNode.RealScalar("N", N);
            yield return new InspectableNode("character", summary: Character.ToString());
            yield return InspectableNode.RealScalar("sector dim", SectorDim);
            yield return InspectableNode.RealScalar("nnz total", NnzTotal);
            yield return InspectableNode.RealScalar("max nnz/row", MaxNnzPerRow);
            yield return InspectableNode.RealScalar("mean nnz/row", MeanNnzPerRow, "F2");
        }
    }
}
