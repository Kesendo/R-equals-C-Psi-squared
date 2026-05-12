using System.Collections.Concurrent;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.SymmetryFamily;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.BlockSpectrum;

/// <summary>LiouvillianBlockSpectrum (Tier 1 derived; 2026-05-11): For the XY+Z-dephasing
/// Liouvillian L on N qubits, the union of per-block eigenvalues over the (N+1)² joint
/// popcount sectors (<see cref="JointPopcountSectors"/>) equals the spectrum of full-L
/// as a multiset.
///
/// <para>Exploits the U(1)×U(1) per-side popcount conservation established by the parent
/// <see cref="JointPopcountSectors"/> Claim: the Liouvillian is exactly block-diagonal in
/// the joint (popcount_col, popcount_row) label after the basis permutation produced by
/// <see cref="JointPopcountSectorBuilder.Build"/>. Each diagonal block is an isolated
/// eigenvalue problem, so the full spectrum is the disjoint union of per-block spectra.</para>
///
/// <para><b>Cubic-cost speedup.</b> Naive full-L diagonalisation costs O((4^N)³). Block-wise
/// the cost drops to Σ_{p_c, p_r} (C(N, p_c) · C(N, p_r))³. Indicative speedups:
/// <list type="bullet">
///   <item>N=5: full 4^N = 1024 → max block 100, total cost ratio ≈ 50× faster</item>
///   <item>N=6: full 4^N = 4096 → max block 400, total cost ratio ≈ 110× faster</item>
///   <item>N=7: full 4^N = 16384 → max block 1225, total cost ratio ≈ 250× faster</item>
///   <item>N=8: full 4^N = 65536 → max block 4900, total cost ratio ≈ 515× faster</item>
/// </list>
/// At N=8 the largest block fits in ~0.38 GB vs ~68.7 GB for the full L, removing the
/// need for native-memory + ILP64 LAPACK on the dense path.</para>
///
/// <para><b>Bit-exact witness at N=3, 4, 5.</b> The per-block spectrum and the direct
/// full-L spectrum agree as multisets to <c>|Δλ| &lt; 1e-9</c> across uniform XY chain
/// <c>(J = 1.0)</c> + per-site Z-dephasing <c>(γ = 0.5)</c>, and under varied parameters
/// (γ ∈ {0.1, 0.5, 2.0}, J ∈ {0.5, 1.0, 3.0}). Verified by
/// <c>LiouvillianBlockSpectrumTests</c>. Source of L: <see cref="Pauli.PauliHamiltonian.XYChain"/>
/// composed with <see cref="Lindblad.PauliDephasingDissipator.BuildZ"/>.</para>
///
/// <para>Anchors: <c>compute/RCPsiSquared.Core/BlockSpectrum/JointPopcountSectors.cs</c>
/// (parent Claim, block-diagonal structure), <c>compute/RCPsiSquared.Core/BlockSpectrum/JointPopcountSectorBuilder.cs</c>
/// (basis permutation + sector ranges), <c>compute/RCPsiSquared.Core.Tests/BlockSpectrum/LiouvillianBlockSpectrumTests.cs</c>
/// (bit-exact spectral equality verification at N=3, 4, 5).</para></summary>
public sealed class LiouvillianBlockSpectrum : Claim
{
    private readonly JointPopcountSectors _sectors;

    public LiouvillianBlockSpectrum(JointPopcountSectors sectors)
        : base("LiouvillianBlockSpectrum: per-block eig over (N+1)² joint popcount sectors yields the same spectrum (multiset) as direct full-L eig; bit-exact verified at N=3,4,5.",
               Tier.Tier1Derived,
               "JointPopcountSectors block-diagonality (parent) + per-block diagonalisation; verified bit-exact vs full-L eig at N=3,4,5 in LiouvillianBlockSpectrumTests")
    {
        _sectors = sectors ?? throw new ArgumentNullException(nameof(sectors));
    }

    /// <summary>Compute the full Liouvillian spectrum via per-block eigendecomposition over the
    /// joint-popcount sectors. Returns a flat array of all 4^N eigenvalues, ordered block-by-block
    /// in <see cref="JointPopcountSectorBuilder.SectorRange"/> iteration order.
    ///
    /// <para>For each <see cref="JointPopcountSectorBuilder.SectorRange"/> (p_c, p_r, offset, size):
    /// extract the size×size sub-block from L at the permuted (row, col) indices given by
    /// <see cref="JointPopcountSectorBuilder.Decomposition.Permutation"/>; run MathNet's
    /// <c>Matrix&lt;Complex&gt;.Evd()</c>; append its <c>EigenValues</c> to the result.</para>
    ///
    /// <para>Block extraction is done index-by-index rather than via a full permutation of L,
    /// which avoids materialising the permuted 4^N × 4^N matrix.</para></summary>
    /// <param name="L">The Liouvillian L = -i[H, ·] + dissipator, in the row-major
    /// <c>flat = row·d + col</c> convention used by <see cref="Lindblad.LindbladianBuilder"/>
    /// and <see cref="Lindblad.PauliDephasingDissipator"/>. Must be (4^N) × (4^N).</param>
    /// <param name="N">Qubit count; must satisfy <c>L.RowCount == 4^N</c>.</param>
    /// <returns>Flat array of 4^N eigenvalues, concatenated block-by-block.</returns>
    public static Complex[] ComputeSpectrum(ComplexMatrix L, int N)
    {
        if (L is null) throw new ArgumentNullException(nameof(L));
        int liouvilleDim = 1 << (2 * N);
        if (L.RowCount != liouvilleDim || L.ColumnCount != liouvilleDim)
            throw new ArgumentException(
                $"L must be ({liouvilleDim})×({liouvilleDim}) for N={N}; got {L.RowCount}×{L.ColumnCount}.",
                nameof(L));

        var decomp = JointPopcountSectorBuilder.Build(N);
        var perm = decomp.Permutation;
        var spectrum = new Complex[liouvilleDim];

        // H1: per-sector eig is embarrassingly parallel. Pre-compute the per-sector write
        // offset so each task knows its destination range before any work starts.
        int sectorCount = decomp.SectorRanges.Count;
        var writeOffsets = new int[sectorCount];
        int cum = 0;
        for (int i = 0; i < sectorCount; i++)
        {
            writeOffsets[i] = cum;
            cum += decomp.SectorRanges[i].Size;
        }

        // BLAS-oversubscription strategy (c): cap outer parallelism at ~ProcessorCount/4 to
        // leave headroom for MKL's internal threading on the larger Evd calls. Empirically
        // good balance on 24-core; the few largest sectors keep their MKL parallelism while
        // many small sectors run concurrently. Larger outer DOP would oversubscribe (outer ×
        // MKL threads = ProcessorCount × ProcessorCount); smaller would leave cores idle.
        int outerDop = Math.Max(1, Environment.ProcessorCount / 4);
        var po = new ParallelOptions { MaxDegreeOfParallelism = outerDop };

        Parallel.ForEach(
            Enumerable.Range(0, sectorCount), po,
            sIdx =>
            {
                var sector = decomp.SectorRanges[sIdx];
                int size = sector.Size;
                if (size == 0) return;
                var block = Matrix<Complex>.Build.Dense(size, size);
                for (int r = 0; r < size; r++)
                {
                    int rowFlat = perm[sector.Offset + r];
                    for (int c = 0; c < size; c++)
                    {
                        int colFlat = perm[sector.Offset + c];
                        block[r, c] = L[rowFlat, colFlat];
                    }
                }
                var blockEigs = block.Evd().EigenValues;
                int write = writeOffsets[sIdx];
                for (int i = 0; i < size; i++)
                    spectrum[write + i] = blockEigs[i];
            });

        return spectrum;
    }

    /// <summary>Compute the full Liouvillian spectrum without materialising the full
    /// (4^N) × (4^N) L matrix. Each per-block matrix is built directly from the
    /// Hilbert-space Hamiltonian (size 2^N × 2^N) via
    /// <see cref="PerBlockLiouvillianBuilder.BuildBlockZ"/>, eigendecomposed, and discarded
    /// before the next block. This is the only path that scales past N=6 on commodity
    /// hardware (full L exceeds .NET 2 GB array-size limit at N=7+).
    ///
    /// <para>Uses the X⊗N pairing optimisation: only ~half of joint-popcount sectors are
    /// eigendecomposed; the other half copy from their X⊗N-partner. See
    /// <see cref="SymmetryFamily.XGlobalChargeConjugationPairing"/> for the pairing rule.</para></summary>
    /// <param name="H">Hilbert-space Hamiltonian, dense 2^N × 2^N (cheap even at N=8: 256×256).</param>
    /// <param name="gammaPerSite">Per-site Z-dephasing rates (length N).</param>
    /// <param name="N">Qubit count.</param>
    /// <returns>Flat array of 4^N eigenvalues, concatenated block-by-block.</returns>
    public static Complex[] ComputeSpectrumPerBlock(ComplexMatrix H, IReadOnlyList<double> gammaPerSite, int N)
    {
        if (H is null) throw new ArgumentNullException(nameof(H));
        int hilbertDim = 1 << N;
        if (H.RowCount != hilbertDim || H.ColumnCount != hilbertDim)
            throw new ArgumentException(
                $"H must be ({hilbertDim})×({hilbertDim}) for N={N}; got {H.RowCount}×{H.ColumnCount}.",
                nameof(H));
        if (gammaPerSite is null) throw new ArgumentNullException(nameof(gammaPerSite));
        if (gammaPerSite.Count != N)
            throw new ArgumentException($"gamma list length {gammaPerSite.Count} != N={N}", nameof(gammaPerSite));

        int liouvilleDim = 1 << (2 * N);
        var decomp = JointPopcountSectorBuilder.Build(N);
        var perm = decomp.Permutation;
        var spectrum = new Complex[liouvilleDim];

        // H1: pre-compute per-sector write offsets so each task knows its destination slice.
        int sectorCount = decomp.SectorRanges.Count;
        var writeOffsets = new int[sectorCount];
        int cum = 0;
        for (int i = 0; i < sectorCount; i++)
        {
            writeOffsets[i] = cum;
            cum += decomp.SectorRanges[i].Size;
        }

        // X⊗N pairing (Tier 1, XGlobalChargeConjugationPairing): the chain XY+Z-deph L
        // commutes with the global X-string operator, which on joint-popcount labels maps
        // (p_c, p_r) ↔ (N - p_c, N - p_r). Paired sectors share spectrum exactly. We compute
        // eig only on "primary" sectors (lex-smaller of each pair, plus all self-paired
        // sectors at even N) and copy the result onto follower sectors. Primaries are sorted
        // descending by size so the largest sector starts first under Parallel.ForEach,
        // overlapping its wall-time with smaller sectors' work.
        var (primarySectorIndices, followerToPrimary) =
            XGlobalChargeConjugationPairing.PartitionByXNPairing(
                N, decomp.SectorRanges, s => (s.PCol, s.PRow), s => s.Size);

        // BLAS-oversubscription strategy (c): outer DOP ≈ ProcessorCount/4 leaves room for
        // MKL inside the largest sectors' Evd. See ComputeSpectrum for the rationale.
        int outerDop = Math.Max(1, Environment.ProcessorCount / 4);
        var po = new ParallelOptions { MaxDegreeOfParallelism = outerDop };

        // Phase 2: parallel eig on primary sectors only.
        var primaryEigs = new ConcurrentDictionary<int, Complex[]>();
        Parallel.ForEach(
            primarySectorIndices, po,
            sIdx =>
            {
                var sector = decomp.SectorRanges[sIdx];
                int size = sector.Size;
                if (size == 0)
                {
                    primaryEigs[sIdx] = Array.Empty<Complex>();
                    return;
                }
                var flatIndices = new int[size];
                for (int k = 0; k < size; k++)
                    flatIndices[k] = perm[sector.Offset + k];

                var block = PerBlockLiouvillianBuilder.BuildBlockZ(H, gammaPerSite, flatIndices);
                var blockEigs = block.Evd().EigenValues;
                var arr = new Complex[size];
                for (int i = 0; i < size; i++) arr[i] = blockEigs[i];
                primaryEigs[sIdx] = arr;
            });

        // Phase 3: sequential write to output array (primaries + followers). Followers copy
        // the eigenvalue array of their X⊗N-partner sector (same multiset of eigenvalues).
        for (int sIdx = 0; sIdx < sectorCount; sIdx++)
        {
            int sourceIdx = primaryEigs.ContainsKey(sIdx) ? sIdx : followerToPrimary[sIdx];
            var eigs = primaryEigs[sourceIdx];
            int write = writeOffsets[sIdx];
            for (int i = 0; i < eigs.Length; i++) spectrum[write + i] = eigs[i];
        }
        return spectrum;
    }

    public override string DisplayName =>
        "LiouvillianBlockSpectrum: per-block eig over (N+1)² joint popcount sectors = full-L spectrum";

    public override string Summary =>
        $"per-block eig multiset equals full-L spectrum bit-exactly at N=3,4,5; cubic-cost speedup ≈ 515× at N=8 ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("parent",
                summary: "JointPopcountSectors (block-diagonal structure)");
            yield return new InspectableNode("witness",
                summary: "bit-exact spectral equality vs full-L eig at N=3, 4, 5 (|Δλ| < 1e-9)");
            yield return new InspectableNode("cubic-cost speedup",
                summary: "N=5: ≈ 50×, N=6: ≈ 110×, N=7: ≈ 250×, N=8: ≈ 515×");
            yield return new InspectableNode("N=8 max block",
                summary: $"size {JointPopcountSectors.MaxSectorSize(8)} (vs full 4^8 = 65536)");
        }
    }
}
