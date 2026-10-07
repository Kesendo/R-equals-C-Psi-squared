using System.Collections.Concurrent;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Core.SymmetryFamily;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.BlockSpectrum;

/// <summary>F71MirrorBlockRefinement (Tier 1 derived; 2026-05-11):
/// the chain-topology Liouvillian L commutes with the spatial-mirror operator
/// <c>P_F71 ⊗ P_F71</c> on Liouville space (P_F71: site b ↔ site N−1−b on each Hilbert
/// side) when H and the dephasing rates are symmetric under that reflection. This Z₂
/// symmetry refines each joint-popcount sector from
/// <see cref="JointPopcountSectors"/> into F71-even and F71-odd sub-blocks (≈ factor-2
/// further block-size reduction for non-fixed-point orbits).
///
/// <para>In each (p_c, p_r) sector, Liouville-space basis pairs (col, row) split into
/// F71-orbits of size 1 (fixed points: col and row are both spatially-palindromic bit
/// patterns) and size 2 (swapped pairs). The F71-even sub-block is spanned by
/// <c>(s + P_F71 s)/√2</c> for size-2 orbit representatives plus all size-1 fixed points;
/// the F71-odd sub-block by <c>(s − P_F71 s)/√2</c> from size-2 orbits only. Fixed points
/// have no F71-odd component.</para>
///
/// <para>The basis change Q is a real orthogonal matrix with entries 0, ±1, ±1/√2;
/// <c>Q^T L Q</c> is then block-diagonal in the refined sector layout. Verified at
/// N=3, 4 chain XY+Z-deph (off-block Frobenius below 1e-10).</para>
///
/// <para><b>N=8 projection.</b> Without F71 refinement, the largest joint-popcount
/// sector at N=8 is C(8,4)² = 4900. With F71 refinement, the F71-even and F71-odd halves
/// roughly halve this; the exact split depends on how many of the 4900 basis pairs are
/// F71-fixed (col, row both palindromic), which biases the even half slightly larger.</para>
///
/// <para>Distinct from <see cref="F86.F71MirrorInvariance"/>, which is a different Claim
/// about Q_peak(b) = Q_peak(N−2−b) for the F86 c=2 observable. This Claim is purely
/// about *spectral block-decomposition* of L.</para>
///
/// <para>Anchors: <c>compute/RCPsiSquared.Core/BlockSpectrum/JointPopcountSectorBuilder.cs</c>
/// (parent decomposition), <c>compute/RCPsiSquared.Core.Tests/BlockSpectrum/F71MirrorBlockRefinementTests.cs</c>
/// (off-block Frobenius below 1e-10 + spectrum match within 1e-9 at N=3, 4).</para></summary>
public sealed class F71MirrorBlockRefinement : Claim
{
    private readonly JointPopcountSectors _sectors;

    public F71MirrorBlockRefinement(JointPopcountSectors sectors)
        : base("F71MirrorBlockRefinement: chain L with reflection-symmetric H and rates commutes with P_F71 ⊗ P_F71 (site b ↔ N−1−b on each Hilbert side); each joint-popcount sector splits into F71-even + F71-odd sub-blocks; verified at N=3, 4 (off-block Frobenius below 1e-10).",
               Tier.Tier1Derived,
               "JointPopcountSectors block-diagonality (parent) + chain spatial-mirror Z₂ symmetry; verified off-block Frobenius < 1e-10 in F71MirrorBlockRefinementTests")
    {
        _sectors = sectors ?? throw new ArgumentNullException(nameof(sectors));
    }

    public override string DisplayName =>
        "F71MirrorBlockRefinement: Z₂ split of joint-popcount sectors into F71-even/odd sub-blocks";

    public override string Summary =>
        $"chain L with reflection-symmetric H and rates commutes with P_F71⊗P_F71; each (p_c,p_r) sector → (even, odd) sub-blocks; ≈ 2× further block-size reduction ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("parent",
                summary: "JointPopcountSectors (joint-popcount block-diagonal structure)");
            yield return new InspectableNode("symmetry",
                summary: "P_F71 ⊗ P_F71 acts on Liouville space; basis change Q is real orthogonal");
            yield return new InspectableNode("sub-blocks",
                summary: "F71-even = (s + P_F71 s)/√2 + fixed-points; F71-odd = (s − P_F71 s)/√2 from size-2 orbits");
            yield return new InspectableNode("witness",
                summary: "off-block Frobenius < 1e-10 in Q^T L Q at N=3, 4 chain (XY+Z-deph)");
        }
    }

    /// <summary>One refined sub-block of L after F71 refinement of a joint-popcount sector.
    /// Spatial-mirror parity is +1 ('+') for the F71-even component or −1 ('−') for F71-odd.
    /// Empty sub-blocks (e.g. the odd component of a sector consisting entirely of F71-fixed
    /// pairs) are still emitted with <c>Size = 0</c> for structural completeness; the
    /// <see cref="F71RefinedDecomposition"/> contract guarantees that all sectors with size
    /// &gt; 0 appear in offset order with no gaps.</summary>
    public sealed record F71RefinedSector(int PCol, int PRow, char Parity, int Offset, int Size);

    /// <summary>Refined block-diagonal decomposition: a single real-orthogonal basis change
    /// Q maps L, when H and the rates are symmetric under the reflection, into a matrix that
    /// is block-diagonal in the (p_c, p_r, parity) labels.
    /// <c>BasisChange.Transpose() * L * BasisChange</c> has the structure described by
    /// <see cref="SectorRanges"/>. Empty sub-blocks (Size = 0) may appear when a sector
    /// has no F71-odd component (all basis pairs are F71-fixed); these are listed but do
    /// not occupy any rows/cols in the basis-change image.</summary>
    public sealed class F71RefinedDecomposition
    {
        public int N { get; }
        public int D { get; }   // 2^N
        public ComplexMatrix BasisChange { get; }   // unitary Q (real orthogonal); Q^T L Q block-diagonal for reflection-symmetric H and rates
        public IReadOnlyList<F71RefinedSector> SectorRanges { get; }
        public F71RefinedDecomposition(int n, int d, ComplexMatrix q, IReadOnlyList<F71RefinedSector> sectors)
        {
            N = n; D = d; BasisChange = q; SectorRanges = sectors;
        }
    }

    /// <summary>Refine a <see cref="JointPopcountSectorBuilder.Decomposition"/> by the
    /// chain spatial-mirror Z₂ symmetry. Produces an orthogonal basis change Q whose columns
    /// are organised so that <c>Q^T L Q</c> is block-diagonal in (p_c, p_r, parity) whenever H
    /// and the rates are reflection-symmetric.
    ///
    /// <para>For each input sector, the (col, row) flat indices in that sector are partitioned
    /// into F71-orbits via the involution <c>flat ↔ mirror(flat) := P_F71(row)·d + P_F71(col)</c>
    /// (since flat = row·d + col). Size-1 orbits (fixed points) contribute one even-parity
    /// basis vector each (the standard basis vector itself). Size-2 orbits contribute one even
    /// vector <c>(e_s + e_{P s})/√2</c> and one odd vector <c>(e_s − e_{P s})/√2</c>.</para>
    ///
    /// <para>Sub-block ordering: for each input sector (in the original
    /// <see cref="JointPopcountSectorBuilder.SectorRange"/> order), emit the F71-even sub-block
    /// first then the F71-odd sub-block. Empty sub-blocks (Size = 0) are still recorded for
    /// structural visibility but consume no columns of Q.</para></summary>
    public static F71RefinedDecomposition RefineWithF71(JointPopcountSectorBuilder.Decomposition baseDecomp)
    {
        if (baseDecomp is null) throw new ArgumentNullException(nameof(baseDecomp));
        int N = baseDecomp.N;
        int d = baseDecomp.D;
        int liouvilleDim = d * d;

        // Build the per-Hilbert-side P_F71 mirror map: a |b_0...b_{N-1}⟩ ↔ |b_{N-1}...b_0⟩.
        var mirrorBits = F71MirrorIndexHelper.BuildHilbertMirrorLookup(N);

        // Liouville-space mirror: flat = row*d + col → mirror_flat = mirror(row)*d + mirror(col).
        int Mirror(int flat) => mirrorBits[flat / d] * d + mirrorBits[flat % d];

        // Per-sector orbit decomposition + Q-column emission.
        var Q = Matrix<Complex>.Build.Sparse(liouvilleDim, liouvilleDim);
        var refinedSectors = new List<F71RefinedSector>();
        int writeCol = 0;

        foreach (var sector in baseDecomp.SectorRanges)
        {
            // Collect flat indices in this sector and find F71-orbits within it.
            // Note: F71 preserves popcount per Hilbert side, so mirror(flat) lands in the same
            // (p_c, p_r) sector. We work with the global flat indices directly (not permuted).
            int size = sector.Size;
            var sectorIndices = new int[size];
            for (int k = 0; k < size; k++)
                sectorIndices[k] = baseDecomp.Permutation[sector.Offset + k];

            // Walk indices in ascending flat order; for each unseen flat, partner with
            // Mirror(flat) (which is also in this sector). Fixed points: Mirror(flat)==flat.
            var (fixedPoints, pairs) = F71MirrorIndexHelper.FindOrbitsInSector(sectorIndices, Mirror);

            // F71-even sub-block: fixed points + (e_s + e_{ps})/√2 for each pair.
            int evenStart = writeCol;
            foreach (int fp in fixedPoints)
            {
                Q[fp, writeCol] = Complex.One;
                writeCol++;
            }
            double invSqrt2 = 1.0 / Math.Sqrt(2.0);
            foreach (var (s, ps) in pairs)
            {
                Q[s, writeCol] = invSqrt2;
                Q[ps, writeCol] = invSqrt2;
                writeCol++;
            }
            int evenSize = writeCol - evenStart;
            refinedSectors.Add(new F71RefinedSector(sector.PCol, sector.PRow, '+', evenStart, evenSize));

            // F71-odd sub-block: (e_s − e_{ps})/√2 for each pair (no fixed-point contribution).
            int oddStart = writeCol;
            foreach (var (s, ps) in pairs)
            {
                Q[s, writeCol] = invSqrt2;
                Q[ps, writeCol] = -invSqrt2;
                writeCol++;
            }
            int oddSize = writeCol - oddStart;
            refinedSectors.Add(new F71RefinedSector(sector.PCol, sector.PRow, '-', oddStart, oddSize));
        }

        if (writeCol != liouvilleDim)
            throw new InvalidOperationException(
                $"F71 refinement produced {writeCol} basis vectors, expected {liouvilleDim} (4^N).");

        return new F71RefinedDecomposition(N, d, Q, refinedSectors);
    }

    /// <summary>Compute the eigenvalues of L's F71-refined sub-blocks, its full spectrum when
    /// L commutes with the site reflection, via per-sub-block eigendecomposition over the
    /// F71-refined sectors. Returns a flat array of all 4^N eigenvalues, ordered
    /// sub-block-by-sub-block in <see cref="F71RefinedSector"/> iteration order.
    ///
    /// <para>Block extraction uses the basis change Q implicitly: the sub-block
    /// <c>(Q^T L Q)[off:off+size, off:off+size]</c> is materialised explicitly via
    /// <c>Q^T_{cols=off..} · L · Q_{cols=off..}</c> on the relevant column slice. Empty
    /// sub-blocks (Size = 0) contribute no eigenvalues.</para></summary>
    /// <param name="L">The Liouvillian L = -i[H, ·] + dissipator, in the row-major
    /// <c>flat = row·d + col</c> convention used by <see cref="Lindblad.LindbladianBuilder"/>
    /// and <see cref="Lindblad.PauliDephasingDissipator"/>. Must be (4^N) × (4^N).</param>
    /// <param name="N">Qubit count; must satisfy <c>L.RowCount == 4^N</c>.</param>
    /// <returns>Flat array of 4^N eigenvalues, concatenated sub-block-by-sub-block.</returns>
    public static Complex[] ComputeSpectrum(ComplexMatrix L, int N)
    {
        if (L is null) throw new ArgumentNullException(nameof(L));
        int liouvilleDim = 1 << (2 * N);
        if (L.RowCount != liouvilleDim || L.ColumnCount != liouvilleDim)
            throw new ArgumentException(
                $"L must be ({liouvilleDim})×({liouvilleDim}) for N={N}; got {L.RowCount}×{L.ColumnCount}.",
                nameof(L));

        var baseDecomp = JointPopcountSectorBuilder.Build(N);
        var refined = RefineWithF71(baseDecomp);

        var spectrum = new Complex[liouvilleDim];
        int write = 0;
        foreach (var sub in refined.SectorRanges)
        {
            if (sub.Size == 0) continue;
            // Materialise the sub-block: B = Q_sub^T · L · Q_sub, where Q_sub is the
            // (liouvilleDim × sub.Size) slice of refined.BasisChange starting at column sub.Offset.
            var qSub = refined.BasisChange.SubMatrix(0, liouvilleDim, sub.Offset, sub.Size);
            var block = qSub.Transpose() * L * qSub;
            var blockEigs = block.Evd().EigenValues;
            for (int i = 0; i < sub.Size; i++)
                spectrum[write++] = blockEigs[i];
        }
        return spectrum;
    }

    /// <summary>Compute the spectrum of the F71-refined sub-blocks, L's full spectrum when H
    /// and the rates are symmetric under the site reflection, via per-sub-block
    /// eigendecomposition WITHOUT materialising the full L matrix. For each joint-popcount
    /// sector, builds the (even, odd) F71 sub-blocks directly from the Hilbert-space H +
    /// γ_per_site by:
    /// <list type="number">
    ///   <item>Identifying the sector's F71-orbits (fixed-points + size-2 pairs).</item>
    ///   <item>Constructing the union of paired flat indices (size 2·n_pairs + n_fix).</item>
    ///   <item>Calling <see cref="PerBlockLiouvillianBuilder.BuildBlockZ"/> to get the
    ///         union block.</item>
    ///   <item>Applying the local 1/√2 sign-walk basis change to project onto F71-even
    ///         and F71-odd sub-blocks (a real-orthogonal transform on the union block).</item>
    /// </list>
    /// <para>Memory footprint is per-block: O(blockSize²) only, never O(4^N · 4^N), so it runs
    /// at N=7, 8 where full-L cannot be allocated.</para>
    ///
    /// <para>The result is the union of the eigenvalues of the F71-even and F71-odd diagonal
    /// sub-blocks of every sector. It is the spectrum of L when H and the rates are symmetric
    /// under the site reflection; when they are not, the off-diagonal F71 blocks are nonzero
    /// and the result is the spectrum of the compressed diagonal blocks, the object F91, F92
    /// and F93 compare.</para>
    ///
    /// <para>Uses the F1 Π-orbit pairing when H commutes with X^⊗N: the F1 palindrome
    /// conjugation Π is order-4 and on joint-popcount labels acts as the whole-sector cycle
    /// (p_c, p_r) ↦ (N − p_r, p_c), grouping the (N+1)² sectors into orbits of 4. The map
    /// that carries a sector to its Π-image (ρ ↦ ρ†·X^⊗N, see
    /// <see cref="LiouvillianBlockSpectrum"/>) and conjugation by X^⊗N both commute with the
    /// site reflection, so the pairing respects F71 parity whether or not H and the rates are
    /// F71-symmetric. Only one "primary" sector per orbit is eigendecomposed; the three
    /// followers are derived from it, the Π²-image (X⊗N partner) by a verbatim copy of both the
    /// F71-even and F71-odd sub-arrays and the Π/Π³-images by the F1 reflection λ ↦ −2·Σγ − λ
    /// applied to both. When H does not commute with X^⊗N every sector is solved
    /// (<see cref="LiouvillianBlockSpectrum.SectorPairing"/>).</para></summary>
    /// <param name="H">Hilbert-space Hamiltonian, dense 2^N × 2^N; must be popcount-conserving
    /// (checked exactly, <see cref="ArgumentException"/> otherwise).</param>
    /// <param name="gammaPerSite">Per-site Z-dephasing rates (length N).</param>
    /// <param name="N">Qubit count.</param>
    /// <returns>Flat array of 4^N eigenvalues.</returns>
    public static Complex[] ComputeSpectrumPerBlock(ComplexMatrix H, IReadOnlyList<double> gammaPerSite, int N) =>
        ComputeSpectrumPerBlock(H, gammaPerSite, N, PauliLetter.Z);

    /// <summary>Overload choosing the sector pairing
    /// (<see cref="LiouvillianBlockSpectrum.SectorPairing"/>).</summary>
    public static Complex[] ComputeSpectrumPerBlock(
        ComplexMatrix H, IReadOnlyList<double> gammaPerSite, int N,
        LiouvillianBlockSpectrum.SectorPairing pairing) =>
        ComputeSpectrumPerBlock(H, gammaPerSite, N, PauliLetter.Z, pairing);

    /// <summary>Overload with an explicit dephase letter, as on
    /// <see cref="LiouvillianBlockSpectrum.ComputeSpectrumPerBlock(ComplexMatrix, IReadOnlyList{double}, int, LiouvillianBlockSpectrum.EigenPath, PauliLetter)"/>:
    /// the letter is validated and does not change the computation.
    ///
    /// <para>Only <see cref="PauliLetter.Z"/> is currently supported by the per-block builder
    /// (<see cref="PerBlockLiouvillianBuilder.BuildBlockZ"/> is hardcoded Z-only and the
    /// joint-popcount sector structure is not preserved under X- or Y-dephasing). X and Y
    /// throw <see cref="NotSupportedException"/> (design-permanent under the current basis);
    /// I throws <see cref="ArgumentException"/> (not a valid dephase letter for a Lindblad
    /// dissipator). The message points to the overload taking a <see cref="PauliHamiltonian"/>, which
    /// serves X and Y by an exact letter turn.</para></summary>
    /// <param name="dephaseLetter">Dephase letter. Only Z is
    /// currently supported by the per-block construction; X and Y throw
    /// <see cref="NotSupportedException"/>, I throws <see cref="ArgumentException"/>.</param>
    public static Complex[] ComputeSpectrumPerBlock(
        ComplexMatrix H, IReadOnlyList<double> gammaPerSite, int N, PauliLetter dephaseLetter) =>
        ComputeSpectrumPerBlock(H, gammaPerSite, N, dephaseLetter, LiouvillianBlockSpectrum.SectorPairing.PiOrbit);

    /// <summary>X-, Y- or Z-dephasing for a Hamiltonian given as Pauli terms, as on
    /// <see cref="LiouvillianBlockSpectrum.ComputeSpectrumPerBlock(PauliHamiltonian, IReadOnlyList{double}, PauliLetter, LiouvillianBlockSpectrum.SectorPairing, LiouvillianBlockSpectrum.EigenPath)"/>:
    /// the letters are turned exactly so that the dephasing letter becomes Z, then the Z path runs.
    /// The turn is the same on every site, so it commutes with the chain reflection and the F71
    /// refinement keeps its meaning.</summary>
    public static Complex[] ComputeSpectrumPerBlock(
        PauliHamiltonian H, IReadOnlyList<double> gammaPerSite, PauliLetter dephaseLetter,
        LiouvillianBlockSpectrum.SectorPairing pairing = LiouvillianBlockSpectrum.SectorPairing.PiOrbit)
    {
        if (H is null) throw new ArgumentNullException(nameof(H));
        if (dephaseLetter == PauliLetter.I)
            throw new ArgumentException(
                "PauliLetter.I is not a valid dephase letter (the Lindblad dissipator requires a non-identity operator).",
                nameof(dephaseLetter));
        var turned = LetterTurn.Turn(H, dephaseLetter);
        return ComputeSpectrumPerBlock(turned.ToMatrix(), gammaPerSite, H.N, PauliLetter.Z, pairing);
    }

    /// <summary>The full overload: dephase letter (Z only) and sector pairing. The pairing is
    /// used only when H commutes exactly with X^⊗N.</summary>
    public static Complex[] ComputeSpectrumPerBlock(
        ComplexMatrix H, IReadOnlyList<double> gammaPerSite, int N, PauliLetter dephaseLetter,
        LiouvillianBlockSpectrum.SectorPairing pairing)
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
        if (dephaseLetter == PauliLetter.I)
            throw new ArgumentException(
                "PauliLetter.I is not a valid dephase letter (the Lindblad dissipator requires a " +
                "non-identity operator); use Z (canonical / F108 Part 1; X and Y, F108 Parts 2 and 3, are " +
                "refused under the current basis).",
                nameof(dephaseLetter));
        if (dephaseLetter != PauliLetter.Z)
            throw new NotSupportedException(
                $"The dense F71MirrorBlockRefinement.ComputeSpectrumPerBlock overloads only support Z-dephasing in the " +
                $"joint-popcount basis (PerBlockLiouvillianBuilder.BuildBlockZ is hardcoded to Z, which is " +
                $"diagonal in the computational basis and popcount-conserving); got {dephaseLetter}. " +
                "X- and Y-dephasing break the joint-popcount sector structure that JointPopcountSectors and " +
                "the F71 mirror refinement rely on. Use the PauliHamiltonian overload, which turns the Pauli letters " +
                "exactly so that the dephasing letter becomes Z.");

        // LiouvillianBlockSpectrum's contract check: H must be popcount-conserving in the 2^N
        // Hilbert basis, otherwise the joint-popcount sector decomposition misses part of L.
        // Exact, on every entry.
        LiouvillianBlockSpectrum.RequirePopcountConservingH(H, N);

        // Reflection constant: the genuine Σ of per-site Z-dephasing rates (NOT N·γ), so
        // non-uniform γ stays exact. When H commutes with X^⊗N a sector's spectrum maps to its
        // Π-image's via λ ↦ −2·Σγ − λ (LiouvillianBlockSpectrum's class summary); F1's
        // identity is the truly case.
        double sumGamma = gammaPerSite.Sum();

        int liouvilleDim = 1 << (2 * N);
        int d = hilbertDim;

        // Build the per-Hilbert-side P_F71 mirror map.
        var mirrorBits = F71MirrorIndexHelper.BuildHilbertMirrorLookup(N);
        int Mirror(int flat) => mirrorBits[flat / d] * d + mirrorBits[flat % d];

        var baseDecomp = JointPopcountSectorBuilder.Build(N);
        var spectrum = new Complex[liouvilleDim];

        // H1: each (joint-popcount sector → even+odd sub-blocks) chunk is independent. Pre-
        // compute the per-sector write offset (= sum of even+odd sizes for sectors before it)
        // so each task writes its slice deterministically. Even-size = nFix + nPairs, odd-size
        // = nPairs, so per-sector write = sectorSize (== nFix + 2·nPairs).
        int sectorCount = baseDecomp.SectorRanges.Count;
        var writeOffsets = new int[sectorCount];
        int cum = 0;
        for (int i = 0; i < sectorCount; i++)
        {
            writeOffsets[i] = cum;
            cum += baseDecomp.SectorRanges[i].Size;
        }

        // Which sectors are solved and how the rest are derived (shared with
        // LiouvillianBlockSpectrum): with [H, X^⊗N] = 0 the requested pairing, otherwise every
        // sector. The pairing respects F71 parity, the F71-even sub-block of a sector mapping to
        // the F71-even sub-block of its image (ρ ↦ X^⊗N·ρ·X^⊗N and ρ ↦ ρ†·X^⊗N, which carry a
        // sector onto its images, commute with the site reflection; the
        // Pi_CommutesWith_F71MirrorPair test checks the same of Π), so Phase 3 copies or reflects BOTH
        // sub-arrays. Primaries are sorted descending by size so the largest starts first
        // under Parallel.ForEach.
        var (primarySectorIndices, followerToPrimary) =
            LiouvillianBlockSpectrum.PartitionSectors(H, N, baseDecomp.SectorRanges, pairing);

        // BLAS-oversubscription strategy (c): outer DOP ≈ ProcessorCount/4. See
        // LiouvillianBlockSpectrum.ComputeSpectrum for rationale.
        int outerDop = Math.Max(1, Environment.ProcessorCount / 4);
        var po = new ParallelOptions { MaxDegreeOfParallelism = outerDop };

        // Phase 2: parallel eig on primary sectors only. Stores (evenEigs, oddEigs) per
        // primary sector. Followers reuse these arrays via copy in Phase 3.
        var primaryEigs = new ConcurrentDictionary<int, (Complex[] Even, Complex[] Odd)>();

        Parallel.ForEach(
            primarySectorIndices, po,
            sIdx =>
            {
                var sector = baseDecomp.SectorRanges[sIdx];
                int size = sector.Size;
                if (size == 0)
                {
                    primaryEigs[sIdx] = (Array.Empty<Complex>(), Array.Empty<Complex>());
                    return;
                }

                // Find F71 orbits in this sector.
                var sectorFlat = new int[size];
                for (int k = 0; k < size; k++) sectorFlat[k] = baseDecomp.Permutation[sector.Offset + k];
                var (fixedPoints, pairs) = F71MirrorIndexHelper.FindOrbitsInSector(sectorFlat, Mirror);

                // Build the union block over (fixed-points ++ pair-firsts ++ pair-seconds).
                int nFix = fixedPoints.Count;
                int nPairs = pairs.Count;
                int unionSize = nFix + 2 * nPairs;
                var unionFlat = new int[unionSize];
                for (int i = 0; i < nFix; i++) unionFlat[i] = fixedPoints[i];
                for (int k = 0; k < nPairs; k++) unionFlat[nFix + k] = pairs[k].S;
                for (int k = 0; k < nPairs; k++) unionFlat[nFix + nPairs + k] = pairs[k].Ps;

                var unionBlock = PerBlockLiouvillianBuilder.BuildBlockZ(H, gammaPerSite, unionFlat);

                // H3: in-place F71 rotation B' = R^T · unionBlock · R is computed entry-by-entry
                // from the F71-orbit Hadamard structure (O(n²) instead of two O(n³) matmuls).
                var rotated = RotateUnionBlockF71InPlace(unionBlock, nFix, nPairs);

                int evenSize = nFix + nPairs;
                int oddSize = nPairs;
                Complex[] evenArr = Array.Empty<Complex>();
                Complex[] oddArr = Array.Empty<Complex>();
                if (evenSize > 0)
                {
                    var evenBlock = rotated.SubMatrix(0, evenSize, 0, evenSize);
                    var evenEigs = evenBlock.Evd().EigenValues;
                    evenArr = new Complex[evenSize];
                    for (int i = 0; i < evenSize; i++) evenArr[i] = evenEigs[i];
                }
                if (oddSize > 0)
                {
                    var oddBlock = rotated.SubMatrix(evenSize, oddSize, evenSize, oddSize);
                    var oddEigs = oddBlock.Evd().EigenValues;
                    oddArr = new Complex[oddSize];
                    for (int i = 0; i < oddSize; i++) oddArr[i] = oddEigs[i];
                }
                primaryEigs[sIdx] = (evenArr, oddArr);
            });

        // Phase 3: sequential write to output array. Each sector writes its even-block
        // eigenvalues then its odd-block eigenvalues in the same layout the original code
        // produced. A primary writes its own (even, odd) arrays verbatim. A follower derives
        // them from its orbit primary's: an X⊗N-image follower (Π²-image) copies both arrays
        // verbatim (exact, since the partition pairs only when H commutes with X^⊗N); a
        // Π/Π³-image follower reflects each λ in
        // BOTH arrays through the F1 palindrome map λ ↦ −2·Σγ − λ. The map that carries a sector
        // onto its Π-image, ρ ↦ ρ†·X^⊗N, commutes with the F71 mirror, so the reflection respects
        // F71 parity: even stays even, odd stays odd.
        for (int sIdx = 0; sIdx < sectorCount; sIdx++)
        {
            int write = writeOffsets[sIdx];
            (Complex[] Even, Complex[] Odd) source;
            bool reflect;
            if (primaryEigs.ContainsKey(sIdx))
            {
                source = primaryEigs[sIdx];
                reflect = false;
            }
            else
            {
                var follower = followerToPrimary[sIdx];
                source = primaryEigs[follower.PrimaryIndex];
                reflect = follower.Kind == F1PalindromeOrbitPairing.F1FollowerKind.F1Reflect;
            }

            var (evenArr, oddArr) = source;
            if (!reflect)
            {
                for (int i = 0; i < evenArr.Length; i++) spectrum[write + i] = evenArr[i];
                write += evenArr.Length;
                for (int i = 0; i < oddArr.Length; i++) spectrum[write + i] = oddArr[i];
            }
            else
            {
                // F1Reflect: λ ↦ −2·Σγ − λ on both the F71-even and F71-odd sub-arrays.
                for (int i = 0; i < evenArr.Length; i++)
                    spectrum[write + i] = new Complex(
                        -2.0 * sumGamma - evenArr[i].Real, -evenArr[i].Imaginary);
                write += evenArr.Length;
                for (int i = 0; i < oddArr.Length; i++)
                    spectrum[write + i] = new Complex(
                        -2.0 * sumGamma - oddArr[i].Real, -oddArr[i].Imaginary);
            }
        }
        return spectrum;
    }

    /// <summary>H3: in-place F71 union-block rotation B' = R^T · B · R, computed entry-by-entry
    /// from the F71-orbit Hadamard structure rather than via two O(n³) MathNet matmuls.
    ///
    /// <para>R is identity on the first <paramref name="nFix"/> rows/cols (F71 fixed-point
    /// pairs) and a 2×2 Hadamard <c>[[1/√2, 1/√2], [1/√2, −1/√2]]</c> per (s, ps) pair acting
    /// across rows/cols (nFix + k, nFix + nPairs + k) for k = 0..nPairs-1. The new basis is
    /// laid out as <c>(fixed-points, even-pairs, odd-pairs)</c> in that order.</para>
    ///
    /// <para>Cost: O(unionSize²) work, no matmul. Equivalent to <c>R^T · B · R</c>
    /// up to floating-point rounding (verified via the existing
    /// <c>F71MirrorBlockRefinement_OffBlockFrobenius_IsZero_ChainXYZDeph</c> tolerance).</para></summary>
    public static ComplexMatrix RotateUnionBlockF71InPlace(ComplexMatrix B, int nFix, int nPairs)
    {
        int unionSize = nFix + 2 * nPairs;
        if (B.RowCount != unionSize || B.ColumnCount != unionSize)
            throw new ArgumentException(
                $"B must be ({unionSize})×({unionSize}); got {B.RowCount}×{B.ColumnCount}.",
                nameof(B));

        var Bp = Matrix<Complex>.Build.Dense(unionSize, unionSize);
        const double half = 0.5;
        double invSqrt2 = 1.0 / Math.Sqrt(2.0);

        // Layout in the new basis (B'):
        //   rows/cols [0 .. nFix)             = fixed-point block (R = I)
        //   rows/cols [nFix .. nFix+nPairs)   = even-pair block (e_s + e_ps)/√2
        //   rows/cols [nFix+nPairs .. union)  = odd-pair block  (e_s − e_ps)/√2
        // Source layout in B:
        //   rows/cols [0 .. nFix)             = fixed points
        //   rows/cols [nFix .. nFix+nPairs)   = s_k indices  (one per pair)
        //   rows/cols [nFix+nPairs .. union)  = ps_k indices (one per pair)

        // Block FF (fixed × fixed): identity on both sides → B'[i, j] = B[i, j].
        for (int i = 0; i < nFix; i++)
            for (int j = 0; j < nFix; j++)
                Bp[i, j] = B[i, j];

        // Block FE / FO (fixed-row × pair-col): identity row, Hadamard col.
        // Even-col k:  B'[i, nFix + k] = (B[i, sk] + B[i, psk]) / √2
        // Odd-col  k:  B'[i, nFix + nPairs + k] = (B[i, sk] − B[i, psk]) / √2
        for (int i = 0; i < nFix; i++)
        {
            for (int k = 0; k < nPairs; k++)
            {
                int sk = nFix + k;
                int psk = nFix + nPairs + k;
                Complex bs = B[i, sk];
                Complex bp = B[i, psk];
                Bp[i, sk] = (bs + bp) * invSqrt2;          // even col
                Bp[i, psk] = (bs - bp) * invSqrt2;          // odd col
            }
        }

        // Block EF / OF (pair-row × fixed-col): Hadamard row, identity col.
        // Even-row k:  B'[nFix + k, j] = (B[sk, j] + B[psk, j]) / √2
        // Odd-row  k:  B'[nFix + nPairs + k, j] = (B[sk, j] − B[psk, j]) / √2
        for (int k = 0; k < nPairs; k++)
        {
            int sk = nFix + k;
            int psk = nFix + nPairs + k;
            for (int j = 0; j < nFix; j++)
            {
                Complex bs = B[sk, j];
                Complex bp = B[psk, j];
                Bp[sk, j] = (bs + bp) * invSqrt2;
                Bp[psk, j] = (bs - bp) * invSqrt2;
            }
        }

        // Block EE / EO / OE / OO (pair-row × pair-col): Hadamard on both sides → 4-term
        // combination with signs from the 2×2 Hadamard pattern.
        //   sk_l, psk_l = source indices for col-pair l
        //   sk_k, psk_k = source indices for row-pair k
        //   B'[evenRow_k, evenCol_l] = (B[sk_k, sk_l] + B[sk_k, psk_l] + B[psk_k, sk_l] + B[psk_k, psk_l]) / 2
        //   B'[evenRow_k, oddCol_l]  = (B[sk_k, sk_l] − B[sk_k, psk_l] + B[psk_k, sk_l] − B[psk_k, psk_l]) / 2
        //   B'[oddRow_k,  evenCol_l] = (B[sk_k, sk_l] + B[sk_k, psk_l] − B[psk_k, sk_l] − B[psk_k, psk_l]) / 2
        //   B'[oddRow_k,  oddCol_l]  = (B[sk_k, sk_l] − B[sk_k, psk_l] − B[psk_k, sk_l] + B[psk_k, psk_l]) / 2
        for (int k = 0; k < nPairs; k++)
        {
            int skR = nFix + k;
            int pskR = nFix + nPairs + k;
            for (int l = 0; l < nPairs; l++)
            {
                int skC = nFix + l;
                int pskC = nFix + nPairs + l;
                Complex a = B[skR, skC];     // (s, s)
                Complex b = B[skR, pskC];    // (s, ps)
                Complex c = B[pskR, skC];    // (ps, s)
                Complex d = B[pskR, pskC];   // (ps, ps)

                Bp[skR, skC] = (a + b + c + d) * half;     // even-row × even-col
                Bp[skR, pskC] = (a - b + c - d) * half;     // even-row × odd-col
                Bp[pskR, skC] = (a + b - c - d) * half;     // odd-row  × even-col
                Bp[pskR, pskC] = (a - b - c + d) * half;     // odd-row  × odd-col
            }
        }

        return Bp;
    }
}
