using System.Collections.Concurrent;
using System.Diagnostics;
using System.Numerics;
using System.Runtime.InteropServices;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Numerics;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.Symmetry;
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
/// <para><b>Contract.</b> Both <see cref="ComputeSpectrum"/> and
/// <see cref="ComputeSpectrumPerBlock"/> require that the input Liouvillian (or its
/// underlying Hamiltonian) be block-diagonal in the joint-popcount basis. This holds
/// for popcount-conserving H (XX+YY, ZZ, XXZ, Heisenberg, any sum of these) and FAILS
/// for H that breaks popcount conservation (XX+YZ, XY+YX, anything with shadow-crossing
/// Pauli pairs like X_iZ_j or Y_iZ_j). Calling with non-popcount-conserving H returns
/// silently wrong spectra: the per-block eigenvalues miss the cross-sector entries of L,
/// so <c>Σ_blocks ‖L_b‖²_F ≠ ‖L_full‖²_F</c> and the spectrum is incomplete. The 2026-05-18
/// F1 general-topology N=7 dogfood discovered this empirically (XX+YZ at N=5 gave block
/// sum 4403 vs dense 16691, factor ~3.8 off). In DEBUG builds a sample-based assertion in
/// each entry point (<see cref="DebugAssertBlockDiagonalL"/> and
/// <see cref="DebugAssertPopcountConservingH"/>) throws an
/// <see cref="InvalidOperationException"/> when the contract is violated; in RELEASE
/// builds the assertion is stripped (no production-perf cost). Callers must validate H
/// structurally before invoking on a hot path. Both methods ensure MKL is initialized
/// via <see cref="MathNetSetup.EnsureInitialized"/> on entry, so no caller needs to
/// pre-initialize; the lazy global guard makes the redundant call free after the
/// assembly-level <c>CoreModuleInitializer</c> has already run.</para>
///
/// <para><b>F108 Π_5bilinear at the Builder layer (2026-05-25).</b> The orbit-pairing
/// primitive <see cref="SymmetryFamily.F1PalindromeOrbitPairing.PartitionByPiOrbit"/> is
/// Π-agnostic at the sector-label level: it consumes only the joint-popcount permutation
/// rule (p_c, p_r) ↦ (N − p_r, p_c), not the matrix form of Π. F108 Part 1
/// (<see cref="Symmetry.F108Part1Pi2EvenAlwaysPalindromic"/> via
/// <see cref="Symmetry.Pi5BilinearOperator"/> at <c>dephaseLetter = Z</c>) extends the F1
/// conjugation identity Π·L·Π⁻¹ = −L − 2σ·I from chain Heisenberg/XY (truly) to every
/// Π²_Z-even 2-site bilinear H (XX, YY, YZ, ZY, ZZ combos). Canonical Π (F1) and
/// Π_5bilinear-Z (F108 Part 1) induce the SAME joint-popcount sector cycle (only the
/// per-letter phases differ at the operator level), so on the Z-dephasing domain F108
/// brings no extra speedup at this layer: the existing orbit-pairing already covers the
/// broader H set the contract admits. The F108-aware <c>ComputeSpectrumPerBlock</c>
/// overload exists for intent-declaration and forward-compatibility, not for new gain.</para>
///
/// <para><b>Z-dephasing structural constraint.</b> The joint-popcount sector basis is
/// tied to Z-dephasing because <see cref="PerBlockLiouvillianBuilder.BuildBlockZ"/> builds
/// the dissipator element-wise from the computational-basis bit-parity disagreement, a
/// representation that requires the dephase letter to be diagonal in that basis. F108 Part 2
/// (X-deph, <see cref="Symmetry.F108Part2Pi2XEvenAlwaysPalindromic"/>) and Part 3 (Y-deph,
/// <see cref="Symmetry.F108Part3Pi2YEvenAlwaysPalindromic"/>) prove operator-level
/// palindromicity for those dephase channels, but the Builder cannot exploit them in its
/// current basis (X- and Y-dephasing break popcount conservation in the computational
/// basis). The F108-aware overload throws <see cref="NotSupportedException"/> for
/// <see cref="Pauli.PauliLetter.X"/> / <see cref="Pauli.PauliLetter.Y"/> (design-permanent
/// under the current basis; CLR convention reserves <see cref="NotImplementedException"/>
/// for stubs awaiting a body, <see cref="NotSupportedException"/> for combinations the
/// design refuses by intent) and <see cref="ArgumentException"/> for
/// <see cref="Pauli.PauliLetter.I"/> (not a valid dephase letter; the Lindblad dissipator
/// requires a non-identity operator) rather than silently producing wrong eigenvalues.
/// Lifting the X/Y restriction is tracked in <see cref="BlockSpectrumOpenQuestions"/> and
/// requires a per-dephase-letter rotated basis (e.g. apply per-site U_X = H for X-deph,
/// conjugate H → H' = U H U†, then run the existing BuildBlockZ; L's eigenvalues are
/// basis-independent so the result transports back).</para>
///
/// <para>Anchors: <c>compute/RCPsiSquared.Core/BlockSpectrum/JointPopcountSectors.cs</c>
/// (parent Claim, block-diagonal structure), <c>compute/RCPsiSquared.Core/BlockSpectrum/JointPopcountSectorBuilder.cs</c>
/// (basis permutation + sector ranges), <c>compute/RCPsiSquared.Core.Tests/BlockSpectrum/LiouvillianBlockSpectrumTests.cs</c>
/// (bit-exact spectral equality verification at N=3, 4, 5).</para></summary>
public sealed class LiouvillianBlockSpectrum : Claim
{
    private readonly JointPopcountSectors _sectors;

    /// <summary>Per-block eigensolver selection knob for
    /// <see cref="ComputeSpectrumPerBlock(ComplexMatrix, IReadOnlyList{double}, int, EigenPath)"/>.
    /// The default <see cref="Auto"/> picks <see cref="MathNet"/> for blocks below
    /// <see cref="Lp64ComplexCeiling"/> and <see cref="MklDirectNative"/> above; the explicit
    /// overrides exist only for the parity-witness test
    /// <c>PerBlockLiouvillianBuilderNativeMemoryParityTests</c> which forces both paths on the
    /// same small block to demonstrate bit-exact agreement.</summary>
    public enum EigenPath
    {
        /// <summary>Auto-select by block size: MathNet for size ≤ <see cref="Lp64ComplexCeiling"/>,
        /// MklDirect + NativeMemory + ILP64 above. The production default.</summary>
        Auto = 0,
        /// <summary>Force the MathNet <c>Matrix&lt;Complex&gt;.Evd()</c> path on every block.
        /// Will throw inside MathNet's marshaller for blocks &gt; <see cref="Lp64ComplexCeiling"/>.
        /// Test-only.</summary>
        MathNet = 1,
        /// <summary>Force the <see cref="MklDirect.EigenvaluesOnlyNative"/> path on every block.
        /// Allocates and frees a native column-major Complex buffer per block. Test-only; the
        /// production code uses <see cref="Auto"/> which only takes this branch above the
        /// LP64 ceiling.</summary>
        MklDirectNative = 2,
    }

    /// <summary>Largest square Complex block (size n × n with n² ≤ 134 217 728) whose total
    /// byte count (n² × 16) fits inside the LP64 2 GB single-native-array marshalling ceiling
    /// enforced by MathNet's <c>MklLinearAlgebraProvider.EigenDecomp</c>. n = 11 585 gives
    /// n² × 16 ≈ 2.147 GB (just inside the marshaller's <c>int.MaxValue</c> byte limit when the
    /// fixed-pinned <c>Complex[]</c> is rounded into a single P/Invoke array). Blocks at or
    /// below this size are routed through MathNet's well-tested managed wrapper; blocks above
    /// are routed through <see cref="MklDirect.EigenvaluesOnlyNative"/> on a
    /// <see cref="NativeMemory.AllocZeroed(nuint)"/>-backed buffer, which also auto-selects
    /// ILP64 when n &gt; 46 340 (n² &gt; <c>int.MaxValue</c>; not currently reachable at
    /// joint-popcount block sizes for N ≤ 12).
    ///
    /// <para>The threshold is the same value the N=9 test
    /// (<c>F1GeneralTopologyN9BlockSpectrumChainTests.Lp64EvdSquareMatrixCeiling</c>) uses for
    /// its pre-flight check; keeping both in sync means a future bump (e.g. if MathNet relaxes
    /// the marshaller) only needs updating one constant per file.</para></summary>
    public const int Lp64ComplexCeiling = 11_585;

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
    /// and <see cref="Lindblad.PauliDephasingDissipator"/>. Must be (4^N) × (4^N).
    /// Must be block-diagonal in the joint-popcount basis (see class Contract).</param>
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

        // Belt-and-braces with CoreModuleInitializer: makes the MKL dependency
        // self-documenting at the API boundary; the lazy guard makes it free after first call.
        MathNetSetup.EnsureInitialized();

        var decomp = JointPopcountSectorBuilder.Build(N);
        DebugAssertBlockDiagonalL(L, decomp);
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
    /// <see cref="PerBlockLiouvillianBuilder.BuildBlockZ"/> (MathNet path) or
    /// <see cref="PerBlockLiouvillianBuilder.BuildBlockZIntoNativeMemory"/> (MklDirect +
    /// NativeMemory + ILP64 path), eigendecomposed, and discarded before the next block.
    /// This is the only path that scales past N=6 on commodity hardware (full L exceeds
    /// .NET 2 GB array-size limit at N=7+).
    ///
    /// <para>Eigensolver selection is automatic by block size: blocks at or below
    /// <see cref="Lp64ComplexCeiling"/> (11 585² ≈ 2 GB Complex matrix) go through MathNet's
    /// well-tested managed wrapper; blocks above the ceiling route through
    /// <see cref="MklDirect.EigenvaluesOnlyNative"/> on a
    /// <see cref="NativeMemory.AllocZeroed(nuint)"/>-backed column-major buffer. The bridge
    /// unlocks N ≥ 9 where the largest joint-popcount sector exceeds the LP64 marshaller's
    /// 2 GB cap (N=9 max block C(9, 4) · C(9, 5) = 15 876² ≈ 4 GB; see
    /// <c>F1GeneralTopologyN9BlockSpectrumChainTests</c>).</para>
    ///
    /// <para>Uses the F1 Π-orbit pairing optimisation: the F1 palindrome conjugation Π is
    /// order-4 and on joint-popcount labels acts as the whole-sector cycle
    /// (p_c, p_r) ↦ (N − p_r, p_c), grouping the (N+1)² sectors into orbits of 4. Only one
    /// "primary" sector per orbit is eigendecomposed; the three followers are derived from
    /// it, the Π²-image (X⊗N partner) by a verbatim spectrum copy and the Π/Π³-images by the
    /// F1 reflection λ ↦ −2·Σγ − λ. This quarters the eigendecomposition count (the X⊗N
    /// pairing of <see cref="SymmetryFamily.XGlobalChargeConjugationPairing"/>, which Π²
    /// equals, only halved it). See <see cref="SymmetryFamily.F1PalindromeOrbitPairing"/> for
    /// the orbit rule.</para></summary>
    /// <param name="H">Hilbert-space Hamiltonian, dense 2^N × 2^N (cheap even at N=9: 512×512).
    /// Must be popcount-conserving (see class Contract); non-conserving H gives silently
    /// wrong spectra.</param>
    /// <param name="gammaPerSite">Per-site Z-dephasing rates (length N).</param>
    /// <param name="N">Qubit count.</param>
    /// <returns>Flat array of 4^N eigenvalues, concatenated block-by-block.</returns>
    public static Complex[] ComputeSpectrumPerBlock(ComplexMatrix H, IReadOnlyList<double> gammaPerSite, int N) =>
        ComputeSpectrumPerBlock(H, gammaPerSite, N, EigenPath.Auto);

    /// <summary>Test-aware overload that lets the parity witness force a specific eigensolver
    /// path on every block. Production callers should use the parameterless overload (or pass
    /// <see cref="EigenPath.Auto"/> explicitly) which routes by block size against
    /// <see cref="Lp64ComplexCeiling"/>.</summary>
    /// <param name="path">Force the MathNet path, force the MklDirect + NativeMemory path,
    /// or let the per-block size decide. <see cref="EigenPath.MathNet"/> will throw inside
    /// MathNet's marshaller for blocks larger than the LP64 ceiling; use
    /// <see cref="EigenPath.Auto"/> in production code.</param>
    public static Complex[] ComputeSpectrumPerBlock(
        ComplexMatrix H, IReadOnlyList<double> gammaPerSite, int N, EigenPath path) =>
        ComputeSpectrumPerBlock(H, gammaPerSite, N, path, PauliLetter.Z);

    /// <summary>F108-aware overload accepting an explicit dephase letter D ∈ {X, Y, Z}
    /// that scopes the orbit-pairing optimisation. The default <see cref="PauliLetter.Z"/>
    /// preserves the historical behaviour (chain XY+Z-deph orbit-pairing via canonical
    /// <see cref="PiOperator"/>); explicit values trigger the F108-family Π_5bilinear
    /// dispatch documented in
    /// <see cref="F108Part1Pi2EvenAlwaysPalindromic"/>,
    /// <see cref="F108Part2Pi2XEvenAlwaysPalindromic"/>, and
    /// <see cref="F108Part3Pi2YEvenAlwaysPalindromic"/>.
    ///
    /// <para><b>Auto-detect.</b> The orbit-pairing's F1 reflection
    /// λ ↦ −2·Σγ − λ is justified iff <c>Π·L·Π⁻¹ = −L − 2σ·I</c> holds for some Π that
    /// permutes the (p_c, p_r) joint-popcount sectors as the whole-sector cycle
    /// (p_c, p_r) ↦ (N − p_r, p_c). F1 establishes this for chain XY+Z-deph
    /// (Π = canonical <see cref="PiOperator"/>); F108 Part 1 extends it to every Π²_Z-even
    /// 2-site bilinear H + Z-dephasing via the Z-deph variant of
    /// <see cref="Pi5BilinearOperator"/>. Both variants share the same per-letter permutation
    /// pattern (I↔X, Y↔Z), so they induce the same joint-popcount sector cycle; the orbit-
    /// pairing primitive <see cref="F1PalindromeOrbitPairing.PartitionByPiOrbit"/> is
    /// already Π-agnostic at the sector-label level and consumes only the sector-permutation
    /// rule, not the matrix Π.</para>
    ///
    /// <para><b>Z-only restriction.</b> The per-block Liouvillian construction in
    /// <see cref="PerBlockLiouvillianBuilder.BuildBlockZ"/> is hardcoded to Z-dephasing
    /// (the dissipator is built element-wise in the computational basis from the Hilbert-
    /// side bit-parity disagreement, not from a general P_l ⊗ P_l⁺ Kronecker construction).
    /// Non-Z-dephasing matchings (F108 Part 2 X-deph, F108 Part 3 Y-deph) are not currently
    /// supported by this entry point and throw <see cref="NotSupportedException"/>
    /// (design-permanent under the current basis); lifting the restriction requires a
    /// parallel <c>BuildBlockX</c> / <c>BuildBlockY</c> path plus a rederived per-letter-D
    /// sector decomposition compatible with the chosen dephase letter (X- and Y-dephasing
    /// break popcount conservation in the computational basis; see the class Contract).
    /// <see cref="Pauli.PauliLetter.I"/> is rejected up-front with
    /// <see cref="ArgumentException"/> (not a valid dephase letter).</para>
    ///
    /// <para>Mismatch handling: <paramref name="dephaseLetter"/> = X or Y throws cleanly
    /// rather than silently producing wrong eigenvalues (the underlying per-block builder
    /// would otherwise stamp Z-deph entries onto an X- or Y-deph problem).</para></summary>
    /// <param name="dephaseLetter">The dephase letter the orbit-pairing is matched to. Only
    /// <see cref="PauliLetter.Z"/> is currently supported by the per-block construction; X and
    /// Y throw <see cref="NotSupportedException"/> (design-permanent under the current basis),
    /// I throws <see cref="ArgumentException"/> (not a valid dephase letter).</param>
    public static Complex[] ComputeSpectrumPerBlock(
        ComplexMatrix H, IReadOnlyList<double> gammaPerSite, int N, EigenPath path,
        PauliLetter dephaseLetter)
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
                "non-identity operator); use Z (canonical / F108 Part 1), X (F108 Part 2), or Y (F108 Part 3).",
                nameof(dephaseLetter));
        if (dephaseLetter != PauliLetter.Z)
            throw new NotSupportedException(
                $"ComputeSpectrumPerBlock only supports Z-dephasing under the current joint-popcount " +
                $"basis (PerBlockLiouvillianBuilder.BuildBlockZ is hardcoded to Z, which is diagonal in " +
                $"the computational basis and popcount-conserving); got {dephaseLetter}. X- and Y-dephasing " +
                "break the joint-popcount sector structure that JointPopcountSectors relies on; see " +
                "BlockSpectrumOpenQuestions for the X/Y basis-rotation extension path (e.g. apply per-site " +
                "U_X = H for X-deph, conjugate H → H' = U H U†, then call this entry point with Z).");
        DebugAssertPopcountConservingH(H, N);

        // F1 palindrome reflection constant: the genuine Σ of per-site Z-dephasing rates
        // (NOT N·γ), so non-uniform γ stays exact. F1 maps a sector's spectrum to its
        // Π-image's via λ ↦ −2·Σγ − λ (docs/proofs/MIRROR_SYMMETRY_PROOF.md).
        double sumGamma = gammaPerSite.Sum();

        // Belt-and-braces with CoreModuleInitializer: makes the MKL dependency
        // self-documenting at the API boundary; the lazy guard makes it free after first call.
        MathNetSetup.EnsureInitialized();

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

        // F1 Π-orbit pairing (Tier 1, F1PalindromeOrbitPairing): the chain XY+Z-deph L
        // satisfies the F1 palindrome Π·L·Π⁻¹ = −L − 2·Σγ·I. Π is order-4 and on joint-
        // popcount labels acts as the whole-sector cycle (p_c, p_r) ↦ (N − p_r, p_c),
        // grouping the (N+1)² sectors into orbits of 4 (plus the single Π-fixed (N/2, N/2)
        // sector at even N). This subsumes the X⊗N pairing (Π² = X⊗N): one eigendecomposition
        // per orbit feeds three followers. We compute eig only on the lex-smallest "primary"
        // of each orbit (plus the Π-fixed sector) and derive the followers in Phase 3 — the
        // Π²-image by a verbatim spectrum copy, the Π/Π³-images by the F1 reflection
        // λ ↦ −2·Σγ − λ. Primaries are sorted descending by size so the largest sector starts
        // first under Parallel.ForEach, overlapping its wall-time with smaller sectors' work.
        //
        // F108 generalisation: F108 Part 1 (Pi5BilinearOperator with dephaseLetter = Z)
        // extends the F1 conjugation identity from chain Heisenberg/XY (truly) to every
        // Π²_Z-even bilinear H (XX, YY, YZ, ZY, ZZ combos). At the BUILDER level the gain
        // is operator-equivalent for the Hamiltonians this entry point can actually accept:
        // popcount conservation in the computational basis (required by JointPopcountSectors)
        // restricts 2-body H to linear combinations of (XX+YY) and ZZ, all of which are also
        // Π²_Z-even. So canonical Π (F1) and Π_5bilinear (F108 Part 1) induce the same
        // sector cycle and the same F1 reflection on this domain. The dephaseLetter parameter
        // declares the user's intent; the partition rule is Π-agnostic at the sector-label
        // level (F1PalindromeOrbitPairing.PartitionByPiOrbit consumes only the (p_c, p_r) ↦
        // (N − p_r, p_c) rule). F108 Part 2 (X-deph) / Part 3 (Y-deph) dispatch is reserved
        // for a future per-block-D builder.
        var (primarySectorIndices, followerToPrimary) =
            F1PalindromeOrbitPairing.PartitionByPiOrbit(
                N, decomp.SectorRanges, s => (s.PCol, s.PRow), s => s.Size);

        // BLAS-oversubscription strategy (c): outer DOP ≈ ProcessorCount/4 leaves room for
        // MKL inside the largest sectors' Evd. See ComputeSpectrum for the rationale.
        //
        // EigenPath.MklDirectNative serialises (DOP=1) because each block can hold ~4 GB
        // of NativeMemory while LAPACK runs; allowing multiple large MklDirect blocks
        // concurrently risks running the working set past 128 GB on N=9 chains where
        // the four largest paired-primary sectors each carry ≥ 1 GB of block data.
        // Smaller MathNet blocks at the same N are unaffected (they live below the
        // ceiling and use ProcessorCount/4 outer DOP).
        int outerDop = path == EigenPath.MklDirectNative
            ? 1
            : Math.Max(1, Environment.ProcessorCount / 4);
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

                primaryEigs[sIdx] = SolveSectorBlock(H, gammaPerSite, flatIndices, size, path);
            });

        // Phase 3: sequential write to output array (primaries + followers). A primary writes
        // its own eigenvalues verbatim. A follower derives its eigenvalues from its orbit
        // primary's: an X⊗N-image follower (Π²-image) copies them verbatim (X⊗N is a genuine
        // symmetry); a Π/Π³-image follower reflects each λ through the F1 palindrome map
        // λ ↦ −2·Σγ − λ. The eigenvalue MULTISET is what each sector needs; the per-block
        // ordering within a follower is irrelevant since the output is a flat union.
        for (int sIdx = 0; sIdx < sectorCount; sIdx++)
        {
            int write = writeOffsets[sIdx];
            if (primaryEigs.ContainsKey(sIdx))
            {
                var eigs = primaryEigs[sIdx];
                for (int i = 0; i < eigs.Length; i++) spectrum[write + i] = eigs[i];
                continue;
            }

            var follower = followerToPrimary[sIdx];
            var primary = primaryEigs[follower.PrimaryIndex];
            if (follower.Kind == F1PalindromeOrbitPairing.F1FollowerKind.XnCopy)
            {
                for (int i = 0; i < primary.Length; i++) spectrum[write + i] = primary[i];
            }
            else
            {
                // F1Reflect: λ ↦ −2·Σγ − λ (real part reflected about −Σγ; imaginary negated).
                for (int i = 0; i < primary.Length; i++)
                    spectrum[write + i] = new Complex(
                        -2.0 * sumGamma - primary[i].Real, -primary[i].Imaginary);
            }
        }
        return spectrum;
    }

    /// <summary>Per-block eigensolver dispatch. <see cref="EigenPath.Auto"/> picks MathNet for
    /// blocks at or below <see cref="Lp64ComplexCeiling"/> (well-tested wrapper, MKL multi-
    /// threaded BLAS-3) and MklDirect + NativeMemory + ILP64-aware for blocks above (bypasses
    /// the LP64 2 GB marshaller cap). Explicit overrides serve the parity-witness test only.
    ///
    /// <para>The MklDirect path allocates a single native column-major Complex buffer per
    /// block, calls <see cref="MklDirect.EigenvaluesOnlyNative"/> which automatically routes
    /// to ILP64 OpenBLAS when n &gt; 46 340, and frees the buffer in a <c>finally</c>. The
    /// LAPACK convention is that <c>zgeev</c> destroys the input matrix, so the buffer has
    /// single-use semantics; we discard it once eigenvalues are out.</para></summary>
    private static unsafe Complex[] SolveSectorBlock(
        ComplexMatrix H, IReadOnlyList<double> gammaPerSite, int[] flatIndices, int size,
        EigenPath path)
    {
        bool useNative = path switch
        {
            EigenPath.MathNet => false,
            EigenPath.MklDirectNative => true,
            // Auto: route by block size against the LP64 ceiling. The decisive cutoff comes
            // from the LP64 MathNet MklLinearAlgebraProvider.EigenDecomp marshaller cap; any
            // block of size > 11 585 would throw "Array size exceeds addressing limitations"
            // out of MngdNativeArrayMarshaler.ConvertSpaceToNative if routed via MathNet.
            _ => size > Lp64ComplexCeiling,
        };

        if (!useNative)
        {
            // MathNet path: identical to the pre-bridge behaviour, exercised at every block
            // size up through N=8 in the SLOW_N8 dogfood. Bit-exact lineage preserved.
            var block = PerBlockLiouvillianBuilder.BuildBlockZ(H, gammaPerSite, flatIndices);
            var blockEigs = block.Evd().EigenValues;
            var arr = new Complex[size];
            for (int i = 0; i < size; i++) arr[i] = blockEigs[i];
            return arr;
        }

        // MklDirect + NativeMemory + ILP64-aware path: the bridge that unlocks N ≥ 9.
        IntPtr ptr = PerBlockLiouvillianBuilder.BuildBlockZIntoNativeMemory(H, gammaPerSite, flatIndices);
        try
        {
            // EigenvaluesOnlyNative auto-selects LP64 vs ILP64 based on size against the 46 340
            // threshold (sqrt(int.MaxValue)). For block sizes 11 586..46 340 it stays on LP64
            // OpenBLAS but reads from NativeMemory rather than a managed Complex[], which is
            // precisely the marshaller bypass we need. Above 46 340 (not reachable until
            // N ≥ 13 joint-popcount sectors) it switches to ILP64 OpenBLAS automatically.
            return MklDirect.EigenvaluesOnlyNative(ptr, size);
        }
        finally
        {
            NativeMemory.Free((void*)ptr);
        }
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

    /// <summary>DEBUG-only sample-based check that L is block-diagonal in the joint-popcount
    /// basis. Picks 20 random index pairs from distinct sectors and asserts each
    /// <c>|L[f1, f2]| &lt; 1e-12</c>. Throws <see cref="InvalidOperationException"/> on the
    /// first violation, with sector labels in the message. Stripped from RELEASE builds.</summary>
    [Conditional("DEBUG")]
    private static void DebugAssertBlockDiagonalL(ComplexMatrix L, JointPopcountSectorBuilder.Decomposition decomp)
    {
        const double Tol = 1e-12;
        const int Samples = 20;
        var ranges = decomp.SectorRanges;
        int sectorCount = ranges.Count;
        if (sectorCount < 2) return;

        var rng = new Random(0);
        for (int sample = 0; sample < Samples; sample++)
        {
            int sA = rng.Next(sectorCount);
            int sB = rng.Next(sectorCount);
            while (sB == sA) sB = rng.Next(sectorCount);
            var rangeA = ranges[sA];
            var rangeB = ranges[sB];
            if (rangeA.Size == 0 || rangeB.Size == 0) continue;
            int f1 = decomp.Permutation[rangeA.Offset + rng.Next(rangeA.Size)];
            int f2 = decomp.Permutation[rangeB.Offset + rng.Next(rangeB.Size)];
            double mag = L[f1, f2].Magnitude;
            if (mag > Tol)
                throw new InvalidOperationException(
                    $"LiouvillianBlockSpectrum.ComputeSpectrum contract violation: L is NOT block-diagonal " +
                    $"in joint-popcount basis. Found |L[{f1}, {f2}]| = {mag:E3} between sectors " +
                    $"(p_c={rangeA.PCol}, p_r={rangeA.PRow}) and (p_c={rangeB.PCol}, p_r={rangeB.PRow}) " +
                    $"(tolerance {Tol:E0}). This routine requires popcount-conserving H; " +
                    $"non-conserving H (e.g., XX+YZ, XY+YX) silently returns wrong spectra. See class XML doc.");
        }
    }

    /// <summary>DEBUG-only sample-based check that H is popcount-conserving in the 2^N
    /// Hilbert basis. Picks 20 random index pairs (i, j) with <c>popcount(i) != popcount(j)</c>
    /// and asserts each <c>|H[i, j]| &lt; 1e-12</c>. Throws on the first violation. Stripped
    /// from RELEASE builds. <c>internal</c> rather than <c>private</c> so the matching guard
    /// in <see cref="F71MirrorBlockRefinement.ComputeSpectrumPerBlock(ComplexMatrix, IReadOnlyList{double}, int, PauliLetter)"/>
    /// can reuse the same routine instead of duplicating it.</summary>
    [Conditional("DEBUG")]
    internal static void DebugAssertPopcountConservingH(ComplexMatrix H, int N)
    {
        const double Tol = 1e-12;
        const int Samples = 20;
        int d = 1 << N;
        if (d < 2) return;

        var rng = new Random(0);
        int found = 0;
        for (int attempt = 0; attempt < Samples * 8 && found < Samples; attempt++)
        {
            int i = rng.Next(d);
            int j = rng.Next(d);
            int pi = BitOperations.PopCount((uint)i);
            int pj = BitOperations.PopCount((uint)j);
            if (pi == pj) continue;
            found++;
            double mag = H[i, j].Magnitude;
            if (mag > Tol)
                throw new InvalidOperationException(
                    $"LiouvillianBlockSpectrum.ComputeSpectrumPerBlock contract violation: H is NOT " +
                    $"popcount-conserving. Found |H[{i}, {j}]| = {mag:E3} between popcount {pi} and {pj} " +
                    $"(tolerance {Tol:E0}). This routine requires popcount-conserving H (XX+YY, ZZ, XXZ, " +
                    $"Heisenberg, sums of these). Non-conserving H (e.g., XX+YZ, XY+YX, X_iZ_j) silently " +
                    $"returns wrong spectra because L is no longer block-diagonal in the joint-popcount " +
                    $"basis. See class XML doc.");
        }
    }
}
