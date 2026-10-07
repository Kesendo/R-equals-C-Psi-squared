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
///   <item>N=5: full 4^N = 1024 → max block 100, total cost ratio ≈ 212× faster</item>
///   <item>N=6: full 4^N = 4096 → max block 400, total cost ratio ≈ 298× faster</item>
///   <item>N=7: full 4^N = 16384 → max block 1225, total cost ratio ≈ 399× faster</item>
///   <item>N=8: full 4^N = 65536 → max block 4900, total cost ratio ≈ 515× faster</item>
/// </list>
/// At N=8 the largest block fits in ~0.38 GB vs ~68.7 GB for the full L, removing the
/// need for native-memory + ILP64 LAPACK on the dense path.</para>
///
/// <para><b>Witness at N=3, 4, 5.</b> The per-block spectrum and the direct
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
/// Pauli pairs like X_iZ_j or Y_iZ_j). For such an H the per-block eigenvalues would miss
/// the cross-sector entries of L, so <c>Σ_blocks ‖L_b‖²_F ≠ ‖L_full‖²_F</c> and the spectrum
/// would be incomplete. The 2026-05-18
/// F1 general-topology N=7 dogfood discovered this empirically (XX+YZ at N=5 gave block
/// sum 4403 vs dense 16691, factor ~3.8 off). <see cref="ComputeSpectrumPerBlock(ComplexMatrix, IReadOnlyList{double}, int)"/>
/// checks every entry of H between basis states of different popcount and throws an
/// <see cref="ArgumentException"/> unless each is exactly zero
/// (<see cref="RequirePopcountConservingH"/>, O(4^N) reads, negligible beside the
/// eigensolves); <see cref="ComputeSpectrum"/>, which receives L itself, keeps a
/// DEBUG-only sampled check (<see cref="DebugAssertBlockDiagonalL"/>). Both methods ensure MKL is initialized
/// via <see cref="MathNetSetup.EnsureInitialized"/> on entry, so no caller needs to
/// pre-initialize; the lazy global guard makes the redundant call free after the
/// assembly-level <c>CoreModuleInitializer</c> has already run.</para>
///
/// <para><b>Sector pairing, and when it is exact.</b> <c>ComputeSpectrumPerBlock</c> solves
/// one sector per orbit of the F1 mirror Π on joint-popcount labels,
/// (p_c, p_r) ↦ (N − p_r, p_c) (<see cref="SymmetryFamily.F1PalindromeOrbitPairing"/>), and
/// fills the other three: the Π²-image, the X⊗N partner, by copying, and the Π- and
/// Π³-images by the reflection λ ↦ −2σ − λ. Both rules hold whenever a unitary W = X^⊗N·D
/// with D diagonal commutes with H. W anticommutes with every Z jump, so ρ ↦ ρ·W carries L to
/// −L† − 2σ (the far-end relation of <c>experiments/THE_PAIRING_CONDITION.md</c>, written
/// there with U on the left); composed with F119's dagger ρ ↦ ρ†, which every Lindbladian
/// commutes with and which sends sector (p_c, p_r) to (p_r, p_c) and conjugates its spectrum,
/// the map ρ ↦ ρ†·W moves sector s onto Π(s), W moving sectors as X^⊗N does, so
/// spec(L|Π(s)) = −2σ − spec(L|s), and twice over the copy. With every site dephased at a
/// positive rate such a W exists exactly when the spectrum is palindromic (F158's far kernel
/// read on the hopping graph, Theorem 2 of
/// <c>docs/proofs/PROOF_PALINDROME_COMPLEMENT_CONNECTION.md</c>, whose flat section is D), so
/// there the fill is exact wherever the spectrum is palindromic, and at odd N, where no
/// sector is Π-fixed, only there; at even N the Π-fixed (N/2, N/2) sector is solved on its
/// own, and at every even N it can break the palindrome under an exact fill. Add the split
/// a·(|0101…⟩⟨0101…| − |1010…⟩⟨1010…|), a ≠ 0, to the uniform chain Σ(XX + YY): it sits at
/// popcount N/2 alone, and ρ ↦ ρ†·X^⊗N lets X^⊗N act on a sector's ket side only, so at any rates
/// the reflection still holds from every sector whose ket popcount is not N/2; at rates
/// symmetric under the chain reflection R, R·X^⊗N fixes both Néel states, commutes with H
/// and carries the copy, and every orbit of four keeps all its relations. No unitary X^⊗N·D commutes
/// with H, as it would move the diagonal entry a onto the −a at the other Néel state, so at
/// positive rates the spectrum is not palindromic (Theorem 2) and, the orbits of four being
/// palindromic, neither is that sector: 1.385 in Hausdorff distance from its own reflection
/// at N = 4 with a = 2.8 and rates 0.5, 0.955 at N = 6. At N = 2 the split with a = 1.4 is
/// the field 0.7·(Z₀ − Z₁).
/// The method checks the
/// constant section D = 1, [H, X^⊗N] = 0 (F63's parity), exactly on the entries of H
/// (<see cref="CommutesWithXN"/>); it suffices at any rates, zero included, and every sector
/// is solved when it fails. The truly Hamiltonians (every Pauli string with #Y and #Z even,
/// H real and X^⊗N-symmetric) are the real ones among them. A nonzero Z field added to an
/// X^⊗N-symmetric H fails the check, X^⊗N flipping its sign; filled anyway, the spectrum of
/// XX + YY + 0.7·Z on one site at N = 2, rates 0.3 and 0.4, lies 0.73 in Hausdorff distance
/// from the full eigensolver's.
/// Flux Φ through a ring, written with complex hopping phases, keeps the copy, X^⊗N
/// carrying such an H to Hᵀ with L(Hᵀ) = L(H)ᵀ, and fails the check. A palindromic H that
/// commutes with
/// X^⊗N only in a diagonal frame, such as an XY − YX bond on a chain (the frame theorem of
/// <c>experiments/THE_PALINDROME_AS_A_COLOURING.md</c>), pairs exactly as well but fails the
/// check and is solved sector by sector. <see cref="SectorPairing"/> chooses the rule: the
/// Π-orbit pairing (the default), the X⊗N copy alone, under which a sector and its Π-image
/// are solved separately so that a palindrome check on the output stays a check, or none.
/// F108 Part 1's Π_5bilinear runs the same sector orbits in the opposite direction
/// (Π_5bilinear = Π_Z ∘ Ad_{Y^⊗N}); the popcount-conserving members of its family,
/// combinations of (XX+YY) and ZZ, commute with X^⊗N, so it adds nothing at this layer.</para>
///
/// <para><b>Dephasing letter.</b> The joint-popcount sector basis is tied to Z-dephasing because
/// <see cref="PerBlockLiouvillianBuilder.BuildBlockZ"/> builds the dissipator element-wise from the
/// computational-basis bit-parity disagreement, which requires the dephase letter to be diagonal in that
/// basis. X- and Y-dephasing reach it through the overload taking a <see cref="PauliHamiltonian"/>: the
/// letters are turned exactly so that the dephasing letter becomes Z (<see cref="LetterTurn"/>, a
/// per-site Clifford acting as a signed letter permutation), and the turned Hamiltonian runs the Z path
/// with its exact popcount check and pairing guard; the spectrum is unchanged by the turn. The dense
/// overloads refuse X and Y with <see cref="NotSupportedException"/>, since a dense H cannot be turned
/// without rounding residue across popcounts, which the exact check refuses, and refuse
/// <see cref="Pauli.PauliLetter.I"/> with <see cref="ArgumentException"/> (not a valid dephase letter).
/// Under X-dephasing the guard reads [H, Z^⊗N] = 0 in the frame of H, under Y-dephasing [H, X^⊗N] = 0.</para>
///
/// <para>Anchors: <c>compute/RCPsiSquared.Core/BlockSpectrum/JointPopcountSectors.cs</c>
/// (parent Claim, block-diagonal structure), <c>compute/RCPsiSquared.Core/BlockSpectrum/JointPopcountSectorBuilder.cs</c>
/// (basis permutation + sector ranges), <c>compute/RCPsiSquared.Core.Tests/BlockSpectrum/LiouvillianBlockSpectrumTests.cs</c>
/// (spectral equality verification at N=3, 4, 5, to 1e-9),
/// <c>compute/RCPsiSquared.Core.Tests/BlockSpectrum/BlockSpectrumPairingGuardTests.cs</c> (the
/// pairing guard, and both engines read against the full eigensolver and the pairing rule
/// under the error model at three scales).</para></summary>
public sealed class LiouvillianBlockSpectrum : Claim
{
    private readonly JointPopcountSectors _sectors;

    /// <summary>Per-block eigensolver selection knob for
    /// <see cref="ComputeSpectrumPerBlock(ComplexMatrix, IReadOnlyList{double}, int, EigenPath)"/>.
    /// The default <see cref="Auto"/> picks <see cref="MathNet"/> for blocks below
    /// <see cref="Lp64ComplexCeiling"/> and <see cref="MklDirectNative"/> above; the explicit
    /// overrides exist only for the parity-witness test
    /// <c>PerBlockLiouvillianBuilderNativeMemoryParityTests</c> which forces both paths on the
    /// same small block to demonstrate their agreement.</summary>
    public enum EigenPath
    {
        /// <summary>Auto-select by block size: MathNet for size ≤ <see cref="Lp64ComplexCeiling"/>,
        /// MklDirect + NativeMemory + ILP64-aware LAPACK above. The production default.</summary>
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

    /// <summary>How <c>ComputeSpectrumPerBlock</c> derives sector spectra from one another.
    /// Both pairings are used only when H commutes exactly with X^⊗N
    /// (<see cref="CommutesWithXN"/>); otherwise every sector is solved, whatever is
    /// requested.</summary>
    public enum SectorPairing
    {
        /// <summary>One eigensolve per orbit of the F1 mirror Π on sector labels: the
        /// Π²-image copied, the Π- and Π³-images reflected by λ ↦ −2σ − λ: about a quarter of
        /// the eigensolves, exactly a quarter at odd N. A palindrome check on its output passes by construction at odd N,
        /// where every sector lies in an orbit of four, and at even N on every sector but the
        /// Π-fixed (N/2, N/2).</summary>
        PiOrbit = 0,
        /// <summary>One eigensolve per X⊗N pair, the partner copied: about half the
        /// eigensolves, exactly half at odd N.
        /// A sector and its Π-image are solved separately, so a palindrome check on the output
        /// is a check.</summary>
        XNCopy = 1,
        /// <summary>Every sector solved.</summary>
        None = 2,
    }

    /// <summary>Largest square Complex block (size n × n with n² ≤ 134 217 728) whose total
    /// byte count (n² × 16) fits inside the LP64 2 GB single-native-array marshalling ceiling
    /// enforced by MathNet's <c>MklLinearAlgebraProvider.EigenDecomp</c>. n = 11 585 gives
    /// n² × 16 ≈ 2.147 GB (just inside the marshaller's <c>int.MaxValue</c> byte limit when the
    /// fixed-pinned <c>Complex[]</c> is rounded into a single P/Invoke array). Blocks at or
    /// below this size are routed through MathNet's well-tested managed wrapper; blocks above
    /// are routed through <see cref="MklDirect.EigenvaluesOnlyNative"/> on a
    /// <see cref="NativeMemory.AllocZeroed(nuint)"/>-backed buffer, which also auto-selects
    /// ILP64 when n &gt; 46 340 (n² &gt; <c>int.MaxValue</c>; first reached at N = 10, whose
    /// central block has C(10,5)² = 63 504 rows).
    ///
    /// <para>The threshold is the same value the N=9 test
    /// (<c>F1GeneralTopologyN9BlockSpectrumChainTests.Lp64EvdSquareMatrixCeiling</c>) uses for
    /// its pre-flight check; keeping both in sync means a future bump (e.g. if MathNet relaxes
    /// the marshaller) only needs updating one constant per file.</para></summary>
    public const int Lp64ComplexCeiling = 11_585;

    public LiouvillianBlockSpectrum(JointPopcountSectors sectors)
        : base("LiouvillianBlockSpectrum: per-block eig over (N+1)² joint popcount sectors yields the same spectrum (multiset) as direct full-L eig; verified at N=3,4,5 to 1e-9.",
               Tier.Tier1Derived,
               "JointPopcountSectors block-diagonality (parent) + per-block diagonalisation; verified vs full-L eig at N=3,4,5 to 1e-9 in LiouvillianBlockSpectrumTests")
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
    /// NativeMemory + ILP64-aware path), eigendecomposed, and discarded before the next block.
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
    /// <para>Uses the F1 Π-orbit pairing when H commutes with X^⊗N: the F1 palindrome
    /// conjugation Π is order-4 and on joint-popcount labels acts as the whole-sector cycle
    /// (p_c, p_r) ↦ (N − p_r, p_c), grouping the (N+1)² sectors into orbits of 4. Only one
    /// "primary" sector per orbit is eigendecomposed; the three followers are derived from
    /// it, the Π²-image (X⊗N partner) by a verbatim spectrum copy and the Π/Π³-images by the
    /// F1 reflection λ ↦ −2·Σγ − λ. This quarters the eigendecomposition count, where the
    /// X⊗N copy alone (<see cref="SectorPairing.XNCopy"/>, the pairing of
    /// <see cref="SymmetryFamily.XGlobalChargeConjugationPairing"/>, which Π² equals) halves
    /// it. When H does not commute with X^⊗N every sector is solved;
    /// the class summary gives the reason, and <see cref="SectorPairing"/> the
    /// alternatives.</para></summary>
    /// <param name="H">Hilbert-space Hamiltonian, dense 2^N × 2^N (cheap even at N=9: 512×512).
    /// Must be popcount-conserving; checked exactly, <see cref="ArgumentException"/>
    /// otherwise.</param>
    /// <param name="gammaPerSite">Per-site Z-dephasing rates (length N).</param>
    /// <param name="N">Qubit count.</param>
    /// <returns>Flat array of 4^N eigenvalues, concatenated block-by-block.</returns>
    public static Complex[] ComputeSpectrumPerBlock(ComplexMatrix H, IReadOnlyList<double> gammaPerSite, int N) =>
        ComputeSpectrumPerBlock(H, gammaPerSite, N, EigenPath.Auto);

    /// <summary>Overload choosing the sector pairing. <see cref="SectorPairing.XNCopy"/> keeps
    /// a palindrome check on the output a check, at about twice the eigensolves of the default,
    /// exactly twice at odd N;
    /// <see cref="SectorPairing.None"/> solves every sector.</summary>
    public static Complex[] ComputeSpectrumPerBlock(
        ComplexMatrix H, IReadOnlyList<double> gammaPerSite, int N, SectorPairing pairing) =>
        ComputeSpectrumPerBlock(H, gammaPerSite, N, EigenPath.Auto, PauliLetter.Z, pairing);

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

    /// <summary>Overload with an explicit dephase letter. Only <see cref="PauliLetter.Z"/>
    /// is supported; the letter is validated and does not change the computation, and no
    /// F108 operator (<see cref="Pi5BilinearOperator"/>,
    /// <see cref="F108Part1Pi2EvenAlwaysPalindromic"/>) is used: the popcount-conserving
    /// members of F108's Z family commute with X^⊗N and pair under the same rule as the
    /// canonical case.
    ///
    /// <para><b>Z-only restriction.</b> The per-block Liouvillian construction in
    /// <see cref="PerBlockLiouvillianBuilder.BuildBlockZ"/> is hardcoded to Z-dephasing
    /// (the dissipator is built element-wise in the computational basis from the Hilbert-
    /// side bit-parity disagreement, not from a general P_l ⊗ P_l⁺ Kronecker construction).
    /// X- and Y-dephasing (F108 Parts 2 and 3) throw <see cref="NotSupportedException"/> on this dense
    /// entry point; the overload taking a <see cref="PauliHamiltonian"/> serves them by an exact letter
    /// turn (X- and Y-dephasing break popcount conservation in the computational basis).
    /// <see cref="Pauli.PauliLetter.I"/> is rejected up-front with
    /// <see cref="ArgumentException"/> (not a valid dephase letter).</para>
    ///
    /// <para>Mismatch handling: <paramref name="dephaseLetter"/> = X or Y throws cleanly
    /// rather than silently producing wrong eigenvalues (the underlying per-block builder
    /// would otherwise stamp Z-deph entries onto an X- or Y-deph problem).</para></summary>
    /// <param name="dephaseLetter">The dephase letter. Only
    /// <see cref="PauliLetter.Z"/> is currently supported by the per-block construction; X and
    /// Y throw <see cref="NotSupportedException"/> (design-permanent under the current basis),
    /// I throws <see cref="ArgumentException"/> (not a valid dephase letter).</param>
    public static Complex[] ComputeSpectrumPerBlock(
        ComplexMatrix H, IReadOnlyList<double> gammaPerSite, int N, EigenPath path,
        PauliLetter dephaseLetter) =>
        ComputeSpectrumPerBlock(H, gammaPerSite, N, path, dephaseLetter, SectorPairing.PiOrbit);

    /// <summary>X-, Y- or Z-dephasing on every site for a Hamiltonian given as Pauli terms. The
    /// letters are turned exactly so that the dephasing letter becomes Z (<see cref="LetterTurn"/>),
    /// and the turned Hamiltonian goes to the Z path, which checks popcount conservation exactly and
    /// pairs sectors under its own guard. The spectrum is the one of the J-dephased Lindbladian of
    /// <paramref name="H"/>, since the turn is a unitary that carries the jump to Z.
    ///
    /// <para>Read in the frame of <paramref name="H"/>, the two conditions become: the turned H
    /// conserves popcount exactly when H commutes with Σ_l J_l, the total magnetisation along the jump
    /// letter (which implies, but is not implied by, [H, J^⊗N] = 0), and the pairing guard,
    /// [H′, X^⊗N] = 0 after the turn, reads [H, Z^⊗N] = 0 under X-dephasing and [H, X^⊗N] = 0 under
    /// Y-dephasing; given [H, J^⊗N] = 0 either reading says that H commutes with all three letter
    /// strings. The turn itself is exact, a relabelling with signs; the dense matrix is then built by
    /// <see cref="PauliHamiltonian.ToMatrix"/>, whose term-by-term sum can leave a residue of one
    /// rounding on a popcount-changing entry where terms cancel in exact arithmetic but not in the
    /// order they are summed. The exact check then refuses rather than return a wrong spectrum, as on
    /// the Z path. The dense overloads refuse X and Y, since a dense H cannot be turned at all without
    /// rounding.</para></summary>
    public static Complex[] ComputeSpectrumPerBlock(
        PauliHamiltonian H, IReadOnlyList<double> gammaPerSite, PauliLetter dephaseLetter,
        SectorPairing pairing = SectorPairing.PiOrbit, EigenPath path = EigenPath.Auto)
    {
        if (H is null) throw new ArgumentNullException(nameof(H));
        if (dephaseLetter == PauliLetter.I)
            throw new ArgumentException(
                "PauliLetter.I is not a valid dephase letter (the Lindblad dissipator requires a non-identity operator).",
                nameof(dephaseLetter));
        var turned = LetterTurn.Turn(H, dephaseLetter);
        return ComputeSpectrumPerBlock(turned.ToMatrix(), gammaPerSite, H.N, path, PauliLetter.Z, pairing);
    }

    /// <summary>The full overload: eigensolver path, dephase letter (Z only) and sector
    /// pairing. The pairing is used only when H commutes exactly with X^⊗N.</summary>
    public static Complex[] ComputeSpectrumPerBlock(
        ComplexMatrix H, IReadOnlyList<double> gammaPerSite, int N, EigenPath path,
        PauliLetter dephaseLetter, SectorPairing pairing)
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
                $"The dense ComputeSpectrumPerBlock overloads only support Z-dephasing in the joint-popcount " +
                $"basis (PerBlockLiouvillianBuilder.BuildBlockZ is hardcoded to Z, which is diagonal in " +
                $"the computational basis and popcount-conserving); got {dephaseLetter}. X- and Y-dephasing " +
                "break the joint-popcount sector structure that JointPopcountSectors relies on; see " +
                "the PauliHamiltonian overload, which turns the Pauli letters exactly (LetterTurn) so that the " +
                "dephasing letter becomes Z; a dense H cannot be turned without rounding residue that the exact " +
                "popcount check refuses.");
        RequirePopcountConservingH(H, N);

        // Reflection constant: the genuine Σ of per-site Z-dephasing rates (NOT N·γ), so
        // non-uniform γ stays exact. When H commutes with X^⊗N a sector's spectrum maps to its
        // Π-image's via λ ↦ −2·Σγ − λ (class summary); F1's identity is the truly case.
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

        // Which sectors are solved and how the rest are derived (PartitionSectors): with
        // [H, X^⊗N] = 0 the requested pairing, otherwise every sector. Under the Π-orbit
        // pairing one sector per orbit of Π, (p_c, p_r) ↦ (N − p_r, p_c), is solved (plus the
        // Π-fixed (N/2, N/2) at even N) and Phase 3 derives the Π²-image by a verbatim copy and
        // the Π/Π³-images by the F1 reflection λ ↦ −2·Σγ − λ; both rules are exact when H
        // commutes with X^⊗N (class summary). Primaries are sorted descending by size so the
        // largest sector starts first under Parallel.ForEach.
        var (primarySectorIndices, followerToPrimary) =
            PartitionSectors(H, N, decomp.SectorRanges, pairing);

        // BLAS-oversubscription strategy (c): outer DOP ≈ ProcessorCount/4 leaves room for
        // MKL inside the largest sectors' Evd. See ComputeSpectrum for the rationale.
        //
        // EigenPath.MklDirectNative serialises (DOP=1) because each block holds its whole
        // NativeMemory buffer while LAPACK runs: N = 9's four largest blocks hold about 4 GB
        // each, and at N = 10, whose central block holds about 64 GB and the next ones about
        // 45 GB, three at once would pass 128 GB. Under EigenPath.Auto the outer DOP stays
        // ProcessorCount/4 for every block, those routed to MklDirect included.
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
        // primary's: an X⊗N-image follower (Π²-image) copies them verbatim (exact, since the
        // partition pairs only when H commutes with X^⊗N); a Π/Π³-image follower reflects
        // each λ through the F1 palindrome map
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
            // precisely the marshaller bypass we need. Above 46 340 (first reached at N = 10,
            // whose central block has 63 504 rows) it switches to ILP64 OpenBLAS automatically.
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
        $"per-block eig multiset equals the full-L spectrum at N=3,4,5 to 1e-9; cubic-cost speedup ≈ 515× at N=8 ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("parent",
                summary: "JointPopcountSectors (block-diagonal structure)");
            yield return new InspectableNode("witness",
                summary: "spectral equality vs full-L eig at N=3, 4, 5 (|Δλ| < 1e-9)");
            yield return new InspectableNode("cubic-cost speedup",
                summary: "N=5: ≈ 212×, N=6: ≈ 298×, N=7: ≈ 399×, N=8: ≈ 515×. These are FLOP ratios from the cost model above, (4^N)³ / Σ(C(N,p)·C(N,q))³; measured wall-clock speedups are a different and much smaller quantity, around 8× to 10× per BlockSpectrumPerformanceWitness");
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

    /// <summary>Checks that H is popcount-conserving in the 2^N Hilbert basis: every entry
    /// between basis states of different popcount must be exactly zero, or the joint-popcount
    /// blocks miss part of L and the spectrum is wrong. Exact, on every entry, always on;
    /// throws <see cref="ArgumentException"/> on the first nonzero one. <c>internal</c> so that
    /// <see cref="F71MirrorBlockRefinement"/> uses the same check.</summary>
    internal static void RequirePopcountConservingH(ComplexMatrix H, int N)
    {
        int d = 1 << N;
        for (int i = 0; i < d; i++)
        {
            int pi = BitOperations.PopCount((uint)i);
            for (int j = 0; j < d; j++)
            {
                if (BitOperations.PopCount((uint)j) == pi) continue;
                var h = H[i, j];
                if (h != Complex.Zero)
                    throw new ArgumentException(
                        $"H is not popcount-conserving: H[{i}, {j}] = {h} between popcount {pi} and " +
                        $"{BitOperations.PopCount((uint)j)}. The joint-popcount blocks would miss this " +
                        $"entry of L and return a wrong spectrum (XX+YY, ZZ, XXZ, Heisenberg, Z fields and " +
                        $"sums of these conserve popcount; XX+YZ, XY+YX, X_iZ_j do not).", nameof(H));
            }
        }
    }

    /// <summary>True when H commutes with X^⊗N, checked exactly on the entries:
    /// H[i, j] == H[ī, j̄] with ī the bitwise complement. Under this condition the sector
    /// pairing of <c>ComputeSpectrumPerBlock</c> is exact (class summary); it holds for the
    /// truly Hamiltonians (XX+YY, ZZ, XXZ, Heisenberg with real couplings) and fails for a
    /// Z field, an XY − YX bond, XZX + YZY or ZZZ. An H whose entries are X^⊗N-symmetric only
    /// up to rounding fails the check and is solved sector by sector.</summary>
    public static bool CommutesWithXN(ComplexMatrix H, int N)
    {
        if (H is null) throw new ArgumentNullException(nameof(H));
        int d = 1 << N;
        if (H.RowCount != d || H.ColumnCount != d)
            throw new ArgumentException($"H must be {d}×{d} for N = {N}, got {H.RowCount}×{H.ColumnCount}.", nameof(H));
        int mask = d - 1;
        for (int i = 0; i < d; i++)
            for (int j = 0; j < d; j++)
                if (H[i, j] != H[i ^ mask, j ^ mask]) return false;
        return true;
    }

    /// <summary>The sectors <c>ComputeSpectrumPerBlock</c> solves (sorted descending by size)
    /// and, for each other sector, the solved sector it is derived from and how. The requested
    /// pairing applies only when H commutes exactly with X^⊗N; otherwise every sector is
    /// solved. Shared with <see cref="F71MirrorBlockRefinement"/>.</summary>
    internal static (List<int> Primaries, Dictionary<int, F1PalindromeOrbitPairing.F1Follower> FollowerToPrimary)
        PartitionSectors(ComplexMatrix H, int N, IReadOnlyList<JointPopcountSectorBuilder.SectorRange> sectors,
            SectorPairing pairing)
    {
        if (!Enum.IsDefined(pairing))
            throw new ArgumentOutOfRangeException(nameof(pairing), pairing, "undefined SectorPairing");
        if (pairing != SectorPairing.None && CommutesWithXN(H, N))
        {
            if (pairing == SectorPairing.PiOrbit)
                return F1PalindromeOrbitPairing.PartitionByPiOrbit(
                    N, sectors, s => (s.PCol, s.PRow), s => s.Size);
            var (xnPrimaries, xnFollowers) = XGlobalChargeConjugationPairing.PartitionByXNPairing(
                N, sectors, s => (s.PCol, s.PRow), s => s.Size);
            var copies = new Dictionary<int, F1PalindromeOrbitPairing.F1Follower>(xnFollowers.Count);
            foreach (var (follower, primary) in xnFollowers)
                copies[follower] = new F1PalindromeOrbitPairing.F1Follower(
                    primary, F1PalindromeOrbitPairing.F1FollowerKind.XnCopy);
            return (xnPrimaries, copies);
        }
        var all = Enumerable.Range(0, sectors.Count).ToList();
        all.Sort((a, b) => sectors[b].Size.CompareTo(sectors[a].Size));
        return (all, new Dictionary<int, F1PalindromeOrbitPairing.F1Follower>());
    }
}
