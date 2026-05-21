using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.BlockSpectrum;

/// <summary>Open computational items for the Liouvillian-spectrum builder: the BlockSpectrum
/// family (<see cref="LiouvillianBlockSpectrum"/> and the per-block, Jordan-Wigner, and
/// Prosen engines around it).
///
/// <para>This is the first typed OpenQuestion surface for <i>computational abbreviation</i>:
/// shortcuts that skip or cheapen building and diagonalising the 4^N Liouvillian. It is
/// deliberately a different kind of list from the physics-characterisation OpenQuestions of
/// F1, F71, F86, and the Π² Klein layer (which describe the shape of an answer, not the cost
/// of computing it). The builder's abbreviation frontier was, until this file, scattered
/// across <c>[Fact(Skip = "...")]</c> markers, "Phase 3" comments, and
/// <c>docs/SYMMETRY_FAMILY_INVENTORY.md</c>; collecting it here makes it a tracked discovery
/// channel, the way the F-family OpenQuestions are for the physics side.</para>
///
/// <para>The recurring shape of these items: a result already proven, sitting as a typed
/// claim or a proof, that is not yet wired into a compute path and so abbreviates nothing
/// yet. The frontier is our own under-exploited theorems (F80, F89 below) plus the one
/// genuinely unbuilt engine (the full Prosen two-sided reduction); it is defined by what this
/// repository has and has not built, not by the published literature.</para>
///
/// <para>The first item, "F1 palindrome as a spectrum-halving builder step", was closed
/// 2026-05-22: Π acts on the joint-popcount sector labels as the whole-sector cycle
/// (p_c, p_r) ↦ (N − p_r, p_c), grouping the (N+1)² sectors into orbits of 4. The builder
/// now eigendecomposes one primary per orbit and derives three followers (the Π²-image by a
/// verbatim copy, the Π/Π³-images by the F1 reflection λ ↦ −2Σγ − λ); see
/// <see cref="SymmetryFamily.F1PalindromeOrbitPairing"/> and the orbit-pairing wiring in
/// <see cref="LiouvillianBlockSpectrum"/> and <see cref="F71MirrorBlockRefinement"/>.</para></summary>
public static class BlockSpectrumOpenQuestions
{
    public static IReadOnlyList<OpenQuestion> Standard { get; } = new[]
    {
        new OpenQuestion(
            "F80 spectral map Spec(M) = ±2i·Spec(H) wired as a builder skip",
            "F80 proves a closed map from the Hamiltonian spectrum to the spectrum of the " +
            "residual operator M: for the chain Π²-odd two-body case, Spec(M) = ±2i·Spec(H), " +
            "verified bit-exact at N=3..7. It currently exists only as the proof " +
            "PROOF_F80_BLOCH_SIGNWALK.md and Python sweep scripts; no Core primitive computes an " +
            "M-spectrum by reading it off the cheap 2^N Hamiltonian. Anywhere the builder or a " +
            "diagnostic would diagonalise or SVD M directly, F80 makes that unnecessary in the " +
            "cases it covers: a proven shortcut that abbreviates nothing yet.",
            "Add a Core primitive returning Spec(M) from PauliHamiltonian eigenvalues for the " +
            "F80-covered class, and route the M-spectrum consumers through it. Pin down which " +
            "(k-body, Π²-parity, topology) cases F80 covers, so the skip is applied only where " +
            "the map is proven.",
            "docs/ANALYTICAL_FORMULAS.md F80; docs/proofs/PROOF_F80_BLOCH_SIGNWALK.md; " +
            "the F78-F80 M-decomposition family"),

        new OpenQuestion(
            "Full Prosen third-quantization: the two-sided (m, m̃) sectors",
            "Prosen/OneSidedSectorClosedForm gives a genuine closed form for the one-sided " +
            "joint-popcount sectors (p_c = 0 or p_r = 0): the sector spectrum is the subset-sums " +
            "of N Prosen rapidities, with no diagonalisation. This is the simplest leaf of the " +
            "third-quantization programme; the two-sided sectors (m, m̃ ≥ 1), which hold the " +
            "bulk of the 4^N spectrum, are not implemented. docs/SYMMETRY_FAMILY_INVENTORY.md " +
            "names the full reduction to a 2N × 2N Nambu spectral problem as 'Phase 3' of the " +
            "N=10 push. This is the one structurally transformative shortcut the repo knows it " +
            "lacks: it would replace dense per-block diagonalisation with a polynomial-cost route.",
            "Extend the one-sided closed form to the two-sided sectors via the third-quantization " +
            "/ Nambu reduction the Prosen leaf and SYMMETRY_FAMILY_INVENTORY point to. Apply the " +
            "obstruction-result discipline: characterise the reduction before assuming it closes, " +
            "the way F86's g_eff was shown not to.",
            "compute/RCPsiSquared.Core/BlockSpectrum/Prosen/OneSidedSectorClosedForm.cs; " +
            "docs/SYMMETRY_FAMILY_INVENTORY.md (Phase 3); " +
            "compute/RCPsiSquared.Core/BlockSpectrum/JordanWigner/JwSlaterPairF1PalindromeProbe.cs " +
            "(cites the 2N × 2N Nambu problem)"),

        new OpenQuestion(
            "F89 topology-orbit equivalence as sweep deduplication",
            "F89 (topology-orbit closure, F89TopologyOrbitClosure) establishes that the closure " +
            "observable depends on a bond configuration only through its topology orbit: " +
            "distinct bond configurations in the same orbit give bit-identical results. " +
            "Parameter sweeps over bond configurations (the F86 Q-scans, the per-orbit tables) " +
            "currently recompute orbit-equivalent configurations independently. F89 " +
            "orbit-equivalence is typed as a claim but is not used to deduplicate a sweep down " +
            "to one representative per orbit.",
            "Before a bond-configuration sweep, partition the configurations by topology orbit " +
            "via the F89 orbit machinery and compute one representative per orbit, copying the " +
            "rest. Quantify the orbit-count versus configuration-count reduction per N.",
            "docs/ANALYTICAL_FORMULAS.md F89; " +
            "compute/RCPsiSquared.Core/Symmetry/F89TopologyOrbitClosure.cs"),
    };
}
