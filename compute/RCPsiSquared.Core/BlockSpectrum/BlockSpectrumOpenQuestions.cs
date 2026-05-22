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
/// yet. After the 2026-05-22 review (below) the frontier is the one genuinely unbuilt
/// engine, the full Prosen two-sided reduction; it is defined by what this repository has
/// and has not built, not by the published literature.</para>
///
/// <para>The first item, "F1 palindrome as a spectrum-halving builder step", was closed
/// 2026-05-22: Π acts on the joint-popcount sector labels as the whole-sector cycle
/// (p_c, p_r) ↦ (N − p_r, p_c), grouping the (N+1)² sectors into orbits of 4. The builder
/// now eigendecomposes one primary per orbit and derives three followers (the Π²-image by a
/// verbatim copy, the Π/Π³-images by the F1 reflection λ ↦ −2Σγ − λ); see
/// <see cref="SymmetryFamily.F1PalindromeOrbitPairing"/> and the orbit-pairing wiring in
/// <see cref="LiouvillianBlockSpectrum"/> and <see cref="F71MirrorBlockRefinement"/>.</para>
///
/// <para>Two further first-pass items were removed 2026-05-22 after verification against the
/// code. "F80 spectral map as a builder skip": F80's Spec(M) prediction is already wired, as
/// <c>RCPsiSquared.Diagnostics.F80.BlochSignWalk</c> (tested, consumed by MemoryAxisChain), so
/// it is not an unwired shortcut, and the builder never builds or diagonalises M, only L.
/// "F89 topology-orbit equivalence as sweep deduplication": F89's orbit dedup applies to
/// bond-configuration sweeps of the S(t) propagation observable, not to the
/// Liouvillian-spectrum builder this list scopes; its ChainTopologyClass orbit-label, though
/// unwired, belongs to the F89 family.</para></summary>
public static class BlockSpectrumOpenQuestions
{
    public static IReadOnlyList<OpenQuestion> Standard { get; } = new[]
    {
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
    };
}
