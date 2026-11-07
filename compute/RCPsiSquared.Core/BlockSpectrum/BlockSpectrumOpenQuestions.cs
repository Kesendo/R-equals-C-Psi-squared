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
/// engine, a Bethe-ansatz route to the two-sided sectors, which the literature offers and the
/// repository has neither built nor checked for completeness; it is defined by what this repository has
/// and has not built, not by the published literature.</para>
///
/// <para>The first item, "F1 palindrome as a spectrum-halving builder step", was closed
/// 2026-05-22: Π acts on the joint-popcount sector labels as the whole-sector cycle
/// (p_c, p_r) ↦ (N − p_r, p_c), grouping the (N+1)² sectors into orbits of 4. When H commutes
/// with X⊗N the builder's default pairing eigendecomposes one primary per orbit and derives
/// three followers (the Π²-image by a
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
            "The two-sided (m, m̃) sectors: the imaginary-U Bethe ansatz",
            "Prosen/OneSidedSectorClosedForm gives a genuine closed form for the one-sided " +
            "joint-popcount sectors (p_c = 0 or p_r = 0): the sector spectrum is the subset-sums " +
            "of N Prosen rapidities, with no diagonalisation. The two-sided sectors (m, m̃ ≥ 1), " +
            "which hold the bulk of the 4^N spectrum, have no closed form here apart from pieces, " +
            "among them the images of the one-sided sectors under the ladders Φ and S⁺ of " +
            "PROOF_FROZEN_BAND_SO4 (see the one-sided closed form of experiments/N10_BLOCK_SPECTRUM_PUSH.md; under equal light they " +
            "are F171's grey rooms, PROOF_GREY_GLASS_LADDER_IMAGES, typed as GreyRoomLadderImagesClaim, C(N, k) modes " +
            "at −2γk in every block whose window holds k). The Z-dephasing " +
            "term is quartic in the Jordan-Wigner fermions, and no Nambu reduction holds at every N: " +
            "at odd N the spectrum is not that of a generator quadratic in them (its even-parity " +
            "sector lacks the point symmetry such a generator forces). Medvedyeva, Essler and Prosen (2016) map the dephased " +
            "tight-binding model on any bipartite lattice, and through Jordan-Wigner the XX + YY " +
            "chain, onto a Hubbard model with imaginary interaction and, at uniform γ and J, solve the periodic fermion chain by the Bethe ansatz; the open " +
            "chain becomes the imaginary-interaction Hubbard chain with open ends. With a ZZ term " +
            "the map no longer yields the Hubbard model.",
            "Take the Bethe equations of the open imaginary-U Hubbard chain to the " +
            "two-sided sectors and cross-check them against the dense per-block spectra at " +
            "N = 3..6. Apply the obstruction-result discipline: establish whether their solutions " +
            "cover a whole sector, defective eigenvalues included, before assuming they close, the " +
            "way F86's g_eff was shown not to.",
            "compute/RCPsiSquared.Core/BlockSpectrum/Prosen/OneSidedSectorClosedForm.cs; " +
            "experiments/N10_BLOCK_SPECTRUM_PUSH.md (the one-sided closed form); " +
            "compute/RCPsiSquared.Core/BlockSpectrum/JordanWigner/JwSlaterPairBasis.cs " +
            "(the quartic dissipator)"),

        new OpenQuestion(
            "Builder-side F100 anti-palindromic-J spectral-invariance regression test",
            "F100 (Tier 1 derived; C1QPeakMirrorJParity) predicts that anti-palindromic per-bond " +
            "J profiles give zero c₁/Q_peak mirror-deviation, while palindromic J profiles " +
            "preserve the full F1 palindrome structure. Task A (commit 12f3edd) wired per-bond " +
            "J support into KleinFourGroupSelfPairedRefinement.BuildSubBlockL and " +
            "LiouvillianSectorSweep (plus PauliHamiltonian.XYChain overloads), so the Builder " +
            "can now run F100-territory J profiles directly. What is missing is a dedicated " +
            "Builder-side regression test that sweeps {palindromic, anti-palindromic, generic} " +
            "J profiles and asserts the predicted spectral invariance / asymmetry pattern, " +
            "exercising the per-bond plumbing end-to-end at N=4..6.",
            "Add a test (under compute/RCPsiSquared.Core.Tests/BlockSpectrum/) that builds " +
            "L via ComputeSpectrumPerBlock with palindromic J profiles (e.g. J = [1, 2, 1] " +
            "for N=4) and anti-palindromic profiles (e.g. J = [1, -1, 1] for N=4), then " +
            "verifies the F100-predicted c₁/Q_peak mirror-deviation pattern against an " +
            "external full-L reference at N=4..6. Effort: S. Payoff: S (mostly regression " +
            "coverage of the per-bond plumbing; F100 itself is Tier1Derived and verified " +
            "independently at N=3,4,5 via the framework package).",
            "compute/RCPsiSquared.Core/BlockSpectrum/KleinFourGroupSelfPairedRefinement.cs " +
            "(per-bond BuildSubBlockL overload; commit 12f3edd); " +
            "compute/RCPsiSquared.Core/Pauli/PauliHamiltonian.cs (XYChain per-bond overload); " +
            "compute/RCPsiSquared.Core/F71/C1QPeakMirrorJParity.cs (F100 typed claim); " +
            "framework Confirmations: F100 (C1QPeakMirrorJParity)"),
    };
}
