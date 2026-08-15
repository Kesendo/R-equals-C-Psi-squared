using RCPsiSquared.Core.BlockSpectrum;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>F153, the PINNING CRITERION: a joint-popcount block (p, q) of the Liouvillian sits ENTIRELY on
/// its Absorption-Theorem floor Re λ = −2γ|p − q| exactly when min(p, q) = 0 or max(p, q) = N, which is
/// 4N of the (N+1)² blocks. Live witness: <see cref="PinnedBlockFloorWitness"/>
/// (<c>inspect --root pinned</c>), which recomputes the whole sweep entry-wise at inspect time.
///
/// <para><b>The word, and it is chosen against two others.</b> The FLOOR is the Absorption Theorem's, the
/// bound every block obeys. FLAT is what three older documents call a block that attains it throughout.
/// PINNED is the cause: one index held at 0 or N, leaving a single ket or a single bra, so every cell of
/// the block carries the SAME disagreement count. This claim is named for the cause because the criterion
/// is about which blocks have it, not about which happen to be flat.</para>
///
/// <para><b>The forward direction is not new here, and the citation is the point.</b>
/// <c>docs/proofs/PROOF_CODIM1_BY_ADDITIVITY.md</c> §6 has held it since 2026-07-03: its Edge lemma is the
/// |p−q| = 1 case ("A = −2γ·I is scalar and L(q) = −2γ·I + q·C is a normal pencil"), and its
/// boundary-block sentence is the general one, "every boundary block (p or q̃ ∈ {0,N}) has constant n_diff,
/// hence a normal pencil whose Re-spectrum sits exactly on its even rung". The registry entry was minted
/// without that citation; this claim carries it.</para>
///
/// <para><b>The converse, and it is exact rather than argued.</b> The registry entry first reached the ⟸
/// by saying that off the pinned blocks "both sides can hop, S is no longer the whole story, and the real
/// parts leave it", a mechanism sketch that forces nothing; the entry now carries the proof instead, and
/// this paragraph is where it was written. The proof is one line of trace, recorded
/// on <see cref="ConverseByTrace"/>: Re Tr L = Σ_cells (−2γ·n_diff) and Re Tr L = Σ_λ Re λ, so if every
/// Re λ equals the floor −2γ·n_min then Σ(n_diff − n_min) = 0 over the cells, every term of which is
/// non-negative, forcing n_diff ≡ n_min. Constant n_diff is then the criterion by the window
/// combinatorics of the same §6: n_diff runs over [|p−q|, min(p+q, 2N−p−q)], and that interval has zero
/// width exactly when min(p, q) = 0 or max(p, q) = N. So the ⟺ closes without an eigensolver and without
/// a tolerance, at uniform γ.</para>
///
/// <para><b>What it locates.</b> F1's pairing distance projects onto rates, and that projection is blind
/// to a break that moves only a frequency (pinned in <c>F1SpectrumStatisticsTests</c> with a constructed
/// spectrum). The criterion turns exhibited into LOCATED: on these 4N blocks a rate reading loses nothing
/// at any coupling, and on the other (N−1)² it does.</para>
///
/// <para><b>Scope, and the first fence is the load-bearing one.</b> UNIFORM γ IS REQUIRED. Pinning fixes
/// the disagreement COUNT, not the disagreeing SITES: under a profile a cell's real part is −2·Σ over the
/// sites where it disagrees, which varies cell to cell, so the mechanism fails and the block leaves its
/// FLOOR. Three precisions, because the loose form of this fence is wrong in two directions.
/// (i) The mechanism's failure is entry-wise and needs nothing but the profile:
/// <see cref="PinnedBlockFloorWitness.ProfileRealDiagonalSpread"/> reads it, and the test gates it.
/// (ii) LEAVING THE FLOOR IS NOT THE SAME AS SPREADING. In the PAULI book H = J·Σ(XX+YY+Δ·ZZ), which is
/// <c>WeightCoherenceBlock</c>'s book and the one these numbers were measured in, at γ = linspace(0.1, 1, N),
/// J = 1, Δ = 1 the (0,1) block's spectral spread is 0.4903 at N=4 and 0.9477 at N=5; at Δ = 0 the same
/// profile gives machine zero at J = 1, opening only as J falls (1.77 at J = 0.05), which is the
/// strong-coupling coalescence F152 records, though see the flat-set law below: an OFF-locus profile does
/// not coalesce there at all. And at N=4, Δ=1 eight of the sixteen pinned blocks stay FLAT under the
/// profile, but only FOUR of those sit at the wrong constant, block (0,2) at −2.200 against a floor of
/// −1.000; the other four are the one-cell blocks, which are flat AT their floor.
/// On a PINNED block that constant is the TRACE's, −2·|p−q|·γ̄, and NOT −2·mean(γ), which is −1.100 for
/// this profile; the two coincide at |p−q| = 1, which is where the older reading of this sentence came
/// from. FLAT is not ON THE FLOOR: a flattened block sits at −2·|p−q|·γ̄, strictly below its floor.
/// WHICH blocks flatten is NOT a |p−q| parity, and two earlier readings of it were wrong. Flattening needs
/// an antilinear involution about that constant, from either of two sources: the block's FREE INDEX being
/// N/2, where the complement of the free config stays inside the block and sends rate ↦ Σγ − rate under
/// ANY profile (four blocks at even N, none at odd N; |p−q| = 2 at N=4 but 3 at N=6, which is what killed
/// the parity reading); or γ ANTI-PALINDROMIC, the R₉₀ locus (F91), where the chain reflection does the
/// same for every pinned block and all 4N flatten, and where both the ramp and the gates' dyadic profile
/// happen to lie. What the involution buys at once is a PALINDROME of the Re-spectrum about that constant,
/// holding at EVERY coupling; flatness is its UNBROKEN case and the J threshold is that onset, which is
/// why the threshold is a phenomenon and not a fitted number. Measured N = 3..6 at Δ = 1. What Δ decides is whether
/// the collapse has an onset at all, and that reading is OFF-LOCUS, a clause missing here until 2026-08-15:
/// OFF the locus at Δ = 0 the palindrome holds while the span saturates, and at Δ = 0.5 it saturates at N = 6 but
/// not at N = 4 (measured on γ_l = (l+1)²/16, spread at J = 10³ and 10⁵: (0,2) at N = 4 reads 0.2236 twice at Δ = 0
/// and 0.0000 at Δ = 0.5; (0,3) at N = 6 reads 0.5325 at Δ = 0 and 0.6456 at Δ = 0.5, each J-stable;
/// measured, not gated, the theory below using a different off-locus profile; WHY Δ and N select in
/// those OFF-locus readings is explained in experiments/THE_SPREAD_IS_A_RESONANCE.md, the saturation
/// being a resonance of the sector difference spectrum and the N clause its resonance count, the
/// Δ = 0.5 half-filling double levels numbering 0 / 1 / 8 at N = 4 / 6 / 8). ON the locus every pinned
/// block of dimension above one collapses, past its onset, at every Δ measured.
/// The criterion is about attaining the FLOOR, and a carrier that tested only
/// flatness would pass on the flat ones.
/// (iii) The spectral half is GATED since 2026-08-14, and it is gated EXACTLY, because it turned out not to
/// need an eigensolver: under Z-dephasing, Herm(L_block) = −2·diag(rate) for any Hermitian H and any
/// profile, and being diagonal it decides the real parts by itself, with the trace supplying the
/// converse; the spectral reading still needs this claim's number-conserving H, without which the index
/// set is not an invariant block and its submatrix eigenvalues are not the Liouvillian's. See
/// <see cref="PinnedBlockFloorWitness.ProfileHermitianPartResidual"/>. The committed Python verifier
/// <c>simulations/offdiag_sector_floor.py</c> still takes a SCALAR γ; what it does not cover is now covered
/// in C# rather than merely measured.
/// Beyond γ, the criterion needs the Absorption Theorem (any Hermitian H, Z-dephasing,
/// any graph, any N) plus a NUMBER-CONSERVING Hamiltonian, which is what makes the blocks invariant at
/// all; a transverse field or amplitude damping destroys the grading and the criterion with it. A
/// longitudinal field keeps both and leaves the criterion alone (on a pinned block at uniform γ, the case
/// where the criterion holds, it moves only S), but NOT because it is diagonal in the computational basis
/// (an imaginary diagonal moves the real parts of a non-normal matrix in general) and NOT because −i·ad_H is
/// anti-Hermitian (that proves too much: at uniform γ = 0.5, Δ = 0, J = 1.3, h = (0.5, 0.25, 0.75, 0.125)
/// a field takes the N = 4 block (2,2) from Re range [−4, 0] to [−3.868, 0]; the value is
/// parameter-dependent, J = 1.0 giving −3.805). It is the SCALAR dissipator, one size class AND uniform γ
/// together, which is exactly a pinned block at uniform γ; one size class alone does not suffice, a profile
/// making the rates differ within the class so that the same field moves even a FLAT pinned block (the
/// N = 4 block (0,2) on locus profile (0.25, 0.75, 0.5, 1.0) at Δ = 0, J = 1.25 goes from −2.5000 flat to
/// [−2.7113, −2.2887], J-dependent too); only there does the anti-Hermitian argument close.
/// Corrected 2026-08-15. Derived
/// graph-blind, GATED on the chain only.</para>
///
/// <para><b>The correction that travels with the carrier.</b> The pinned block is c·Id + i·S with S the
/// hops PLUS the ZZ diagonal, because the ZZ term puts a configuration-dependent IMAGINARY part on the
/// diagonal. MirrorWorld's gate certifies the XY-only form c·Id + i·J·A and passes bit-exactly because its
/// <c>Pair.Rate</c> carries no ZZ; a Heisenberg carrier reading c·Id would be reading the wrong chain. The
/// witness reads the real and imaginary diagonal spreads separately for exactly this reason.</para>
///
/// <para><b>Two fences on the FORM, which are not fences on the criterion.</b> REAL q: at complex q the
/// hops put ∓2·Im(q) on the real part off the diagonal and the ZZ term puts +Im(q)·Δ·Δzz on it, so the
/// pinned block leaves its floor and can even gain (block (0,1) at N=4, Δ=1, q = 1+0.3i reaches
/// Re λ = +0.0485 against a floor of −2). PROOF_CODIM1 §6 states its Edge lemma "at real q" for this
/// reason, and the criterion inherits the qualifier. REAL H: "S real symmetric" is the chain's form, not
/// the criterion's. Give the same number-conserving Hermitian H a complex hopping amplitude (a Peierls
/// phase) and the criterion still holds exactly, flat to 5e-15 on the pinned blocks at N=4, but S is
/// HERMITIAN rather than real symmetric. The general form is <b>the Hermitian part of L is scalar</b>;
/// real symmetric is what a REAL H buys, and it is what the entry-wise certificate below tests. So the
/// criterion is graph-blind and the certificate, in the form gated here, is not.</para>
///
/// <para>Typed parents: <see cref="AbsorptionTheoremClaim"/> (the floor itself, and the uniform-γ fence's
/// origin: the per-site law is −Re λ = 2Σ_l γ_l⟨Δ_l⟩, and it is the collapse of that sum to a count that
/// pinning buys) and <see cref="JointPopcountSectors"/> (the (N+1)² grading the 4N count is taken in;
/// without it there is no block to be pinned).</para>
///
/// <para>Anchors: <c>docs/ANALYTICAL_FORMULAS.md</c> F153 + <c>docs/proofs/PROOF_CODIM1_BY_ADDITIVITY.md</c>
/// §6 + <c>docs/proofs/PROOF_ABSORPTION_THEOREM.md</c> §4.3 (which names the saturation this criterion
/// supplies as what its own argument does not force) + <c>simulations/offdiag_sector_floor.py</c> (the
/// spectral gate) + <c>compute/RCPsiSquared.Diagnostics/Foundation/PinnedBlockFloorWitness.cs</c>
/// (PinnedBlockFloorWitness, inspect --root pinned).</para></summary>
public sealed class PinnedBlockFloorClaim : Claim
{
    /// <summary>Typed parent: the Absorption Theorem, whose floor this criterion says when a block attains
    /// throughout.</summary>
    public AbsorptionTheoremClaim Absorption { get; }

    /// <summary>Typed parent: the (N+1)² joint-popcount grading the count 4N is taken in.</summary>
    public JointPopcountSectors Sectors { get; }

    /// <summary>The criterion itself: one index held at an endpoint. Symmetric in (p, q), so the repo's
    /// two index-order conventions (ket-first here, bra-first elsewhere) cannot disagree about it.</summary>
    public static bool IsPinned(int n, int p, int q)
    {
        RequireBlock(n, p, q);
        return Math.Min(p, q) == 0 || Math.Max(p, q) == n;
    }

    /// <summary>How many of the (N+1)² blocks are pinned: 4(N+1) − 4 = 4N, two pinned rows and two pinned
    /// columns with the four one-cell blocks counted twice.</summary>
    public static int PinnedBlockCount(int n) =>
        n >= 1 ? 4 * n
            : throw new ArgumentOutOfRangeException(nameof(n), n, "N ≥ 1 (at N = 0 there is one block and no grading)");

    /// <summary>The block's disagreement window: n_diff = popcount(a ⊕ b) runs over
    /// [|p − q|, min(p + q, 2N − p − q)] in steps of 2. The minimum is forced by the popcount difference,
    /// the maximum by disjoint supports (p + q) or by the complement bound (2N − p − q). This is
    /// PROOF_CODIM1_BY_ADDITIVITY §6's window-combinatorics shell, and the criterion is its zero-width
    /// case: |p − q| = min(p + q, 2N − p − q) exactly when min(p, q) = 0 or max(p, q) = N.
    ///
    /// <para>The same closed form is stated and gated one folder over, by
    /// <c>WindowShellLemmaTests.NDiffRange_IsPopcountCombinatorics_AllBlocks_N5toN8</c>. This is the third
    /// in-repo statement of it (§6, that test, here); the cross-reference is so the three cannot drift.</para></summary>
    public static (int Min, int Max) DisagreementWindow(int n, int p, int q)
    {
        RequireBlock(n, p, q);
        return (Math.Abs(p - q), Math.Min(p + q, 2 * n - p - q));
    }

    /// <summary>The block's floor at uniform γ: Re λ ≥ −2γ·|p − q|, attained throughout exactly on the
    /// pinned blocks.</summary>
    public static double Floor(int n, int p, int q, double gamma)
    {
        RequireBlock(n, p, q);
        return -2.0 * gamma * Math.Abs(p - q);
    }

    private static void RequireBlock(int n, int p, int q)
    {
        if (n < 1) throw new ArgumentOutOfRangeException(nameof(n), n, "N ≥ 1");
        if (p < 0 || p > n) throw new ArgumentOutOfRangeException(nameof(p), p, $"p ∈ [0, {n}]");
        if (q < 0 || q > n) throw new ArgumentOutOfRangeException(nameof(q), q, $"q ∈ [0, {n}]");
    }

    public string ConverseByTrace =>
        "Why the criterion is an ⟺ and not an ⟸-by-mechanism. Suppose every eigenvalue of the (p,q) block " +
        "has Re λ = −2γ·n_min, the floor. The trace is basis-free: Re Tr L = Σ_λ Re λ = −2γ·n_min·dim. " +
        "Read the same trace off the cells: the diagonal entry of the cell |a⟩⟨b| has real part " +
        "−2γ·n_diff(a,b) (the hops are off-diagonal and the ZZ/field terms are purely imaginary), so " +
        "Re Tr L = −2γ·Σ_cells n_diff. Equating, Σ_cells (n_diff − n_min) = 0, and every term is " +
        "non-negative because n_min is the minimum over the block. Hence n_diff ≡ n_min on every cell. " +
        "By the window combinatorics that is min(p,q) = 0 or max(p,q) = N. No eigensolver, no tolerance, " +
        "and the argument needs uniform γ exactly where the criterion does: under a profile the cell's " +
        "real part is a SUM over disagreeing sites and n_diff is no longer what the trace counts.";

    public string TheThreeWords =>
        "FLOOR = the Absorption Theorem's bound, obeyed by every block. It is a floor on the DECAY RATE and " +
        "therefore a CEILING on Re λ: §6 confines Re spec to [−2·n_max, −2·|p−q|], so a pinned block sits at " +
        "the top of its own window and every other block runs downward from there. FLAT = attaining it throughout, " +
        "which is what three older documents assert of the two end blocks alone. PINNED = the cause, one " +
        "index at 0 or N. Also: the object is a BLOCK here and in the typed layer, a SECTOR in the Python " +
        "gate, while 'sector' elsewhere in the registry means the Pauli w = k sector and, in F151, an " +
        "R-parity half. Same word, three grids.";

    public string WhatItSubsumes =>
        "Three tracked documents carry 'only the two end sectors (0,1) and (N−1,N) are flat'. That is this " +
        "criterion read on one rung: among the ADJACENT blocks, |p − q| = 1, the criterion selects exactly " +
        "those two and their transposes. Read as a statement about all off-diagonal blocks it is FALSE: " +
        "(0,2), (0,5) and every block touching the vacuum or the full state are flat too, each at its own " +
        "floor −2γ|p − q|. Carry the |p − q| = 1 qualifier or carry the criterion, never the sentence alone.";

    public PinnedBlockFloorClaim(AbsorptionTheoremClaim absorption, JointPopcountSectors sectors)
        : base("F153 pinning criterion: a joint-popcount block (p,q) sits ENTIRELY on its Absorption floor " +
               "Re λ = −2γ|p−q| exactly when min(p,q) = 0 or max(p,q) = N, which is 4(N+1) − 4 = 4N of the " +
               "(N+1)² blocks. Mechanism: pinning one index leaves a single ket or a single bra, so every cell " +
               "carries the same disagreement count and at UNIFORM γ and REAL coupling the block reads " +
               "L = c·Id + i·S, S = the hops PLUS the ZZ diagonal (which contributes a configuration-dependent " +
               "IMAGINARY diagonal and not a real one), whose eigenvalues are c + i·μ, μ real. S is HERMITIAN in " +
               "general and REAL SYMMETRIC for a real H such as the chain: give the hopping a Peierls phase and " +
               "the criterion survives while the symmetry does not, so the general form is 'the Hermitian part " +
               "of L is scalar'. The ⟸ is a trace identity, " +
               "not a mechanism sketch: Re Tr L = Σ_cells(−2γ·n_diff) = Σ_λ Re λ forces n_diff ≡ n_min, which by " +
               "the §6 window [|p−q|, min(p+q, 2N−p−q)] is the criterion. FORWARD DIRECTION ALREADY OWNED by " +
               "PROOF_CODIM1_BY_ADDITIVITY §6 (Edge lemma = the |p−q|=1 case; the boundary-block sentence = the " +
               "general one; the window shell = the combinatorial half, gated by WindowShellLemmaTests), which " +
               "the registry entry was minted without citing. SCOPE: uniform γ REQUIRED (under a profile the " +
               "block leaves its FLOOR, which is not the same as spreading: at Δ=1, J=1, γ = linspace(0.1,1,N) " +
               "the (0,1) spread is 0.4903 at N=4 and 0.9477 at N=5 in the PAULI book, but at Δ=0 the same " +
               "profile gives MACHINE ZERO at J=1, a coalescence needing the R₉₀ locus the ramp happens to lie on, " +
               "and at J=1 FOUR of the sixteen N=4 " +
               "pinned blocks stay flat at the WRONG constant, (0,2) at −2·|p−q|·γ̄ = −2.200 against a floor of " +
               "−1.000, while the four one-cell blocks are flat AT their floor; " +
               "WHICH blocks flatten is not a |p−q| parity but an antilinear involution about that constant, " +
               "from the free index being N/2 or from γ anti-palindromic on the R₉₀ locus, F91), " +
               "the four ONE-CELL blocks (0,0),(0,N),(N,0),(N,N) EXEMPT from the fence for a dimensional reason " +
               "(not to be confused with F140's (1,1)-type CORNER blocks, a different object on the same locus), " +
               "real coupling REQUIRED, number-conserving H required, " +
               "criterion derived graph-blind and GATED on the chain only. It LOCATES F1's rate-reading blindness: " +
               "on these 4N blocks a rate reading loses nothing, on the other (N−1)² it does.",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F153 + " +
               "docs/proofs/PROOF_CODIM1_BY_ADDITIVITY.md §6 (Edge lemma, boundary blocks, window shell) + " +
               "docs/proofs/PROOF_ABSORPTION_THEOREM.md + " +
               "simulations/offdiag_sector_floor.py + " +
               "compute/RCPsiSquared.Diagnostics/Foundation/PinnedBlockFloorWitness.cs (PinnedBlockFloorWitness, inspect --root pinned)")
    {
        Absorption = absorption ?? throw new ArgumentNullException(nameof(absorption));
        Sectors = sectors ?? throw new ArgumentNullException(nameof(sectors));
    }

    public override string DisplayName =>
        "F153: the pinning criterion — 4N of the (N+1)² blocks sit entirely on their floor";

    public override string Summary =>
        "a joint-popcount block sits entirely on Re λ = −2γ|p−q| iff min(p,q) = 0 or max(p,q) = N, that is 4N " +
        "of (N+1)²; pinned ⇒ L = c·Id + i·S, S real symmetric INCLUDING the ZZ diagonal; the converse by trace " +
        "(Σ_cells n_diff = n_min·dim forces n_diff constant); uniform γ REQUIRED, number-conserving H required; " +
        $"locates F1's rate-reading blindness; live at inspect --root pinned ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("the criterion and the count",
                summary: $"min(p,q) = 0 or max(p,q) = N; 4(N+1) − 4 = 4N blocks, e.g. {PinnedBlockCount(4)} at N=4, " +
                         $"{PinnedBlockCount(5)} at N=5, {PinnedBlockCount(6)} at N=6. Symmetric in (p,q), so the " +
                         "repo's two index-order conventions cannot disagree about which blocks these are.");
            yield return new InspectableNode("the window, and the criterion as its zero-width case",
                summary: "n_diff ∈ [|p−q|, min(p+q, 2N−p−q)] in steps of 2 (PROOF_CODIM1 §6). Zero width exactly " +
                         "when min(p,q) = 0 or max(p,q) = N: |p−q| < p+q iff min(p,q) > 0, and " +
                         "|p−q| < 2N−p−q iff max(p,q) < N.");
            yield return new InspectableNode("the converse, by trace (this is the ⟺'s missing half)", summary: ConverseByTrace);
            yield return new InspectableNode("what it subsumes (and the sentence that is false unqualified)", summary: WhatItSubsumes);
            yield return new InspectableNode("three words for three things", summary: TheThreeWords);
            yield return new InspectableNode("scope: uniform γ is load-bearing",
                summary: "pinning fixes the disagreement COUNT, not the disagreeing SITES; under a profile the cell " +
                         "real part is −2Σ over disagreeing sites and varies, so the block leaves its floor " +
                         "(0.490 at N=4, 0.948 at N=5, γ = [0.1, 0.4, 0.7, 1.0], J = 1). The live witness bakes " +
                         "γ = 1 in and therefore CANNOT see this failure mode; the Python gate carries it.");
            yield return new InspectableNode("the live sweep at N=4, J=1, Δ=1 (Heisenberg)",
                summary: $"{PinnedBlockFloorWitness.WholeCount(4, 1.0, 1.0)} of 25 blocks read pinned entry-wise, " +
                         $"expected 4N = {PinnedBlockCount(4)}; residuals compared to 0.0 exactly, no eigensolver");
            yield return Absorption;   // typed parent edge
            yield return Sectors;      // typed parent edge
        }
    }

    /// <summary>Shared instance with its parents built from the Pi2-Foundation root, mirroring
    /// <see cref="VacuumBlockReductionClaim.Shared"/>: lets metadata be read without the full registry.</summary>
    public static PinnedBlockFloorClaim Shared { get; } =
        new PinnedBlockFloorClaim(new AbsorptionTheoremClaim(new Pi2DyadicLadderClaim()), new JointPopcountSectors());
}
