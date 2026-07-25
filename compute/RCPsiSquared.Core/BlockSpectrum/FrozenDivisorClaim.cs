using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.BlockSpectrum;

/// <summary>F140, typed (Tier 1 derived): <b>the R90 frozen divisor</b>. On the anti-palindromic
/// watching locus of F91 (every reflection pair of site rates carrying the same total,
/// γ_l + γ_{R(l)} = 2γ̄ for every l), the single-excitation corner block of the Z-dephasing
/// Liouvillian, on either chain of the proof's Section 1 (Heisenberg or XY), carries
///
/// <code>
///     λ = −4γ̄   with multiplicity ≥ ⌊N/2⌋,   for EVERY coupling J
/// </code>
///
/// one frozen mode per balanced pair, and exactly ⌊N/2⌋ for all but finitely many J. The value
/// hears only the MEAN watching: not the coupling, not the individual rates, nothing of the
/// Hamiltonian.
///
/// <para><b>No symmetry is behind it.</b> The eigenvectors move with J, the block spectrum is not
/// palindromic about the root, and no invariant subspace carries the modes. What pins the value is
/// a room shortage of the cell mirror τQ : (a,b) ↦ (R(b), R(a)), which fixes exactly the 2⌊N/2⌋
/// anti-diagonal coherences (a, R(a)), one per site. An involution's even rooms outnumber its odd
/// rooms by exactly its fixed count, so an operator odd under τQ must send the bigger half into
/// the smaller and the surplus must freeze; the diagonal cells' mirror pairs (dim D₋ = ⌊N/2⌋) are
/// EVEN under τQ, not odd, and tax away half. The estimate is built from three subspace
/// dimensions, so no deformation that keeps the generator odd can move it. Same lesson as F139 on
/// the other axis: <b>a wall can be a divisor instead of a symmetry.</b>
///
/// <para><b>Where it sits.</b> Four joint-popcount blocks carry, the corners
/// p, q ∈ {1, N−1}: (1,1) and (N−1,N−1) at −4γ̄, (1,N−1) and (N−1,1) at 4γ̄ − 2σ, the root chosen
/// by the parity of how many one-sided gamma folds separate the block from (1,1), each fold
/// sending r ↦ −r − 2σ. (At N = 4 the two roots coincide, (4 − 2N)γ̄ = −4γ̄, and all four corners
/// carry the one root.) Away from the locus the structure disappears entirely: partial balance
/// yields nothing, so the modes do not fade in, they snap in.</para>
///
/// <para><b>The confinement is the ZZ term's, not the divisor's, and it is the quarticity rather
/// than the diagonal on h.</b> That the remaining
/// (N+1)² − 4 blocks carry nothing at either root is verified on the HEISENBERG chain (full
/// census, N = 4, 5, 6). Drop the ZZ diagonal and it fails: on the XY chain the two roots are
/// carried at the same multiplicity ⌊N/2⌋ by many more blocks, counting a block that carries
/// either: nine of twenty-five at N = 4, twenty-four of thirty-six at N = 5, twenty-one of
/// forty-nine at N = 6, of which nine, twelve and fifteen carry the unfolded root itself. Never
/// partially, and −4γ̄ appears only on blocks with p + q even (necessary, not sufficient). The bound and the
/// fold-parity split survive the change of chain; "only the corners" does not. Found in the typing
/// of this claim, 2026-07-25, and pinned by gate G2c. What spreads the root on XY is an exact
/// ladder, ρ ↦ Σ_l d†_l ρ d_l, which commutes with the whole Liouvillian for any γ profile and
/// carries the corner's frozen modes up the diagonal; the ZZ term is quartic in the fermions and
/// removes it, while an R-invariant DIAGONAL on h does not (that one empties only the off-diagonal
/// half of the band, through its seed). See experiments/XY_FROZEN_BAND.md.</para>
///
/// <para><b>The cofactor and the ladder.</b> det(εI − M̃) = ε^⌊N/2⌋·q(ε) with
/// q(0) = (−1)^N·(4γ̄)^⌈N/2⌉·det((X P_{O₊} X)|_{V₋}), one N(N−1)/2 determinant free of γ̄; its
/// nonvanishing IS tightness, and it gives semisimplicity. As J → 0 that determinant vanishes to
/// order 2⌊N²/4⌋ = 2 Σ_c d_c with d_c = N + 1 − 2c the site distance of the balanced pair
/// (c, R(c)): a frozen mode cannot move until the coupling has walked the excitation across its
/// own pair, so it departs at order J^{2d_c} and <b>the far pair is the most protected</b>. At the
/// uniform endpoint the cofactor collapses to J^{N(N−1)}·D_N with
/// D_N = (−1)^{N(N−1)/2}·disc(h)·M^{−⌊(N−1)/2⌋}, M the chain's boundary clock: M = N for
/// Heisenberg (DCT-II Neumann) and M = N+1 for XY (DST-I Dirichlet).</para>
///
/// <para><b>Typed parents.</b> <see cref="F71AntiPalindromicGammaSpectralInvariance"/>: its
/// anti-palindromic γ-locus is exactly where this divisor lives, and the two readings meet, since
/// the whole R90 orbit shares one diagonal-block spectrum with the uniform profile.
/// <see cref="JointPopcountSectors"/>: the corner block is a block at all only because L is
/// exactly block-diagonal in the joint (popcount_row, popcount_col) label, which is also what lets
/// the divisor walk past the N = 8 spectral wall (the corner is N × N, not 4^N).</para>
///
/// <para><b>Scope, and what stays open.</b> The bound and tightness are both theorems: the bound
/// by the τQ pencil argument of Section 3, tightness by Section 7's closed form for the cofactor's
/// leading coefficient, which makes it a nonzero polynomial at every N on the two chains. What
/// tightness does NOT say is which couplings the finitely many exceptions are, and there the story
/// is unfinished: at the real exceptional couplings the root goes DEFECTIVE (one 2×2 Jordan block,
/// exact at N = 3, 4) while its kernel dimension does not move, so the criterion cannot tell that
/// failure from the harmless one at J = 0, where the multiplicity merely doubles and stays
/// semisimple; and how many exceptional couplings are real is not a function of N. Open in a
/// different direction: the upper half of the VALUATION law, that the J → 0 order is exactly
/// 2⌊N²/4⌋ rather than at least it, reduced by the pointed grading χ_x(a,b) = |a−x| + |b−R(x)| to
/// a single nonvanishing (each pair reaching its OWN outer anti-diagonal cell). That one concerns
/// the ladder, not the multiplicity.</para>
///
/// <para>Gate: <c>simulations/r90_frozen_divisor_gate.py</c> (212 checks, G0..G15). Live:
/// <c>inspect --root divisor</c> (<c>FrozenDivisorWitness</c>, the counts recomputed by exact
/// GF(p) ranks at inspect time). Adopted as a MirrorWorld object: run mode
/// <c>divisor N</c>.</para></summary>
public sealed class FrozenDivisorClaim : Claim
{
    // Parent-edge marker: the locus this divisor lives on IS F91's anti-palindromic γ-locus.
    public F71AntiPalindromicGammaSpectralInvariance AntiPalindromicLocus { get; }

    // Parent-edge marker: "the corner block" presupposes the joint-popcount block decomposition.
    public JointPopcountSectors Sectors { get; }

    public FrozenDivisorClaim(
        F71AntiPalindromicGammaSpectralInvariance antiPalindromicLocus,
        JointPopcountSectors sectors)
        : base("The R90 frozen divisor (F140): on F91's anti-palindromic watching locus " +
               "(gamma_l + gamma_{R(l)} = 2*gbar for every l) the single-excitation corner block of the " +
               "Z-dephasing Liouvillian, on either the Heisenberg or the XY chain, carries " +
               "lambda = -4*gbar with multiplicity at least " +
               "floor(N/2) for EVERY coupling J, one frozen mode per balanced pair, and exactly " +
               "floor(N/2) for all but finitely many J; no symmetry is behind it, the value is forced " +
               "by a room shortage of the cell mirror tauQ (2*floor(N/2) fixed anti-diagonal cells, " +
               "floor(N/2) taxed away by the even diagonal pairs), so the modes hear only the mean " +
               "watching and nothing of the Hamiltonian; the four corner blocks p,q in {1, N-1} " +
               "carry, the root picked by gamma-fold parity, and on the Heisenberg chain nothing " +
               "else does (the confinement is the ZZ term's quarticity, not a diagonal on h: on " +
               "the XY chain the same root spreads at the same multiplicity along an exact " +
               "ladder), while partial balance yields nothing on either " +
               "chain; as J -> 0 the cofactor vanishes to order 2*floor(N^2/4), each pair departing " +
               "at J^(2 d_c) with d_c its site distance, so the outermost pair holds longest",
               Tier.Tier1Derived,
               "docs/proofs/PROOF_R90_FROZEN_DIVISOR.md")
    {
        AntiPalindromicLocus = antiPalindromicLocus ?? throw new ArgumentNullException(nameof(antiPalindromicLocus));
        Sectors = sectors ?? throw new ArgumentNullException(nameof(sectors));
    }

    /// <summary>The frozen multiplicity ⌊N/2⌋: one mode per balanced reflection pair.</summary>
    public static int FrozenMultiplicity(int n) =>
        n >= 1 ? n / 2 : throw new ArgumentOutOfRangeException(nameof(n), n, "N must be ≥ 1.");

    /// <summary>τQ's fixed cells among the off-diagonal coherences, 2⌊N/2⌋: the surplus of even
    /// rooms over odd ones, before the diagonal pairs tax half of it away.</summary>
    public static int MirrorFixedCells(int n) => 2 * FrozenMultiplicity(n);

    /// <summary>The frozen value −4γ̄ = −4σ/N.</summary>
    public static double FrozenRoot(double gammaBar) => -4.0 * gammaBar;

    /// <summary>Its gamma-fold partner 4γ̄ − 2σ = (4 − 2N)γ̄, the root the odd-fold corners carry.
    /// </summary>
    public static double FoldedRoot(int n, double gammaBar) => (4.0 - 2.0 * n) * gammaBar;

    /// <summary>The site distance d_c = N + 1 − 2c of the balanced pair (c, R(c)), c = 1..⌊N/2⌋
    /// one-based; its mode departs at order J^{2 d_c}.</summary>
    public static int PairDistance(int n, int c) =>
        c >= 1 && c <= FrozenMultiplicity(n)
            ? n + 1 - 2 * c
            : throw new ArgumentOutOfRangeException(nameof(c), c, $"pair index must lie in 1..{FrozenMultiplicity(n)}");

    /// <summary>The total J → 0 valuation of the cofactor determinant, 2⌊N²/4⌋ = 2 Σ_c d_c.
    /// </summary>
    public static int TotalValuation(int n) => 2 * (n * n / 4);

    public override string DisplayName =>
        "The R90 frozen divisor (F140): a watching locus that pins the decay rate −4γ̄ at every coupling, " +
        "at least ⌊N/2⌋ modes deep";

    public override string Summary =>
        "on F91's anti-palindromic γ-locus the corner block holds λ = −4γ̄ at least ⌊N/2⌋ times for " +
        "every J, and exactly ⌊N/2⌋ for all but finitely many J (at J = 0 it doubles), one " +
        "frozen mode per balanced pair; no symmetry behind it, only a room shortage of the cell mirror τQ " +
        "(2⌊N/2⌋ fixed anti-diagonal cells minus the ⌊N/2⌋ even diagonal pairs), so the modes hear the mean " +
        "watching and nothing else; the four corners carry, the root picked by gamma-fold parity, and on the " +
        "Heisenberg chain nothing else does; partial balance yields nothing; the ladder J^{2d_c} makes the " +
        "outermost pair the last to let go " +
        $"({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("the divisor bound (the proved half)",
                summary: "λ = −4γ̄ with multiplicity ≥ ⌊N/2⌋ at EVERY coupling, proved twice: once by the " +
                         "τQ pencil argument, once as an index (the mirror's own orbit count). Both are " +
                         "dimension estimates, so no deformation that keeps the generator odd can move the " +
                         "number, and the divisor walks past the N = 8 spectral wall because the corner " +
                         "block is N × N.");
            yield return new InspectableNode("the room shortage: why it is not a symmetry",
                summary: "τQ : (a,b) ↦ (R(b), R(a)) fixes exactly the 2⌊N/2⌋ anti-diagonal coherences " +
                         "(a, R(a)), and an involution's even rooms exceed its odd rooms by exactly its " +
                         "fixed count. An odd operator must overflow. The diagonal cells' mirror pairs are " +
                         "EVEN, not odd, and eat half: ⌊N/2⌋ = 2⌊N/2⌋ − ⌊N/2⌋. The eigenvectors still move " +
                         "with J and the block spectrum is not palindromic about the root; what is fixed is " +
                         "a count, not a subspace.");
            yield return new InspectableNode("where it sits: the four corners, and on which chain",
                summary: "p, q ∈ {1, N−1}: (1,1) and (N−1,N−1) at −4γ̄, (1,N−1) and (N−1,1) at 4γ̄ − 2σ, the " +
                         "root chosen by how many one-sided gamma folds separate the block from (1,1), each " +
                         "fold sending r ↦ −r − 2σ; at N = 4 the two roots coincide. That the remaining " +
                         "(N+1)² − 4 blocks carry nothing is a HEISENBERG statement (full census, " +
                         "N = 4, 5, 6): drop the ZZ diagonal and the same root spreads, at the same " +
                         "multiplicity, over many more blocks. The confinement belongs to the ZZ term as a " +
                         "QUARTIC term, not as the diagonal it puts on h: an R-invariant diagonal on h leaves " +
                         "the diagonal band full, and what ZZ removes is the ladder ρ ↦ Σ_l d†_l ρ d_l that " +
                         "carries the corner up it. The bound and the fold parity belong to the divisor.");
            yield return new InspectableNode("the ladder: distance buys immunity",
                summary: "as J → 0 the cofactor vanishes to order 2⌊N²/4⌋, and per pair the mode departs at " +
                         "J^{2 d_c}, d_c = N + 1 − 2c the site distance: the coupling must walk the " +
                         "excitation from one site of the pair to the other before the mode can move. The " +
                         "two ends of the chain are the last to let go. At J = 0 the root carries twice the " +
                         "multiplicity, all semisimple, and exactly ⌊N/2⌋ modes depart as the coupling " +
                         "turns on.");
            yield return new InspectableNode("typed parents",
                summary: $"F71AntiPalindromicGammaSpectralInvariance ({AntiPalindromicLocus.Tier.Label()}): " +
                         "the anti-palindromic γ-locus this divisor lives on, and the reading that an " +
                         "on-locus profile shares its diagonal-block spectrum with the uniform one. " +
                         $"JointPopcountSectors ({Sectors.Tier.Label()}): \"the corner block\" presupposes " +
                         "the exact (N+1)² block decomposition, which is also what puts the divisor past " +
                         "the spectral wall.");
            yield return new InspectableNode("what stays open",
                summary: "not tightness, which Section 7 turned into a theorem at every N, but WHICH " +
                         "couplings the finitely many exceptions are: at the real ones the root goes " +
                         "defective, one 2×2 Jordan block, its kernel dimension unmoved, so the criterion " +
                         "cannot tell that failure from the harmless one at J = 0; and how many are real " +
                         "is not a function of N. Open in another direction: the upper half of the " +
                         "VALUATION law (the J → 0 order exactly 2⌊N²/4⌋ rather than at least it), " +
                         "reduced by the pointed grading to a single nonvanishing. That one is about the " +
                         "ladder, not the multiplicity.");
            yield return new InspectableNode("live witness (inspect --root divisor)",
                summary: "FrozenDivisorWitness recomputes the counts at inspect time through the repo's own " +
                         "committed builders (PauliHamiltonian, Heisenberg by default and XY on request, " +
                         "plus PerBlockLiouvillianBuilder), by " +
                         "exact GF(p) ranks rather than an eigensolver: the departing modes sit at spacing " +
                         "J^{2d} from the root, so a floating-point count would miscount exactly where the " +
                         "theorem is sharpest. Beside each ⌊N/2⌋ read sits a zero (wrong root, off-locus at " +
                         "unchanged σ, a non-corner block) and the J = 0 doubling.");
        }
    }
}
