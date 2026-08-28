using RCPsiSquared.Core.F1;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F158, the palindrome as a count of the two ends: for
///
/// <code>
///     L(rho) = -i[H, rho] + sum_l gamma_l (A_l rho A_l - rho),
///     H Hermitian, every A_l Hermitian with A_l^2 = 1, every gamma_l > 0, sigma = sum_l gamma_l,
///
///     spec(L) is closed under lambda -> -lambda - 2*sigma   <=>   dim ker L == dim ker(L + 2*sigma)
/// </code>
///
/// <para>There is no operator to find and no subspace to sample: the criterion is two nullities
/// compared, each an exact rank. It replaces the existence predicate of
/// <c>experiments/THE_PAIRING_CONDITION.md</c> ("there is an invertible U with [U,H] = 0 and
/// U A_l U^-1 = -A_l"), which is the same statement read one notch less sharply, and it supplies
/// the NECESSITY direction that page could only measure over 140861 rows.</para>
///
/// <para><b>The chain, in four steps, each of which can fail alone.</b>
/// (1) ker L is the COMMUTANT {X : [H,X] = 0, A_l X A_l = +X} and ker(L + 2 sigma) is the
/// ANTICOMMUTANT {W : [H,W] = 0, A_l W A_l = -W}. One inclusion each way is pure algebra; the
/// other is the EQUALITY CASE of Cauchy-Schwarz in the Hilbert-Schmidt norm, since
/// &lt;W, -i[H,W]&gt; is purely imaginary and every Re&lt;W, A_l W A_l&gt; is bounded below by
/// -||W||^2, so a sum of them can reach -sigma||W||^2 only if every term attains its bound.
/// (2) BOTH eigenvalues are semisimple, so geometric multiplicity equals algebraic multiplicity at
/// both ends. This is NOT new here: PROOF_CODIM1_BY_ADDITIVITY's window-edge lemma owns it in a
/// stronger form, at every edge of a rate window rather than at the two extremes.
/// (3) The anticommutant holds an INVERTIBLE element exactly when the two dimensions agree
/// (forward, W = N.U; backward, by Wedderburn multiplicities and AM-GM, with the ungraded case
/// giving a strict inequality instead).
/// (4) Invertible U =&gt; palindrome is the pairing-condition page's one-sided calculation
/// L_U L L_U^-1 = -L^dagger - 2 sigma, collapsed onto the palindrome by hermiticity preservation;
/// palindrome =&gt; the two algebraic multiplicities agree =&gt; by (2) the two nullities agree.</para>
///
/// <para><b>What the fences are for, since all three are load bearing.</b>
/// <c>A_l^2 = 1</c>: F137 records that non-unitary (T1) jumps keep a palindrome about a DIFFERENT
/// centre, which F137 states exactly as trace(L)/dim, so a criterion phrased about -2 sigma is
/// answering a different question there. <c>gamma_l &gt; 0</c>: a zero rate drops that site's
/// condition from both spaces and MOVES the verdict, so the criterion is discontinuous on the
/// boundary of the positive orthant and "the palindrome does not depend on gamma" means inside
/// the open orthant only. EVEN d: an invertible U with U A U^-1 = -A makes A and -A similar, so
/// A needs balanced +1/-1 multiplicities; at odd d no such U exists at any H, and the palindrome
/// is impossible rather than merely absent.</para>
///
/// <para><b>What this class is NOT.</b> It decides the eigenvalue MULTISET, which is what the
/// committed GF(p) kernel tests, and not the Jordan structure away from the two ends;
/// <see cref="F1PalindromeIdentity"/> is the operator identity and sees more. And the char-poly
/// FORM of the multiset statement carries a sign: with D = d^2 the degree,
/// p(x) = p(-x - 2 sigma) is the multiset statement for even d and flips sign at odd d, where by
/// the paragraph above both sides are false anyway.</para>
///
/// <para>Proof <c>docs/proofs/PROOF_PALINDROME_TWO_END_COUNT.md</c>; gate
/// <c>simulations/f138_rank_criterion.py</c> (93 gates, 15415 rows scored in both directions with
/// FP = 0 and FN = 0, including F1's own canonical break under depolarizing, off-axis
/// n.sigma jumps, multi-site Pauli strings, and one float route that shares no code with the
/// rest). Live lab: <c>inspect --root twoend</c>.</para></summary>
public sealed class PalindromeTwoEndCountClaim : Claim
{
    /// <summary>The typed parent: F1 is the OPERATOR identity and this claim is the SPECTRAL one,
    /// and F1's own validity node already draws that distinction and scopes it to F138. The edge
    /// runs this way rather than the other because a row where the spectrum pairs says nothing
    /// about whether Pi conjugates L there, so this claim is the weaker statement of the two and
    /// F1 implies it wherever F1 holds.</summary>
    public F1PalindromeIdentity PalindromeIdentity { get; }

    /// <summary>The canonical Heisenberg chain under Z-dephasing on every site: both counts are
    /// N + 1, the magnetization-sector projectors at the near end and the XOR-sector modes at the
    /// far one. Both numbers are already in the repo, together, on one page
    /// (<c>experiments/DEGENERACY_PALINDROME.md</c> Result 2), which states the bijection between
    /// them and reads it as a consequence of Pi. The theorem inverts that reading: the bijection
    /// is not downstream of the palindrome, it IS the palindrome, and it decides cases where Pi is
    /// not available.</summary>
    public static int CanonicalChainCount(int n) =>
        n >= 1 ? n + 1 : throw new ArgumentOutOfRangeException(nameof(n), n, "N must be at least 1");

    /// <summary>Whether the palindrome is possible at all at this local dimension. False at odd d
    /// for every H and every jump set, by the balanced-multiplicity argument above; this is a
    /// corollary, not a measurement, and is gated at d = 3 and d = 5.</summary>
    public static bool PalindromePossibleAtDimension(int d) =>
        d >= 1 ? d % 2 == 0 : throw new ArgumentOutOfRangeException(nameof(d), d, "d must be at least 1");

    public PalindromeTwoEndCountClaim(F1PalindromeIdentity palindromeIdentity)
        : base("F158 the palindrome is a count of the two ends: for L(rho) = -i[H,rho] + sum_l gamma_l (A_l rho A_l " +
               "- rho) with H Hermitian, every A_l Hermitian and squaring to 1, and every gamma_l > 0, the spectrum " +
               "multiset is closed under lambda -> -lambda - 2*sigma IFF dim ker L = dim ker(L + 2*sigma). The near " +
               "kernel is the commutant of the algebra generated by H and the jumps and the far one is the same " +
               "space with the jump sign flipped, both by a Cauchy-Schwarz equality case; both eigenvalues are " +
               "semisimple, so the char-poly identity forces the two nullities to agree, and equal nullities force " +
               "an invertible U with [U,H] = 0 and U A_l U^-1 = -A_l, which reflects the spectrum. No operator to " +
               "exhibit and no subspace to sample: two ranks compared. Fences, all three load bearing: A_l^2 = 1 " +
               "(F137 recentres for T1 jumps), gamma_l > 0 (a zero rate drops a condition and moves the verdict), " +
               "and even d (at odd d no invertible U exists at any H, so the palindrome is impossible)",
               Tier.Tier1Derived,
               "docs/proofs/PROOF_PALINDROME_TWO_END_COUNT.md (primary: Lemmas 1-3 and the theorem) + " +
               "experiments/THE_PAIRING_CONDITION.md (the criterion this sharpens, and the sufficiency " +
               "calculation consumed unchanged) + " +
               "docs/ANALYTICAL_FORMULAS.md (F158; F1 for the operator identity, F138 for the law whose converse " +
               "was withdrawn 2026-08-03, F137 for the non-unitary fence, F4 for the near count) + " +
               "docs/proofs/PROOF_CODIM1_BY_ADDITIVITY.md (the window-edge lemma, which owns the semisimplicity " +
               "step in a stronger form and is cited rather than re-derived) + " +
               "docs/proofs/PROOF_F103_F87_Z2_CUBED_REFINEMENT.md (7.5 the forcing step in a narrower setting, " +
               "7.12 the three rows gated as an agreement test) + " +
               "experiments/DEGENERACY_PALINDROME.md (Result 2: both counts, on one page)")
    {
        PalindromeIdentity = palindromeIdentity ?? throw new ArgumentNullException(nameof(palindromeIdentity));
    }

    public override string DisplayName =>
        "F158: the palindrome is a count of the two ends (dim ker L = dim ker(L + 2 sigma))";

    public override string Summary =>
        "the spectrum pairs about -sigma exactly when the commutant and the anticommutant have the same " +
        $"dimension; on the canonical Heisenberg chain under Z-dephasing both are N + 1 ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("the two ends are one object with one sign flipped (Lemma 1)",
                summary: "ker L = {X : [H,X] = 0, A_l X A_l = +X} and ker(L + 2 sigma) = {W : [H,W] = 0, " +
                         "A_l W A_l = -W}. The easy inclusion each way is pure algebra and holds over any field, " +
                         "which is why the modular layer of the gate can check it at every prime. The reverse " +
                         "needs C: <W, -i[H,W]> is purely imaginary because the trace of a product of two " +
                         "Hermitian operators is real, so taking real parts leaves a sum of terms each bounded " +
                         "below by -||W||^2 with positive weights summing to sigma, which can reach -sigma||W||^2 " +
                         "only if every term attains its bound. Equality in Cauchy-Schwarz between vectors of " +
                         "equal norm then gives A_l W A_l = -W per site, and substituting back leaves [H,W] = 0. " +
                         "The kernel end is the same four lines with one sign changed.");

            yield return new InspectableNode("both ends are semisimple (Lemma 2), and this is NOT new",
                summary: "Re<X, L X> lies in [-2 sigma ||X||^2, 0], so L is HS-dissipative, L + 2 sigma is " +
                         "HS-accretive, and an accretive operator has no Jordan block at 0 (if M W = V with " +
                         "M V = 0 then Re<W + tV, M(W + tV)> = Re<W, M W> + t||V||^2 >= 0 for every real t, so " +
                         "V = 0). Geometric = algebraic at both ends, which is what turns the char-poly identity " +
                         "into a statement about two nullities. PROOF_CODIM1_BY_ADDITIVITY's window-edge lemma " +
                         "owns this in a stronger form, at EVERY edge of a rate window; the arc ledger carries a " +
                         "standing instruction to cite it rather than re-derive it, and the first draft of the " +
                         "proof file re-derived it a third time.");

            yield return new InspectableNode("invertibility IS the rank equality (Lemma 3)",
                summary: "Forward: an invertible U in the anticommutant gives W = N.U, so the dimensions agree. " +
                         "Backward, by representation theory of the algebra A generated by H and the jumps: if " +
                         "the sign flip alpha (H -> H, A_l -> -A_l) is well defined on A then the anticommutant is " +
                         "Hom_A(rho, rho o alpha), of dimension sum_k m_k m_tau(k) <= sum_k m_k^2 = dim commutant " +
                         "by AM-GM, with equality exactly when rho and rho o alpha are isomorphic, i.e. exactly " +
                         "when an invertible intertwiner exists. If alpha is NOT well defined, some c != 0 is " +
                         "both an even and an odd word, so cW = -cW = 0 kills every W on both sides, no element " +
                         "is invertible, and the inequality is STRICT. This is what turns the pairing-condition " +
                         "page's sampled predicate into a decided one: on every row that page reports as " +
                         "nonempty-but-singular, dim W < dim N, and the strict inequality PROVES no invertible " +
                         "element is hiding there.");

            yield return new InspectableNode("the kernel always dominates, and the far space's name",
                summary: "dim ker(L + 2 sigma) <= dim ker L for every H, every unitary-Hermitian jump set and " +
                         "every positive profile, with no palindrome in sight. The tempting name for the far " +
                         "space, 'the fastest-decaying modes', is wrong twice: when it is 0 the spectrum still " +
                         "has a leftmost point, and even when the value -2 sigma is attained the LINE " +
                         "Re = -2 sigma can carry oscillating modes the kernel at the POINT -2 sigma does not " +
                         "contain (one qubit, H = Z, A = Z: the kernel is 0 while two eigenvalues sit at " +
                         "-2 gamma +- 2i). The honest name is the non-oscillating modes at the left edge.");

            yield return new InspectableNode("the canonical chain, where the repo held both numbers",
                summary: $"Heisenberg chain, Z-dephasing on every site: both counts are N + 1 " +
                         $"(N = 3 gives {CanonicalChainCount(3)}, N = 5 gives {CanonicalChainCount(5)}), the " +
                         "magnetization-sector projectors at the near end and the XOR-sector modes at the far " +
                         "one. experiments/DEGENERACY_PALINDROME.md Result 2 carries BOTH counts, on one page, " +
                         "and states the bijection between them, reading it as a consequence of Pi. What the " +
                         "theorem adds is that their equality is not downstream of the palindrome but the same " +
                         "statement, and decides cases where Pi is not available. Note the two gamma books: that " +
                         "page writes the far end as -N gamma at gamma = 0.1 per site, this claim as -2 sigma in " +
                         "the canonical Lindblad book at gamma = 0.05; the same point, and the repo fences the " +
                         "trap in docs/GLOSSARY.md and docs/CAUGHT_ERRORS.md.");

            yield return new InspectableNode("gamma-blindness, and where it stops",
                summary: "Neither defining condition mentions gamma_l, so the criterion PREDICTS that the " +
                         "palindrome cannot depend on the rate profile at all, before any run. What that does not " +
                         "say: the spectrum depends on gamma throughout, the centre -sigma moves with it, and at " +
                         "the boundary gamma_l = 0 the VERDICT moves, because a zero rate drops a condition from " +
                         "both spaces. Gated: a ZZ bond with an X field on site 0 and X-dephasing on both sites " +
                         "is broken while both rates are on and palindromic the moment site 0 stops being watched.");

            yield return new InspectableNode("odd local dimension: impossible, not merely absent",
                summary: $"An invertible U with U A U^-1 = -A makes A and -A similar, so A needs balanced " +
                         $"+1/-1 eigenvalue multiplicities and d must be even. d = 2 possible: " +
                         $"{PalindromePossibleAtDimension(2)}; d = 3 possible: {PalindromePossibleAtDimension(3)}; " +
                         $"d = 4: {PalindromePossibleAtDimension(4)}. So at odd d the anticommutant is strictly " +
                         "smaller than the commutant at every H, and the theorem holds there with both sides " +
                         "false. Consistently, the char-poly FORM p(x) = p(-x - 2 sigma) is unsatisfiable at odd " +
                         "d, a monic polynomial of odd degree d^2 picking up a sign under the reflection.");

            yield return new InspectableNode("scope, and the fences that do not lift",
                summary: "Proved for any Hermitian H, any finite set of Hermitian jumps squaring to 1 (single " +
                         "letters, off-axis n.sigma at unit directions, multi-site Pauli strings and full " +
                         "depolarizing sites alike), any strictly positive profile, any topology, any N, any " +
                         "finite dimension. Gated at d = 2^N with N <= 5, plus d = 3 and 5 for the odd-d " +
                         "corollary. OUTSIDE: jumps with A^2 != 1 (F137 recentres), rates that are zero or " +
                         "negative, and the Jordan structure away from the two ends. What is decided is the " +
                         "MULTISET; F1PalindromeIdentity is the operator identity and sees the Jordan structure " +
                         "too, which is why a row where the spectrum pairs anyway falsifies a spectral converse " +
                         "and says nothing about whether Pi conjugates L there.");

            yield return new InspectableNode("live lab (the witness)",
                summary: "PalindromeTwoEndCountWitness recomputes both nullities at inspect time as exact GF(p) " +
                         "eliminations and, independently, decides the palindrome by the characteristic-polynomial " +
                         "identity over the same field, then compares the two verdicts: inspect --root twoend. " +
                         "Two independent computations meeting, in the house pattern. Gate: " +
                         "simulations/f138_rank_criterion.py (93 gates), whose companion " +
                         "simulations/f138_pairing_condition.py carries the 140861-row census the necessity " +
                         "direction no longer needs.");
        }
    }
}
