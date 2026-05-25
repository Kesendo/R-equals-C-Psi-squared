using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Pauli;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F108 Part 1 (Tier1Derived): every Π²_Z-even Hamiltonian H built from the
/// five Π²-even 2-site bilinears {XX, YY, YZ, ZY, ZZ} with arbitrary real bond
/// coefficients admits an EXACT operator-level palindrome under Z-dephasing:
///
/// <para>  Π_5bilinear · L · Π_5bilinear⁻¹ = −L − 2σ·I exactly, where σ = Σ_l γ_l.</para>
///
/// <para>Consequence: no Π²-even Pauli pair can be F87-hard. The empirical fact that
/// 5346+ Π²-even pairs across F103/F105/F106 anchors were observed soft (zero hard)
/// is now closed-form, not just empirical.</para>
///
/// <para>The proof switches from the canonical Heisenberg Π (where the palindrome
/// residual M = Π·L·Π⁻¹ + L + 2σ·I is generically nonzero for Π²-even non-truly H)
/// to <see cref="Pi5BilinearOperator"/>, a phase variant of the same P1 permutation
/// I↔X, Y↔Z with the X→I and Z→Y arrows sign-flipped. The phase-flipped per-site
/// map anti-commutes with the commutator superoperator [B, ·] of every Π²-even
/// 2-body bilinear B ∈ {XX, YY, YZ, ZY, ZZ}, while the per-site dissipator identity
/// M·D[Z]·M⁻¹ = −D[Z] − 2γ·I holds for the Z-dephasing Lindblad term. Combining
/// the two yields the operator-level palindrome with zero residual.</para>
///
/// <para>Derivation chain (see PROOF_F108_PART1):</para>
/// <list type="number">
///   <item>{Q, [B, ·]} = 0 for every Π²-even 2-body bilinear B (verified bit-exact
///         at the 2-qubit level; 4 Π²-odd bilinears produce residual = 8.00, clean
///         separation).</item>
///   <item>Q · L_H · Q⁻¹ = −L_H for any linear combination of Π²-even bilinears.</item>
///   <item>Q · D[Z_l] · Q⁻¹ = −D[Z_l] − 2γ_l·I per site, via diagonal permutation in
///         the Pauli basis: D[Z]_pauli = γ·diag(0, −2, −2, 0) on {I, X, Y, Z}; M's
///         (I↔X, Y↔Z) swap permutes the diagonal to γ·diag(−2, 0, 0, −2) =
///         −D[Z]_pauli − 2γ·I_4 (phase factors cancel pairwise on each 2-cycle).</item>
///   <item>Combining 2+3: Q · L · Q⁻¹ = −L − 2σ·I exactly, σ = Σγ.</item>
///   <item>spec(L) = spec(Q·L·Q⁻¹) = {−λ − 2σ : λ ∈ spec(L)}: palindromic around −σ.
///         Hence no Π²-even pair can be F87-hard.</item>
/// </list>
///
/// <para>Empirical confirmation: bit-exact residual = 0 across all 9 pure-Π²-even
/// non-truly pairs (single-bilinear YZ/ZY plus 7 two-term combinations) at N=3, 4, 5,
/// plus 15 random non-uniform-J instances + 9 asymmetric J_YZ ≠ J_ZY instances. Pure
/// dissipator D[Z]^⊗N also bit-exact at N=3, 4, 5. Reproduction:
/// <c>simulations/_f108_part1_pi_family_scan.py</c>.</para>
///
/// <para>Implements <see cref="IZ2AxisClaim"/> with <see cref="Z2Axis.BitB"/>; the
/// BitA twin (F108 Part 2, X-dephasing + Π²_X-even Hamiltonians via the analogous
/// P4-family phase-variant operator) requires its own derivation of the per-site
/// permutation + phases + bilinear set + D[X] identity; classified as
/// <see cref="BitATwinClassification.NeedsDerivation"/> (the bit_a-axis statement
/// involves different algebra than the bit_b-axis one). Proof:
/// <c>docs/proofs/PROOF_F108_PART1_PI2_EVEN_ALWAYS_PALINDROMIC.md</c>.</para></summary>
public sealed class F108Part1Pi2EvenAlwaysPalindromic : Claim, IZ2AxisClaim
{
    /// <summary>BitB axis (Π²_Z = X⊗N, bit_b parity = #Y + #Z mod 2). F108 Part 1
    /// constrains the Π²-even (bit_b = 0) sub-cell to always admit a palindrome
    /// operator.</summary>
    public Z2Axis Z2Axis => Z2Axis.BitB;

    /// <summary>The typed BitA-twin (X-dephasing + Π²_X-even Hamiltonians) is not
    /// yet implemented as a separate Claim; see <see cref="BitATwinStatus"/>.</summary>
    public Claim? BitATwin => null;

    /// <summary>The BitA-twin (F108 Part 2, F108 Part 1 under X-dephasing) is NOT a
    /// mechanical letter-swap mirror. The lift requires four pieces of new structural
    /// content: (a) a different per-site Π operator (analogous Π_5bilinear-X-variant,
    /// likely in the P4 family with its own phase choices); (b) a different Π²_X-even
    /// bilinear set (the bit_a parity = 0 pairs at the bilinear level, i.e.,
    /// {XX, XY, YX, YY, ZZ}); (c) a re-derived per-bond anti-commutation algebra
    /// {Q, [B, ·]} = 0 against that new bilinear set; (d) a new per-site D[X]
    /// dissipator identity Q·D[X]·Q⁻¹ = −D[X] − 2γ·I. Each piece needs verification.
    /// Classified <see cref="BitATwinClassification.NeedsDerivation"/> rather than
    /// TrivialNotYetTyped: the bit_a-axis statement involves different algebra than
    /// the bit_b-axis one even if the structural pattern transfers.</summary>
    public BitATwinClassification BitATwinStatus =>
        BitATwinClassification.NeedsDerivation;

    /// <summary>The five Π²_Z-even 2-site bilinears that Π_5bilinear handles
    /// with exact palindrome residual. {XX, YY, YZ, ZY, ZZ} = the bit_b=0 subset
    /// of all 16 letter pairs; equivalently the (BitA, BitB) Klein cells (0,0)
    /// {XX, YY, ZZ truly} and (1,0) {YZ, ZY non-truly}.</summary>
    public IReadOnlyList<(PauliLetter, PauliLetter)> Pi2EvenBilinears { get; } =
        new (PauliLetter, PauliLetter)[]
        {
            (PauliLetter.X, PauliLetter.X),
            (PauliLetter.Y, PauliLetter.Y),
            (PauliLetter.Y, PauliLetter.Z),
            (PauliLetter.Z, PauliLetter.Y),
            (PauliLetter.Z, PauliLetter.Z),
        };

    /// <summary>Returns true iff the (letter1, letter2) ordered pair belongs to the
    /// Π²_Z-even bilinear set {XX, YY, YZ, ZY, ZZ} that Π_5bilinear handles with
    /// exact palindrome residual. I-containing pairs (II, IX, XI, etc.) are
    /// considered trivially Π²-even but carry no 2-body Hamiltonian and are
    /// excluded.</summary>
    public static bool IsPi2EvenBilinear(PauliLetter letter1, PauliLetter letter2)
    {
        if (letter1 == PauliLetter.I || letter2 == PauliLetter.I) return false;
        // bit_b for each letter: I=0, X=0, Y=1, Z=1. bit_b of pair = sum mod 2.
        int bitB = letter1.BitB() + letter2.BitB();
        return (bitB & 1) == 0;
    }

    /// <summary>Returns true iff the given Pauli term is a 2-body bilinear whose
    /// two non-I letters form a Π²-even pair from
    /// <see cref="Pi2EvenBilinears"/>. F108 Part 1's exact-palindrome guarantee
    /// holds for any H assembled as a linear combination of such terms.</summary>
    public static bool IsPi2EvenBilinearTerm(PauliTerm term)
    {
        if (term is null) throw new ArgumentNullException(nameof(term));
        if (term.KBody != 2) return false;
        var nonI = new List<PauliLetter>();
        foreach (var L in term.Letters)
            if (L != PauliLetter.I) nonI.Add(L);
        if (nonI.Count != 2) return false;
        return IsPi2EvenBilinear(nonI[0], nonI[1]);
    }

    /// <summary>Returns true iff every term in the given Hamiltonian's term list is
    /// a Π²-even 2-body bilinear (per <see cref="IsPi2EvenBilinearTerm"/>).
    /// F108 Part 1's exact-palindrome guarantee holds iff this returns true.</summary>
    public static bool IsPi2EvenBilinearHamiltonian(IEnumerable<PauliTerm> terms)
    {
        if (terms is null) throw new ArgumentNullException(nameof(terms));
        foreach (var t in terms)
            if (!IsPi2EvenBilinearTerm(t)) return false;
        return true;
    }

    /// <summary>The theorem statement in one line.</summary>
    public string Theorem =>
        "For H = Σ_b α_b · B_b with B_b ∈ {XX, YY, YZ, ZY, ZZ} bilinears and α_b ∈ ℝ, " +
        "and Z-dephasing on every site: Π_5bilinear · L · Π_5bilinear⁻¹ = −L − 2σ·I exactly.";

    /// <summary>The F87 corollary, scoped to Z-dephasing (the F108 Part 1 proof's
    /// scope). The X- and Y-dephasing analogs are tracked on the BitATwin slot
    /// (see <see cref="BitATwinStatus"/>) and are not covered by this Claim.</summary>
    public string F87Corollary =>
        "Under Z-dephasing: no Π²_Z-even Pauli pair (truly or non-truly) is F87-hard; every such pair has palindromic spec(L).";

    public F108Part1Pi2EvenAlwaysPalindromic()
        : base("F108 Part 1: Π²-even H + Z-dephasing always admits exact operator-level palindrome via Π_5bilinear",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F108 + " +
               "docs/proofs/PROOF_F108_PART1_PI2_EVEN_ALWAYS_PALINDROMIC.md + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi5BilinearOperator.cs + " +
               "experiments/NON_HEISENBERG_PALINDROME.md")
    {
    }

    public override string DisplayName =>
        "F108 Part 1 Π²-even always palindromic via Π_5bilinear (closed-form)";

    public override string Summary =>
        $"{Theorem} {F87Corollary} ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("Theorem", summary: Theorem);
            yield return new InspectableNode("F87 corollary", summary: F87Corollary);
            yield return new InspectableNode("Π_5bilinear per-site map",
                summary: "I → +X, X → −I, Y → +iZ, Z → −iY (P1 permutation with X→I and Z→Y phase flips); M² = diag(−1, −1, +1, +1) on {I, X, Y, Z}");
            yield return new InspectableNode("Π²-even bilinears handled",
                summary: "{XX, YY, YZ, ZY, ZZ} = 5 of 16 2-site Pauli pairs; covers Klein (0,0) truly + (1,0) non-truly");
            yield return new InspectableNode("Empirical verification",
                summary: "9 pure-Π²-even non-truly pairs × N=3,4,5 = 27 instances, residual = 0; " +
                         "15 random non-uniform-J + 9 asymmetric J_YZ≠J_ZY instances, residual = 0; " +
                         "pure D[Z]^⊗N dissipator N=3,4,5, residual = 0")
            ;
            yield return new InspectableNode("Closes F109 dependency",
                summary: "F109 (mother sector Klein (0,0) soft y_par=1 purity) was Tier1Derived modulo F108 Part 1; " +
                         "F108 Part 1 closure promotes F109 to fully unconditional Tier1Derived");
            yield return new InspectableNode("Open siblings",
                summary: "F108 Part 2 (BitA twin via X-deph + P4-family phase variant): NeedsDerivation. " +
                         "F108 Part 3 (Y-deph analog, no covering Claim yet). " +
                         "F110 (hard cells y_par-pure with Y-inversion): higher difficulty, per-dephase-letter algebra");
        }
    }
}
