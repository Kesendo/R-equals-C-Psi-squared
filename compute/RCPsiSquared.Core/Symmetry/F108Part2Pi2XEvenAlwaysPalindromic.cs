using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Pauli;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F108 Part 2 (Tier1Derived, BitA twin of F108 Part 1): every Π²_X-even
/// Hamiltonian H built from the five Π²_X-even 2-site bilinears {ZZ, XX, XY, YX, YY}
/// with arbitrary real bond coefficients admits an EXACT operator-level palindrome
/// under X-dephasing:
///
/// <para>  Π_5bilinear · L · Π_5bilinear⁻¹ = −L − 2σ·I exactly (with the X-deph
/// variant of <see cref="Pi5BilinearOperator"/>), where σ = Σ_l γ_l.</para>
///
/// <para>Consequence: no Π²_X-even Pauli pair can be F87-hard under X-dephasing.
/// Closes the X-dephasing branch of F109's Step 5 (mother sector soft = y_par=1
/// purity), leaving only the Y-dephasing branch (F108 Part 3) open.</para>
///
/// <para>Mirrors F108 Part 1's mechanism with the X-dephasing-specific per-site
/// permutation and bilinear set. Per-site M:
/// I → +Z, Z → −I, X → −i·Y, Y → +i·X. Same I↔Z, X↔Y permutation as the canonical
/// X-dephasing Π (P1 family on bit_b axis); the back-arrows Z → I and Y → X carry
/// the analogous sign flips that distinguish the 5-bilinear variant from the
/// canonical X-dephasing Π. M² = diag(−1, +1, +1, −1) on {I, X, Y, Z}; M⁴ = I.</para>
///
/// <para>Derivation chain (see PROOF_F108_PART2):</para>
/// <list type="number">
///   <item>{Q, [B, ·]} = 0 for every Π²_X-even 2-body bilinear B ∈ {ZZ, XX, XY,
///         YX, YY} (verified bit-exact at the 2-qubit level; 4 Π²_X-odd
///         bilinears produce residual = 8.00, clean separation).</item>
///   <item>Q · L_H · Q⁻¹ = −L_H for any linear combination of Π²_X-even
///         bilinears.</item>
///   <item>Q · D[X_l] · Q⁻¹ = −D[X_l] − 2γ_l·I per site, via diagonal permutation
///         in the Pauli basis: D[X]_pauli = γ·diag(0, 0, −2, −2) on {I, X, Y, Z};
///         M's (I↔Z, X↔Y) swap permutes the diagonal entries to
///         γ·diag(−2, −2, 0, 0) = −D[X]_pauli − 2γ·I_4 (phase factors cancel
///         pairwise on each 2-cycle).</item>
///   <item>Combining 2+3: Q · L · Q⁻¹ = −L − 2σ·I exactly, σ = Σγ.</item>
///   <item>spec(L) = spec(Q·L·Q⁻¹) = {−λ − 2σ : λ ∈ spec(L)}: palindromic around
///         −σ. Hence no Π²_X-even pair is F87-hard under X-dephasing.</item>
/// </list>
///
/// <para>Empirical confirmation: bit-exact residual = 0 across all 9 pure-Π²_X-even
/// non-truly pairs (single-bilinear XY/YX plus 7 two-term combinations) at N=3, 4, 5,
/// plus random non-uniform-J instances. Pure dissipator D[X]^⊗N also bit-exact at
/// N=3, 4, 5. Reproduction:
/// <c>simulations/_f108_part2_x_dephasing_scan.py</c>.</para>
///
/// <para>Implements <see cref="IZ2AxisClaim"/> with <see cref="Z2Axis.BitA"/>; the
/// BitB twin (F108 Part 1) is <see cref="F108Part1Pi2EvenAlwaysPalindromic"/>. F108
/// Part 1 references this Claim as its BitATwin (closing the Stage-2a slot from
/// NeedsDerivation to Filled). Proof:
/// <c>docs/proofs/PROOF_F108_PART2_PI2X_EVEN_ALWAYS_PALINDROMIC.md</c>.</para></summary>
public sealed class F108Part2Pi2XEvenAlwaysPalindromic : Claim, IZ2AxisClaim
{
    /// <summary>BitA axis (Π²_X = Z⊗N, bit_a parity = #X + #Y mod 2). F108 Part 2
    /// constrains the Π²_X-even (bit_a = 0) sub-cell under X-dephasing to always
    /// admit a palindrome operator.</summary>
    public Z2Axis Z2Axis => Z2Axis.BitA;

    /// <summary>BitA-axis claims do not have BitATwin slots (the twin concept lives
    /// on BitB-axis claims pointing at BitA-axis siblings). Always null per the
    /// IZ2AxisClaim contract.</summary>
    public Claim? BitATwin => null;

    /// <summary>Explicit override returning <see cref="BitATwinClassification.NotApplicableForThisAxis"/>:
    /// the BitA-twin slot semantics does not apply to BitA-axis Claims (it is for
    /// BitB-axis Claims pointing at BitA-axis siblings). Override exposes the
    /// status on the concrete class for tests and the inspector.</summary>
    public BitATwinClassification BitATwinStatus => BitATwinClassification.NotApplicableForThisAxis;

    /// <summary>The five Π²_X-even 2-site bilinears that Π_5bilinear (X-deph
    /// variant) handles with exact palindrome residual. {ZZ, XX, XY, YX, YY} =
    /// the bit_a=0 subset of all 16 letter pairs at the bilinear level.</summary>
    public IReadOnlyList<(PauliLetter, PauliLetter)> Pi2XEvenBilinears { get; } =
        new (PauliLetter, PauliLetter)[]
        {
            (PauliLetter.Z, PauliLetter.Z),
            (PauliLetter.X, PauliLetter.X),
            (PauliLetter.X, PauliLetter.Y),
            (PauliLetter.Y, PauliLetter.X),
            (PauliLetter.Y, PauliLetter.Y),
        };

    /// <summary>Returns true iff the (letter1, letter2) ordered pair belongs to the
    /// Π²_X-even bilinear set {ZZ, XX, XY, YX, YY} that Π_5bilinear (X-deph
    /// variant) handles with exact palindrome residual. I-containing pairs (II,
    /// IZ, ZI, etc.) are considered trivially Π²_X-even but carry no 2-body
    /// Hamiltonian and are excluded.</summary>
    public static bool IsPi2XEvenBilinear(PauliLetter letter1, PauliLetter letter2)
    {
        if (letter1 == PauliLetter.I || letter2 == PauliLetter.I) return false;
        int bitA = letter1.BitA() + letter2.BitA();
        return (bitA & 1) == 0;
    }

    /// <summary>Returns true iff the given Pauli term is a 2-body bilinear whose
    /// two non-I letters form a Π²_X-even pair from
    /// <see cref="Pi2XEvenBilinears"/>. F108 Part 2's exact-palindrome guarantee
    /// holds for any H assembled as a linear combination of such terms under
    /// X-dephasing.</summary>
    public static bool IsPi2XEvenBilinearTerm(PauliTerm term)
    {
        if (term is null) throw new ArgumentNullException(nameof(term));
        if (term.KBody != 2) return false;
        var nonI = new List<PauliLetter>();
        foreach (var L in term.Letters)
            if (L != PauliLetter.I) nonI.Add(L);
        return IsPi2XEvenBilinear(nonI[0], nonI[1]);
    }

    /// <summary>Returns true iff every term in the given Hamiltonian's term list is
    /// a Π²_X-even 2-body bilinear (per <see cref="IsPi2XEvenBilinearTerm"/>).
    /// F108 Part 2's exact-palindrome guarantee under X-dephasing holds iff this
    /// returns true.</summary>
    public static bool IsPi2XEvenBilinearHamiltonian(IEnumerable<PauliTerm> terms)
    {
        if (terms is null) throw new ArgumentNullException(nameof(terms));
        foreach (var t in terms)
            if (!IsPi2XEvenBilinearTerm(t)) return false;
        return true;
    }

    /// <summary>The theorem statement in one line.</summary>
    public string Theorem =>
        "For H = Σ_b α_b · B_b with B_b ∈ {ZZ, XX, XY, YX, YY} bilinears and α_b ∈ ℝ, " +
        "and X-dephasing on every site: Π_5bilinear · L · Π_5bilinear⁻¹ = −L − 2σ·I exactly " +
        "(with the X-deph variant of Pi5BilinearOperator).";

    /// <summary>The F87 corollary, scoped to X-dephasing.</summary>
    public string F87Corollary =>
        "Under X-dephasing: no Π²_X-even Pauli pair (truly or non-truly) is F87-hard; every such pair has palindromic spec(L).";

    public F108Part2Pi2XEvenAlwaysPalindromic()
        : base("F108 Part 2: Π²_X-even H + X-dephasing always admits exact operator-level palindrome via Π_5bilinear (X-deph variant); BitA twin of F108 Part 1",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F108 + " +
               "docs/proofs/PROOF_F108_PART2_PI2X_EVEN_ALWAYS_PALINDROMIC.md + " +
               "docs/proofs/PROOF_F108_PART1_PI2_EVEN_ALWAYS_PALINDROMIC.md + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi5BilinearOperator.cs + " +
               "experiments/NON_HEISENBERG_PALINDROME.md")
    {
    }

    public override string DisplayName =>
        "F108 Part 2 Π²_X-even always palindromic via Π_5bilinear (X-deph variant, closed-form)";

    public override string Summary =>
        $"{Theorem} {F87Corollary} ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("Theorem", summary: Theorem);
            yield return new InspectableNode("F87 corollary", summary: F87Corollary);
            yield return new InspectableNode("Π_5bilinear per-site map (X-deph variant)",
                summary: "I → +Z, Z → −I, X → −iY, Y → +iX (canonical X-deph permutation I↔Z, X↔Y " +
                         "with Z→I and Y→X back-arrow phase flips); M² = diag(−1, +1, +1, −1) on {I, X, Y, Z}");
            yield return new InspectableNode("Π²_X-even bilinears handled",
                summary: "{ZZ, XX, XY, YX, YY} = 5 of 16 2-site Pauli pairs with bit_a=0");
            yield return new InspectableNode("Empirical verification",
                summary: "9 pure-Π²_X-even non-truly pairs × N=3,4,5 = 27 instances, residual = 0; " +
                         "random non-uniform-J trials at N=3,4,5, residual = 0; " +
                         "pure D[X]^⊗N dissipator N=3,4,5, residual = 0");
            yield return new InspectableNode("Closes F109 X-dephasing branch",
                summary: "F109 (mother sector Klein (0,0) soft y_par=1 purity) had X-dephasing branch " +
                         "empirically anchored only after F108 Part 1; this Claim closes the X branch closed-form. " +
                         "Only the Y-dephasing branch remains empirically anchored (awaits F108 Part 3).");
            yield return new InspectableNode("Sibling on BitB axis",
                summary: "F108 Part 1 (Pi2EvenAlwaysPalindromic, Z-dephasing). Together they cover Z- and X-dephasing " +
                         "branches of the F108 Π²-even palindrome family. Y-dephasing analog (F108 Part 3) is open.");
        }
    }
}
